/*------------------------------------------------------------------------------
* ppp_ar.c : ppp ambiguity resolution
*
* reference :
*    [1] H.Okumura, C-gengo niyoru saishin algorithm jiten (in Japanese),
*        Software Technology, 1991
*
*          Copyright (C) 2012-2015 by T.TAKASU, All rights reserved.
*
* version : $Revision:$ $Date:$
* history : 2013/03/11  1.0  new
*           2016/05/10  1.1  delete codes
            2025/07/29  2.0  AR
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

#define SQR(x)      ((x)*(x))
#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))
#define MAX(x,y)    ((x)>(y)?(x):(y))
#define MIN(x,y)    ((x)<(y)?(x):(y))
#define ROUND(x)    (int)floor((x)+0.5)
#define LOG_PI      1.14472988584940017
#define SQRT2       1.41421356237309510
#define SWAP_I(x,y) do {int _tmp=x; x=y; y=_tmp;} while (0)
#define SWAP_D(x,y) do {double _tmp=x; x=y; y=_tmp;} while (0)
#define MIN_AMB_RES 4         /* min number of ambiguities for ILS-AR */
#define MIN_LOCK_AR 15

/* number and index of ekf states */
#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
#define NP(opt)     ((opt)->dynamics?9:3)
#define NC(opt)     ((opt)->sdopt?0:NSYS)
#define NT(opt)     (((opt)->tropopt<TROPOPT_EST||(opt)->tropopt==TROPOPT_GPT3)?0:((opt)->tropopt==TROPOPT_ESTG?3:1))
#define NI(opt)     ((opt)->ionoopt==IONOOPT_EST?MAXSAT:0)
#define ND(opt)     ((opt)->nf>=3?1:0)
#define NR(opt)     (NP(opt)+NC(opt)+NT(opt)+NI(opt)+ND(opt))
#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1)

static int __attribute__((unused)) test_sys(int sys, int m){
    switch(sys){
        case SYS_GPS: return m==0;
        case SYS_SBS: return m==0;
        case SYS_GLO: return m==1;
        case SYS_GAL: return m==2;
        case SYS_CMP: return m==3;
        case SYS_QZS: return m==4;
    }
    return 0;
}

/* GLONASS frequency channel number (returns INVALID_FCN if not found) */
#define INVALID_FCN 99
static int get_glo_fcn(int sat, const nav_t *nav)
{
    int i, prn;
    if (satsys(sat, &prn) != SYS_GLO || !nav) return INVALID_FCN;
    for (i = 0; i < nav->ng; i++) {
        if (nav->geph[i].sat == sat) return nav->geph[i].frq;
    }
    if (nav->glo_fcn[prn-1] > 0) return nav->glo_fcn[prn-1] - 8;
    return INVALID_FCN;
}

/* Least-squares estimate of GLO IFB rate (cycles/channel) from (Δk, frac) pairs.
 * Model: frac_i = δ * Δk_i + noise  →  δ = Σ(Δk * frac) / Σ(Δk²).
 * Returns 0.0 if there are fewer than 2 distinct Δk values or sum_dk2 == 0. */
static double est_glo_ifb(const int *dk, const double *frac, int n)
{
    double sum_dk2 = 0.0, sum_dk_frac = 0.0;
    int i;
    for (i = 0; i < n; i++) {
        sum_dk2     += (double)dk[i] * dk[i];
        sum_dk_frac += (double)dk[i] * frac[i];
    }
    return (sum_dk2 > 0.0) ? sum_dk_frac / sum_dk2 : 0.0;
}

/* complemetaty error function [1] */
static double q_gamma(double a, double x, double log_gamma_a);
static double p_gamma(double a, double x, double log_gamma_a){
    double y,w;
    int i;

    if (x==0.0) return 0.0;
    if (x>=a+1.0) return 1.0-q_gamma(a,x,log_gamma_a);
    y=w=exp(a*log(x)-x-log_gamma_a)/a;
    for (i=1;i<100;i++) {
        w*=x/(a+i);
        y+=w;
        if (fabs(w)<1E-15) break;
    }
    return y;
}

static double q_gamma(double a, double x, double log_gamma_a){
    double y,w,la=1.0,lb=x+1.0-a,lc;
    int i;

    if (x<a+1.0) return 1.0-p_gamma(a,x,log_gamma_a);
    w=exp(-x+a*log(x)-log_gamma_a);
    y=w/lb;
    for (i=2;i<100;i++) {
        lc=((i-1-a)*(lb-la)+(i+x)*lb)/i;
        la=lb; lb=lc;
        w*=(i-1-a)/i;
        y+=w/la/lb;
        if (fabs(w/la/lb)<1E-15) break;
    }
    return y;
}

static double f_erfc(double x){
    return x>=0.0?q_gamma(0.5,x*x,LOG_PI/2.0):1.0+p_gamma(0.5,x*x,LOG_PI/2.0);
}

/* confidence function of integer ambiguity ----------------------------------*/
static double conf_func(int N, double B, double var){
    double x,p=1.0,sig=sqrt(var);
    int i;

    x=fabs(B-N);
    for (i=1;i<10;i++) {
        p-=f_erfc((i-x)/(SQRT2*sig))-f_erfc((i+x)/(SQRT2*sig));
    }
    return p;
}

static int matchnlupd(const gtime_t obst, int sat1, int sat2, 
                      double *nl_upd1, double *nl_upd2, const nav_t *nav){
    double upd1=0.0,upd2=0.0;
    int i,stat=0;

    for(i=0;i<nav->upds->nls.n;i++){
        if(timediff(nav->upds->nls.data[i].ts,obst)<=0&&timediff(nav->upds->nls.data[i].te,obst)>=0){
            upd1=nav->upds->nls.data[i].nl[sat1-1];
            upd2=nav->upds->nls.data[i].nl[sat2-1];
            stat=1;
            break;
        }
    }

    if(nl_upd1) *nl_upd1=upd1;
    if(nl_upd2) *nl_upd2=upd2;

    return stat;
}

static int __attribute__((unused)) matchnlfcb(const gtime_t obst, int sat1, int sat2,
                      double *nl_fcb1, double *nl_fcb2, const nav_t *nav){
    double fcb1=0.0,fcb2=0.0;
    int i,stat=0;

    for(i=0;i<nav->fcbs->n;i++){
        if(timediff(nav->fcbs->data[i].ts,obst)<=0&&timediff(nav->fcbs->data[i].te,obst)>=0){
            fcb1=nav->fcbs->data[i].bias[sat1-1];
            fcb2=nav->fcbs->data[i].bias[sat2-1];
            stat=1;
            break;
        }
    }

    if(nl_fcb1) *nl_fcb1=fcb1;
    if(nl_fcb2) *nl_fcb2=fcb2;

    return stat;
}

/* Generete SD between sat2 and sat1. Sat2 is ref, sat1 is supportive (unit) */
static int gen_sat_sd(rtk_t *rtk,const nav_t *nav, const obsd_t *obs,
                      int n, const int *exc, int *sat1, int *sat2, int *iu,
                      int *ir, int f, double *el){
    const int sat_sys[]={SYS_GPS,SYS_GLO,SYS_GAL,SYS_CMP,0};
    double elmask,el_temp[MAXOBS]={0};
    int i,j,k,m=0,ns=0,prn,idxs[MAXOBS]={0},sat_no[MAXOBS]={0},frq=f==-1?0:f;
    int ref_sat_idx=0,sys;

    elmask=rtk->opt.elmin;

    for(i=0;sat_sys[i];i++){
        if((sat_sys[i]&SYS_GPS)&&rtk->opt.gpsmodear==ARMODE_OFF) continue;

        m=0; /* reset per-system satellite count */

        for(j=0;j<n;j++){
            sys=satsys(obs[j].sat,&prn);
            if(!rtk->ssat[obs[j].sat-1].vs) continue;
        

            if(sys!=sat_sys[i]) continue;

            if(rtk->ssat[obs[j].sat-1].slip[0]||rtk->ssat[obs[j].sat-1].slip[1]){
                rtk->ssat[obs[j].sat-1].lock[f]=-rtk->opt.minlock;
                continue;
            }

            if(exc[j]||!rtk->ssat[obs[j].sat-1].vsat[frq]||rtk->ssat[obs[j].sat-1].azel[1]<elmask
                ||rtk->ssat[obs[j].sat-1].lock[f]<0){continue;}

            if (rtk->ssat[obs[j].sat-1].lock[f] < MIN_LOCK_AR) {
                continue;
            }
#if 0
        int iamb=0;
        iamb=rtk->tc?xiAmb(&rtk->opt.insopt,obs[j].sat,0):IB(obs[j].sat,0,&rtk->opt);
        if(SQRT(rtk->P[iamb+iamb*rtk->nx])>0.30){continue;}
#endif
            sat_no[m]=obs[j].sat;
            el_temp[m]=rtk->ssat[obs[j].sat-1].azel[1];
            idxs[m]=j;
            m++;
        }

        if(m<=0) continue;

        for(j=0;j<m-1;j++) for(k=j+1;k<m;k++){
            if(el_temp[j]>=el_temp[k]) continue;
            SWAP_I(sat_no[j],sat_no[k]);
            SWAP_I(idxs[j],idxs[k]);
            SWAP_D(el_temp[j],el_temp[k]);
        }

        /* generate SD referenced to max elecation angel */
        ref_sat_idx=0;
        for(j=0;j<m;j++){
            if(j==0) continue;
            sat1[ns]=sat_no[j];
            iu[ns]=idxs[j];
            el[ns]=el_temp[j];
            sat2[ns]=sat_no[ref_sat_idx];
            ir[ns]=idxs[ref_sat_idx];
            ns++;
        }
    }
    return ns;
}

static void resetamb(rtk_t *rtk, double *xa, const double *Bc, const int *sat1, const int *sat2, int nb){
    int i,iamb,jamb;
    for(i=0;i<nb;i++){
        iamb=IB(sat1[i],0,&rtk->opt);
        jamb=IB(sat2[i],0,&rtk->opt);
        xa[iamb]=Bc[i]+rtk->x[jamb];
    }
}

/* Build design matrix for ambiguity WL,NL,IF */
static int SDmat(rtk_t *rtk,const obsd_t *obs,int ns,const nav_t *nav,double *H_nl,double *H_if,int *sat1,
        int *sat2,int *iu,int *ir,double *el,double *Nw,double *Bw,double *Nl,double *Nc,double *sd_nl_fcb){
    prcopt_t opt=rtk->opt;
    int i,j,f2,sat,ref_sat,prn,sys_idx=-1,iamb,jamb,nb=0;
    double frq1,frq2,lam1,lam2,lam_nl,sd_if;

    /* clear fix flag for all sats (1=float, 2=fix) */
    for(i=0;i<MAXSAT;i++) for(j=0;j<NFREQ;j++){
        rtk->ssat[i].fix[j]=0;
    }

    for(i=0;i<ns;i++){
        sat=sat1[i];
        ref_sat=sat2[i];
        satsys(sat,&prn);
        sys_idx=satsysidx(sat);
        if(sys_idx==-1) continue;
        if(rtk->sdamb[sat-1].fix_nl_flag!=1) continue;

        f2=obs[iu[i]].L[1]==0.0?2:1;

        frq1 = sat2freq(sat, obs[iu[i]].code[0], nav);
        frq2 = sat2freq(sat, obs[iu[i]].code[f2], nav);

        lam1=CLIGHT/frq1;
        lam2=CLIGHT/frq2;
        lam_nl=lam1*lam2/(lam1+lam2);

        iamb=IB(sat,0,&opt);
        jamb=IB(ref_sat,0,&opt);
        sd_if=rtk->x[iamb]-rtk->x[jamb];

        H_nl[iamb+nb*rtk->nx]=1.0/lam_nl;
        H_nl[jamb+nb*rtk->nx]=-1.0/lam_nl;
        H_if[iamb+(rtk->na+nb)*rtk->nx]=1.0;
        H_if[jamb+(rtk->na+nb)*rtk->nx]=-1.0;

        sat1[nb]=sat;
        sat2[nb]=ref_sat;
        iu[nb]=iu[i];
        ir[nb]=ir[i];
        el[nb]=el[i];
        Nw[nb]=rtk->sdamb[sat-1].wl;
        Bw[nb]=newround(rtk->sdamb[sat-1].wl);
        Nl[nb]=rtk->sdamb[sat-1].nl;
        Nc[nb]=sd_if;
        
        if(opt.arprod == AR_PROD_UPD){
            double nl_fcb1=0.0,nl_fcb2=0.0;
            matchnlupd(obs[i].time,sat,ref_sat,&nl_fcb1,&nl_fcb2,nav);

            sd_nl_fcb[nb]=nl_fcb1-nl_fcb2;
        }
        trace(2, "SAT: %d | IF: %f, WL: %f, NL: %f\n\r",sat-1, sd_if, rtk->sdamb[sat-1].wl, rtk->sdamb[sat-1].nl);

        rtk->ssat[sat1[nb]-1].fix[0]=2;
        rtk->ssat[sat2[nb]-1].fix[0]=2;
        nb++;
    }
    return nb;
}

/* Narrow-lane ambiguity resolution */
static int resamb_nl(rtk_t *rtk,double *H_nl,double *nl_amb,int num_nl)
{
    int nb=num_nl,stat=0;
    double *Qnl,*HP,s[2]={0},*b,*pb;

    Qnl=mat(nb,nb);HP=mat(nb,rtk->nx);b=mat(nb,2);pb=zeros(nb,2);
    matmul("TN",nb,rtk->nx,rtk->nx,1.0,H_nl,rtk->P,0.0,HP);
    matmul("NN",nb,nb,rtk->nx,1.0,HP,H_nl,0.0,Qnl);

    if(!lambda(nb,2,nl_amb,Qnl,b,s)){

        rtk->sol.ratio=s[0]>0?(float)(s[1]/s[0]):0.0f;
        if (rtk->sol.ratio>999.9) rtk->sol.ratio=999.9f;
        /* NL validation uses the LAMBDA ratio (s1/s0, always >=1), so the
         * threshold must be the ratio threshold thresar[0] (e.g. 3.0), NOT
         * thresar[1] which is a rounding-confidence probability (~0.9999). */
        rtk->sol.thres=(float)rtk->opt.thresar[0];

        if(rtk->sol.ratio<rtk->sol.thres&&nb>MIN_AMB_RES){
            stat=0;
        }
        else if(rtk->sol.ratio>rtk->sol.thres){
            matcpy(nl_amb,b,nb,1);
            stat=1;
        }
        else stat=0;
    }
    else stat=0;

    free(Qnl);free(HP);free(b);free(pb);

    return stat?nb:0;
}

static int fix_sol(rtk_t *rtk,const obsd_t *obs,const nav_t *nav,const double *sd_nl_fcb,double *H_if,
                   const double *Bl,const double *Bw,int nb,const int *sat1,const int *sat2,const int *iu,double *xa){
    prcopt_t opt=rtk->opt;
    int i,j,f2,ny,na=rtk->na,sat,sys_idx=-1,prn,stat=1;
    double *y,*db,*Qb_if,*Qab,*QQ,*Qy,*DP,*Bc,*dx;
    double frq1=0.0,frq2=0.0,lam1,lam2,lam_nl,gamma;

    ny=na+nb;y=mat(ny,1);
    db=mat(nb,1);Qb_if=mat(nb,nb);Qab=mat(na,nb);
    QQ=mat(na,nb);Qy=mat(ny,ny);DP=mat(ny,rtk->nx);
    Bc=mat(nb,1);
    dx=mat(na,1);
    matmul("TN",ny,1,rtk->nx,1.0,H_if,rtk->x,0.0,y);
    matmul("TN",ny,rtk->nx,rtk->nx,1.0,H_if,rtk->P,0.0,DP);
    matmul("NN",ny,ny,rtk->nx,1.0,DP,H_if,0.0,Qy);

    for(i=0;i<nb;i++) for(j=0;j<nb;j++){
        Qb_if[i+j*nb]=Qy[na+i+(na+j)*ny]; /* SD IF ambiguity covariance */
    }

    for(i=0;i<na;i++) for(j=0;j<nb;j++) Qab[i+j*na]=Qy[i+(na+j)*ny];

    for(i=0;i<na;i++){
        rtk->xa[i]=xa[i];
        for(j=0;j<na;j++) rtk->Pa[i+j*na]=rtk->P[i+j*rtk->nx];
    }

    /* fix IF ambiguity */
    for(i=0;i<nb;i++){
        sat=sat1[i];
        satsys(sat,&prn);
        sys_idx=satsysidx(sat);
        if(sys_idx==-1) continue;
        f2=obs[iu[i]].L[1]==0.0?2:1;

        frq1 = sat2freq(sat, obs[iu[i]].code[0], nav);
        frq2 = sat2freq(sat, obs[iu[i]].code[f2], nav);
        lam1=CLIGHT/frq1;
        lam2=CLIGHT/frq2;
        lam_nl=lam1*lam2/(lam2+lam1);
        gamma=CLIGHT*frq2/(SQR(frq1)-SQR(frq2));

        if(opt.arprod == AR_PROD_OSB_COD){
            /* Bl[i] holds the LAMBDA-fixed SD-NL integer for every entry
             * 0..nb-1; a fixed value of 0 is a valid integer, so reconstruct
             * the IF bias unconditionally (do not fall back to the float). */
            Bc[i]=lam_nl*Bl[i]+gamma*Bw[i];

            rtk->sdamb[sat1[i]-1].lc_fix=Bc[i];
            rtk->sdamb[sat1[i]-1].lc_res=Bc[i]-rtk->sdamb[sat1[i]-1].lc;
        } else if(opt.arprod == AR_PROD_UPD || opt.arprod == AR_PROD_FCB){
            Bc[i]=lam_nl*(Bl[i]+sd_nl_fcb[i])+gamma*Bw[i];
            rtk->sdamb[sat1[i]-1].lc_fix=Bc[i];
        }
    }
    /* y=differences between float and fixed ionospheric-free bias */
    for(i=0;i<nb;i++){
        y[na+i]-=Bc[i];
    }

    /* adjuset non phase-bias states and covariance using fixed slution */
    if(!matinv(Qb_if,nb)){
        /* rtk->Pa=rtk->P-Qab*Qb^-1*Qab') */
        matmul("NN",na,nb,nb, 1.0,Qab,Qb_if ,0.0,QQ);  /* QQ = Qab*Qb^-1 */
        matmul("NT",na,na,nb,-1.0,QQ ,Qab,1.0,rtk->Pa); /* rtk->Pa = rtk->P-QQ*Qab' */

        /* rtk->xa = rtk->x-Qab*Qb^-1*(b0-b) */
        matmul("NN",nb,1,nb, 1.0,Qb_if ,y+na,0.0,db); /* db = Qb^-1*(b0-b) */
        matmul("NN",na,1,nb,-1.0,Qab,db  ,0.0,dx); /* rtk->xa = rtk->x-Qab*db */

        for(i=0;i<na;i++){
            rtk->xa[i]+=dx[i];
        }
        for(i=0;i<na;i++){
            xa[i]=rtk->xa[i];
        }

        /* reset ambiguity */
        resetamb(rtk,xa,Bc,sat1,sat2,nb);

        if(nb>4){
            rtk->holdamb=1;
        }
    } else stat=0;
    free(y);free(db);free(Qb_if);free(Qab);free(QQ);free(Qy);free(DP);free(Bc);free(dx);

    return stat;
}

static int pppar_IF_ILS(rtk_t *rtk,double *xa,double *bias, const obsd_t *obs,
                        int n, const int *exc,const nav_t *nav){
    prcopt_t opt = rtk->opt;
    int ns=0,nb=0;
    int i,f2,sat,prn,ref_sat,sys,sys_idx=-1,sat1[MAXOBS]={0},sat2[MAXOBS]={0},iu[MAXOBS]={0},ir[MAXOBS]={0},na=rtk->na,stat=0;
    double frq1=0.0,frq2=0.0,lam1,lam2,lam_nl,gamma,el[MAXOBS]={0};
    double wl_amb,sd_wl,wl_fcb,wl_var,sd_if,nl_amb,nl_fcb1,nl_fcb2;
    int iamb,jamb,k1,k2;
    /* GLONASS IFB rates estimated below (cycles per channel unit) */
    double glo_ifb_wl = 0.0;  /* WL IFB rate */
    double glo_ifb_nl = 0.0;  /* NL IFB rate */
    double *H_nl,*H_if;
    double *Nw,*Nl,*Nc; /* float ambiguity */
    double *Bw,*Bl,*Bc; /* float ambiguity */
    double *sd_nl_fcb,*Qnl;

    /* generate satellite SD */
    if(!(ns=gen_sat_sd(rtk,nav,obs,n,exc,sat1,sat2,iu,ir,0,el))) return 0;

    for(i=0;i<MAXSAT;i++){
        rtk->sdamb[i].nl=0.0;
        rtk->sdamb[i].wl=rtk->sdamb[i].nl_res=rtk->sdamb[i].wl_res=0.0;
        rtk->sdamb[i].fix_nl_flag=rtk->sdamb[i].fix_wl_flag=0;
        rtk->sdamb[i].nl_fix=0.0;
        rtk->sdamb[i].wl_fix=0;
        rtk->sdamb[i].ref_sat_no=0;
    }
    
    H_nl=zeros(rtk->nx,ns);
    H_if=zeros(rtk->nx,rtk->na+ns);

    for(i=0;i<na;i++) H_if[i+i*rtk->nx]=1.0;

    Nw=zeros(ns,1);
    Nl=zeros(ns,1);
    Nc=zeros(ns,1);
    Bw=zeros(ns,1);
    Bl=zeros(ns,1);
    Bc=zeros(ns,1);
    sd_nl_fcb=zeros(ns,1);
    Qnl=zeros(ns,ns);

    /* ================================================================
     * GLONASS IFB pre-estimation (two-step least-squares).
     *
     * GLONASS uses FDMA: each satellite has a unique frequency channel
     * k ∈ [-7,+6].  For a between-satellite SD, the receiver hardware
     * inter-frequency bias (IFB) does NOT cancel but leaves a residual
     * proportional to the channel difference Δk = k_i − k_j:
     *   WL_ij = N_WL_ij + Δk * δ_WL           (WL IFB)
     *   NL_ij = N_NL_ij + Δk * δ_NL           (NL IFB)
     *
     * Step 1: estimate δ_WL from fractional parts of raw MW-smoothed WL.
     * Step 2: estimate δ_NL from fractional parts of NL computed with
     *         the corrected (integer-rounded) WL from step 1.
     * Both use simple ridge-less LS:  δ = Σ(Δk · frac) / Σ(Δk²).
     * ================================================================ */
    if (opt.navsys & SYS_GLO) {
        int    glo_dk[MAXOBS];
        double glo_wl_frac[MAXOBS], glo_nl_frac[MAXOBS];
        int    n_glo_wl = 0, n_glo_nl = 0;

        /* Use separate dk arrays for WL and NL (NL is a subset of WL pairs) */
        int glo_dk_nl[MAXOBS];

        for (i = 0; i < ns; i++) {
            double wl_a, wl_fcb_g, sd_wl_g, frq1_g, frq2_g, lam1_g, lam2_g, gamma_g, lam_nl_g;
            double sd_if_g, raw_nl, nl_frac;
            int k1, k2, dk, iamb, jamb;

            sat     = sat1[i];
            ref_sat = sat2[i];
            if (satsys(sat, NULL) != SYS_GLO) continue;

            k1 = get_glo_fcn(sat,     nav);
            k2 = get_glo_fcn(ref_sat, nav);
            if (k1 == INVALID_FCN || k2 == INVALID_FCN) continue;
            dk = k1 - k2;

            /* require enough MW-smoothed epochs for reliable WL */
            if (rtk->ssat[sat-1].mw[2] < 10 || rtk->ssat[ref_sat-1].mw[2] < 10) continue;

            f2 = obs[iu[i]].L[1] == 0.0 ? 2 : 1;
            frq1_g = sat2freq(sat, obs[iu[i]].code[0],  nav);
            frq2_g = sat2freq(sat, obs[iu[i]].code[f2], nav);
            if (frq1_g == 0.0 || frq2_g == 0.0) continue;

            lam1_g  = CLIGHT / frq1_g;
            lam2_g  = CLIGHT / frq2_g;
            lam_nl_g = lam1_g * lam2_g / (lam1_g + lam2_g);
            gamma_g  = CLIGHT * frq2_g / (SQR(frq1_g) - SQR(frq2_g));

            sd_wl_g = rtk->ssat[sat-1].mw[1] - rtk->ssat[ref_sat-1].mw[1];
            if (opt.arprod == AR_PROD_UPD && nav->upds) {
                wl_fcb_g = nav->upds->wls.wl[sat-1] - nav->upds->wls.wl[ref_sat-1];
            } else {
                wl_fcb_g = nav->wlbias[sat-1] - nav->wlbias[ref_sat-1];
            }
            wl_a = (opt.arprod == AR_PROD_OSB_COD) ? sd_wl_g : sd_wl_g - wl_fcb_g;

            /* collect WL fractional part for Step-1 IFB estimation */
            glo_dk[n_glo_wl]      = dk;
            glo_wl_frac[n_glo_wl] = wl_a - floor(wl_a + 0.5);
            n_glo_wl++;

            /* Step-2: NL IFB — only usable if WL is close enough to an integer */
            if (fabs(wl_a - floor(wl_a + 0.5)) > 0.35) continue;

            iamb    = IB(sat,     0, &opt);
            jamb    = IB(ref_sat, 0, &opt);
            sd_if_g = rtk->x[iamb] - rtk->x[jamb];
            raw_nl  = (sd_if_g - gamma_g * floor(wl_a + 0.5)) / lam_nl_g;
            nl_frac = raw_nl - floor(raw_nl + 0.5);

            /* store dk separately for NL pairs (subset of WL pairs) */
            glo_dk_nl[n_glo_nl]    = dk;
            glo_nl_frac[n_glo_nl]  = nl_frac;
            n_glo_nl++;
        }

        if (n_glo_wl >= 2) {
            glo_ifb_wl = est_glo_ifb(glo_dk,    glo_wl_frac, n_glo_wl);
            trace(2, "GLO IFB WL: %.4f cyc/ch  (n=%d)\n\r", glo_ifb_wl, n_glo_wl);
        }
        if (n_glo_nl >= 2) {
            glo_ifb_nl = est_glo_ifb(glo_dk_nl, glo_nl_frac, n_glo_nl);
            trace(2, "GLO IFB NL: %.4f cyc/ch  (n=%d)\n\r", glo_ifb_nl, n_glo_nl);
        }
    }

    /* generate WL-IF-NL */
    for(i=0;i<ns;i++){
        sat=sat1[i];
        ref_sat=sat2[i];
        sys=satsys(sat,&prn);
        sys_idx=satsysidx(sat);
        if(sys_idx==-1) continue;

        rtk->sdamb[ref_sat-1].ref_sat_no=0;

        f2=obs[iu[i]].L[1]==0.0?2:1;

        frq1 = sat2freq(sat, obs[iu[i]].code[0], nav); 
        frq2 = sat2freq(sat, obs[iu[i]].code[f2], nav);
        lam1=CLIGHT/frq1;
        lam2=CLIGHT/frq2;
        lam_nl=lam1*lam2/(lam2+lam1);
        gamma=CLIGHT*frq2/(SQR(frq1)-SQR(frq2));

        /* round WL ambiguity */
        wl_amb=0.0;
        sd_wl=rtk->ssat[sat-1].mw[1]-rtk->ssat[ref_sat-1].mw[1];
        wl_fcb=0.0;
        if(opt.arprod == AR_PROD_UPD){
            wl_fcb=nav->upds->wls.wl[sat-1]-nav->upds->wls.wl[ref_sat-1];
        } else {
            wl_fcb=nav->wlbias[sat-1]-nav->wlbias[ref_sat-1];
        }

        wl_var=rtk->ssat[sat-1].mw[3] + rtk->ssat[ref_sat-1].mw[3];

        if(opt.arprod == AR_PROD_FCB || opt.arprod == AR_PROD_UPD){
            /* FCB/UPD products contain WL biases computed for the L1+L2
             * combination.  When the satellite is tracked on L1+L5 (f2==2),
             * those biases are wrong — skip them to avoid corrupting WL. */
            wl_amb = (f2==2) ? sd_wl : sd_wl-wl_fcb;
        } else if(opt.arprod == AR_PROD_OSB_COD){
            wl_amb = sd_wl;  /* OSBs already applied per-signal in mwmeas() */
        }

        /* GLONASS FDMA: apply receiver WL IFB correction (cyc/channel) */
        sys = satsys(sat, NULL);
        if (sys == SYS_GLO && glo_ifb_wl != 0.0) {
            k1 = get_glo_fcn(sat,     nav);
            k2 = get_glo_fcn(ref_sat, nav);
            if (k1 != INVALID_FCN && k2 != INVALID_FCN)
                wl_amb -= (k1 - k2) * glo_ifb_wl;
        }

        rtk->sdamb[sat-1].wl=wl_amb;
        rtk->sdamb[sat-1].wl_fix=newround(wl_amb);
        rtk->sdamb[sat-1].wl_res=wl_amb-newround(wl_amb);

        if(conf_func(newround(wl_amb), wl_amb, wl_var) < rtk->opt.thresar[1] || fabs(newround(wl_amb)-wl_amb) > 0.25){
            rtk->sdamb[sat-1].fix_wl_flag=0;
            continue;
        }
        trace(2,"SAT: %d; WL: %f < TRE: %f\n\r", i, conf_func(newround(wl_amb), wl_amb, wl_var), rtk->opt.thresar[1]);
        rtk->sdamb[sat-1].fix_wl_flag=1;
        rtk->sdamb[sat-1].ref_sat_no=ref_sat;

        trace(2,"SAT: %d; FIX: %d REF: %d\n\r", i, rtk->sdamb[sat-1].fix_wl_flag, rtk->sdamb[sat-1].ref_sat_no);
        /* float if ambiguity */
        iamb=IB(sat,0,&opt);
        jamb=IB(ref_sat,0,&opt);
        sd_if = rtk->x[iamb]-rtk->x[jamb];
        trace(2, "SAT: %d; SD IF: %f\n\r", i, sd_if);

        /* nl ambiguity */
        nl_amb=0.0;
        nl_fcb1=0.0; nl_fcb2=0.0;

        if(opt.arprod == AR_PROD_OSB_COD){
            nl_amb = (sd_if-gamma*newround(wl_amb))/lam_nl;
            trace(2, "SAT: %d; NL: %f\n\r", i, nl_amb);
        } else if(opt.arprod == AR_PROD_UPD&&f2!=2&&
                  !matchnlupd(obs[i].time,sat,ref_sat,&nl_fcb1,&nl_fcb2,nav)){
            /* UPD NL FCBs are computed for L1+L2; skip for L1+L5 (f2==2) */
            nl_amb=(sd_if-gamma*newround(wl_amb))/lam_nl-(nl_fcb1-nl_fcb2);
        } else {
            /* L1+L5 with FCB/UPD: OSBs applied in corr_meas(), no separate FCB */
            nl_amb=(sd_if-gamma*newround(wl_amb))/lam_nl;
        }
        if(isnan(nl_amb)) nl_amb=0.0;

        /* GLONASS FDMA: apply receiver NL IFB correction (cyc/channel) */
        if (sys == SYS_GLO && glo_ifb_nl != 0.0) {
            k1 = get_glo_fcn(sat,     nav);
            k2 = get_glo_fcn(ref_sat, nav);
            if (k1 != INVALID_FCN && k2 != INVALID_FCN)
                nl_amb -= (k1 - k2) * glo_ifb_nl;
        }

        rtk->sdamb[sat-1].nl=nl_amb;
        rtk->sdamb[sat-1].nl_fix=newround(nl_amb);
        rtk->sdamb[sat-1].nl_res=nl_amb-newround(nl_amb);
        rtk->sdamb[sat-1].lc=sd_if;
        rtk->sdamb[sat-1].fix_nl_flag=1;
        trace(2, "SAT: %d; NL: %f; NL_FIX: %d; NL_RES: %f; SD_IF: %f\n\r", i, nl_amb, newround(nl_amb), nl_amb-newround(nl_amb), sd_if);
        
    }

    nb=SDmat(rtk,obs,ns,nav,H_nl,H_if,sat1,sat2,iu,ir,el,Nw,Bw,Nl,Nc,sd_nl_fcb);
    rtk->nb_ar=nb;

    if(nb >= MIN_AMB_RES){
        trace(2, "NB: %d || MIN_AMB %d\n\r", nb, MIN_AMB_RES);
        nb=resamb_nl(rtk,H_nl,Nl,nb);
        trace(2, "POST NB: %d || MIN_AMB %d\n\r", nb, MIN_AMB_RES);
        if(nb&&fix_sol(rtk,obs,nav,sd_nl_fcb,H_if,Nl,Bw,nb,sat1,sat2,iu,xa)){
            stat=1;
        }
    }
    trace(2,"STAT: %d\n\r", stat);
    free(H_nl);free(H_if);
    free(Nw);free(Nl);free(Nc);
    free(Bw);free(Bl);free(Bc);
    free(sd_nl_fcb);free(Qnl);

    return stat?nb:0;
}

static int arfilter(rtk_t *rtk,const obsd_t *obs,int ns,int nf)
{
    int rerun=0,dly=0,i,f;
    dly=2;

    for(i=0;i<ns;i++){
        for(f=0;f<nf;f++){
            if (rtk->ssat[obs[i].sat-1].fix[f]!=2) continue;
            if(rtk->ssat[obs[i].sat-1].lock[f]==0){
                rtk->ssat[obs[i].sat-1].lock[f]=-rtk->opt.minlock-dly;
                dly+=2;
                rerun=1;
            }
        }
    }
    return rerun;
}

/* hold integer ambiguity (fix-and-hold for PPP) -------------------------------
 * Applies SD-IF ambiguity pseudo-measurements to the float Kalman filter.
 * For each fixed satellite pair (sat, ref), the constraint forces the float
 * SD-IF ambiguity to stay near the previously resolved integer value.
 * Variance is taken from opt.varholdamb (cycle^2).
 *-----------------------------------------------------------------------------*/
extern void holdamb_ppp(rtk_t *rtk, const double *xa)
{
    double *v, *H, *R;
    int i, iamb, jamb, nv=0, nb=rtk->nx-rtk->na, ref, info;
    double var_hold;

    trace(3,"holdamb_ppp: na=%d nx=%d\n",rtk->na,rtk->nx);

    var_hold = rtk->opt.varholdamb>0.0 ? rtk->opt.varholdamb : 0.001;

    v=mat(nb,1);
    H=zeros(nb,rtk->nx);

    for (i=0;i<MAXSAT;i++) {
        if (rtk->ssat[i].fix[0]!=2) continue;
        if (rtk->sdamb[i].ref_sat_no<=0) continue;
        if (rtk->ssat[i].azel[1]<rtk->opt.elmaskhold) continue;

        ref=rtk->sdamb[i].ref_sat_no-1;
        iamb=IB(i+1,0,&rtk->opt);
        jamb=IB(ref+1,0,&rtk->opt);

        /* residual: fixed SD minus float SD */
        v[nv] = rtk->sdamb[i].lc_fix - (rtk->x[iamb] - rtk->x[jamb]);
        H[iamb+nv*rtk->nx]= 1.0;
        H[jamb+nv*rtk->nx]=-1.0;
        rtk->ssat[i].fix[0]=3; /* mark as held */
        nv++;
    }
    if (nv>0) {
        R=zeros(nv,nv);
        for (i=0;i<nv;i++) R[i+i*nv]=var_hold;

        if ((info=filter(rtk->x,rtk->P,H,v,R,rtk->nx,nv))) {
            trace(2,"holdamb_ppp: filter error info=%d\n",info);
        }
        rtk->holdamb=1;
        free(R);
    }
    free(v); free(H);
}

/* ambiguity resolution in ppp -----------------------------------------------*/
extern int ppp_ar(rtk_t *rtk,double *bias, double *xa,double *Pa,int nf, const obsd_t *obs,int ns,const nav_t *nav, int *exc)
{
    int i,ipos=0,npos=3,nb=0;
    float ratio_post=0.0;
    double var=0.0;
    prcopt_t opt=rtk->opt;
    rtk->sol.ratio=0.0;

    if(opt.thresar[0]<1.0){
        rtk->nb_ar=0;
        return 0;
    }

    /* skip AR if position variance too high to avoid false fix */
    for(i=ipos;i<ipos+npos;i++) var+=SQRT(rtk->P[i+i*rtk->nx]);
    var=var/3.0; /* maintain compatibility with previous code */
    if(var>opt.thresar[2]){
        trace(2,"position variance too large, var=%7.3f thres=%7.3f\n", var,opt.thresar[2]);
        return 0;
    }

    if(rtk->opt.ionoopt==IONOOPT_IFLC){
        nb=pppar_IF_ILS(rtk,xa,bias,obs,ns,exc,nav);
    }

    ratio_post=rtk->sol.ratio;
    if(nb>=0){
        if(arfilter(rtk,obs,ns,nf)) nb=pppar_IF_ILS(rtk,xa,bias,obs,ns,exc,nav);
    }
    rtk->sol.prev_ratio=ratio_post>0?ratio_post:rtk->sol.ratio;
    rtk->sol.prevf_ratio=rtk->sol.ratio;

    return nb;
}

extern int manage_ppp_ar(rtk_t *rtk,double *bias,double *xa,double *Pa,int nf,const obsd_t *obs,int ns,const nav_t *nav,int *exc){
    int i,f,lockc[NFREQ],ar=0,excflag=0,arsats[MAXOBS]={0},sat,nb=0;

    /* if no FIX on previous sample and enough sats, exclude next sat in list */
    if(rtk->sol.prevf_ratio<rtk->sol.thres&&rtk->nb_ar>=rtk->opt.mindropsats){
        for(f=0;f<nf;f++) for(i=0;i<ns;i++){
            sat=obs[i].sat;
            if(rtk->ssat[sat-1].vsat[f]&&rtk->ssat[sat-1].lock[f]>=0&&rtk->ssat[sat-1].azel[1]>=rtk->opt.elmaskar){
                arsats[ar++]=i;
            }
        }
        if(rtk->excsats<ar){
            sat=obs[arsats[rtk->excsats]].sat;
            for(f=0;f<nf;f++){
                lockc[f]=rtk->ssat[sat-1].lock[f];
                rtk->ssat[sat-1].lock[f]=-rtk->nb_ar;
            }
            excflag=1;
        }
        else rtk->excsats=0; /* exclude none and reset to beginning of list */
    }
    
    nb=ppp_ar(rtk,bias,xa,Pa,nf,obs,ns,nav,exc);


    /* restore exclude sat if still no fix or significant increase in ar ration */
    if(excflag&&(rtk->sol.ratio<rtk->sol.thres)&&(rtk->sol.ratio<(1.5*rtk->sol.prevf_ratio))){
        sat=obs[arsats[rtk->excsats++]].sat;
        for(f=0;f<nf;f++) rtk->ssat[sat-1].lock[f]=lockc[f];
    }

    return nb;
}