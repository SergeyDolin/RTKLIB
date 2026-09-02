/*------------------------------------------------------------------------------
* postls.c : smartphone-oriented relative post-processing by batch LS
*
*          Copyright (C) 2026, all rights reserved.
*
* notes    : This is the first implementation slice for a Bernese-like
*            post-processor. It estimates a static float double-difference
*            solution from a whole RINEX session. Integer ambiguity resolution,
*            arc repair, and kinematic knots are intentionally layered later.
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

#define SQR(x)      ((x)*(x))
#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))
#define MIN_(x,y)   ((x)<(y)?(x):(y))
#define MAX_(x,y)   ((x)>(y)?(x):(y))

#define POSTLS_NSYS 6
#define POSTLS_MAXITR 8
#define POSTLS_MINOBS 8
#define POSTLS_MAXROB 4
#define POSTLS_REJSTD 10.0
#define POSTLS_MAXAMB 4096
#define POSTLS_MINARCOBS 5
#define POSTLS_MINPHASEOBS 4
#define POSTLS_MAXFIXAMB 12
#define POSTLS_ARC_REJSTD 6.0
#define POSTLS_MAXBADARC 4

#define SMARTPHONE(opt) ((opt)->nf==6)

typedef struct {
    int sat,f,aur,abr,arr,arb,aur2,abr2,arr2,arb2,index,n,use,bad;
} postls_amb_t;
typedef struct {
    int iu,nu,ir,nr,index,nv,ns;
    gtime_t time;
    double ru[3];
} postls_epoch_t;

static int postls_nf(const prcopt_t *opt)
{
    int nf=SMARTPHONE(opt)?2:opt->nf;
    return MIN_(MAX_(nf,1),NFREQ);
}
static int postls_obsidx(const prcopt_t *opt, int f)
{
    return SMARTPHONE(opt)&&f==1?2:f;
}
static int postls_sysidx(int sys)
{
    switch (sys) {
        case SYS_GPS: return 0;
        case SYS_SBS: return 0;
        case SYS_GLO: return 1;
        case SYS_GAL: return 2;
        case SYS_CMP: return 3;
        case SYS_QZS: return 4;
        case SYS_IRN: return 5;
    }
    return -1;
}
static int postls_iflc(const prcopt_t *opt)
{
    return SMARTPHONE(opt)&&opt->ionoopt==IONOOPT_IFLC&&postls_nf(opt)>=2;
}
static int nextobsf_postls(const obs_t *obs, int *i, int rcv)
{
    double tt;
    int n;

    for (;*i<obs->n;(*i)++) if (obs->data[*i].rcv==rcv) break;
    for (n=0;*i+n<obs->n;n++) {
        tt=timediff(obs->data[*i+n].time,obs->data[*i].time);
        if (obs->data[*i+n].rcv!=rcv||tt>DTTOL) break;
    }
    return n;
}
static int find_sat_obs(const obsd_t *obs, int n, int sat)
{
    int i;

    for (i=0;i<n;i++) if (obs[i].sat==sat) return i;
    return -1;
}
static int valid_sig(const obsd_t *obs, int p, int phase)
{
    if (p<0||p>=NFREQ+NEXOBS||obs->code[p]==CODE_NONE) return 0;
    if (phase&&(obs->LLI[p]&LLI_HALFC)) return 0;
    return phase?obs->L[p]!=0.0:obs->P[p]!=0.0;
}
static int build_phase_arcs(const obs_t *obs, const prcopt_t *opt, int nf,
                            int *arc)
{
    gtime_t last_time[2][MAXSAT][NFREQ];
    double last_L[2][MAXSAT][NFREQ],gapmax;
    int active[2][MAXSAT][NFREQ]={{{0}}},cur[2][MAXSAT][NFREQ]={{{0}}};
    int i,j,f,p,rcv,sat,narc=0;

    for (i=0;i<obs->n*nf;i++) arc[i]=0;
    memset(last_time,0,sizeof(last_time));
    memset(last_L,0,sizeof(last_L));
    gapmax=MAX_(opt->maxtdiff>0.0?opt->maxtdiff*2.0:60.0,60.0);

    for (i=0;i<obs->n;i++) {
        rcv=obs->data[i].rcv;
        if (rcv<1||rcv>2) continue;
        j=rcv-1;
        sat=obs->data[i].sat;
        if (sat<=0||sat>MAXSAT) continue;
        for (f=0;f<nf;f++) {
            p=postls_obsidx(opt,f);
            if (!valid_sig(obs->data+i,p,1)) {
                arc[i*nf+f]=0;
                continue;
            }
            if (!active[j][sat-1][f]||
                (obs->data[i].LLI[p]&LLI_SLIP)||
                fabs(timediff(obs->data[i].time,last_time[j][sat-1][f]))>gapmax||
                fabs(obs->data[i].L[p]-last_L[j][sat-1][f])>1E5) {
                cur[j][sat-1][f]=++narc;
                active[j][sat-1][f]=1;
            }
            arc[i*nf+f]=cur[j][sat-1][f];
            last_time[j][sat-1][f]=obs->data[i].time;
            last_L[j][sat-1][f]=obs->data[i].L[p];
        }
    }
    return narc;
}
static double obs_sigma(const obsd_t *obs, const nav_t *nav,
                        const prcopt_t *opt, int f, int phase, double el)
{
    double sig,cn0,cnfact=1.0,reported=0.0,freq,sinel;
    int p=postls_obsidx(opt,f);

    sinel=MAX_(sin(el),0.10);
    sig=opt->err[1]+opt->err[2]/sinel;
    if (!phase) sig*=opt->eratio[f]>0.0?opt->eratio[f]:opt->eratio[0];
    if (sig<=0.0) sig=phase?0.01:3.0;

    cn0=obs->SNR[p]*SNR_UNIT;
    if (SMARTPHONE(opt)&&cn0>0.0) {
        cnfact=pow(10.0,MAX_(0.0,(f==1?35.0:40.0)-cn0)/20.0);
        cnfact=MIN_(cnfact,10.0);
        sig*=cnfact;
    }
    if (phase&&obs->Lstd[p]>0.0&&
        (freq=sat2freq(obs->sat,obs->code[p],nav))>0.0) {
        reported=obs->Lstd[p]*CLIGHT/freq;
    }
    else if (!phase&&obs->Pstd[p]>0.0) {
        reported=obs->Pstd[p];
    }
    return MAX_(sig,reported);
}
static int modeled_obs(const obsd_t *obs, const double *rs, const double *dts,
                       double vare, int svh, const nav_t *nav,
                       const prcopt_t *opt, const double *rr, int rcv,
                       int f, int phase, double *y, double *e,
                       double *azel, double *sig)
{
    double r,pos[3],trp=0.0,vtrp=0.0,ion=0.0,vion=0.0,freq;
    double dant[NFREQ]={0},dants[NFREQ]={0},ru[3],rz[3],eu[3],ez[3],cosa,nadir;
    int p=postls_obsidx(opt,f),sys=satsys(obs->sat,NULL),ionoopt,tropopt;

    if (!(sys&opt->navsys)||opt->exsats[obs->sat-1]==1) return 0;
    if (!valid_sig(obs,p,phase)) return 0;
    if ((freq=sat2freq(obs->sat,obs->code[p],nav))<=0.0) return 0;
    if (satexclude(obs->sat,vare,svh,opt)) return 0;
    if ((r=geodist(rs,rr,e))<=0.0) return 0;

    ecef2pos(rr,pos);
    if (satazel(pos,e,azel)<opt->elmin) return 0;

    r+=-CLIGHT*dts[0];

    tropopt=opt->tropopt>=TROPOPT_EST?TROPOPT_SAAS:opt->tropopt;
    if (!tropcorr(obs->time,nav,pos,azel,tropopt,&trp,&vtrp)) return 0;
    r+=trp;

    ionoopt=opt->ionoopt;
    if (ionoopt==IONOOPT_EST) ionoopt=IONOOPT_BRDC;
    else if (ionoopt==IONOOPT_IFLC) ionoopt=IONOOPT_OFF;
    if (ionoopt!=IONOOPT_OFF) {
        if (!ionocorr(obs->time,nav,obs->sat,pos,azel,ionoopt,&ion,&vion)) {
            ion=0.0;
        }
        ion*=SQR(FREQ1/freq);
        r+=phase?-ion:ion;
    }
    if (opt->posopt[1]&&(rcv==1||rcv==2)) {
        antmodel(opt->pcvr+rcv-1,opt->antdel[rcv-1],azel,opt->posopt[1],dant);
        r+=dant[p<NFREQ?p:0];
    }
    if (opt->posopt[0]&&*nav->pcvs[obs->sat-1].type) {
        for (p=0;p<3;p++) {
            ru[p]=rr[p]-rs[p];
            rz[p]=-rs[p];
        }
        if (normv3(ru,eu)&&normv3(rz,ez)) {
            cosa=dot(eu,ez,3);
            cosa=cosa<-1.0?-1.0:(cosa>1.0?1.0:cosa);
            nadir=acos(cosa);
            antmodel_s(nav->pcvs+obs->sat-1,nadir,dants);
            p=postls_obsidx(opt,f);
            r+=dants[p<NFREQ?p:0];
        }
    }
    p=postls_obsidx(opt,f);
    *y=(phase?obs->L[p]*CLIGHT/freq:obs->P[p])-corr_code_bias(opt,nav,obs,p)-r;
    *sig=obs_sigma(obs,nav,opt,f,phase,azel[1]);
    return 1;
}
static int iflc_coeff(const obsd_t *obs, const nav_t *nav,
                      const prcopt_t *opt, double *c1, double *c2)
{
    double f1,f2;
    int p1=postls_obsidx(opt,0),p2=postls_obsidx(opt,1);

    if ((f1=sat2freq(obs->sat,obs->code[p1],nav))<=0.0) return 0;
    if ((f2=sat2freq(obs->sat,obs->code[p2],nav))<=0.0) return 0;
    if (fabs(f1-f2)<1.0) return 0;
    *c1=SQR(f1)/(SQR(f1)-SQR(f2));
    *c2=-SQR(f2)/(SQR(f1)-SQR(f2));
    return 1;
}
static int modeled_iflc(const obsd_t *obs, const double *rs, const double *dts,
                        double vare, int svh, const nav_t *nav,
                        const prcopt_t *opt, const double *rr, int rcv,
                        int phase, double *y, double *e, double *azel,
                        double *sig)
{
    double y1,y2,e1[3],e2[3],az1[2],az2[2],s1,s2,c1,c2;

    if (!iflc_coeff(obs,nav,opt,&c1,&c2)) return 0;
    if (!modeled_obs(obs,rs,dts,vare,svh,nav,opt,rr,rcv,0,phase,
                     &y1,e1,az1,&s1)) return 0;
    if (!modeled_obs(obs,rs,dts,vare,svh,nav,opt,rr,rcv,1,phase,
                     &y2,e2,az2,&s2)) return 0;
    *y=c1*y1+c2*y2;
    *sig=SQRT(SQR(c1*s1)+SQR(c2*s2));
    e[0]=e1[0]; e[1]=e1[1]; e[2]=e1[2];
    azel[0]=az1[0]; azel[1]=az1[1];
    return 1;
}
static int avepos_postls(const obs_t *obs, const nav_t *nav,
                         const prcopt_t *opt, double *rr)
{
    obsd_t data[MAXOBS];
    sol_t sol={{0}};
    char msg[128];
    int i,j,n,m,iobs=0,count=0;

    for (i=0;i<3;i++) rr[i]=0.0;
    while ((m=nextobsf_postls(obs,&iobs,1))>0) {
        for (i=j=0;i<m&&j<MAXOBS;i++) {
            data[j]=obs->data[iobs+i];
            if ((satsys(data[j].sat,NULL)&opt->navsys)&&
                opt->exsats[data[j].sat-1]!=1) j++;
        }
        if (j>0&&pntpos(data,j,nav,opt,&sol,NULL,NULL,msg)) {
            for (n=0;n<3;n++) rr[n]+=sol.rr[n];
            count++;
        }
        iobs+=m;
    }
    if (count<=0) return 0;
    for (i=0;i<3;i++) rr[i]/=count;
    return 1;
}
static int pair_epochs_tol(const obs_t *obs, int *iu, int *nu, int *ir, int *nr,
                           double maxdt)
{
    double dt,bestdt=1E99;
    int j,nb,prev=-1,nprev=0,best=-1,nbest=0;

    if ((*nu=nextobsf_postls(obs,iu,1))<=0) return 0;
    for (j=*ir;(nb=nextobsf_postls(obs,&j,2))>0;) {
        dt=timediff(obs->data[j].time,obs->data[*iu].time);
        if (dt>=0.0) break;
        prev=j;
        nprev=nb;
        j+=nb;
    }
    if (prev>=0) {
        bestdt=fabs(timediff(obs->data[prev].time,obs->data[*iu].time));
        best=prev; nbest=nprev;
    }
    if (nb>0) {
        dt=fabs(timediff(obs->data[j].time,obs->data[*iu].time));
        if (dt<bestdt) {
            bestdt=dt;
            best=j; nbest=nb;
        }
    }
    if (best<0) return 0;
    if (bestdt>maxdt) {
        *iu+=*nu;
        return -1;
    }
    *ir=best;
    *nr=nbest;
    return 1;
}
static int pair_epochs(const obs_t *obs, int *iu, int *nu, int *ir, int *nr)
{
    return pair_epochs_tol(obs,iu,nu,ir,nr,DTTOL);
}
static void select_refs(const obs_t *obs, const nav_t *nav,
                        const prcopt_t *opt, const double *ru,
                        int refsat[POSTLS_NSYS][NFREQ], double maxdt)
{
    obsd_t data[MAXOBS*2];
    double *rs,*dts,*vare,e[3],azel[2],pos[3];
    double score[POSTLS_NSYS][NFREQ][MAXSAT];
    int *svh;
    int iu=0,ir=0,nu,nr,ret,i,j,f,sat,bi,sys,si,n,nf=postls_nf(opt);

    for (si=0;si<POSTLS_NSYS;si++) for (f=0;f<NFREQ;f++) {
        refsat[si][f]=0;
        for (sat=0;sat<MAXSAT;sat++) score[si][f][sat]=0.0;
    }
    rs=mat(6,MAXOBS*2); dts=mat(2,MAXOBS*2);
    vare=mat(1,MAXOBS*2); svh=imat(1,MAXOBS*2);

    while ((ret=pair_epochs_tol(obs,&iu,&nu,&ir,&nr,maxdt))!=0) {
        if (ret<0) continue;
        for (i=0;i<nu;i++) data[i]=obs->data[iu+i];
        for (i=0;i<nr;i++) data[nu+i]=obs->data[ir+i];
        n=nu+nr;
        satposs(data[0].time,data,n,nav,opt->sateph,rs,dts,vare,svh);
        for (i=0;i<nu;i++) {
            sat=data[i].sat;
            if ((bi=find_sat_obs(data+nu,nr,sat))<0) continue;
            if ((sys=satsys(sat,NULL))==0||(si=postls_sysidx(sys))<0) continue;
            if (!(sys&opt->navsys)) continue;
            if (geodist(rs+i*6,ru,e)<=0.0) continue;
            ecef2pos(ru,pos);
            if (satazel(pos,e,azel)<opt->elmin) continue;
            if (postls_iflc(opt)) {
                if (valid_sig(data+i,postls_obsidx(opt,0),0)&&
                    valid_sig(data+i,postls_obsidx(opt,1),0)&&
                    valid_sig(data+nu+bi,postls_obsidx(opt,0),0)&&
                    valid_sig(data+nu+bi,postls_obsidx(opt,1),0)&&
                    iflc_coeff(data+i,nav,opt,&e[0],&e[1])) {
                    score[si][0][sat-1]+=1.0+azel[1]*R2D/90.0;
                }
                continue;
            }
            for (f=0;f<nf;f++) {
                j=postls_obsidx(opt,f);
                if (!valid_sig(data+i,j,0)||!valid_sig(data+nu+bi,j,0)) continue;
                if (sat2freq(sat,data[i].code[j],nav)<=0.0) continue;
                score[si][f][sat-1]+=1.0+azel[1]*R2D/90.0;
            }
        }
        iu+=nu;
    }
    for (si=0;si<POSTLS_NSYS;si++) for (f=0;f<nf;f++) {
        double best=0.0;
        for (sat=1;sat<=MAXSAT;sat++) {
            if (score[si][f][sat-1]>best) {
                best=score[si][f][sat-1];
                refsat[si][f]=sat;
            }
        }
    }
    free(rs); free(dts); free(vare); free(svh);
}
static int find_amb(postls_amb_t *ambs, int namb, int sat, int f,
                    int aur, int abr, int arr, int arb,
                    int aur2, int abr2, int arr2, int arb2)
{
    int i;

    for (i=0;i<namb;i++) {
        if (ambs[i].sat==sat&&ambs[i].f==f&&ambs[i].aur==aur&&
            ambs[i].abr==abr&&ambs[i].arr==arr&&ambs[i].arb==arb&&
            ambs[i].aur2==aur2&&ambs[i].abr2==abr2&&
            ambs[i].arr2==arr2&&ambs[i].arb2==arb2) return i;
    }
    return -1;
}
static int add_amb(postls_amb_t *ambs, int *namb, int sat, int f,
                   int aur, int abr, int arr, int arb,
                   int aur2, int abr2, int arr2, int arb2)
{
    int i;

    if ((i=find_amb(ambs,*namb,sat,f,aur,abr,arr,arb,
                    aur2,abr2,arr2,arb2))>=0) return i;
    if (*namb>=POSTLS_MAXAMB) return -1;
    i=*namb;
    ambs[i].sat=sat;
    ambs[i].f=f;
    ambs[i].aur=aur;
    ambs[i].abr=abr;
    ambs[i].arr=arr;
    ambs[i].arb=arb;
    ambs[i].aur2=aur2;
    ambs[i].abr2=abr2;
    ambs[i].arr2=arr2;
    ambs[i].arb2=arb2;
    ambs[i].index=3+i;
    ambs[i].n=0;
    ambs[i].use=1;
    ambs[i].bad=0;
    (*namb)++;
    return i;
}
static int prepare_ambs(postls_amb_t *ambs, int namb)
{
    int i,n=0;

    for (i=0;i<namb;i++) {
        ambs[i].use=!ambs[i].bad&&ambs[i].n>=POSTLS_MINPHASEOBS;
        ambs[i].index=ambs[i].use?3+n++:-1;
    }
    return n;
}
static int scan_rows(const obs_t *obs, const nav_t *nav, const prcopt_t *opt,
                     const double *ru, const int refsat[POSTLS_NSYS][NFREQ],
                     const int *arc, postls_amb_t *ambs, int *namb)
{
    obsd_t data[MAXOBS*2];
    double *rs,*dts,*vare,yu,yb,yr,ybr,ebr[3],e[3],az[2],sig;
    int *svh;
    int iu=0,ir=0,nu,nr,ret,i,f,sat,bi,ri,rbi,sys,si,n,nv=0,nf=postls_nf(opt);
    int aur,abr,arr,arb,ia;

    rs=mat(6,MAXOBS*2); dts=mat(2,MAXOBS*2);
    vare=mat(1,MAXOBS*2); svh=imat(1,MAXOBS*2);

    while ((ret=pair_epochs(obs,&iu,&nu,&ir,&nr))!=0) {
        if (ret<0) continue;
        for (i=0;i<nu;i++) data[i]=obs->data[iu+i];
        for (i=0;i<nr;i++) data[nu+i]=obs->data[ir+i];
        n=nu+nr;
        satposs(data[0].time,data,n,nav,opt->sateph,rs,dts,vare,svh);
        for (i=0;i<nu;i++) {
            sat=data[i].sat;
            if ((bi=find_sat_obs(data+nu,nr,sat))<0) continue;
            if ((sys=satsys(sat,NULL))==0||(si=postls_sysidx(sys))<0) continue;
            if (postls_iflc(opt)) {
                f=0;
                if (refsat[si][f]<=0||sat==refsat[si][f]) continue;
                if ((ri=find_sat_obs(data,nu,refsat[si][f]))<0) continue;
                if ((rbi=find_sat_obs(data+nu,nr,refsat[si][f]))<0) continue;
                rbi+=nu;
                if (!modeled_iflc(data+i,rs+i*6,dts+i*2,vare[i],svh[i],nav,opt,
                                  ru,1,0,&yu,e,az,&sig)) continue;
                if (!modeled_iflc(data+nu+bi,rs+(nu+bi)*6,dts+(nu+bi)*2,
                                  vare[nu+bi],svh[nu+bi],nav,opt,opt->rb,2,0,
                                  &yb,ebr,az,&sig)) continue;
                if (!modeled_iflc(data+ri,rs+ri*6,dts+ri*2,vare[ri],svh[ri],
                                  nav,opt,ru,1,0,&yr,e,az,&sig)) continue;
                if (!modeled_iflc(data+rbi,rs+rbi*6,dts+rbi*2,vare[rbi],svh[rbi],
                                  nav,opt,opt->rb,2,0,&ybr,ebr,az,&sig)) continue;
                nv++;
                if (modeled_iflc(data+i,rs+i*6,dts+i*2,vare[i],svh[i],nav,opt,
                                 ru,1,1,&yu,e,az,&sig)&&
                    modeled_iflc(data+nu+bi,rs+(nu+bi)*6,dts+(nu+bi)*2,
                                 vare[nu+bi],svh[nu+bi],nav,opt,opt->rb,2,1,
                                 &yb,ebr,az,&sig)&&
                    modeled_iflc(data+ri,rs+ri*6,dts+ri*2,vare[ri],svh[ri],
                                 nav,opt,ru,1,1,&yr,e,az,&sig)&&
                    modeled_iflc(data+rbi,rs+rbi*6,dts+rbi*2,vare[rbi],svh[rbi],
                                 nav,opt,opt->rb,2,1,&ybr,ebr,az,&sig)) {
                    aur=arc[(iu+i)*nf];
                    abr=arc[(ir+bi)*nf];
                    arr=arc[(iu+ri)*nf];
                    arb=arc[(ir+rbi-nu)*nf];
                    if (aur>0&&abr>0&&arr>0&&arb>0&&
                        arc[(iu+i)*nf+1]>0&&arc[(ir+bi)*nf+1]>0&&
                        arc[(iu+ri)*nf+1]>0&&arc[(ir+rbi-nu)*nf+1]>0) {
                        if (!ambs) {
                            nv++;
                        }
                        else if ((ia=add_amb(ambs,namb,sat,-1,aur,abr,arr,arb,
                                             arc[(iu+i)*nf+1],arc[(ir+bi)*nf+1],
                                             arc[(iu+ri)*nf+1],
                                             arc[(ir+rbi-nu)*nf+1]))>=0) {
                            ambs[ia].n++;
                            nv++;
                        }
                    }
                }
                continue;
            }
            for (f=0;f<nf;f++) {
                if (refsat[si][f]<=0||sat==refsat[si][f]) continue;
                if ((ri=find_sat_obs(data,nu,refsat[si][f]))<0) continue;
                if ((rbi=find_sat_obs(data+nu,nr,refsat[si][f]))<0) continue;
                rbi+=nu;
                if (!modeled_obs(data+i,rs+i*6,dts+i*2,vare[i],svh[i],nav,opt,
                                 ru,1,f,0,&yu,e,az,&sig)) continue;
                if (!modeled_obs(data+nu+bi,rs+(nu+bi)*6,dts+(nu+bi)*2,
                                 vare[nu+bi],svh[nu+bi],nav,opt,opt->rb,2,f,0,
                                 &yb,ebr,az,&sig)) continue;
                if (!modeled_obs(data+ri,rs+ri*6,dts+ri*2,vare[ri],svh[ri],
                                 nav,opt,ru,1,f,0,&yr,e,az,&sig)) continue;
                if (!modeled_obs(data+rbi,rs+rbi*6,dts+rbi*2,vare[rbi],svh[rbi],
                                 nav,opt,opt->rb,2,f,0,&ybr,ebr,az,&sig)) continue;
                nv++;
                if (modeled_obs(data+i,rs+i*6,dts+i*2,vare[i],svh[i],nav,opt,
                                ru,1,f,1,&yu,e,az,&sig)&&
                    modeled_obs(data+nu+bi,rs+(nu+bi)*6,dts+(nu+bi)*2,
                                vare[nu+bi],svh[nu+bi],nav,opt,opt->rb,2,f,1,
                                &yb,ebr,az,&sig)&&
                    modeled_obs(data+ri,rs+ri*6,dts+ri*2,vare[ri],svh[ri],
                                nav,opt,ru,1,f,1,&yr,e,az,&sig)&&
                    modeled_obs(data+rbi,rs+rbi*6,dts+rbi*2,vare[rbi],svh[rbi],
                                nav,opt,opt->rb,2,f,1,&ybr,ebr,az,&sig)) {
                    aur=arc[(iu+i)*nf+f];
                    abr=arc[(ir+bi)*nf+f];
                    arr=arc[(iu+ri)*nf+f];
                    arb=arc[(ir+rbi-nu)*nf+f];
                    if (aur>0&&abr>0&&arr>0&&arb>0) {
                        if (!ambs) {
                            nv++;
                        }
                        else if ((ia=add_amb(ambs,namb,sat,f,aur,abr,arr,arb,
                                             0,0,0,0))>=0) {
                            ambs[ia].n++;
                            nv++;
                        }
                    }
                }
            }
        }
        iu+=nu;
    }
    free(rs); free(dts); free(vare); free(svh);
    return nv;
}
static int fill_system(const obs_t *obs, const nav_t *nav, const prcopt_t *opt,
                       const double *ru,
                       const int refsat[POSTLS_NSYS][NFREQ],
                       const int *arc, const postls_amb_t *ambs, int namb,
                       int nx, double *A, double *v, int *used_sat,
                       int *row_amb)
{
    obsd_t data[MAXOBS*2];
    double *rs,*dts,*vare,yus,ybs,yur,ybr,sus,sbs,sur,sbr;
    double eus[3],ebs[3],eur[3],ebr[3],az[2],sig,dd;
    int *svh;
    int iu=0,ir=0,nu,nr,ret,i,k,f,sat,bi,ri,rbi,sys,si,n,nv=0,nf=postls_nf(opt);
    int aur,abr,arr,arb,ia;

    rs=mat(6,MAXOBS*2); dts=mat(2,MAXOBS*2);
    vare=mat(1,MAXOBS*2); svh=imat(1,MAXOBS*2);

    while ((ret=pair_epochs(obs,&iu,&nu,&ir,&nr))!=0) {
        if (ret<0) continue;
        for (i=0;i<nu;i++) data[i]=obs->data[iu+i];
        for (i=0;i<nr;i++) data[nu+i]=obs->data[ir+i];
        n=nu+nr;
        satposs(data[0].time,data,n,nav,opt->sateph,rs,dts,vare,svh);
        for (i=0;i<nu;i++) {
            sat=data[i].sat;
            if ((bi=find_sat_obs(data+nu,nr,sat))<0) continue;
            if ((sys=satsys(sat,NULL))==0||(si=postls_sysidx(sys))<0) continue;
            if (postls_iflc(opt)) {
                f=0;
                if (refsat[si][f]<=0||sat==refsat[si][f]) continue;
                if ((ri=find_sat_obs(data,nu,refsat[si][f]))<0) continue;
                if ((rbi=find_sat_obs(data+nu,nr,refsat[si][f]))<0) continue;
                rbi+=nu;
                if (!modeled_iflc(data+i,rs+i*6,dts+i*2,vare[i],svh[i],nav,opt,
                                  ru,1,0,&yus,eus,az,&sus)) continue;
                if (!modeled_iflc(data+nu+bi,rs+(nu+bi)*6,dts+(nu+bi)*2,
                                  vare[nu+bi],svh[nu+bi],nav,opt,opt->rb,2,0,
                                  &ybs,ebs,az,&sbs)) continue;
                if (!modeled_iflc(data+ri,rs+ri*6,dts+ri*2,vare[ri],svh[ri],
                                  nav,opt,ru,1,0,&yur,eur,az,&sur)) continue;
                if (!modeled_iflc(data+rbi,rs+rbi*6,dts+rbi*2,vare[rbi],svh[rbi],
                                  nav,opt,opt->rb,2,0,&ybr,ebr,az,&sbr)) continue;
                sig=SQRT(SQR(sus)+SQR(sbs)+SQR(sur)+SQR(sbr));
                if (sig<=0.0) sig=1.0;
                dd=(yus-ybs)-(yur-ybr);
                v[nv]=dd/sig;
                for (k=0;k<nx;k++) A[k+nv*nx]=0.0;
                for (k=0;k<3;k++) A[k+nv*nx]=(-eus[k]+eur[k])/sig;
                if (row_amb) row_amb[nv]=-1;
                used_sat[sat-1]=used_sat[refsat[si][f]-1]=1;
                nv++;
                if (modeled_iflc(data+i,rs+i*6,dts+i*2,vare[i],svh[i],nav,opt,
                                 ru,1,1,&yus,eus,az,&sus)&&
                    modeled_iflc(data+nu+bi,rs+(nu+bi)*6,dts+(nu+bi)*2,
                                 vare[nu+bi],svh[nu+bi],nav,opt,opt->rb,2,1,
                                 &ybs,ebs,az,&sbs)&&
                    modeled_iflc(data+ri,rs+ri*6,dts+ri*2,vare[ri],svh[ri],
                                 nav,opt,ru,1,1,&yur,eur,az,&sur)&&
                    modeled_iflc(data+rbi,rs+rbi*6,dts+rbi*2,vare[rbi],svh[rbi],
                                 nav,opt,opt->rb,2,1,&ybr,ebr,az,&sbr)) {
                    aur=arc[(iu+i)*nf];
                    abr=arc[(ir+bi)*nf];
                    arr=arc[(iu+ri)*nf];
                    arb=arc[(ir+rbi-nu)*nf];
                    if (aur>0&&abr>0&&arr>0&&arb>0&&
                        arc[(iu+i)*nf+1]>0&&arc[(ir+bi)*nf+1]>0&&
                        arc[(iu+ri)*nf+1]>0&&arc[(ir+rbi-nu)*nf+1]>0) {
                        if ((ia=find_amb((postls_amb_t *)ambs,namb,sat,-1,
                                         aur,abr,arr,arb,arc[(iu+i)*nf+1],
                                         arc[(ir+bi)*nf+1],arc[(iu+ri)*nf+1],
                                         arc[(ir+rbi-nu)*nf+1]))<0) continue;
                        if (!ambs[ia].use) continue;
                        sig=SQRT(SQR(sus)+SQR(sbs)+SQR(sur)+SQR(sbr));
                        if (sig<=0.0) sig=0.02;
                        dd=(yus-ybs)-(yur-ybr);
                        v[nv]=dd/sig;
                        for (k=0;k<nx;k++) A[k+nv*nx]=0.0;
                        for (k=0;k<3;k++) A[k+nv*nx]=(-eus[k]+eur[k])/sig;
                        A[ambs[ia].index+nv*nx]=1.0/sig;
                        if (row_amb) row_amb[nv]=ia;
                        used_sat[sat-1]=used_sat[refsat[si][f]-1]=1;
                        nv++;
                    }
                }
                continue;
            }
            for (f=0;f<nf;f++) {
                if (refsat[si][f]<=0||sat==refsat[si][f]) continue;
                if ((ri=find_sat_obs(data,nu,refsat[si][f]))<0) continue;
                if ((rbi=find_sat_obs(data+nu,nr,refsat[si][f]))<0) continue;
                rbi+=nu;

                if (!modeled_obs(data+i,rs+i*6,dts+i*2,vare[i],svh[i],nav,opt,
                                 ru,1,f,0,&yus,eus,az,&sus)) continue;
                if (!modeled_obs(data+nu+bi,rs+(nu+bi)*6,dts+(nu+bi)*2,
                                 vare[nu+bi],svh[nu+bi],nav,opt,opt->rb,2,f,0,
                                 &ybs,ebs,az,&sbs)) continue;
                if (!modeled_obs(data+ri,rs+ri*6,dts+ri*2,vare[ri],svh[ri],
                                 nav,opt,ru,1,f,0,&yur,eur,az,&sur)) continue;
                if (!modeled_obs(data+rbi,rs+rbi*6,dts+rbi*2,vare[rbi],svh[rbi],
                                 nav,opt,opt->rb,2,f,0,&ybr,ebr,az,&sbr)) continue;

                sig=SQRT(SQR(sus)+SQR(sbs)+SQR(sur)+SQR(sbr));
                if (sig<=0.0) sig=1.0;
                dd=(yus-ybs)-(yur-ybr);
                v[nv]=dd/sig;
                for (k=0;k<nx;k++) A[k+nv*nx]=0.0;
                for (k=0;k<3;k++) A[k+nv*nx]=(-eus[k]+eur[k])/sig;
                if (row_amb) row_amb[nv]=-1;
                used_sat[sat-1]=used_sat[refsat[si][f]-1]=1;
                nv++;

                if (!modeled_obs(data+i,rs+i*6,dts+i*2,vare[i],svh[i],nav,opt,
                                 ru,1,f,1,&yus,eus,az,&sus)) continue;
                if (!modeled_obs(data+nu+bi,rs+(nu+bi)*6,dts+(nu+bi)*2,
                                 vare[nu+bi],svh[nu+bi],nav,opt,opt->rb,2,f,1,
                                 &ybs,ebs,az,&sbs)) continue;
                if (!modeled_obs(data+ri,rs+ri*6,dts+ri*2,vare[ri],svh[ri],
                                 nav,opt,ru,1,f,1,&yur,eur,az,&sur)) continue;
                if (!modeled_obs(data+rbi,rs+rbi*6,dts+rbi*2,vare[rbi],svh[rbi],
                                 nav,opt,opt->rb,2,f,1,&ybr,ebr,az,&sbr)) continue;
                aur=arc[(iu+i)*nf+f];
                abr=arc[(ir+bi)*nf+f];
                arr=arc[(iu+ri)*nf+f];
                arb=arc[(ir+rbi-nu)*nf+f];
                if (aur<=0||abr<=0||arr<=0||arb<=0) continue;
                if ((ia=find_amb((postls_amb_t *)ambs,namb,sat,f,aur,abr,arr,arb,
                                 0,0,0,0))<0) continue;
                if (!ambs[ia].use) continue;

                sig=SQRT(SQR(sus)+SQR(sbs)+SQR(sur)+SQR(sbr));
                if (sig<=0.0) sig=0.02;
                dd=(yus-ybs)-(yur-ybr);
                v[nv]=dd/sig;
                for (k=0;k<nx;k++) A[k+nv*nx]=0.0;
                for (k=0;k<3;k++) A[k+nv*nx]=(-eus[k]+eur[k])/sig;
                A[ambs[ia].index+nv*nx]=1.0/sig;
                if (row_amb) row_amb[nv]=ia;
                used_sat[sat-1]=used_sat[refsat[si][f]-1]=1;
                nv++;
            }
        }
        iu+=nu;
    }
    free(rs); free(dts); free(vare); free(svh);
    return nv;
}
static int reject_bad_ambs(const double *A, const double *v, const double *x,
                           int nx, int nv, const int *row_amb,
                           postls_amb_t *ambs, int namb)
{
    double *sum,*maxr,r,rms,worst=0.0;
    int *cnt,i,j,ia,worstia=-1;

    sum=mat(namb,1); maxr=mat(namb,1); cnt=imat(namb,1);
    for (j=0;j<nv;j++) {
        if (!row_amb||(ia=row_amb[j])<0||ia>=namb||!ambs[ia].use) continue;
        r=v[j];
        for (i=0;i<nx;i++) r-=A[i+j*nx]*x[i];
        sum[ia]+=r*r;
        cnt[ia]++;
        if (fabs(r)>maxr[ia]) maxr[ia]=fabs(r);
    }
    for (ia=0;ia<namb;ia++) {
        if (!ambs[ia].use||cnt[ia]<POSTLS_MINPHASEOBS) continue;
        rms=sqrt(sum[ia]/cnt[ia]);
        if ((rms>POSTLS_ARC_REJSTD||maxr[ia]>POSTLS_REJSTD)&&rms>worst) {
            worst=rms;
            worstia=ia;
        }
    }
    if (worstia>=0) {
        ambs[worstia].bad=1;
        ambs[worstia].use=0;
    }
    free(sum); free(maxr); free(cnt);
    return worstia>=0;
}
static int compact_rows(const double *A, const double *v, const int *use,
                        int nx, int nv, double *Ac, double *vc)
{
    int i,j,nc=0;

    for (j=0;j<nv;j++) {
        if (!use[j]) continue;
        vc[nc]=v[j];
        for (i=0;i<nx;i++) Ac[i+nc*nx]=A[i+j*nx];
        nc++;
    }
    return nc;
}
static int robust_lsq(const double *A, const double *v, int nx, int nv,
                      double *x, double *Q, int *nout)
{
    double *Ac,*vc,*xp,*Qp,r,worst;
    int *use,i,j,k,nc,info,out,worstj,have=0;

    if (nout) *nout=0;
    Ac=mat(nx,nv); vc=mat(nv,1); xp=mat(nx,1); Qp=mat(nx,nx);
    use=imat(nv,1);
    for (j=0;j<nv;j++) use[j]=1;

    if ((info=lsq(A,v,nx,nv,x,Q))) {
        free(Ac); free(vc); free(xp); free(Qp); free(use);
        return info;
    }
    for (k=0;k<POSTLS_MAXROB;k++) {
        nc=compact_rows(A,v,use,nx,nv,Ac,vc);
        if (nc<=nx+POSTLS_MINOBS) break;
        if ((info=lsq(Ac,vc,nx,nc,xp,Qp))) {
            break;
        }
        have=1;
        worst=0.0; worstj=-1;
        for (j=0;j<nv;j++) {
            if (!use[j]) continue;
            r=v[j];
            for (i=0;i<nx;i++) r-=A[i+j*nx]*xp[i];
            if (fabs(r)>worst) {
                worst=fabs(r);
                worstj=j;
            }
        }
        if (worstj<0||worst<=POSTLS_REJSTD) {
            for (i=0;i<nx;i++) x[i]=xp[i];
            for (i=0;i<nx*nx;i++) Q[i]=Qp[i];
            break;
        }
        use[worstj]=0;
        for (i=0;i<nx;i++) x[i]=xp[i];
        for (i=0;i<nx*nx;i++) Q[i]=Qp[i];
    }
    out=0;
    for (j=0;j<nv;j++) {
        if (!use[j]) out++;
    }
    if (have) {
        for (i=0;i<nx;i++) x[i]=xp[i];
        for (i=0;i<nx*nx;i++) Q[i]=Qp[i];
    }
    if (nout) *nout=out;

    free(Ac); free(vc); free(xp); free(Qp); free(use);
    return 0;
}
static double amb_lambda_postls(const obs_t *obs, const nav_t *nav,
                                const prcopt_t *opt, int sat, int f)
{
    double freq;
    int i,p=postls_obsidx(opt,f),sys=satsys(sat,NULL);

    if (f<0) return 0.0;
    if (sys==SYS_GLO||sys==SYS_SBS) return 0.0;
    for (i=0;i<obs->n;i++) {
        if (obs->data[i].rcv!=1||obs->data[i].sat!=sat) continue;
        if (!valid_sig(obs->data+i,p,1)) continue;
        if ((freq=sat2freq(sat,obs->data[i].code[p],nav))<=0.0) continue;
        return CLIGHT/freq;
    }
    return 0.0;
}
static int fix_ambiguities(const obs_t *obs, const nav_t *nav,
                           const prcopt_t *opt, const postls_amb_t *ambs,
                           int namb, int nx,
                           const double *x, const double *Q,
                           double *ru, float *ratio)
{
    double *af,*Qa,*F,*s,*lam,*invQa,*d,*z,*score,corr[3]={0};
    int *idx,*pool,*used,i,j,k,nb=0,np=0,ntry,info,best;

    if (ratio) *ratio=0.0f;
    if (opt->modear==ARMODE_OFF||namb<4) return 0;
    af=mat(namb,1); Qa=mat(namb,namb); F=mat(namb,2); s=mat(2,1);
    lam=mat(namb,1); invQa=mat(namb,namb); d=mat(namb,1);
    z=mat(namb,1); score=mat(namb,1);
    idx=imat(namb,1); pool=imat(namb,1); used=imat(namb,1);

    for (i=0;i<namb;i++) {
        double qaa;

        if (!ambs[i].use) continue;
        if (ambs[i].n<POSTLS_MINARCOBS) continue;
        lam[np]=amb_lambda_postls(obs,nav,opt,ambs[i].sat,ambs[i].f);
        if (lam[np]<=0.0) continue;
        qaa=Q[ambs[i].index+ambs[i].index*nx]/SQR(lam[np]);
        if (qaa<=0.0||qaa>1.0) continue;
        pool[np]=i;
        score[np]=ambs[i].n/(sqrt(qaa)+1E-6);
        np++;
    }
    if (np<4) {
        free(af); free(Qa); free(F); free(s); free(lam); free(invQa); free(d); free(z);
        free(score); free(idx); free(pool); free(used);
        return 0;
    }
    for (i=0;i<np;i++) used[i]=0;
    for (i=0;i<MIN_(POSTLS_MAXFIXAMB,np);i++) {
        best=-1;
        for (j=0;j<np;j++) {
            if (!used[j]&&(best<0||score[j]>score[best])) best=j;
        }
        if (best<0) break;
        used[best]=1;
        idx[i]=pool[best];
        lam[i]=amb_lambda_postls(obs,nav,opt,ambs[idx[i]].sat,ambs[idx[i]].f);
        af[i]=x[ambs[idx[i]].index]/lam[i];
    }
    nb=MIN_(POSTLS_MAXFIXAMB,np);
    if (nb<4) {
        free(af); free(Qa); free(F); free(s); free(lam); free(invQa); free(d); free(z);
        free(score); free(idx); free(pool); free(used);
        return 0;
    }
    for (ntry=nb;ntry>=4;ntry--) {
        for (i=0;i<ntry;i++) for (j=0;j<ntry;j++) {
            Qa[i+j*ntry]=Q[ambs[idx[i]].index+ambs[idx[j]].index*nx]/
                         (lam[i]*lam[j]);
        }
        if ((info=lambda(ntry,2,af,Qa,F,s))||s[0]<=0.0) continue;
        if (ratio) *ratio=(float)(s[1]/s[0]);
        if (s[1]/s[0]>=MAX_(SMARTPHONE(opt)?10.0:3.0,
                             opt->thresar[0]>0.0?opt->thresar[0]:3.0)) break;
    }
    if (ntry<4) {
        free(af); free(Qa); free(F); free(s); free(lam); free(invQa); free(d); free(z);
        free(score); free(idx); free(pool); free(used);
        return 0;
    }
    for (i=0;i<ntry;i++) d[i]=af[i]-F[i];
    matcpy(invQa,Qa,ntry,ntry);
    if (matinv(invQa,ntry)) {
        free(af); free(Qa); free(F); free(s); free(lam); free(invQa); free(d); free(z);
        free(score); free(idx); free(pool); free(used);
        return 0;
    }
    matmul("NN",ntry,1,ntry,1.0,invQa,d,0.0,z);
    for (k=0;k<3;k++) {
        for (i=0;i<ntry;i++) {
            corr[k]+=Q[k+ambs[idx[i]].index*nx]/lam[i]*z[i];
        }
        ru[k]-=corr[k];
    }
    free(af); free(Qa); free(F); free(s); free(lam); free(invQa); free(d); free(z);
    free(score); free(idx); free(pool); free(used);
    return ntry;
}
extern int postls_relpos(const obs_t *obs, const nav_t *nav,
                         const prcopt_t *opt, sol_t *sol, char *msg)
{
    prcopt_t opt_=*opt;
    double ru[3],*A,*v,*dx,*Q,dpos;
    float ratio=0.0f;
    int refsat[POSTLS_NSYS][NFREQ];
    postls_amb_t *ambs;
    int *arc,*row_amb;
    int used_sat[MAXSAT],i,namb=0,nuseamb,nx,nv,maxnv,info,iter,nused=0,nf,nout=0,nout_tot=0,nfix=0;
    int nbadamb=0;

    trace(3,"postls_relpos: nobs=%d\n",obs?obs->n:0);
    if (msg) *msg='\0';
    if (!obs||!nav||!opt||!sol) return 0;
    if (obs->n<=0) {
        if (msg) strcpy(msg,"no observation data");
        return 0;
    }
    if (norm(opt_.rb,3)<=0.0) {
        if (msg) strcpy(msg,"base position is required");
        return 0;
    }
    if (!opt_.navsys) opt_.navsys=SYS_GPS|SYS_GLO|SYS_GAL|SYS_QZS|SYS_CMP|SYS_IRN;
    nf=postls_nf(&opt_);

    if (norm(opt_.ru,3)>0.0) for (i=0;i<3;i++) ru[i]=opt_.ru[i];
    else if (!avepos_postls(obs,nav,&opt_,ru)) {
        if (msg) strcpy(msg,"could not estimate initial rover position by SPP");
        return 0;
    }
    arc=imat(obs->n,nf);
    ambs=(postls_amb_t *)calloc(POSTLS_MAXAMB,sizeof(postls_amb_t));
    if (!arc||!ambs) {
        if (msg) strcpy(msg,"memory allocation error");
        free(arc); free(ambs);
        return 0;
    }
    build_phase_arcs(obs,&opt_,nf,arc);
    select_refs(obs,nav,&opt_,ru,refsat,DTTOL);
    maxnv=scan_rows(obs,nav,&opt_,ru,refsat,arc,ambs,&namb);
    nuseamb=prepare_ambs(ambs,namb);
    nx=3+nuseamb;
    if (maxnv<nx||maxnv<POSTLS_MINOBS) {
        if (msg) sprintf(msg,"not enough DD observations (nv=%d nx=%d)",maxnv,nx);
        free(arc); free(ambs);
        return 0;
    }
    A=mat(nx,maxnv); v=mat(maxnv,1); dx=mat(nx,1); Q=mat(nx,nx);
    row_amb=imat(maxnv,1);
    if (!A||!v||!dx||!Q||!row_amb) {
        if (msg) strcpy(msg,"memory allocation error");
        free(A); free(v); free(dx); free(Q); free(row_amb); free(arc); free(ambs);
        return 0;
    }

    for (iter=0;iter<POSTLS_MAXITR;iter++) {
        prcopt_t iteropt=opt_;

        if (SMARTPHONE(&iteropt)&&!postls_iflc(&iteropt)&&iter>0) {
            for (i=0;i<NFREQ;i++) if (iteropt.eratio[i]<1000.0) {
                iteropt.eratio[i]=1000.0;
            }
        }
        nuseamb=prepare_ambs(ambs,namb);
        nx=3+nuseamb;
        for (i=0;i<MAXSAT;i++) used_sat[i]=0;
        nv=fill_system(obs,nav,&iteropt,ru,refsat,arc,ambs,namb,nx,A,v,used_sat,
                       row_amb);
        if (nv<nx) {
            if (msg) sprintf(msg,"rank-deficient DD system (nv=%d nx=%d)",nv,nx);
            free(A); free(v); free(dx); free(Q); free(row_amb); free(arc); free(ambs);
            return 0;
        }
        if ((info=robust_lsq(A,v,nx,nv,dx,Q,&nout))) {
            if (msg) sprintf(msg,"least-squares error info=%d nv=%d nx=%d",info,nv,nx);
            free(A); free(v); free(dx); free(Q); free(row_amb); free(arc); free(ambs);
            return 0;
        }
        nout_tot+=nout;
        if (nbadamb<POSTLS_MAXBADARC&&
            reject_bad_ambs(A,v,dx,nx,nv,row_amb,ambs,namb)) {
            nbadamb++;
            continue;
        }
        for (i=0;i<3;i++) ru[i]+=dx[i];
        dpos=norm(dx,3);
        if (dpos<1E-4) break;
    }
    for (i=0;i<MAXSAT;i++) if (used_sat[i]) nused++;
    nfix=fix_ambiguities(obs,nav,&opt_,ambs,namb,nx,dx,Q,ru,&ratio);
    memset(sol,0,sizeof(*sol));
    for (i=0;i<obs->n;i++) if (obs->data[i].rcv==1) {
        sol->time=obs->data[i].time;
        break;
    }
    for (i=0;i<3;i++) sol->rr[i]=ru[i];
    sol->qr[0]=(float)Q[0];
    sol->qr[1]=(float)Q[1+nx];
    sol->qr[2]=(float)Q[2+2*nx];
    sol->qr[3]=(float)Q[1];
    sol->qr[4]=(float)Q[1+2*nx];
    sol->qr[5]=(float)Q[2];
    sol->type=0;
    sol->stat=nfix>0?SOLQ_FIX:SOLQ_FLOAT;
    sol->ns=(uint8_t)MIN_(nused,255);
    sol->age=0.0f;
    sol->ratio=nfix>0?ratio:0.0f;
    if (msg) sprintf(msg,"postls %s static: iter=%d nv=%d nx=%d namb=%d/%d badarc=%d nfix=%d out=%d ratio=%.2f",
                     nfix>0?"fix":"float",iter+1,nv,nx,nuseamb,namb,nbadamb,nfix,nout_tot,ratio);

    free(A); free(v); free(dx); free(Q); free(row_amb); free(arc); free(ambs);
    return 1;
}
static int build_epochs(const obs_t *obs, const nav_t *nav, const prcopt_t *opt,
                        postls_epoch_t **epochs, int *nep, double *ru0,
                        double maxdt)
{
    postls_epoch_t *eps=NULL,*tmp;
    obsd_t data[MAXOBS];
    sol_t sol={{0}};
    char msg[128];
    int iu=0,ir=0,nu,nr,ret,i,j,n=0,nmax=0,have=0;

    *epochs=NULL; *nep=0;
    for (i=0;i<3;i++) ru0[i]=0.0;
    while ((ret=pair_epochs_tol(obs,&iu,&nu,&ir,&nr,maxdt))!=0) {
        if (ret<0) continue;
        if (n>=nmax) {
            nmax=nmax==0?256:nmax*2;
            if (!(tmp=(postls_epoch_t *)realloc(eps,sizeof(postls_epoch_t)*nmax))) {
                free(eps);
                return 0;
            }
            eps=tmp;
        }
        eps[n].iu=iu;
        eps[n].nu=nu;
        eps[n].ir=ir;
        eps[n].nr=nr;
        eps[n].index=3*n;
        eps[n].nv=eps[n].ns=0;
        eps[n].time=obs->data[iu].time;
        for (i=j=0;i<nu&&j<MAXOBS;i++) {
            data[j]=obs->data[iu+i];
            if ((satsys(data[j].sat,NULL)&opt->navsys)&&
                opt->exsats[data[j].sat-1]!=1) j++;
        }
        if (j>0&&pntpos(data,j,nav,opt,&sol,NULL,NULL,msg)) {
            for (i=0;i<3;i++) eps[n].ru[i]=sol.rr[i];
            for (i=0;i<3;i++) ru0[i]+=eps[n].ru[i];
            have++;
        }
        else if (norm(ru0,3)>0.0&&have>0) {
            for (i=0;i<3;i++) eps[n].ru[i]=ru0[i]/have;
        }
        else {
            for (i=0;i<3;i++) eps[n].ru[i]=opt->ru[i];
        }
        iu+=nu;
        n++;
    }
    if (have<=0&&norm(opt->ru,3)<=0.0) {
        free(eps);
        return 0;
    }
    if (have>0) for (i=0;i<3;i++) ru0[i]/=have;
    else for (i=0;i<3;i++) ru0[i]=opt->ru[i];
    for (i=0;i<n;i++) if (norm(eps[i].ru,3)<=0.0) {
        for (j=0;j<3;j++) eps[i].ru[j]=ru0[j];
    }
    *epochs=eps; *nep=n;
    return n>0;
}
static int scan_rows_kin(const obs_t *obs, const nav_t *nav,
                         const prcopt_t *opt, const postls_epoch_t *eps,
                         int nep, const int refsat[POSTLS_NSYS][NFREQ],
                         const int *arc, postls_amb_t *ambs, int *namb)
{
    obsd_t data[MAXOBS*2];
    double *rs,*dts,*vare,yu,yb,yr,ybr,e[3],az[2],sig;
    int *svh;
    int ie,i,f,sat,bi,ri,rbi,sys,si,n,nv=0,nf=postls_nf(opt);
    int aur,abr,arr,arb,ia;

    rs=mat(6,MAXOBS*2); dts=mat(2,MAXOBS*2);
    vare=mat(1,MAXOBS*2); svh=imat(1,MAXOBS*2);
    for (ie=0;ie<nep;ie++) {
        for (i=0;i<eps[ie].nu;i++) data[i]=obs->data[eps[ie].iu+i];
        for (i=0;i<eps[ie].nr;i++) data[eps[ie].nu+i]=obs->data[eps[ie].ir+i];
        n=eps[ie].nu+eps[ie].nr;
        satposs(data[0].time,data,n,nav,opt->sateph,rs,dts,vare,svh);
        for (i=0;i<eps[ie].nu;i++) {
            sat=data[i].sat;
            if ((bi=find_sat_obs(data+eps[ie].nu,eps[ie].nr,sat))<0) continue;
            if ((sys=satsys(sat,NULL))==0||(si=postls_sysidx(sys))<0) continue;
            for (f=0;f<nf;f++) {
                if (refsat[si][f]<=0||sat==refsat[si][f]) continue;
                if ((ri=find_sat_obs(data,eps[ie].nu,refsat[si][f]))<0) continue;
                if ((rbi=find_sat_obs(data+eps[ie].nu,eps[ie].nr,refsat[si][f]))<0) continue;
                rbi+=eps[ie].nu;
                if (!modeled_obs(data+i,rs+i*6,dts+i*2,vare[i],svh[i],nav,opt,
                                 eps[ie].ru,1,f,0,&yu,e,az,&sig)) continue;
                if (!modeled_obs(data+eps[ie].nu+bi,rs+(eps[ie].nu+bi)*6,
                                 dts+(eps[ie].nu+bi)*2,vare[eps[ie].nu+bi],
                                 svh[eps[ie].nu+bi],nav,opt,opt->rb,2,f,0,
                                 &yb,e,az,&sig)) continue;
                if (!modeled_obs(data+ri,rs+ri*6,dts+ri*2,vare[ri],svh[ri],
                                 nav,opt,eps[ie].ru,1,f,0,&yr,e,az,&sig)) continue;
                if (!modeled_obs(data+rbi,rs+rbi*6,dts+rbi*2,vare[rbi],svh[rbi],
                                 nav,opt,opt->rb,2,f,0,&ybr,e,az,&sig)) continue;
                nv++;
                if (modeled_obs(data+i,rs+i*6,dts+i*2,vare[i],svh[i],nav,opt,
                                eps[ie].ru,1,f,1,&yu,e,az,&sig)&&
                    modeled_obs(data+eps[ie].nu+bi,rs+(eps[ie].nu+bi)*6,
                                dts+(eps[ie].nu+bi)*2,vare[eps[ie].nu+bi],
                                svh[eps[ie].nu+bi],nav,opt,opt->rb,2,f,1,
                                &yb,e,az,&sig)&&
                    modeled_obs(data+ri,rs+ri*6,dts+ri*2,vare[ri],svh[ri],
                                nav,opt,eps[ie].ru,1,f,1,&yr,e,az,&sig)&&
                    modeled_obs(data+rbi,rs+rbi*6,dts+rbi*2,vare[rbi],svh[rbi],
                                nav,opt,opt->rb,2,f,1,&ybr,e,az,&sig)) {
                    aur=arc[(eps[ie].iu+i)*nf+f];
                    abr=arc[(eps[ie].ir+bi)*nf+f];
                    arr=arc[(eps[ie].iu+ri)*nf+f];
                    arb=arc[(eps[ie].ir+rbi-eps[ie].nu)*nf+f];
                    if (aur>0&&abr>0&&arr>0&&arb>0&&
                        (ia=add_amb(ambs,namb,sat,f,aur,abr,arr,arb,
                                    0,0,0,0))>=0) {
                        ambs[ia].n++;
                        nv++;
                    }
                }
            }
        }
    }
    free(rs); free(dts); free(vare); free(svh);
    return nv;
}
static int fill_system_kin(const obs_t *obs, const nav_t *nav,
                           const prcopt_t *opt, const postls_epoch_t *eps,
                           int nep, const int refsat[POSTLS_NSYS][NFREQ],
                           const int *arc, const postls_amb_t *ambs, int namb,
                           int nx, double *A, double *v, int *row_amb,
                           int *row_epoch, int *usedsat)
{
    obsd_t data[MAXOBS*2];
    double *rs,*dts,*vare,yus,ybs,yur,ybr,sus,sbs,sur,sbr;
    double eus[3],ebs[3],eur[3],ebr[3],az[2],sig,dd;
    int *svh;
    int ie,i,k,f,sat,bi,ri,rbi,sys,si,n,nv=0,nf=postls_nf(opt);
    int aur,abr,arr,arb,ia;

    rs=mat(6,MAXOBS*2); dts=mat(2,MAXOBS*2);
    vare=mat(1,MAXOBS*2); svh=imat(1,MAXOBS*2);
    for (ie=0;ie<nep;ie++) {
        for (i=0;i<eps[ie].nu;i++) data[i]=obs->data[eps[ie].iu+i];
        for (i=0;i<eps[ie].nr;i++) data[eps[ie].nu+i]=obs->data[eps[ie].ir+i];
        n=eps[ie].nu+eps[ie].nr;
        satposs(data[0].time,data,n,nav,opt->sateph,rs,dts,vare,svh);
        for (i=0;i<eps[ie].nu;i++) {
            sat=data[i].sat;
            if ((bi=find_sat_obs(data+eps[ie].nu,eps[ie].nr,sat))<0) continue;
            if ((sys=satsys(sat,NULL))==0||(si=postls_sysidx(sys))<0) continue;
            for (f=0;f<nf;f++) {
                if (refsat[si][f]<=0||sat==refsat[si][f]) continue;
                if ((ri=find_sat_obs(data,eps[ie].nu,refsat[si][f]))<0) continue;
                if ((rbi=find_sat_obs(data+eps[ie].nu,eps[ie].nr,refsat[si][f]))<0) continue;
                rbi+=eps[ie].nu;
                if (!modeled_obs(data+i,rs+i*6,dts+i*2,vare[i],svh[i],nav,opt,
                                 eps[ie].ru,1,f,0,&yus,eus,az,&sus)) continue;
                if (!modeled_obs(data+eps[ie].nu+bi,rs+(eps[ie].nu+bi)*6,
                                 dts+(eps[ie].nu+bi)*2,vare[eps[ie].nu+bi],
                                 svh[eps[ie].nu+bi],nav,opt,opt->rb,2,f,0,
                                 &ybs,ebs,az,&sbs)) continue;
                if (!modeled_obs(data+ri,rs+ri*6,dts+ri*2,vare[ri],svh[ri],
                                 nav,opt,eps[ie].ru,1,f,0,&yur,eur,az,&sur)) continue;
                if (!modeled_obs(data+rbi,rs+rbi*6,dts+rbi*2,vare[rbi],svh[rbi],
                                 nav,opt,opt->rb,2,f,0,&ybr,ebr,az,&sbr)) continue;
                sig=SQRT(SQR(sus)+SQR(sbs)+SQR(sur)+SQR(sbr));
                if (sig<=0.0) sig=1.0;
                dd=(yus-ybs)-(yur-ybr);
                v[nv]=dd/sig;
                for (k=0;k<nx;k++) A[k+nv*nx]=0.0;
                for (k=0;k<3;k++) A[eps[ie].index+k+nv*nx]=(-eus[k]+eur[k])/sig;
                if (row_amb) row_amb[nv]=-1;
                if (row_epoch) row_epoch[nv]=ie;
                if (usedsat) {
                    usedsat[ie*MAXSAT+sat-1]=1;
                    usedsat[ie*MAXSAT+refsat[si][f]-1]=1;
                }
                nv++;
                if (!modeled_obs(data+i,rs+i*6,dts+i*2,vare[i],svh[i],nav,opt,
                                 eps[ie].ru,1,f,1,&yus,eus,az,&sus)) continue;
                if (!modeled_obs(data+eps[ie].nu+bi,rs+(eps[ie].nu+bi)*6,
                                 dts+(eps[ie].nu+bi)*2,vare[eps[ie].nu+bi],
                                 svh[eps[ie].nu+bi],nav,opt,opt->rb,2,f,1,
                                 &ybs,ebs,az,&sbs)) continue;
                if (!modeled_obs(data+ri,rs+ri*6,dts+ri*2,vare[ri],svh[ri],
                                 nav,opt,eps[ie].ru,1,f,1,&yur,eur,az,&sur)) continue;
                if (!modeled_obs(data+rbi,rs+rbi*6,dts+rbi*2,vare[rbi],svh[rbi],
                                 nav,opt,opt->rb,2,f,1,&ybr,ebr,az,&sbr)) continue;
                aur=arc[(eps[ie].iu+i)*nf+f];
                abr=arc[(eps[ie].ir+bi)*nf+f];
                arr=arc[(eps[ie].iu+ri)*nf+f];
                arb=arc[(eps[ie].ir+rbi-eps[ie].nu)*nf+f];
                if (aur<=0||abr<=0||arr<=0||arb<=0) continue;
                if ((ia=find_amb((postls_amb_t *)ambs,namb,sat,f,aur,abr,arr,arb,
                                 0,0,0,0))<0) continue;
                if (!ambs[ia].use) continue;
                sig=SQRT(SQR(sus)+SQR(sbs)+SQR(sur)+SQR(sbr));
                if (sig<=0.0) sig=0.02;
                dd=(yus-ybs)-(yur-ybr);
                v[nv]=dd/sig;
                for (k=0;k<nx;k++) A[k+nv*nx]=0.0;
                for (k=0;k<3;k++) A[eps[ie].index+k+nv*nx]=(-eus[k]+eur[k])/sig;
                A[ambs[ia].index+nv*nx]=1.0/sig;
                if (row_amb) row_amb[nv]=ia;
                if (row_epoch) row_epoch[nv]=ie;
                if (usedsat) {
                    usedsat[ie*MAXSAT+sat-1]=1;
                    usedsat[ie*MAXSAT+refsat[si][f]-1]=1;
                }
                nv++;
            }
        }
    }
    free(rs); free(dts); free(vare); free(svh);
    return nv;
}
#define POSTLS_KINWIN 60.0

static int postls_relpos_kin_all(const obs_t *obs, const nav_t *nav,
                                 const prcopt_t *opt, solbuf_t *solbuf,
                                 char *msg)
{
    prcopt_t opt_=*opt,iteropt;
    postls_epoch_t *eps=NULL;
    postls_amb_t *ambs=NULL;
    double ru0[3],*A=NULL,*v=NULL,*dx=NULL,*Q=NULL,dmax;
    int refsat[POSTLS_NSYS][NFREQ];
    int *arc=NULL,*row_amb=NULL,*row_epoch=NULL,*usedsat=NULL;
    int i,j,ie,nf,nep=0,namb=0,nuseamb,nx,nv,maxnv,iter,info,nout,nout_tot=0;

    if (msg) *msg='\0';
    if (!obs||!nav||!opt||!solbuf) return 0;
    if (norm(opt_.rb,3)<=0.0) {
        if (msg) strcpy(msg,"base position is required");
        return 0;
    }
    if (!opt_.navsys) opt_.navsys=SYS_GPS|SYS_GLO|SYS_GAL|SYS_QZS|SYS_CMP|SYS_IRN;
    nf=postls_nf(&opt_);
    if (!build_epochs(obs,nav,&opt_,&eps,&nep,ru0,
                      opt_.maxtdiff>0.0?opt_.maxtdiff:DTTOL)) {
        if (msg) strcpy(msg,"could not initialize kinematic epochs");
        return 0;
    }
    arc=imat(obs->n,nf);
    ambs=(postls_amb_t *)calloc(POSTLS_MAXAMB,sizeof(postls_amb_t));
    usedsat=imat(MAXSAT,nep);
    if (!arc||!ambs||!usedsat) {
        if (msg) strcpy(msg,"memory allocation error");
        free(eps); free(arc); free(ambs); free(usedsat);
        return 0;
    }
    build_phase_arcs(obs,&opt_,nf,arc);
    select_refs(obs,nav,&opt_,ru0,refsat,opt_.maxtdiff>0.0?opt_.maxtdiff:DTTOL);
    maxnv=scan_rows_kin(obs,nav,&opt_,eps,nep,refsat,arc,ambs,&namb);
    nuseamb=prepare_ambs(ambs,namb);
    nx=3*nep+nuseamb;
    if (maxnv<nx||maxnv<POSTLS_MINOBS) {
        if (msg) sprintf(msg,"not enough kinematic DD observations (nv=%d nx=%d epochs=%d)",
                         maxnv,nx,nep);
        free(eps); free(arc); free(ambs); free(usedsat);
        return 0;
    }
    for (j=0,i=0;i<namb;i++) if (ambs[i].use) ambs[i].index=3*nep+j++;
    A=mat(nx,maxnv); v=mat(maxnv,1); dx=mat(nx,1); Q=mat(nx,nx);
    row_amb=imat(maxnv,1); row_epoch=imat(maxnv,1);
    if (!A||!v||!dx||!Q||!row_amb||!row_epoch) {
        if (msg) strcpy(msg,"memory allocation error");
        free(A); free(v); free(dx); free(Q); free(row_amb); free(row_epoch);
        free(eps); free(arc); free(ambs); free(usedsat);
        return 0;
    }
    for (iter=0;iter<POSTLS_MAXITR;iter++) {
        iteropt=opt_;
        if (SMARTPHONE(&iteropt)&&iter>0) for (i=0;i<NFREQ;i++) {
            if (iteropt.eratio[i]<1000.0) iteropt.eratio[i]=1000.0;
        }
        for (i=0;i<nep*MAXSAT;i++) usedsat[i]=0;
        nuseamb=prepare_ambs(ambs,namb);
        for (j=0,i=0;i<namb;i++) if (ambs[i].use) ambs[i].index=3*nep+j++;
        nx=3*nep+nuseamb;
        nv=fill_system_kin(obs,nav,&iteropt,eps,nep,refsat,arc,ambs,namb,nx,
                           A,v,row_amb,row_epoch,usedsat);
        if (nv<nx) {
            if (msg) sprintf(msg,"rank-deficient kinematic DD system (nv=%d nx=%d epochs=%d)",
                             nv,nx,nep);
            free(A); free(v); free(dx); free(Q); free(row_amb); free(row_epoch);
            free(eps); free(arc); free(ambs); free(usedsat);
            return 0;
        }
        if ((info=robust_lsq(A,v,nx,nv,dx,Q,&nout))) {
            if (msg) sprintf(msg,"kinematic least-squares error info=%d nv=%d nx=%d",
                             info,nv,nx);
            free(A); free(v); free(dx); free(Q); free(row_amb); free(row_epoch);
            free(eps); free(arc); free(ambs); free(usedsat);
            return 0;
        }
        nout_tot+=nout;
        dmax=0.0;
        for (ie=0;ie<nep;ie++) {
            for (i=0;i<3;i++) eps[ie].ru[i]+=dx[eps[ie].index+i];
            dmax=MAX_(dmax,norm(dx+eps[ie].index,3));
        }
        if (dmax<1E-4) break;
    }
    for (ie=0;ie<nep;ie++) for (i=0;i<MAXSAT;i++) {
        if (usedsat[ie*MAXSAT+i]) eps[ie].ns++;
    }
    for (ie=0;ie<nep;ie++) {
        sol_t sol={{0}};
        sol.time=eps[ie].time;
        for (i=0;i<3;i++) sol.rr[i]=eps[ie].ru[i];
        sol.qr[0]=(float)Q[eps[ie].index+eps[ie].index*nx];
        sol.qr[1]=(float)Q[eps[ie].index+1+(eps[ie].index+1)*nx];
        sol.qr[2]=(float)Q[eps[ie].index+2+(eps[ie].index+2)*nx];
        sol.qr[3]=(float)Q[eps[ie].index+(eps[ie].index+1)*nx];
        sol.qr[4]=(float)Q[eps[ie].index+1+(eps[ie].index+2)*nx];
        sol.qr[5]=(float)Q[eps[ie].index+2+eps[ie].index*nx];
        sol.type=0;
        sol.stat=SOLQ_FLOAT;
        sol.ns=(uint8_t)MIN_(eps[ie].ns,255);
        if (!addsol(solbuf,&sol)) {
            if (msg) strcpy(msg,"solution buffer allocation error");
            free(A); free(v); free(dx); free(Q); free(row_amb); free(row_epoch);
            free(eps); free(arc); free(ambs); free(usedsat);
            return 0;
        }
    }
    if (msg) sprintf(msg,"postls float kinematic: epochs=%d nv=%d nx=%d namb=%d/%d out=%d",
                     nep,nv,nx,nuseamb,namb,nout_tot);
    free(A); free(v); free(dx); free(Q); free(row_amb); free(row_epoch);
    free(eps); free(arc); free(ambs); free(usedsat);
    return 1;
}
static int slice_obs_time(const obs_t *obs, gtime_t ts, gtime_t te,
                          double pad, obs_t *sub)
{
    obsd_t *data;
    double dt1,dt2;
    int i,n=0;

    for (i=0;i<obs->n;i++) {
        dt1=timediff(obs->data[i].time,ts);
        dt2=timediff(obs->data[i].time,te);
        if ((obs->data[i].rcv==1&&dt1>=-DTTOL&&dt2<-DTTOL)||
            (obs->data[i].rcv==2&&dt1>=-pad-DTTOL&&dt2<=pad+DTTOL)) {
            n++;
        }
    }
    memset(sub,0,sizeof(*sub));
    if (n<=0) return 0;
    if (!(data=(obsd_t *)malloc(sizeof(obsd_t)*n))) return 0;
    for (i=n=0;i<obs->n;i++) {
        dt1=timediff(obs->data[i].time,ts);
        dt2=timediff(obs->data[i].time,te);
        if ((obs->data[i].rcv==1&&dt1>=-DTTOL&&dt2<-DTTOL)||
            (obs->data[i].rcv==2&&dt1>=-pad-DTTOL&&dt2<=pad+DTTOL)) {
            data[n++]=obs->data[i];
        }
    }
    sub->data=data;
    sub->n=sub->nmax=n;
    return 1;
}
extern int postls_relpos_kin(const obs_t *obs, const nav_t *nav,
                             const prcopt_t *opt, solbuf_t *solbuf, char *msg)
{
    obs_t sub={0};
    solbuf_t winbuf={0};
    gtime_t ts={0},te={0},tend={0};
    double pad;
    int i,j,nwin=0,nfail=0,nsol0;
    char wmsg[MAXERRMSG];

    if (msg) *msg='\0';
    if (!obs||obs->n<=0||!nav||!opt||!solbuf) return 0;
    pad=opt->maxtdiff>0.0?opt->maxtdiff:DTTOL;
    for (i=0;i<obs->n;i++) if (obs->data[i].rcv==1) {
        if (ts.time==0) ts=obs->data[i].time;
        tend=obs->data[i].time;
    }
    if (ts.time==0) {
        if (msg) strcpy(msg,"no rover observations");
        return 0;
    }
    while (timediff(tend,ts)>=-DTTOL) {
        te=timeadd(ts,POSTLS_KINWIN);
        if (!slice_obs_time(obs,ts,te,pad,&sub)) {
            ts=te;
            continue;
        }
        memset(&winbuf,0,sizeof(winbuf));
        nsol0=winbuf.n;
        if (postls_relpos_kin_all(&sub,nav,opt,&winbuf,wmsg)) {
            for (j=0;j<winbuf.n;j++) addsol(solbuf,winbuf.data+j);
            nwin++;
        }
        else {
            nfail++;
        }
        if (winbuf.n==nsol0&&sub.n>0) {
            trace(3,"postls kinematic window skipped: %s\n",wmsg);
        }
        freesolbuf(&winbuf);
        free(sub.data); sub.data=NULL; sub.n=sub.nmax=0;
        ts=te;
    }
    if (solbuf->n<=0) {
        if (msg) sprintf(msg,"no kinematic solutions (windows failed=%d)",nfail);
        return 0;
    }
    if (msg) sprintf(msg,"postls float kinematic-window: epochs=%d windows=%d failed=%d dt=%.1f",
                     solbuf->n,nwin,nfail,pad);
    return 1;
}
