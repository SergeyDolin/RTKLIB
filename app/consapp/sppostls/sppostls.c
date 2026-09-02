/*------------------------------------------------------------------------------
* sppostls.c : command line front-end for smartphone batch LS post-processing
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

#define PROGNAME "sppostls"
#define MAXINFILE 128

static const char *help[]={
"",
" usage: sppostls [options] rover.obs base.obs navfile...",
"",
" options:",
" -k file       input options file",
" -o file       output solution file",
" -r x y z      base station ECEF position (m)",
" -l lat lon h  base station latitude/longitude/height (deg,m)",
" -f freq       frequency mode (1:L1,2:L1+L2,3:L1+L2+L5,6:phone L1+L5)",
" -m mask       elevation mask (deg)",
" -sys systems  navigation systems, comma separated letters G,R,E,J,C,I",
" -precise      use SP3/CLK precise ephemerides if supplied",
" -e           output ECEF XYZ",
" -t           output calendar time",
" -x level      trace level",
""
};
static void printhelp(void)
{
    int i;

    for (i=0;i<(int)(sizeof(help)/sizeof(*help));i++) fprintf(stderr,"%s\n",help[i]);
    exit(0);
}
static void freeobsnav(obs_t *obs, nav_t *nav)
{
    free(obs->data); obs->data=NULL; obs->n=obs->nmax=0;
    freenav(nav,0);
}
static int extmatch(const char *file, const char *ext)
{
    const char *p=strrchr(file,'.');

    if (!p) return 0;
    return !strcmp(p,ext);
}
static int readnavfile(const char *file, const prcopt_t *opt, nav_t *nav)
{
    obs_t obs={0};
    sta_t sta={0};
    gtime_t ts={0},te={0};
    int ret;

    if (extmatch(file,".sp3")||extmatch(file,".SP3")||
        extmatch(file,".eph")||extmatch(file,".EPH")) {
        readsp3(file,nav,0);
        return 1;
    }
    if (extmatch(file,".clk")||extmatch(file,".CLK")) {
        readrnxc(file,nav);
        return 1;
    }
    ret=readrnxt(file,0,ts,te,0.0,opt->rnxopt[0],&obs,nav,&sta);
    free(obs.data);
    return ret>=0;
}
static int readobsnav_simple(char **infile, int n, const prcopt_t *opt,
                             obs_t *obs, nav_t *nav, sta_t *sta)
{
    gtime_t ts={0},te={0};
    int i,rcv,hasnav=0;

    memset(obs,0,sizeof(*obs));
    memset(nav,0,sizeof(*nav));
    memset(sta,0,sizeof(sta_t)*MAXRCV);

    for (i=0;i<2;i++) {
        rcv=i==0?1:(i==1?2:0);
        if (readrnxt(infile[i],rcv,ts,te,0.0,
                     opt->rnxopt[rcv<=1?0:1],obs,nav,rcv?sta+rcv-1:NULL)<0) {
            return 0;
        }
    }
    for (i=2;i<n;i++) {
        if (!readnavfile(infile[i],opt,nav)) return 0;
        hasnav=1;
    }
    if (obs->n<=0||sortobs(obs)<=0) return 0;
    if (!hasnav||(nav->n<=0&&nav->ng<=0&&nav->ns<=0&&nav->ne<=0)) return 0;
    uniqnav(nav);
    return 1;
}
static void set_ant_delta(prcopt_t *opt, const sta_t *sta, int rcv)
{
    double pos[3],del[3],ref[3];
    int i=rcv-1,j;

    if (sta[i].deltype==1) {
        for (j=0;j<3;j++) ref[j]=norm(sta[i].pos,3)>0.0?sta[i].pos[j]:
                                (i==0?opt->ru[j]:opt->rb[j]);
        if (norm(ref,3)>0.0) {
            ecef2pos(ref,pos);
            ecef2enu(pos,sta[i].del,del);
            for (j=0;j<3;j++) opt->antdel[i][j]=del[j];
        }
    }
    else {
        for (j=0;j<3;j++) opt->antdel[i][j]=sta[i].del[j];
    }
}
static int setup_pcv(gtime_t time, const filopt_t *fopt, prcopt_t *opt,
                     nav_t *nav, const sta_t *sta, pcvs_t *satpcvs,
                     pcvs_t *rcvpcvs, char *msg)
{
    pcv_t *pcv,pcv0={0};
    char id[64];
    int i;

    if (*fopt->satantp&&!readpcv(fopt->satantp,satpcvs)) {
        if (msg) sprintf(msg,"satellite antenna PCV read error: %s",fopt->satantp);
        return 0;
    }
    if (*fopt->rcvantp&&!readpcv(fopt->rcvantp,rcvpcvs)) {
        if (msg) sprintf(msg,"receiver antenna PCV read error: %s",fopt->rcvantp);
        return 0;
    }
    for (i=0;i<MAXSAT;i++) {
        nav->pcvs[i]=pcv0;
        if (!(satsys(i+1,NULL)&opt->navsys)) continue;
        if (!(pcv=searchpcv(i+1,"",time,satpcvs))) {
            satno2id(i+1,id);
            trace(3,"no satellite antenna pcv: %s\n",id);
            continue;
        }
        nav->pcvs[i]=*pcv;
    }
    for (i=0;i<2;i++) {
        opt->pcvr[i]=pcv0;
        if ((!*opt->anttype[i]||!strcmp(opt->anttype[i],"*"))&&
            *sta[i].antdes&&strcmp(sta[i].antdes,"NA")) {
            strcpy(opt->anttype[i],sta[i].antdes);
        }
        set_ant_delta(opt,sta,i+1);
        if (!*opt->anttype[i]||!*fopt->rcvantp) continue;
        if (!(pcv=searchpcv(0,opt->anttype[i],time,rcvpcvs))) {
            trace(2,"no receiver antenna pcv: %s\n",opt->anttype[i]);
            continue;
        }
        strcpy(opt->anttype[i],pcv->type);
        opt->pcvr[i]=*pcv;
    }
    return 1;
}
int main(int argc, char **argv)
{
    prcopt_t prcopt=prcopt_default;
    solopt_t solopt=solopt_default;
    filopt_t filopt={0};
    sta_t sta[MAXRCV];
    pcvs_t satpcvs={0},rcvpcvs={0};
    obs_t *obs;
    nav_t *nav;
    solbuf_t solbuf={0};
    sol_t sol={{0}};
    FILE *fp=stdout;
    double pos[3];
    char *infile[MAXINFILE],*outfile="",*p,msg[MAXERRMSG]="";
    int i,j,n=0,ret=0,precise=0;

    if (!(obs=(obs_t *)calloc(1,sizeof(obs_t)))||
        !(nav=(nav_t *)calloc(1,sizeof(nav_t)))) {
        fprintf(stderr,"error: memory allocation\n");
        return -1;
    }

    prcopt.mode=PMODE_STATIC;
    prcopt.navsys=SYS_GPS|SYS_GLO|SYS_GAL|SYS_QZS|SYS_CMP|SYS_IRN;
    prcopt.sateph=EPHOPT_BRDC;
    prcopt.ionoopt=IONOOPT_BRDC;
    prcopt.tropopt=TROPOPT_SAAS;
    prcopt.modear=ARMODE_OFF;
    prcopt.nf=6;
    solopt.timef=0;
    solopt.posf=SOLF_LLH;
    sprintf(solopt.prog,"%s ver.%s %s",PROGNAME,VER_RTKLIB,PATCH_LEVEL);

    for (i=1;i<argc;i++) {
        if (!strcmp(argv[i],"-k")&&i+1<argc) {
            resetsysopts();
            if (!loadopts(argv[++i],sysopts)) return -1;
            getsysopts(&prcopt,&solopt,&filopt);
        }
    }
    for (i=1;i<argc;i++) {
        if      (!strcmp(argv[i],"-o")&&i+1<argc) outfile=argv[++i];
        else if (!strcmp(argv[i],"-k")&&i+1<argc) {i++; continue;}
        else if (!strcmp(argv[i],"-f")&&i+1<argc) prcopt.nf=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-m")&&i+1<argc) prcopt.elmin=atof(argv[++i])*D2R;
        else if (!strcmp(argv[i],"-e")) solopt.posf=SOLF_XYZ;
        else if (!strcmp(argv[i],"-t")) solopt.timef=1;
        else if (!strcmp(argv[i],"-precise")) precise=1;
        else if (!strcmp(argv[i],"-x")&&i+1<argc) solopt.trace=atoi(argv[++i]);
        else if (!strcmp(argv[i],"-r")&&i+3<argc) {
            for (j=0;j<3;j++) prcopt.rb[j]=atof(argv[++i]);
        }
        else if (!strcmp(argv[i],"-l")&&i+3<argc) {
            for (j=0;j<3;j++) pos[j]=atof(argv[++i]);
            pos[0]*=D2R; pos[1]*=D2R;
            pos2ecef(pos,prcopt.rb);
        }
        else if (!strcmp(argv[i],"-sys")&&i+1<argc) {
            prcopt.navsys=0;
            for (p=argv[++i];*p;p++) {
                if      (*p=='G') prcopt.navsys|=SYS_GPS;
                else if (*p=='R') prcopt.navsys|=SYS_GLO;
                else if (*p=='E') prcopt.navsys|=SYS_GAL;
                else if (*p=='J') prcopt.navsys|=SYS_QZS;
                else if (*p=='C') prcopt.navsys|=SYS_CMP;
                else if (*p=='I') prcopt.navsys|=SYS_IRN;
            }
        }
        else if (*argv[i]=='-') printhelp();
        else if (n<MAXINFILE) infile[n++]=argv[i];
    }
    if (n<3) printhelp();
    if (solopt.trace>0) {
        traceopen(PROGNAME ".trace");
        tracelevel(solopt.trace);
    }
    if (!readobsnav_simple(infile,n,&prcopt,obs,nav,sta)) {
        fprintf(stderr,"error: failed to read observations/navigation\n");
        ret=-1;
    }
    else {
        if (precise&&nav->ne>0) prcopt.sateph=EPHOPT_PREC;
        if (!setup_pcv(obs->data[0].time,&filopt,&prcopt,nav,sta,&satpcvs,
                       &rcvpcvs,msg)) {
            fprintf(stderr,"error: %s\n",msg);
            ret=-1;
        }
        else if (prcopt.mode==PMODE_KINEMA) {
            if (!postls_relpos_kin(obs,nav,&prcopt,&solbuf,msg)) {
                if (prcopt.sateph==EPHOPT_PREC) {
                    prcopt.sateph=EPHOPT_BRDC;
                    if (!postls_relpos_kin(obs,nav,&prcopt,&solbuf,msg)) {
                        fprintf(stderr,"error: %s\n",msg);
                        ret=-1;
                    }
                    else strcat(msg," (precise unavailable, used broadcast)");
                }
                else {
                    fprintf(stderr,"error: %s\n",msg);
                    ret=-1;
                }
            }
        }
        else if (!postls_relpos(obs,nav,&prcopt,&sol,msg)) {
            if (prcopt.sateph==EPHOPT_PREC) {
                prcopt.sateph=EPHOPT_BRDC;
                if (!postls_relpos(obs,nav,&prcopt,&sol,msg)) {
                    fprintf(stderr,"error: %s\n",msg);
                    ret=-1;
                }
                else strcat(msg," (precise unavailable, used broadcast)");
            }
            else {
                fprintf(stderr,"error: %s\n",msg);
                ret=-1;
            }
        }
    }
    if (!ret) {
        if (*outfile&&!(fp=fopen(outfile,"w"))) {
            fprintf(stderr,"error: output open error %s\n",outfile);
            ret=-1;
        }
        else {
            solopt.pmode=prcopt.mode;
            outsolhead(fp,&solopt);
            if (prcopt.mode==PMODE_KINEMA) {
                for (i=0;i<solbuf.n;i++) outsol(fp,solbuf.data+i,prcopt.rb,&solopt);
            }
            else {
                outsol(fp,&sol,prcopt.rb,&solopt);
            }
            fprintf(fp,"%s postls : %s\n",COMMENTH,msg);
            if (fp!=stdout) fclose(fp);
        }
    }
    if (solopt.trace>0) traceclose();
    freesolbuf(&solbuf);
    free(satpcvs.pcv); free(rcvpcvs.pcv);
    freeobsnav(obs,nav);
    free(obs); free(nav);
    return ret;
}
