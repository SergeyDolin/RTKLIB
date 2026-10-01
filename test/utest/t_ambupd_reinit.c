/* Production PPP -> AMBUPD boundary regression. Run in an empty directory.
 * gcc -O1 -g -std=gnu99 -ffunction-sections -fdata-sections -Isrc \
 * -DNFREQ=5 -DNEXOBS=2 -DENAGLO -DENAQZS -DENAGAL -DENACMP -DENAIRN \
 * test/utest/t_ambupd_reinit.c src/rtkcmn.c -Wl,--gc-sections -lm -o TEST
 */
#include <assert.h>
#include "../../src/ppp.c"
#include "../../src/ambupd.c"
static rtk_t tr;
static nav_t tn;
static gtime_t start;

static void fixture(int short_fixture)
{
    double ep[]={2020,1,1,0,0,30};
    ambupd_atexit_hook();
    ambupd_reset();
    memset(&g_ambflag,0,sizeof(g_ambflag));
    free(tr.x); free(tr.P);
    memset(&tr,0,sizeof(tr)); memset(&tn,0,sizeof(tn));
    tr.opt.ionoopt=IONOOPT_IFLC; tr.opt.nf=2; tr.opt.maxout=1000;
    tr.opt.elmin=10*D2R; tr.tt=30; tr.sol.rr[0]=RE_WGS84;
    tr.nx=pppnx(&tr.opt); tr.x=zeros(tr.nx,1); tr.P=zeros(tr.nx,tr.nx);
    tr.ssat[0].azel[1]=45*D2R;
    strcpy(tr.opt.station_name,"TEST");
    /* Small A/B examples use the existing test-only gate override.
     * Spike and consecutive-reinit tests exercise the unchanged default 16. */
    if(short_fixture) strcpy(tr.opt.pppopt,"-AMBUPD_MIN_ARC_EPOCHS=2");
    start=epoch2time(ep);
}
static void epoch(int k,double value,int reinit,int gap,int reject)
{
    obsd_t o={0}; int j=IB(1,0,&tr.opt);
    o.time=timeadd(start,30*(k-1)+gap);o.sat=1;
    o.code[0]=CODE_L1W;o.code[1]=CODE_L2W;
    o.P[0]=o.P[1]=CLIGHT;
    o.L[0]=(CLIGHT+value)*FREQ1/CLIGHT;
    o.L[1]=(CLIGHT+value)*FREQ2/CLIGHT;
    tr.ssat[0].slip[0]=reinit;
    udbias_ppp(&tr,&o,1,&tn);
    assert(tr.ssat[0].slip[0]==0&&tr.ssat[0].slip[1]==0);
    assert(fabs(tr.x[j]-value)<1e-6);
    /* The ordinary numerical gate may reject the first new epoch. */
    tr.P[j+j*tr.nx]=reject?NAN:1.0; /* accepted float covariance */
    ambupd_write_epoch(&tr,&o,1,&tn);
    if(reject) tr.P[j+j*tr.nx]=1;
}
static void rows(int first,int n,double value,int first_epoch)
{
    int i;
    for(i=0;i<n;i++) {
        amb_pending_row_t *p=g_pending+first+i;
        assert(fabs(p->row.if_out-value)<1e-6);
        assert(p->row.sod==30*(first_epoch+i));
    }
}
static void finish(int expected_rows,int expected_arcs)
{
    char line[256],station[32],sat[8];double mjd,sod,bc,wl,sd,prev=-1;int n=0;
    ambupd_flush_sat_arc(1,"test-end");
    assert(g_pending_n==expected_rows&&g_arc_kept==(unsigned long)expected_arcs);
    pending_flush_to_file();fflush(g_fp_amb);
    /* Output is opened write-only by production: reopen it for verification. */
    { FILE *fp=fopen(g_amb_file,"r");assert(fp);
      while(fgets(line,sizeof(line),fp)) {
        assert(sscanf(line,"%lf %lf %31s %7s %lf %lf %lf",&mjd,&sod,station,sat,&bc,&wl,&sd)==7);
        assert(sod>prev);prev=sod;n++;
      } fclose(fp); }
    assert(n==expected_rows);
}
int main(void)
{
    int k; ambflag_rec_t flag={1,AMBFLAG_AMB,4,6};
    fixture(1);initx(&tr,10,1,IB(1,0,&tr.opt));
    for(k=1;k<=5;k++)epoch(k,10,0,0,0);
    assert(g_arc_rows_n[0]==5);finish(5,1);puts("PASS A continuous: one arc / five samples");

    fixture(1);for(k=1;k<=3;k++)epoch(k,10,0,0,0);
    epoch(4,50,1,0,0);
    assert(g_pending_n==3&&g_arc_rows_n[0]==1&&wl_arc[0].n==1);
    rows(0,3,10,1);
    for(k=5;k<=6;k++)epoch(k,50,0,0,0);
    ambupd_flush_sat_arc(1,"test-end");rows(3,3,50,4);
    assert(fabs(g_pending[0].wl_mean-g_pending[3].wl_mean)>1);
    finish(6,2);puts("PASS B reinit: separate IF means and WL accumulators");

    fixture(0);for(k=1;k<=16;k++)epoch(k,-3,0,0,0);
    epoch(17,-49,1,0,0);epoch(18,-3,1,0,0);
    for(k=19;k<=33;k++)epoch(k,-3,0,0,0);
    assert(g_drop_arc_sh==1);rows(0,16,-3,1);finish(32,2);
    puts("PASS C spike: old mean unchanged; one-sample arc rejected by default gate");

    fixture(1);epoch(1,10,1,0,0);assert(g_pending_n==0&&g_drop_arc_sh==0);
    epoch(2,10,0,0,0);finish(2,1);puts("PASS D first-epoch reinit: no empty arc");

    fixture(0);for(k=1;k<=5;k++)epoch(k,10*k,1,0,0);
    assert(g_drop_arc_sh==4&&g_arc_rows_n[0]==1&&wl_arc[0].n==1);
    finish(0,0);assert(g_drop_arc_sh==5);puts("PASS E consecutive reinit / EOF");

    fixture(1);for(k=1;k<=3;k++)epoch(k,10,0,0,0);
    epoch(4,50,1,300,0);epoch(5,50,0,300,0);
    assert(g_pending_n==3&&g_arc_rows_n[0]==2);finish(5,2);
    puts("PASS F reinit plus gap: no duplicate flush or epoch");

    fixture(1);for(k=1;k<=3;k++)epoch(k,10,0,0,0);
    g_ambflag.recs=&flag;g_ambflag.n=1;g_ambflag.interval=30;g_ambflag.t0=start;
    epoch(4,50,1,0,0);epoch(5,50,0,0,0);finish(5,2);
    memset(&g_ambflag,0,sizeof(g_ambflag));
    puts("PASS reinit plus existing reset: old arc emitted once");

    fixture(1);for(k=1;k<=3;k++)epoch(k,10,0,0,0);
    epoch(4,50,1,0,1);assert(g_pending_n==3&&g_arc_rows_n[0]==0&&wl_arc[0].n==0);
    epoch(5,50,0,0,0);epoch(6,50,0,0,0);finish(5,2);
    puts("PASS rejected first new epoch cannot erase the boundary");

    fixture(1);for(k=1;k<=3;k++)epoch(k,10,0,0,0);
    tr.opt.maxout=0;epoch(4,50,0,0,0);assert(g_pending_n==3&&g_arc_rows_n[0]==1);
    tr.opt.maxout=1000;epoch(5,50,0,0,0);finish(5,2);
    puts("PASS actual outage reinitialization without slip flag");
    fixture(1);for(k=1;k<=3;k++)epoch(k,10,0,0,0);
    epoch(4,0,1,0,0);
    assert(g_pending_n==3&&g_arc_rows_n[0]==0&&wl_arc[0].n==0);
    finish(3,1);puts("PASS zero-valued initialized bias still emits the boundary");
    ambupd_if_reinit(0);ambupd_if_reinit(MAXSAT+1);
    ambupd_atexit_hook();free(tr.x);free(tr.P);tr.x=tr.P=NULL;
    return 0;
}
