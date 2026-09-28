/* Regression: partial PPP-AR cursor bounds and solution initialization.
 * Link with the production RTKLIB sources (see run_memory_regressions.py).
 * ppp_ar() takes its existing thresar[0]<1 early return, so the test exercises
 * the real manager without needing an integer-ambiguity solution fixture.
 */
#include <assert.h>
#include <limits.h>
#include "rtklib.h"

static rtk_t rtk;
static obsd_t obs[2];

static void init_solution(void)
{
    prcopt_t opt=prcopt_default;
    opt.mode=PMODE_PPP_STATIC;
    opt.ionoopt=IONOOPT_IFLC;
    opt.nf=2;
    opt.thresar[0]=0.5;
    opt.mindropsats=1;
    opt.elmaskar=0.0;
    rtkinit(&rtk,&opt);
    assert(rtk.excsats==0);
}

static void attempt(void)
{
    int exc[2]={0};
    rtk.nb_ar=10;
    rtk.sol.prevf_ratio=1.0;
    rtk.sol.thres=3.0;
    assert(manage_ppp_ar(&rtk,NULL,NULL,NULL,1,obs,2,NULL,exc)==0);
    assert(rtk.ssat[0].lock[0]==7);
    assert(rtk.ssat[1].lock[0]==9);
}

static void bounds(int cursor,int expected)
{
    rtk.excsats=cursor;
    attempt();
    assert(rtk.excsats==expected);
}

int main(void)
{
    /* A caller need not pre-zero rtk_t; initialization owns the cursor. */
    memset(&rtk,0xA5,sizeof(rtk));
    init_solution();
    obs[0].sat=1; obs[1].sat=2;
    rtk.ssat[0].vsat[0]=rtk.ssat[1].vsat[0]=1;
    rtk.ssat[0].lock[0]=7; rtk.ssat[1].lock[0]=9;
    bounds(0,1);
    puts("PASS excsats=0: exclude, restore, advance");
    /* Cursor must persist between attempts, not reset at every epoch. */
    attempt(); assert(rtk.excsats==2);
    attempt(); assert(rtk.excsats==0);
    puts("PASS repeated AR attempts: 0 -> 1 -> 2 -> 0");
    bounds(-1,0); bounds(INT_MIN,0);
    puts("PASS negative excsats: no indexing, reset");
    bounds(2,0);
    puts("PASS excsats==ar: no indexing, reset");
    bounds(3,0); bounds(INT_MAX,0);
    puts("PASS oversized excsats: no indexing, reset");
    rtk.ssat[0].vsat[0]=rtk.ssat[1].vsat[0]=0;
    bounds(0,0);
    puts("PASS empty candidate list");
    rtk.excsats=INT_MAX;
    rtkfree(&rtk);
    init_solution();
    rtkfree(&rtk);
    puts("PASS free/reinitialize rtk_t resets cursor");
    return 0;
}
