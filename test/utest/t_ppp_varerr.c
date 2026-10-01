/* PPP observation-variance regression: tests the production private varerr().
 * Run from the repository root (optionally add -fsanitize=address,undefined):
 * gcc -O1 -g -ansi -ffunction-sections -fdata-sections -Isrc \
 *   -DNFREQ=5 -DNEXOBS=2 -DENAGLO -DENAQZS -DENAGAL -DENACMP -DENAIRN \
 *   test/utest/t_ppp_varerr.c -Wl,--gc-sections -lm -o /tmp/t_ppp_varerr
 * /tmp/t_ppp_varerr
 */
#include <assert.h>
#include "../../src/ppp.c"

static void close_to(double actual,double expected)
{
    assert(isfinite(actual));
    assert(fabs(actual-expected)<=1e-12*fabs(expected)+1e-20);
}

static prcopt_t options(int iono)
{
    prcopt_t opt={0};
    opt.ionoopt=iono;
    opt.eratio[0]=opt.eratio[1]=100.0;
    opt.err[1]=opt.err[2]=0.003;
    return opt;
}

static void scaling(void)
{
    double elevations[]={10.0,30.0,90.0};
    int modes[]={IONOOPT_OFF,IONOOPT_IFLC};
    prcopt_t opt;
    double v1,v2,el,floor_var;
    int i,j,k;
    for(i=0;i<2;i++) for(j=0;j<3;j++) for(k=0;k<3;k++) {
        opt=options(modes[i]);
        if(k==0) opt.err[1]=0.0; /* zero additive floor: exact total scaling */
        if(k==1) opt.err[2]=0.0; /* constant component only */
        el=elevations[j]*D2R;
        floor_var=SQR(opt.err[1]*0.01*0.0025);
        v1=varerr(1,SYS_GPS,el,0.0,1,&opt,NULL);
        opt.eratio[0]*=2.0; /* GPS code fact: 200 -> 400 */
        v2=varerr(1,SYS_GPS,el,0.0,1,&opt,NULL);
        /* The existing independent additive floor is not scaled by fact. */
        close_to(v2-floor_var,4.0*(v1-floor_var));
        close_to(sqrt(v2-floor_var),2.0*sqrt(v1-floor_var));
    }
    puts("PASS fact x2 -> variance x4, sigma x2: raw/IF, 10/30/90 degrees, both noise components");
}

static void absolute_values(void)
{
    prcopt_t opt=options(IONOOPT_IFLC);
    double v,phase,raw,floor_code=SQR(0.003*0.01*0.0025);
    double floor_phase=SQR(0.003*0.0025*0.004*0.2);
    v=varerr(1,SYS_GPS,PI/2.0,0.0,1,&opt,NULL);
    close_to(v,6.48+floor_code);
    close_to(sqrt(v),sqrt(6.48+floor_code));
    phase=varerr(1,SYS_GPS,PI/2.0,0.0,0,&opt,NULL);
    close_to(phase,0.000162+floor_phase);
    opt.eratio[0]=200.0;
    close_to(varerr(1,SYS_GPS,PI/2.0,0.0,0,&opt,NULL),phase);
    opt=options(IONOOPT_OFF);
    raw=varerr(1,SYS_GPS,PI/2.0,0.0,1,&opt,NULL);
    close_to(raw,0.72+floor_code);
    close_to(v-floor_code,9.0*(raw-floor_code));
    printf("PASS GPS zenith: IF code variance=%.12g m^2 sigma=%.12g m; phase and IF x9 preserved\n",v,sqrt(v));
}

static void snr_scaling(void)
{
    prcopt_t opt=options(IONOOPT_IFLC);
    obsd_t obs={0};
    double base,v1,v2,floor_var=SQR(opt.err[1]*0.01*0.0025);
    base=varerr(1,SYS_GPS,PI/6.0,0.0,1,&opt,NULL);
    obs.SNR[0]=(unsigned short)(40.0/SNR_UNIT);
    close_to(varerr(1,SYS_GPS,PI/6.0,0.0,1,&opt,&obs),base);
    obs.SNR[0]=(unsigned short)(20.0/SNR_UNIT);
    v1=varerr(1,SYS_GPS,PI/6.0,0.0,1,&opt,&obs);
    close_to(v1,10.0*base);
    opt.eratio[0]*=2.0;
    v2=varerr(1,SYS_GPS,PI/6.0,0.0,1,&opt,&obs);
    close_to(v2-10.0*floor_var,4.0*(v1-10.0*floor_var));
    puts("PASS existing SNR rule preserved; low-SNR variance still scales quadratically with fact");
}

int main(void)
{
    scaling();
    absolute_values();
    snr_scaling();
    return 0;
}
