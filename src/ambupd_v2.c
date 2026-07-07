/*------------------------------------------------------------------------------
 * ambupd.c : write ambupd file (GREAT-UPD style) from the PPP state
 * The file contains:
 *     1) Modified Julian Day
 *     2) Time of day, (unit: second)
 *     3) Station name
 *     4) PRN number of GNSS satellite
 *     5) IF ambiguity derived from PPP float solution, (unit: meter)
 *     6) WL ambiguity derived from Melbourne-Wübbena combination, (unit: cycle)
 *     7) Standard deviation of WL ambiguity.
 * The ambupd_write_epoch() function is intended to be called from ppp.c -> pppos() once per epoch. 
 * -----------------------------------------------------------------------------*/
#include "rtklib.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#define SQR(x)      ((x)*(x))
#define SQRT(x)     ((x)<=0.0||(x)!=(x)?0.0:sqrt(x))
#define MAX(x,y)    ((x)>(y)?(x):(y))
#define MIN(x,y)    ((x)<(y)?(x):(y))
#define ROUND(x)    (int)floor((x)+0.5)

/* number and index of states */
#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
#define NP(opt)     ((opt)->dynamics?9:3)
#define NC(opt)     (NSYS)
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt==TROPOPT_EST?1:3))
#define NI(opt)     ((opt)->ionoopt==IONOOPT_EST?MAXSAT:0)
#define ND(opt)     ((opt)->nf>=3?1:0)
#define NR(opt)     (NP(opt)+NC(opt)+NT(opt)+NI(opt)+ND(opt))
#define NB(opt)     (NF(opt)*MAXSAT)
#define NX(opt)     (NR(opt)+NB(opt))
#define IC(s,opt)   (NP(opt)+(s))
#define IT(opt)     (NP(opt)+NC(opt))
#define II(s,opt)   (NP(opt)+NC(opt)+NT(opt)+(s)-1)
#define ID(opt)     (NP(opt)+NC(opt)+NT(opt)+NI(opt))
#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1)
#define IRDCB(s,opt) (NP(opt)+NC(opt)+(s))
#define IIFCB(s,opt) (NP(opt)+NC(opt)+(s))

#define dbgprint(...)  do { printf(__VA_ARGS__); } while(0)

#define GPS_L1_HZ 1575.42e6
#define GPS_L2_HZ 1227.60e6
#define GPS_L5_HZ 1176.45e6
#define F_TOL_HZ  2.0e6

/* get station name and convert to uppercase */
static void get_station_name(const rtk_t *rtk, char *out, int outsz)
{
    const char *src;
    char *p;

    if (rtk && rtk->opt.station_name[0])
        src = rtk->opt.station_name;
    else
        src = "UNKN";

    snprintf(out, outsz, "%s", src);

    p = out;
    while (*p) {
        *p = (char)toupper((unsigned char)*p);
        p++;
    }
}

/* local finite check (avoid isfinite() link issues) */
static int finite_d(double x)
{
    return (x == x) && (x <= DBL_MAX) && (x >= -DBL_MAX);
}

/* find obs of satellite -----------------------------------------------------*/
static const obsd_t *find_obs(const obsd_t *obs, int n, int sat)
{
    int i;
    for (i = 0; i < n; i++) if (obs[i].sat == sat) return obs + i;
    return NULL;
}

/* Return L1/L5 indices (in obs->data arrays) + frequencies ------------------*/
static int get_gps_l1_l5_idx(int sat, const obsd_t *obs, const nav_t *nav,
                            int *kL1, int *kL5, double *f1, double *f5)
{
    if (!obs || !nav || !kL1 || !kL5 || !f1 || !f5) return 0;
    if (satsys(sat, NULL) != SYS_GPS) return 0;

    /* try to detect L1/L5 by code2freq() */
    int i, idx1 = -1, idx5 = -1;

    for (i = 0; i < NFREQ + NEXOBS; i++) {
        unsigned char code = obs->code[i];
        if (!code) continue;

        double fi = code2freq(SYS_GPS, code, nav);
        if (fi <= 0.0) continue;

        if (fabs(fi - FREQ1) < 1e5 && idx1 < 0) idx1 = i;
#ifdef FREQ5
        if (fabs(fi - FREQ5) < 1e5 && idx5 < 0) idx5 = i;
#else
        if (fabs(fi - 1176.45e6) < 1e5 && idx5 < 0) idx5 = i;
#endif
    }

    /* fallback positions if detection failed */
    if (idx1 < 0) idx1 = 0;
    if (idx5 < 0) idx5 = 2;

    *kL1 = idx1;
    *kL5 = idx5;

    *f1 = code2freq(SYS_GPS, obs->code[idx1], nav);
    *f5 = code2freq(SYS_GPS, obs->code[idx5], nav);

    if (*f1 <= 0.0) *f1 = FREQ1;
#ifdef FREQ5
    if (*f5 <= 0.0) *f5 = FREQ5;
#else
    if (*f5 <= 0.0) *f5 = 1176.45e6;
#endif

    return (*f1 > 0.0 && *f5 > 0.0);
}

/* write ambiguity diagnostics per epoch ------------------------------------*/
void ambupd_write_epoch(const rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    if (!rtk || !obs || n <= 0 || !nav) return;

    /* open file per day */
    static FILE *fp = NULL;
    static char curfile[1024] = "";

    double ep[6];
    time2epoch(rtk->sol.time, ep);
    int year = (int)ep[0];
    int doy  = time2doy(rtk->sol.time);

    char stname[32];
    get_station_name(rtk, stname, (int)sizeof(stname));

    char fname[1024];
    snprintf(fname, sizeof(fname), "%s_ambupd_%4d%03d", stname, year, doy);

    if (!fp || strcmp(fname, curfile) != 0) {
        if (fp) fclose(fp);
        fp = fopen(fname, "w");
        if (!fp) {
            trace(1, "ambupd: cannot open %s\n", fname);
            curfile[0] = '\0';
            return;
        }
        strncpy(curfile, fname, sizeof(curfile) - 1);
        curfile[sizeof(curfile) - 1] = '\0';

        trace(1, "ambupd: writing %s\n", fname);
        fprintf(fp, "%%MJD     SOD  STN SAT         B_IF[m]            N_WL[cyc]  SIG_WL\n");
    }

    /* time tags */
    int week;
    double tow = time2gpst(rtk->sol.time, &week);
    double mjd = 44244.0 + (week * 7.0) + (tow / 86400.0);
    double sod = ep[3] * 3600.0 + ep[4] * 60.0 + ep[5];

    const int NFx = NF(&rtk->opt); /* 1 if ionoopt == IONOOPT_IFLC */

    int sat;
    for (sat = 1; sat <= MAXSAT; sat++) {

        const ssat_t *ss = &rtk->ssat[sat - 1];

        /* only satellites used in current solution */
        if (!ss->vsat[0]) continue;

        const obsd_t *pobs = find_obs(obs, n, sat);
        if (!pobs) continue;

        /* only GPS here */
        if (satsys(sat, NULL) != SYS_GPS) continue;

        int kL1 = 0, kL5 = 2;
        double f1 = 0.0, f5 = 0.0;

        if (!get_gps_l1_l5_idx(sat, pobs, nav, &kL1, &kL5, &f1, &f5)) continue;
        if (f1 <= 0.0 || f5 <= 0.0) continue;

        /* WL from MW/SMW (meters -> cycles) */
        double smw_m = ss->mw[1];
        double mw_var_m2 = ss->mw[3];

        /* fallback for forks that store MW in other slots */
        if (!finite_d(smw_m) || fabs(smw_m) < 1e-12) smw_m = ss->mw[0];
        if (!finite_d(mw_var_m2) || mw_var_m2 < 0.0) mw_var_m2 = ss->mw[2];

        double N_wl = 0.0, sigma_wl = 0.0;
        const double denom_wl_hz = (f1 - f5);
        if (finite_d(smw_m) && fabs(denom_wl_hz) > 1e-3) {
            const double lambda_wl = CLIGHT / denom_wl_hz; /* m/cycle */
            if (finite_d(lambda_wl) && fabs(lambda_wl) > 0.0) {
                N_wl = smw_m / lambda_wl;
                if (mw_var_m2 > 0.0) sigma_wl = sqrt(mw_var_m2) / lambda_wl;
            }
        }

        /* IF ambiguity (meters) */
        double B_if = 0.0;

        if (NFx <= 1) {
            /* IONOOPT_IFLC: only ONE ambiguity state per satellite */
            int idx = IB(sat, 0, &rtk->opt);
            if (idx < 0 || idx >= rtk->nx) continue;
            B_if = rtk->x[idx];
        }
        else {
            /* multi-freq (uncombined): combine L1 and L5 ambiguity states */
            if (kL1 < 0 || kL5 < 0 || kL1 >= NFx || kL5 >= NFx) continue;

            int idx1 = IB(sat, kL1, &rtk->opt);
            int idx5 = IB(sat, kL5, &rtk->opt);
            if (idx1 < 0 || idx5 < 0 || idx1 >= rtk->nx || idx5 >= rtk->nx) continue;

            double b1 = rtk->x[idx1];
            double b5 = rtk->x[idx5];

            const double denom_if = f1 * f1 - f5 * f5;
            if (fabs(denom_if) < 1e-6) continue;

            B_if = (f1 * f1 * b1 - f5 * f5 * b5) / denom_if;
        }

        if (!finite_d(B_if)) continue;
        if (!finite_d(N_wl)) N_wl = 0.0;
        if (!finite_d(sigma_wl) || sigma_wl < 0.0) sigma_wl = 0.0;

        char satid[8];
        satno2id(sat, satid);

        fprintf(fp, "%8.0f%10.1f%5s%4s%19.3f%19.3f%10.3f\n",
                floor(mjd), sod, stname, satid, B_if, N_wl, sigma_wl);
    }

    fflush(fp);
}