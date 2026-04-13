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

#include <ctype.h>
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
#define MIN_ARC_EPOCHS  10

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

/* --- per-sat arc stats for WL(MW) -------------------------------- */
typedef struct { int n; double mean, M2; } arcstat_t;
static arcstat_t wl_arc[MAXSAT];
static int wl_arc_init[MAXSAT]; /* 0 = not init, 1 = init */
static gtime_t wl_last_obs[MAXSAT];
static int wl_last_init[MAXSAT];

static void arc_reset(int sat)
{
    if (sat <= 0 || sat > MAXSAT) return;
    wl_arc[sat-1].n = 0;
    wl_arc[sat-1].mean = 0.0;
    wl_arc[sat-1].M2 = 0.0;
    wl_arc_init[sat-1] = 1;
}

static void arc_update_stat(arcstat_t *a, double x)
{
    double d, d2;

    if (!a) return;

    a->n++;
    d = x - a->mean;
    a->mean += d / (double)a->n;
    d2 = x - a->mean;
    a->M2 += d * d2;
}

static void arc_update(int sat, double x)
{
    if (sat <= 0 || sat > MAXSAT) return;
    arc_update_stat(&wl_arc[sat-1], x);
}

static double arc_std(const arcstat_t *a)
{
    if (!a || a->n < 2) return 999.0;
    return sqrt(a->M2 / (double)(a->n - 1)); /* sample STD */
}

static double arc_sigma_mean(const arcstat_t *a)
{
    double std;

    if (!a) return 0.05;
    if (a->n < 2) return 0.05;

    std = arc_std(a);
    return std <= 0.0 ? 0.05 : std / sqrt((double)a->n);
}

/* ---------------------------------------------------------------- */
extern void ambupd_write_epoch(const rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    int week, year, doy;
    double tow, mjd;
    double ep[6], time_of_day;

    char station_name[32];
    const char *src_name;
    char *p;

    char temp_amb_file[1024];
    char satid[8];

    static FILE *fp_amb = NULL;
    static char amb_file[1024] = "";

    int iobs;
    int sat, sys;
    int k1, k2;
    double f1, f2;
    double gap_reset;

    double L1, L2, P1, P2;
    double lam_w, Pn, Nwl_mw;
    double wl_out, sigma_wl, if_out;

    double N_if;
    int idx_if;

    /* Need IFLC PPP (like your original condition) */
    if (!(rtk->opt.ionoopt == IONOOPT_IFLC && rtk->opt.nf >= 2)) return;
    if (n <= 0) return;

    /* station name */
    src_name = rtk->opt.station_name[0] ? rtk->opt.station_name : "UNKN";
    strncpy(station_name, src_name, sizeof(station_name)-1);
    station_name[sizeof(station_name)-1] = '\0';
    for (p = station_name; *p; p++) *p = (char)toupper((unsigned char)*p);

    /* time tags + filename */
    tow = time2gpst(rtk->sol.time, &week);
    time2epoch(rtk->sol.time, ep);
    year = (int)ep[0];
    doy  = time2doy(rtk->sol.time);

    sprintf(temp_amb_file, "%s_ambupd_%4d%03d", station_name, year, doy);

    if (fp_amb == NULL || strcmp(temp_amb_file, amb_file) != 0) {
        if (fp_amb) { fclose(fp_amb); fp_amb = NULL; }
        strcpy(amb_file, temp_amb_file);
        fp_amb = fopen(amb_file, "w");
        if (fp_amb == NULL) {
            trace(1, "Cannot create ambiguity file: %s\n", amb_file);
            return;
        }
        trace(1, "Creating ambiguity file: %s\n", amb_file);
    }

    /* mjd + time-of-day */
    mjd = 44244.0 + (week * 7.0) + (tow / 86400.0);
    time_of_day = ep[3] * 3600.0 + ep[4] * 60.0 + ep[5];
    gap_reset = 1.5 * MAX(fabs(rtk->tt), 1.0);

    /* loop by observations (not 1..MAXSAT) */
    for (iobs = 0; iobs < n; iobs++) {

        sat = obs[iobs].sat;
        if (sat <= 0 || sat > MAXSAT) continue;
        if (!rtk->ssat[sat-1].vsat[0]) continue;

        sys = satsys(sat, NULL);

        /* Current target: GPS L1/L2 only */
        if (sys != SYS_GPS) {
            continue;
        }

        k1 = 0;
        k2 = 1;
        if (obs[iobs].code[k1] == 0 || obs[iobs].code[k2] == 0) continue;

        {
            double L_corr[NFREQ]={0}, P_corr[NFREQ]={0};
            double dantr0[NFREQ]={0}, dants0[NFREQ]={0};
            double Lc_dummy=0.0, Pc_dummy=0.0;

            corr_meas(obs+iobs, nav, rtk->ssat[sat-1].azel, &rtk->opt,
                      dantr0, dants0, rtk->ssat[sat-1].phw,
                      L_corr, P_corr, &Lc_dummy, &Pc_dummy);

            f1 = sat2freq(sat, obs[iobs].code[k1], nav);
            f2 = sat2freq(sat, obs[iobs].code[k2], nav);
            L1 = L_corr[k1];
            L2 = L_corr[k2];
            P1 = P_corr[k1];
            P2 = P_corr[k2];
        }


        if (f1 <= 0.0 || f2 <= 0.0) continue;
        if (fabs(f1 - FREQ1) > 1e6 || fabs(f2 - FREQ2) > 1e6) {
            trace(4, "ambupd: skip sat=%d, non-L1/L2 pair f1=%.3f f2=%.3f\n",
                  sat, f1, f2);
            continue;
        }
        if (L1 == 0.0 || L2 == 0.0 || P1 == 0.0 || P2 == 0.0) continue;

        idx_if = IB(sat, 0, &rtk->opt);
        if (idx_if < 0 || idx_if >= rtk->nx) continue;
        N_if = rtk->x[idx_if]; /* meters */

        if (fabs(N_if) < 1e-12) continue;

        if (!wl_arc_init[sat-1]) arc_reset(sat);
        if (rtk->ssat[sat-1].slip[k1] || rtk->ssat[sat-1].slip[k2]) {
            arc_reset(sat);
            trace(3, "ambupd: reset WL arc on slip sat=%d slip1=%u slip2=%u\n",
                  sat, rtk->ssat[sat-1].slip[k1], rtk->ssat[sat-1].slip[k2]);
        }
        if (wl_last_init[sat-1] &&
            fabs(timediff(obs[iobs].time, wl_last_obs[sat-1])) > gap_reset) {
            arc_reset(sat);
            trace(3, "ambupd: reset WL arc on gap sat=%d gap=%.1f tt=%.1f\n",
                  sat, fabs(timediff(obs[iobs].time, wl_last_obs[sat-1])), fabs(rtk->tt));
        }

        /* WL from Melbourne-Wubbena (cycles) */
        if (fabs(f1 - f2) < 1.0) continue;
        lam_w = CLIGHT / (f1 - f2);               /* meters */
        Pn    = (f1*P1 + f2*P2) / (f1 + f2);      /* meters */
        Nwl_mw = (((f1 * L1) - (f2 * L2)) / (f1 - f2) - Pn) / lam_w; /* cycles */

        arc_update(sat, Nwl_mw);
        wl_last_obs[sat-1] = obs[iobs].time;
        wl_last_init[sat-1] = 1;
        if_out   = N_if;  /* current KF state is the optimal estimate */
        wl_out   = wl_arc[sat-1].mean;
        sigma_wl = arc_sigma_mean(&wl_arc[sat-1]);

        if (wl_arc[sat-1].n < MIN_ARC_EPOCHS) continue;

        satno2id(sat, satid);

        trace(4,
              "ambupd: sat=%s idx_if=%d IF=%.4f WLraw=%.4f WL=%.4f SIG=%.4f n=%d\n",
              satid, idx_if, if_out, Nwl_mw, wl_out, sigma_wl, wl_arc[sat-1].n);

        fprintf(fp_amb, "%8.0f%10.1f%5s%4s%19.3f%19.3f%10.3f\n",
                floor(mjd), time_of_day, station_name, satid,
                if_out, wl_out, sigma_wl);
    }

    fflush(fp_amb);
}
