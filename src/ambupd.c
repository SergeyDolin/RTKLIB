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


/* --- new: select band-pair (example) ---------------------------- */
/* you can put it into rtk->opt as an int, here keep as macro/const */
#define PAIR_L1L2 0
#define PAIR_L1L5 1

static int g_pair = PAIR_L1L5; /* e.g. set from option later */

/* map desired band -> obs index (0..NFREQ-1 in RTKLIB obs array) */
static int obs_index_for_band_gps(const obsd_t *o, int band) {
    /* band: 1,2,5 -> want something like C1*, L1*, etc. in o->code[k]
       In RTKLIB, code2idx()/satsys-dependent helpers may exist in your fork.
       If not, simplest: test sat2freq() > 0 and compare to nominal FREQ*. */
    int k;
    for (k = 0; k < NFREQ; k++) {
        if (o->code[k] == 0) continue;
        /* We decide by frequency value (robust vs RINEX code variants) */
        /* sat2freq needs nav + sat; so do selection in caller where nav/sat known */
    }
    return -1;
}

/* --- per-sat arc stats for WL(MW) (C90) ------------------------- */
typedef struct { int n; double mean, M2; } arcstat_t;
static arcstat_t wl_arc[MAXSAT];
static int wl_arc_init[MAXSAT]; /* 0 = not init, 1 = init */

static void arc_reset(int sat)
{
    if (sat <= 0 || sat > MAXSAT) return;
    wl_arc[sat-1].n = 0;
    wl_arc[sat-1].mean = 0.0;
    wl_arc[sat-1].M2 = 0.0;
    wl_arc_init[sat-1] = 1;
}

static void arc_update(int sat, double x)
{
    arcstat_t *a;
    double d, d2;

    if (sat <= 0 || sat > MAXSAT) return;
    a = &wl_arc[sat-1];

    a->n++;
    d = x - a->mean;
    a->mean += d / (double)a->n;
    d2 = x - a->mean;
    a->M2 += d * d2;
}

static double arc_std(const arcstat_t *a)
{
    if (!a || a->n < 2) return 999.0;
    return sqrt(a->M2 / (double)(a->n - 1)); /* sample STD */
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

    int iobs, k;
    int sat, sys;
    int k1, k2;
    double f1, f2;

    double L1, L2, P1, P2;
    double lam_w, Pn, Nwl_mw;
    double wl_out, sigma_wl;

    double N_if;

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

    /* loop by observations (not 1..MAXSAT) */
    for (iobs = 0; iobs < n; iobs++) {

        sat = obs[iobs].sat;
        if (sat <= 0 || sat > MAXSAT) continue;
        if (!rtk->ssat[sat-1].vs) continue;

        sys = satsys(sat, NULL);

        /* --- choose pair: GPS L1-L5 example ---------------------- */
        /* If you want L1-L2 instead, replace FREQ5 with FREQ2 below */
        if (sys != SYS_GPS) {
            /* for now: skip other systems; you can extend later */
            continue;
        }

        k1 = -1; k2 = -1;
        f1 = 0.0; f2 = 0.0;

        /* find obs indices that match L1 and L5 by frequency */
        for (k = 0; k < NFREQ; k++) {
            double fk;
            if (obs[iobs].code[k] == 0) continue;
            fk = sat2freq(sat, obs[iobs].code[k], nav);
            if (fk <= 0.0) continue;

            if (fabs(fk - FREQ1) < 1e6) { k1 = k; f1 = fk; }
            if (fabs(fk - FREQ5) < 1e6) { k2 = k; f2 = fk; } /* L5 */
        }

        /* фильтрация: нет L5 -> пропускаем */
        if (k1 < 0 || k2 < 0) {
            trace(0,"No L5 in this observation\n");
            continue;
        }
        /* MW requires phase+code on both bands */
        L1 = obs[iobs].L[k1];
        L2 = obs[iobs].L[k2];
        P1 = obs[iobs].P[k1];
        P2 = obs[iobs].P[k2];
        if (L1 == 0.0 || L2 == 0.0 || P1 == 0.0 || P2 == 0.0) continue;

        /* reset arc on slip (you need to confirm slip field exists in your ssat_t)
           Most RTKLIB have rtk->ssat[sat-1].slip[k] */
        if (!wl_arc_init[sat-1]) arc_reset(sat);
        /* if (rtk->ssat[sat-1].slip[k1] || rtk->ssat[sat-1].slip[k2]) arc_reset(sat); */

        /* WL from Melbourne-Wubbena (cycles) */
        if (fabs(f1 - f2) < 1.0) continue;
        lam_w = CLIGHT / (f1 - f2);               /* meters */
        Pn    = (f1*P1 + f2*P2) / (f1 + f2);      /* meters */
        Nwl_mw = (L1 - L2) - (Pn / lam_w);        /* cycles (WL) */

        arc_update(sat, Nwl_mw);
        wl_out   = wl_arc[sat-1].mean;
        sigma_wl = arc_std(&wl_arc[sat-1]);

        /* IF ambiguity from PPP float (meters)
           In IFLC PPP there is typically one ambiguity/bias state per sat.
           Your original code used IB(sat,0) and then recombined with freq^2;
           here we keep it simple: take the IF state directly. */
        {
            int idx_if;
            idx_if = IB(sat, 0, &rtk->opt);
            if (idx_if < 0 || idx_if >= rtk->nx) continue;
            N_if = rtk->x[idx_if]; /* meters */
        }

        satno2id(sat, satid);

        fprintf(fp_amb, "%8.0f%10.1f%5s%4s%19.3f%19.3f%10.3f\n",
                floor(mjd), time_of_day, station_name, satid,
                N_if, wl_out, sigma_wl);
    }

    fflush(fp_amb);
}