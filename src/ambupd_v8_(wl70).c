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

#define MIN_GAP_RESET_S  90.0    /* Minimum gap that triggers an arc reset (seconds) */

/* --- per-sat arc stats for WL(MW) -------------------------------- */
typedef struct { int n; double mean, M2; } arcstat_t;
static arcstat_t wl_arc[MAXSAT];
static int wl_arc_init[MAXSAT]; /* 0 = not init, 1 = init */
static gtime_t wl_last_obs[MAXSAT];
static int wl_last_init[MAXSAT];

 /* AMBFLAG --------------------------------------------------------------------- */
#define MAX_AMB_ARCS  8192

typedef enum { AMBFLAG_AMB=0, AMBFLAG_DEL=1, AMBFLAG_BAD=2 } ambflag_type_t;

typedef struct {
    int           sat;        /* satellite number (satid2no) */
    ambflag_type_t type;      /* AMB / DEL / BAD */
    int           ep_beg;     /* begin epoch index (1-based, from ambflag) */
    int           ep_end;     /* end epoch index */
} ambflag_rec_t;

typedef struct {
    ambflag_rec_t *recs;
    int            n;
    int            nmax;
    int            interval;  /* sampling interval in seconds */
    gtime_t        t0;        /* begin time of the session */
} ambflag_t;

static ambflag_t g_ambflag = {0};  /* global, loaded once per station */
static int time2epno(gtime_t t, const ambflag_t *af);

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
    d       = x - a->mean;
    a->mean += d / (double)a->n;
    d2      = x - a->mean;
    a->M2  += d * d2;
}

static double arc_std(const arcstat_t *a)
{
    if (!a || a->n < 2) return 999.0;
    return sqrt(a->M2 / (double)(a->n - 1));
}
 
/* Standard error of the arc mean: std / sqrt(n) */
static double arc_sigma_mean(const arcstat_t *a)
{
    double std;
    if (!a || a->n < 2) return 999.0;
    std = arc_std(a);
    return std <= 0.0 ? 999.0 : std / sqrt((double)a->n);
}

static int arc_mw_outlier(int sat, double x)
{
    const arcstat_t *a;
    double std, thr;
    if (sat <= 0 || sat > MAXSAT) return 0;
    a = &wl_arc[sat-1];
    if (a->n < 10) return 0; /* no robust stats yet */
    std = arc_std(a);
    if (std <= 0.0) return 0;
    thr = MAX(0.25, 4.0 * std);
    return fabs(x - a->mean) > thr;
}

extern void ambupd_reset(void)
{
    int i;
    for (i = 0; i < MAXSAT; i++) {
        wl_arc[i].n      = 0;
        wl_arc[i].mean   = 0.0;
        wl_arc[i].M2     = 0.0;
        wl_arc_init[i]   = 0;
        wl_last_init[i]  = 0;
    }
    trace(3, "ambupd_reset: all arc statistics cleared\n");
}

/* Convert gtime_t to epoch index (1-based, matches ambflag numbering) */
static int time2epno(gtime_t t, const ambflag_t *af)
{
    if (af->interval <= 0) return -1;
    double dt = timediff(t, af->t0);
    if (dt < -0.5) return -1;
    return (int)(dt / af->interval + 0.5) + 1;
}

/* ambflag_is_valid() — mirrors GREAT-UPD t_gambflag::isValid()
+ *
+ * Returns 0 (skip epoch) or 1 (use epoch), and sets *arc_reset=1
+ * if this is the first epoch of a new AMB arc (KF state must reset).
+ *
+ * Priority: DEL/BAD checked first across ALL records, then AMB.
+ * This matches GREAT-UPD exactly: a DEL overlapping an AMB wins.   */
extern int ambflag_is_valid(int sat, gtime_t t,
                            const ambflag_t *af, int *arc_reset)
{
    if (arc_reset) *arc_reset = 0;
    if (!af || !af->recs || af->n == 0) return 1; /* no ambflag → allow all */

    int ep = time2epno(t, af);
    if (ep < 0) return 0;

    /* --- Step 1: DEL/BAD have priority (scan all records) -------------- */
    /* Mirror of GREAT isValid() first loop */
    int i;
    for (i = 0; i < af->n; i++) {
        ambflag_rec_t *r = &af->recs[i];
        if (r->sat != sat) continue;
        if (ep < r->ep_beg || ep > r->ep_end) continue;
        if (r->type == AMBFLAG_DEL || r->type == AMBFLAG_BAD)
            return 0; /* excluded */
    }

    /* --- Step 2: find AMB arc containing this epoch -------------------- */
    /* Mirror of GREAT isValid() second loop */
    for (i = 0; i < af->n; i++) {
        ambflag_rec_t *r = &af->recs[i];
        if (r->sat != sat) continue;
        if (r->type != AMBFLAG_AMB) continue;
        if (ep < r->ep_beg || ep > r->ep_end) continue;
        /* epoch is inside a valid AMB arc */
        if (arc_reset && ep == r->ep_beg)
            *arc_reset = 1; /* first epoch of arc → reset KF */
        return 1;
    }

    /* epoch not in any AMB arc → exclude */
    return 0;
}

/* Parse one ambflag file. Call once per station before pppos loop. */
extern int ambflag_load(const char *file, ambflag_t *af)
{
    FILE *fp;
    char buff[256], flag[8], prn[8], reason[64];
    int ep_beg, ep_end;
    int in_header = 1;
    double interval = 30.0;
    gtime_t t0 = {0};

    if (!file || !*file) return 0;
    if (!(fp = fopen(file, "r"))) {
        trace(1, "ambflag_load: cannot open %s\n", file);
        return 0;
    }

    /* free previous */
    if (af->recs) { free(af->recs); af->recs = NULL; }
    af->n = af->nmax = 0;

    while (fgets(buff, sizeof(buff), fp)) {
        /* parse header */
        if (in_header) {
            if (strstr(buff, "INTERVAL"))
                interval = atof(buff);
            if (strstr(buff, "BEGIN TIME")) {
                /* format: YYYY MM DD HH MM SS.ss  GPST  BEGIN TIME */
                double ep[6];
                if (sscanf(buff, "%lf %lf %lf %lf %lf %lf",
                           ep, ep+1, ep+2, ep+3, ep+4, ep+5) == 6)
                    t0 = epoch2time(ep);
            }
            if (strstr(buff, "END OF HEADER")) { in_header = 0; continue; }
            continue;
        }

        /* data record: FLAG  PRN  EP_BEG  EP_END  REASON */
        if (sscanf(buff, "%7s %7s %d %d %63s",
                   flag, prn, &ep_beg, &ep_end, reason) < 4) continue;

        ambflag_type_t type;
        if      (!strcmp(flag, "AMB")) type = AMBFLAG_AMB;
        else if (!strcmp(flag, "DEL")) type = AMBFLAG_DEL;
        else if (!strcmp(flag, "BAD")) type = AMBFLAG_BAD;
        else continue;

        int sat = satid2no(prn);
        if (sat <= 0 || sat > MAXSAT) continue;

        /* grow buffer */
        if (af->n >= af->nmax) {
            af->nmax += 1024;
            ambflag_rec_t *tmp = realloc(af->recs,
                                         sizeof(ambflag_rec_t) * af->nmax);
            if (!tmp) { fclose(fp); return 0; }
            af->recs = tmp;
        }
        af->recs[af->n].sat    = sat;
        af->recs[af->n].type   = type;
        af->recs[af->n].ep_beg = ep_beg;
        af->recs[af->n].ep_end = ep_end;
        af->n++;
    }
    fclose(fp);

    af->interval = ROUND(interval);
    af->t0       = t0;
    trace(2, "ambflag_load: %s → %d records, interval=%ds\n",
          file, af->n, af->interval);
    return af->n;
}

extern void ambflag_free(ambflag_t *af)
{
    if (af->recs) { free(af->recs); af->recs = NULL; }
    af->n = af->nmax = 0;
}

static int get_day_key(gtime_t t)
{
    double ep[6];
    int doy;
    time2epoch(t, ep);
    doy = time2doy(t);
    return ((int)ep[0]) * 1000 + doy;
}

static void get_ambflag_path(const prcopt_t *opt, char *path, size_t npath)
{
    const char *p;
    if (!path || npath == 0) return;
    (void)npath;
    path[0] = '\0';
    if (!opt || !opt->pppopt[0]) return;
    if ((p = strstr(opt->pppopt, "-AMBFLAG="))) {
        sscanf(p + 9, "%1023s", path);
    }
}

extern void ambflag_load_from_opt(const prcopt_t *opt, gtime_t t)
{
    static char loaded_for_station[64] = "";
    static char loaded_afpath[1024] = "";
    static int loaded_day_key = -1;
    char afpath[1024] = "";
    const char *sta;
    int day_key;
    int changed_ctx;

    if (!opt) return;
    sta = opt->station_name[0] ? opt->station_name : "UNKN";
    get_ambflag_path(opt, afpath, sizeof(afpath));
    day_key = get_day_key(t);

    changed_ctx = strcmp(sta, loaded_for_station) ||
                  strcmp(afpath, loaded_afpath) ||
                  day_key != loaded_day_key;
    if (!changed_ctx) return;

    ambflag_free(&g_ambflag);
    ambupd_reset();

    strncpy(loaded_for_station, sta, sizeof(loaded_for_station)-1);
    loaded_for_station[sizeof(loaded_for_station)-1] = '\0';
    strncpy(loaded_afpath, afpath, sizeof(loaded_afpath)-1);
    loaded_afpath[sizeof(loaded_afpath)-1] = '\0';
    loaded_day_key = day_key;

    if (afpath[0]) {
        ambflag_load(afpath, &g_ambflag);
        trace(2, "ambupd: loaded ambflag for %s day=%d path=%s n=%d\n",
              sta, day_key, afpath, g_ambflag.n);
    } else {
        trace(3, "ambupd: no -AMBFLAG= for %s day=%d, gap-only mode\n",
              sta, day_key);
    }
}

static int code_rank_l1(uint8_t code)
{
    switch (code) {
    case CODE_L1W: return 0;
    case CODE_L1P: return 1;
    case CODE_L1Y: return 2;
    case CODE_L1C: return 3;
    case CODE_L1S:
    case CODE_L1L:
    case CODE_L1X: return 4;
    default: return 99;
    }
}

static int code_rank_l2(uint8_t code)
{
    switch (code) {
    case CODE_L2W: return 0;
    case CODE_L2P: return 1;
    case CODE_L2Y: return 2;
    case CODE_L2C: return 3;
    case CODE_L2X: return 4;
    case CODE_L2D: return 5;
    case CODE_L2S:
    case CODE_L2L: return 90; /* avoid unless no better L2 available */
    default: return 99;
    }
}

static int select_l1_l2_indices(const obsd_t *obs, const nav_t *nav, int sat,
                                int *k1, int *k2, double *f1, double *f2)
{
    int i, best1 = -1, best2 = -1, rank1 = 999, rank2 = 999;
    double fi;

    if (!obs || !nav || !k1 || !k2 || !f1 || !f2) return 0;

    /* Keep consistency with PPP IF state if canonical slots are valid. */
    if (obs->code[0] && obs->code[1] && obs->L[0] != 0.0 && obs->L[1] != 0.0 &&
        obs->P[0] != 0.0 && obs->P[1] != 0.0) {
        *f1 = sat2freq(sat, obs->code[0], nav);
        *f2 = sat2freq(sat, obs->code[1], nav);
        if (*f1 > 0.0 && *f2 > 0.0 &&
            fabs(*f1 - FREQ1) <= 1e6 && fabs(*f2 - FREQ2) <= 1e6) {
            *k1 = 0;
            *k2 = 1;
            return 1;
        }
    }

    for (i = 0; i < NFREQ + NEXOBS; i++) {
        int r;
        if (!obs->code[i] || obs->L[i] == 0.0 || obs->P[i] == 0.0) continue;
        fi = sat2freq(sat, obs->code[i], nav);
        if (fi <= 0.0) continue;

        if (fabs(fi - FREQ1) <= 1e6) {
            r = code_rank_l1(obs->code[i]);
            if (r < rank1) {
                rank1 = r;
                best1 = i;
                *f1 = fi;
            }
        }
        else if (fabs(fi - FREQ2) <= 1e6) {
            r = code_rank_l2(obs->code[i]);
            if (r < rank2) {
                rank2 = r;
                best2 = i;
                *f2 = fi;
            }
        }
    }

    if (best1 < 0 || best2 < 0) return 0;
    *k1 = best1;
    *k2 = best2;
    return 1;
}

static int get_dsb(const nav_t *nav, int sat, int code1, int code2, double *bias)
{
    double v;
    if (!nav || !bias || sat <= 0 || sat > MAXSAT) return 0;
    v = nav->cbias[sat-1][code1][code2];
    if (fabs(v) > 1E-12) {
        *bias = v;
        return 1;
    }
    v = nav->cbias[sat-1][code2][code1];
    if (fabs(v) > 1E-12) {
        *bias = -v;
        return 1;
    }
    *bias = 0.0;
    return 0;
}

static int map_code_to_corr(uint8_t code, const double cor[5], double *corr)
{
    if (!corr) return 0;
    switch (code) {
    case CODE_L1C: *corr = cor[2]; return 1;
    case CODE_L2C: *corr = cor[3]; return 1;
    case CODE_L1P:
    case CODE_L1Y:
    case CODE_L1W: *corr = cor[0]; return 1;
    case CODE_L2P:
    case CODE_L2Y:
    case CODE_L2W: *corr = cor[1]; return 1;
    case CODE_L2X: *corr = cor[4]; return 1;
    case CODE_L2D: *corr = 0.0;    return 1;
    default:       *corr = 0.0;    return 0;
    }
}

static int mw_dcb_corr_cycles(const obsd_t *obs, const nav_t *nav, int sat,
                              int k1, int k2, double lam1, double lam2,
                              double *corr_cyc)
{
    double p1p2 = 0.0, p1c1 = 0.0, p2c2 = 0.0;
    double cor[5], corr1 = 0.0, corr2 = 0.0;
    double fact, c1, c2;
    int ok1, ok2;

    if (!obs || !nav || !corr_cyc || lam1 <= 0.0 || lam2 <= 0.0) return 0;
    *corr_cyc = 0.0;

    (void)get_dsb(nav, sat, CODE_L1W, CODE_L2W, &p1p2); /* P1-P2 */
    (void)get_dsb(nav, sat, CODE_L1W, CODE_L1C, &p1c1); /* P1-C1 */
    (void)get_dsb(nav, sat, CODE_L2W, CODE_L2C, &p2c2); /* P2-C2 */

    fact = lam1 / lam2;
    if (fabs(1.0 - fact*fact) < 1E-12) return 0;
    c1 = 1.0 / (1.0 - fact*fact);
    c2 = fact*fact * c1;

    cor[0] = c2 * p1p2;
    cor[1] = c1 * p1p2;
    cor[2] = cor[0] + p1c1;
    cor[3] = cor[1] + p2c2;
    cor[4] = cor[1] + p1c1;

    ok1 = map_code_to_corr(obs->code[k1], cor, &corr1);
    ok2 = map_code_to_corr(obs->code[k2], cor, &corr2);
    if (!ok1 || !ok2) return 0;

    *corr_cyc = -(corr1 / lam1 + corr2 / lam2) * (1.0 - fact) / (1.0 + fact);
    return 1;
}

/* ---------------------------------------------------------------------------
 * mw_raw() — compute Melbourne-Wübbena combination from RAW observations.
 *
 * Implements the same structure as GREAT-UPD:
 *   1) MW from raw code/phase
 *   2) additional DCB correction term based on P1P2/P1C1/P2C2 model.
 *
 * Returns MW in cycles, or 0.0 on failure (sets *ok=0).
 * --------------------------------------------------------------------------*/
static double mw_raw(const obsd_t *obs, const nav_t *nav,
                     int sat, int k1, int k2,
                     double f1, double f2, int *ok)
{
    double lam1, lam2, fact;
    double C1, C2, L1, L2;
    double mw, dcb_corr = 0.0;

    *ok = 0;

    if (f1 <= 0.0 || f2 <= 0.0 || fabs(f1 - f2) < 1.0) return 0.0;
    lam1 = CLIGHT / f1;
    lam2 = CLIGHT / f2;
    if (lam1 <= 0.0 || lam2 <= 0.0) return 0.0;

    C1 = obs->P[k1];
    C2 = obs->P[k2];
    L1 = obs->L[k1];
    L2 = obs->L[k2];
    if (C1 == 0.0 || C2 == 0.0 || L1 == 0.0 || L2 == 0.0) return 0.0;

    fact = lam1 / lam2;
    mw = L1 - L2 - (C1 / lam1 + C2 / lam2) * (1.0 - fact) / (1.0 + fact);

    /* GREAT-like DCB term (if required pairs are available). */
    (void)mw_dcb_corr_cycles(obs, nav, sat, k1, k2, lam1, lam2, &dcb_corr);
    mw += dcb_corr;

    *ok = 1;
    return mw;
}

/* ---------------------------------------------------------------- */
extern void ambupd_write_epoch(const rtk_t *rtk, const obsd_t *obs, int n,
                               const nav_t *nav)
{
    int    week, year, doy;
    double tow, mjd;
    double ep[6], time_of_day;
 
    char        station_name[32];
    const char *src_name;
    char       *p;
 
    char temp_amb_file[1024];
    char satid[8];
 
    static FILE *fp_amb   = NULL;
    static char  amb_file[1024] = "";
 
    int    iobs;
    int    sat, sys;
    int    k1, k2;
    double f1, f2;
    double gap_reset;
    double Nwl_mw;
    double wl_out, sigma_wl, if_out;
    double N_if;
    int    idx_if;
    int    mw_ok;
 
    /* Requires IFLC dual-frequency PPP */
    if (!(rtk->opt.ionoopt == IONOOPT_IFLC && rtk->opt.nf >= 2)) return;
    if (n <= 0) return;
 
    /* ---- station name (upper-case, 4 chars for GREAT-UPD format) ---------- */
    src_name = rtk->opt.station_name[0] ? rtk->opt.station_name : "UNKN";
    strncpy(station_name, src_name, sizeof(station_name)-1);
    station_name[sizeof(station_name)-1] = '\0';
    for (p = station_name; *p; p++) *p = (char)toupper((unsigned char)*p);
 
    /* ---- output filename -------------------------------------------------- */
    tow = time2gpst(rtk->sol.time, &week);
    time2epoch(rtk->sol.time, ep);
    year = (int)ep[0];
    doy  = time2doy(rtk->sol.time);
    sprintf(temp_amb_file, "%s_ambupd_%4d%03d", station_name, year, doy);
 
    if (fp_amb == NULL || strcmp(temp_amb_file, amb_file) != 0) {
        if (fp_amb) {
            fclose(fp_amb);
            fp_amb = NULL;
            /* New day/session file: always start with fresh WL arc stats. */
            ambupd_reset();
        }
        strcpy(amb_file, temp_amb_file);
        fp_amb = fopen(amb_file, "w");
        if (!fp_amb) {
            trace(1, "ambupd: cannot create file: %s\n", amb_file);
            return;
        }
        trace(1, "ambupd: creating file: %s\n", amb_file);
    }
 
    /* ---- time tags -------------------------------------------------------- */
    mjd         = 44244.0 + (week * 7.0) + (tow / 86400.0);
    time_of_day = ep[3]*3600.0 + ep[4]*60.0 + ep[5];
 
    /* ---- gap threshold: 3× nominal interval, but at least MIN_GAP_RESET_S --
     * rtk->tt is the time difference to the previous epoch in seconds.
     * At startup it may be 0 or very small, so we clamp from below.        */
    gap_reset = MAX(3.0 * fabs(rtk->tt), MIN_GAP_RESET_S);
 
    /* ---- per-satellite loop ----------------------------------------------- */
    for (iobs = 0; iobs < n; iobs++) {
 
        sat = obs[iobs].sat;
        if (sat <= 0 || sat > MAXSAT) continue;
 
        sys = satsys(sat, NULL);
        if (sys != SYS_GPS) continue;   /* GPS L1/L2 only for now */
 
        if (!select_l1_l2_indices(obs + iobs, nav, sat, &k1, &k2, &f1, &f2)) {
            trace(4, "ambupd: skip sat=%d no valid L1/L2 pair in obs record\n", sat);
            continue;
        }
        if (fabs(f1 - f2) < 1.0) continue; /* degenerate */
 
        /* ---- IF ambiguity from KF state (metres) --------------------------
         * IB() indexes the IFLC phase-bias state, which in RTKLIB ppp.c
         * absorbs: lambda_IF * N_IF + receiver_bias - satellite_bias.
         * GREAT-UPD handles the bias datum itself; we just pass the raw
         * KF float value.  Reject uninitialized (zero) states.            */
        idx_if = IB(sat, 0, &rtk->opt);
        if (idx_if < 0 || idx_if >= rtk->nx) continue;
        N_if = rtk->x[idx_if];
        if (fabs(N_if) < 1e-12) continue;
        
        /* ---- ambflag gate (mirrors GREAT-UPD isValid()) ---------------
         * DEL/BAD: skip epoch entirely.
         * AMB at arc start: reset WL accumulator so MW arc boundaries
         *   match what GREAT-UPD NL will see in ambflag.              */
        {
            int arc_reset_flag = 0;
            if (!ambflag_is_valid(sat, obs[iobs].time,
                                  &g_ambflag, &arc_reset_flag)) {
                trace(4, "ambupd: sat=%d excluded by ambflag DEL/BAD\n", sat);
                continue;
            }
            if (arc_reset_flag) {
                arc_reset(sat);
                trace(3, "ambupd: WL arc reset by ambflag AMB sat=%d\n", sat);
            }
        }
        
        /* ---- fallback: gap-based reset if no ambflag loaded ----------- */
        if (!wl_arc_init[sat-1]) arc_reset(sat);
        if (wl_last_init[sat-1] &&
            fabs(timediff(obs[iobs].time, wl_last_obs[sat-1])) > gap_reset) {
            if (g_ambflag.n == 0) { /* only if ambflag not loaded */
                arc_reset(sat);
            }
        }
 
        /* ---- Melbourne-Wübbena from RAW observations ----------------------
         * Key change v2: we do NOT call corr_meas() here.
         * Antenna PCO/PCV corrections are frequency-dependent; applying them
         * before MW would introduce an elevation-varying bias into the
         * geometry-free combination and corrupt the integer WL.            */
        Nwl_mw = mw_raw(obs+iobs, nav, sat, k1, k2, f1, f2, &mw_ok);
        if (!mw_ok) continue;

        if (arc_mw_outlier(sat, Nwl_mw)) {
            trace(4, "ambupd: sat=%d MW outlier rejected raw=%.3f mean=%.3f n=%d\n",
                  sat, Nwl_mw, wl_arc[sat-1].mean, wl_arc[sat-1].n);
            continue;
        }
 
        arc_update(sat, Nwl_mw);
        wl_last_obs[sat-1]  = obs[iobs].time;
        wl_last_init[sat-1] = 1;
 
        if_out   = N_if;
        wl_out   = wl_arc[sat-1].mean;
        sigma_wl = arc_sigma_mean(&wl_arc[sat-1]);
 
        satno2id(sat, satid);
        trace(3, "ambupd: sat=%-4s idx_if=%3d IF=%10.4f m  "
              "WL_raw=%7.4f cy  WL_mean=%7.4f cy  sig=%6.4f cy  n=%d\n",
              satid, idx_if, if_out, Nwl_mw, wl_out, sigma_wl,
              wl_arc[sat-1].n);
 
        /* GREAT-UPD ambupd format (A.3):
         *   I8   F10.1   A5   A4   F19.3   F19.3   F10.3
         *   MJD  tod(s)  site prn  IF(m)   WL(cy)  sig(cy) */
        fprintf(fp_amb, "%8.0f%10.1f%5s%4s%19.3f%19.3f%10.3f\n",
                floor(mjd), time_of_day, station_name, satid,
                if_out, wl_out, sigma_wl);
    }
 
    fflush(fp_amb);
}
