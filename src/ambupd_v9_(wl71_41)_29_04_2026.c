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

/* --- per-sat WEIGHTED arc stats for WL(MW), Ge et al. 2008 Eq. 31 -- */
typedef struct { double n, w_sum, mean, M2; } arcstat_t;
static arcstat_t wl_arc[MAXSAT];
static int wl_arc_init[MAXSAT]; /* 0 = not init, 1 = init */
static gtime_t wl_last_obs[MAXSAT];
static int wl_last_init[MAXSAT];

/* Robust outlier rejection — short ring buffer of recent MW samples */
#define MW_RBUF_N   32
typedef struct {
    double buf[MW_RBUF_N];
    int    n;          /* how many filled (capped at MW_RBUF_N) */
    int    head;       /* next write index */
} mw_rbuf_t;
static mw_rbuf_t wl_rbuf[MAXSAT];

/* Quality gates. MAX_WL_STD_DEF is compared against arc_sigma_mean (SE of
 * weighted mean), the same field upd_NL gates at 0.130 cy. Default 0.20 cy
 * gives a small headroom over the typical reference (max ~0.114 cy). */
#define MIN_ARC_SEC_DEF   60.0     /* 1 min — minimum arc length */
#define MAX_WL_STD_DEF    0.20     /* cycles — SE of mean gate (upd_NL: 0.130) */
#define MAX_IF_SD_DEF     50.0     /* m — IF state convergence ceiling */
#define ROBUST_K          4.0      /* MAD multiplier */
#define ELEV_FULL_W       (30.0 * D2R) /* full weight above this elevation */

static double g_min_arc_sec = MIN_ARC_SEC_DEF;
static double g_max_wl_std  = MAX_WL_STD_DEF;
static double g_max_if_sd   = MAX_IF_SD_DEF;
static double g_sample_dt   = 30.0;   /* updated from rtk->tt each epoch */

/* Diagnostic counters (cleared by ambupd_reset) */
static unsigned long g_drop_if_sd  = 0UL;
static unsigned long g_drop_arc_sh = 0UL;
static unsigned long g_drop_arc_sd = 0UL;
static unsigned long g_drop_outlr  = 0UL;
static unsigned long g_arc_kept    = 0UL;

typedef struct {
    unsigned long total;
    unsigned long full;
    unsigned long fallback;
    unsigned long map_fail;
} dcb_diag_t;
static dcb_diag_t g_dcb_diag = {0};

typedef enum {
    DCB_MODE_NONE = 0,
    DCB_MODE_FULL,
    DCB_MODE_FALLBACK,
    DCB_MODE_MAP_FAIL
} dcb_mode_t;

typedef struct {
    double mjd;
    double sod;
    char station[32];
    char satid[8];
    double if_out;
    double mw_raw;
    double mw_dcb;
    double mw_dcb_sat;
    double mw_dcb_rcv;
    dcb_mode_t dcb_mode;
    int bsx_rcv_rows;
    char code1[8];
    char code2[8];
    int k1, k2;
    double f1, f2;
} amb_epoch_row_t;

typedef struct {
    int code1;
    int code2;
    double bias_m;
} rcv_dsb_t;

typedef struct {
    rcv_dsb_t *rows;
    int n;
    int nmax;
    char station[32];
    char path[1024];
    int day_key;
} bsx_rcv_cache_t;

static amb_epoch_row_t *g_arc_rows[MAXSAT];
static int g_arc_rows_n[MAXSAT];
static int g_arc_rows_nmax[MAXSAT];

/* Global pending buffer: arc-flushed rows wait here until file close, then are
 * sorted by (mjd, sod, satid) and written. Keeps file in epoch order without
 * changing the per-arc final wl_mean / sigma_wl already computed at flush. */
typedef struct {
    amb_epoch_row_t row;
    double wl_mean;
    double sigma_wl;
} amb_pending_row_t;

static amb_pending_row_t *g_pending = NULL;
static int g_pending_n = 0;
static int g_pending_nmax = 0;

static FILE *g_fp_amb = NULL;
static char g_amb_file[1024] = "";
static FILE *g_fp_diag = NULL;
static char g_diag_file[1024] = "";
static int g_ambupd_atexit_registered = 0;

static bsx_rcv_cache_t g_bsx_rcv = {0};
static char g_bsx_path[1024] = "";
static int g_bsx_diag_notice = 0;

static void ambupd_flush_all_arcs(const char *reason);
static void ambupd_close_outputs(void);
static void bsx_cache_reset(bsx_rcv_cache_t *cache);
static void pending_flush_to_file(void);

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
    wl_arc[sat-1].n     = 0.0;
    wl_arc[sat-1].w_sum = 0.0;
    wl_arc[sat-1].mean  = 0.0;
    wl_arc[sat-1].M2    = 0.0;
    wl_rbuf[sat-1].n    = 0;
    wl_rbuf[sat-1].head = 0;
    wl_arc_init[sat-1]  = 1;
}

/* Ge et al. 2008, Eq. 31 — elevation-dependent weight */
static double elev_weight(double el_rad)
{
    double s;
    if (el_rad >= ELEV_FULL_W) return 1.0;
    if (el_rad <= 0.0)         return 0.0;
    s = 2.0 * sin(el_rad);
    return s > 0.0 ? s : 0.0;
}

/* Weighted online mean & variance (West 1979) */
static void arc_update(int sat, double x, double w)
{
    arcstat_t *a;
    double d;

    if (sat <= 0 || sat > MAXSAT || w <= 0.0) return;
    a = &wl_arc[sat-1];
    a->w_sum += w;
    a->n     += 1.0;
    d         = x - a->mean;
    a->mean  += (w / a->w_sum) * d;
    a->M2    += w * d * (x - a->mean);

    /* Maintain ring buffer for robust MAD test */
    {
        mw_rbuf_t *r = &wl_rbuf[sat-1];
        r->buf[r->head] = x;
        r->head = (r->head + 1) % MW_RBUF_N;
        if (r->n < MW_RBUF_N) r->n++;
    }
}

/* GREAT-UPD getMeanWgt() formulas (gambcommon.cpp:148-149):
 *   sigma    = sqrt( Sum_i w_i*(x_i-xbar)^2 / (n-1) )   // unbiased by count
 *   mean_sig = sigma / sqrt(w_sum)                       // SE of weighted mean
 * GREAT writes mean_sig (NOT sample STD) to the *_ambupd_* file as the σ field;
 * upd_NL then thresholds it against 0.130 cy. Empirically matches reference
 * sample data (AIRA σ range 0.009..0.114). */
static double arc_std(const arcstat_t *a)
{
    double var;
    if (!a || a->n < 2.0 || a->w_sum <= 0.0) return 999.0;
    var = a->M2 / (a->n - 1.0);          /* GREAT divides by count-1 */
    return var > 0.0 ? sqrt(var) : 999.0;
}

/* Standard error of the weighted arc mean — what GREAT writes to ambupd. */
static double arc_sigma_mean(const arcstat_t *a)
{
    double s;
    if (!a || a->n < 2.0 || a->w_sum <= 0.0) return 999.0;
    s = arc_std(a);
    if (s >= 999.0 || s <= 0.0) return 999.0;
    return s / sqrt(a->w_sum);
}

static int cmp_double(const void *a, const void *b)
{
    double da = *(const double*)a, db = *(const double*)b;
    return (da > db) - (da < db);
}

/* MAD-based robust outlier test — replaces 6-sigma gate.
 * Active from the very first epochs (Ge 2008 / Xu 2005 robust spirit). */
static int arc_mw_outlier(int sat, double x)
{
    const mw_rbuf_t *r;
    double tmp[MW_RBUF_N], med, mad, sigma_r, thr;
    int k, i;
    if (sat <= 0 || sat > MAXSAT) return 0;
    if (x != x || fabs(x) > 200.0) return 1;       /* hard sanity */

    r = &wl_rbuf[sat-1];
    if (r->n < 5) return 0;                        /* not enough seed */

    k = r->n;
    memcpy(tmp, r->buf, k * sizeof(double));
    qsort(tmp, k, sizeof(double), cmp_double);
    med = tmp[k/2];

    for (i = 0; i < k; i++) tmp[i] = fabs(r->buf[i] - med);
    qsort(tmp, k, sizeof(double), cmp_double);
    mad = tmp[k/2];

    sigma_r = 1.4826 * mad;                        /* Gaussian-consistent sigma */
    thr     = MAX(0.30, ROBUST_K * sigma_r);       /* floor 0.30 cy */
    return fabs(x - med) > thr;
}

extern void ambupd_reset(void)
{
    int i;
    ambupd_flush_all_arcs("reset");
    for (i = 0; i < MAXSAT; i++) {
        wl_arc[i].n      = 0.0;
        wl_arc[i].w_sum  = 0.0;
        wl_arc[i].mean   = 0.0;
        wl_arc[i].M2     = 0.0;
        wl_rbuf[i].n     = 0;
        wl_rbuf[i].head  = 0;
        wl_arc_init[i]   = 0;
        wl_last_init[i]  = 0;
        g_arc_rows_n[i]  = 0;
    }
    g_dcb_diag.total = g_dcb_diag.full = g_dcb_diag.fallback = g_dcb_diag.map_fail = 0UL;
    g_drop_if_sd = g_drop_arc_sh = g_drop_arc_sd = g_drop_outlr = g_arc_kept = 0UL;
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
 * Returns 0 (skip epoch) or 1 (use epoch), and sets *arc_reset=1
 * if this is the first epoch of a new AMB arc (KF state must reset).
 * Priority: DEL/BAD checked first across ALL records, then AMB.
 * This matches GREAT-UPD exactly: a DEL overlapping an AMB wins. */
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

static void get_pppopt_path(const prcopt_t *opt, const char *token,
                            char *path, size_t npath)
{
    const char *p;
    size_t n;
    if (!path || npath == 0) return;
    path[0] = '\0';
    if (!opt || !token || !*token || !opt->pppopt[0]) return;
    p = strstr(opt->pppopt, token);
    if (!p) return;
    p += strlen(token);
    for (; *p && isspace((unsigned char)*p); ++p) {
        /* skip spaces after token */
    }
    if (!*p) return;
    n = 0;
    while (*p && !isspace((unsigned char)*p) && n + 1 < npath) {
        path[n++] = *p++;
    }
    path[n] = '\0';
}

static void get_ambflag_path(const prcopt_t *opt, char *path, size_t npath)
{
    get_pppopt_path(opt, "-AMBFLAG=", path, npath);
}

static void get_bsx_path(const prcopt_t *opt, char *path, size_t npath)
{
    get_pppopt_path(opt, "-BSX_DCB=", path, npath);
}

/* Parse -AMBUPD_MAX_IF_SD=, -AMBUPD_MIN_ARC=, -AMBUPD_MAX_WL_STD= from pppopt.
 * Re-parsed each call (cheap, runs once per epoch). */
static void apply_gate_overrides(const prcopt_t *opt)
{
    char buf[64];
    double v;
    if (!opt) return;
    get_pppopt_path(opt, "-AMBUPD_MAX_IF_SD=", buf, sizeof(buf));
    if (buf[0] && (v = atof(buf)) > 0.0) g_max_if_sd = v;
    get_pppopt_path(opt, "-AMBUPD_MIN_ARC=", buf, sizeof(buf));
    if (buf[0] && (v = atof(buf)) > 0.0) g_min_arc_sec = v;
    get_pppopt_path(opt, "-AMBUPD_MAX_WL_STD=", buf, sizeof(buf));
    if (buf[0] && (v = atof(buf)) > 0.0) g_max_wl_std = v;
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

static void upper_copy(char *dst, size_t ndst, const char *src)
{
    size_t i;
    if (!dst || ndst == 0) return;
    dst[0] = '\0';
    if (!src) return;
    for (i = 0; src[i] && i + 1 < ndst; i++) {
        dst[i] = (char)toupper((unsigned char)src[i]);
    }
    dst[i] = '\0';
}

static void bsx_cache_reset(bsx_rcv_cache_t *cache)
{
    if (!cache) return;
    if (cache->rows) {
        free(cache->rows);
        cache->rows = NULL;
    }
    cache->n = 0;
    cache->nmax = 0;
    cache->station[0] = '\0';
    cache->path[0] = '\0';
    cache->day_key = -1;
}

static int bsx_cache_add(bsx_rcv_cache_t *cache, int code1, int code2, double bias_m)
{
    int i;
    if (!cache || code1 <= 0 || code2 <= 0) return 0;
    for (i = 0; i < cache->n; i++) {
        if (cache->rows[i].code1 == code1 && cache->rows[i].code2 == code2) {
            cache->rows[i].bias_m = bias_m;
            return 1;
        }
    }
    if (cache->n >= cache->nmax) {
        int nmax = cache->nmax + 64;
        rcv_dsb_t *tmp = (rcv_dsb_t *)realloc(cache->rows, sizeof(rcv_dsb_t) * nmax);
        if (!tmp) return 0;
        cache->rows = tmp;
        cache->nmax = nmax;
    }
    cache->rows[cache->n].code1 = code1;
    cache->rows[cache->n].code2 = code2;
    cache->rows[cache->n].bias_m = bias_m;
    cache->n++;
    return 1;
}

static int bsx_obs_to_code(const char *obs)
{
    char id[4] = "";
    size_t n;
    if (!obs || !*obs) return CODE_NONE;
    n = strlen(obs);
    if (n >= 3 && (obs[0] == 'C' || obs[0] == 'c')) {
        id[0] = obs[1];
        id[1] = (char)toupper((unsigned char)obs[2]);
        id[2] = '\0';
        return obs2code(id);
    }
    if (n >= 2 && isdigit((unsigned char)obs[0])) {
        id[0] = obs[0];
        id[1] = (char)toupper((unsigned char)obs[1]);
        id[2] = '\0';
        return obs2code(id);
    }
    return CODE_NONE;
}

static int bsx_time_to_daysec(const char *s, int *day_key, int *sec_of_day)
{
    int y = 0, d = 0, sod = 0;
    if (!s || sscanf(s, "%d:%d:%d", &y, &d, &sod) < 3) return 0;
    if (day_key) *day_key = y * 1000 + d;
    if (sec_of_day) *sec_of_day = sod;
    return 1;
}

static int bsx_day_in_range(int day_key, const char *start, const char *end)
{
    int ds = 0, de = 0, ss = 0, se = 0;
    if (!bsx_time_to_daysec(start, &ds, &ss) || !bsx_time_to_daysec(end, &de, &se)) {
        return 1; /* no strict time gate if parsing fails */
    }
    if (day_key < ds) return 0;
    if (day_key > de) return 0;
    if (day_key == de && se == 0) return 0; /* BSX end-time is exclusive */
    return 1;
}

static int bsx_unit_to_m(const char *unit, double value, double *out_m)
{
    char up[8] = "";
    size_t i;
    if (!unit || !out_m) return 0;
    for (i = 0; unit[i] && i + 1 < sizeof(up); i++) {
        up[i] = (char)toupper((unsigned char)unit[i]);
    }
    up[i] = '\0';

    if (!strcmp(up, "NS")) {
        *out_m = value * 1E-9 * CLIGHT;
        return 1;
    }
    if (!strcmp(up, "S")) {
        *out_m = value * CLIGHT;
        return 1;
    }
    if (!strcmp(up, "M")) {
        *out_m = value;
        return 1;
    }
    return 0;
}

static int bsx_load_receiver_dsb(const char *path, const char *station, int day_key,
                                 bsx_rcv_cache_t *cache)
{
    FILE *fp;
    char line[2048];
    int in_solution = 0;
    int nline = 0, nmatch = 0;
    char station_up[32] = "";

    if (!cache) return 0;
    bsx_cache_reset(cache);
    if (!path || !*path || !station || !*station) return 0;

    if (!(fp = fopen(path, "r"))) {
        trace(1, "ambupd: cannot open BSX DCB file: %s\n", path);
        return 0;
    }
    upper_copy(station_up, sizeof(station_up), station);

    while (fgets(line, sizeof(line), fp)) {
        char type[8] = "", c2[8] = "", c3[8] = "", sta[16] = "", sta_up[16] = "";
        char obs1[8] = "", obs2[8] = "", tstart[32] = "", tend[32] = "", unit[8] = "";
        double value = 0.0, bias_m = 0.0;
        int code1, code2;
        int ntok;

        if (!in_solution) {
            if (strstr(line, "+BIAS/SOLUTION")) in_solution = 1;
            continue;
        }
        if (strstr(line, "-BIAS/SOLUTION")) break;
        if (line[0] == '*') continue;

        ntok = sscanf(line, "%7s %7s %7s %15s %7s %7s %31s %31s %7s %lf",
                      type, c2, c3, sta, obs1, obs2, tstart, tend, unit, &value);
        if (ntok < 10) continue;
        if (strcmp(type, "DSB")) continue;
        if (strcmp(c2, "G") || strcmp(c3, "G")) continue; /* receiver GPS block */
        nline++;
        upper_copy(sta_up, sizeof(sta_up), sta);
        if (strcmp(sta_up, station_up)) continue;
        if (!bsx_day_in_range(day_key, tstart, tend)) continue;

        code1 = bsx_obs_to_code(obs1);
        code2 = bsx_obs_to_code(obs2);
        if (code1 <= 0 || code2 <= 0) continue;
        if (!bsx_unit_to_m(unit, value, &bias_m)) continue;
        if (!bsx_cache_add(cache, code1, code2, bias_m)) continue;
        nmatch++;
    }
    fclose(fp);

    strncpy(cache->path, path, sizeof(cache->path) - 1);
    cache->path[sizeof(cache->path) - 1] = '\0';
    strncpy(cache->station, station_up, sizeof(cache->station) - 1);
    cache->station[sizeof(cache->station) - 1] = '\0';
    cache->day_key = day_key;

    trace(2,
          "ambupd: BSX receiver DSB station=%s day=%d path=%s rows=%d matches=%d\n",
          cache->station, day_key, path, nline, nmatch);
    return cache->n;
}

static void bsx_prepare_from_opt(const prcopt_t *opt, gtime_t t, const char *station)
{
    static char loaded_station[32] = "";
    static char loaded_path[1024] = "";
    static int loaded_day_key = -1;
    char path[1024] = "";
    char station_up[32] = "";
    int day_key;

    if (!opt || !station) return;
    day_key = get_day_key(t);
    upper_copy(station_up, sizeof(station_up), station);
    get_bsx_path(opt, path, sizeof(path));

    if (!strcmp(loaded_station, station_up) &&
        !strcmp(loaded_path, path) &&
        loaded_day_key == day_key) {
        return;
    }

    bsx_cache_reset(&g_bsx_rcv);
    g_bsx_path[0] = '\0';
    if (path[0]) {
        (void)bsx_load_receiver_dsb(path, station_up, day_key, &g_bsx_rcv);
        strncpy(g_bsx_path, path, sizeof(g_bsx_path) - 1);
        g_bsx_path[sizeof(g_bsx_path) - 1] = '\0';
        if (!g_bsx_diag_notice) {
            trace(2,
                  "ambupd: WL DCB uses sat-only mode; receiver BSX is diagnostics-only\n");
            g_bsx_diag_notice = 1;
        }
    } else {
        trace(3, "ambupd: no -BSX_DCB= for station=%s day=%d\n", station_up, day_key);
    }

    strncpy(loaded_station, station_up, sizeof(loaded_station) - 1);
    loaded_station[sizeof(loaded_station) - 1] = '\0';
    strncpy(loaded_path, path, sizeof(loaded_path) - 1);
    loaded_path[sizeof(loaded_path) - 1] = '\0';
    loaded_day_key = day_key;
}

static char code_attr(uint8_t code)
{
    char *obsid = code2obs(code);  /* ex: "1C", "2W" */
    if (!obsid || !obsid[0] || !obsid[1]) return '\0';
    return (char)toupper((unsigned char)obsid[1]);
}

static int code_band(uint8_t code)
{
    char *obsid = code2obs(code);  /* ex: "1C", "2W" */
    if (!obsid || !obsid[0]) return 0;
    if (obsid[0] == '1') return 1;
    if (obsid[0] == '2') return 2;
    return 0;
}

/* Match GREAT select_range/select_phase priority for GPS.
 * Note: GREAT picks the highest index in these strings. */
static int gps_attr_rank(int band, char attr)
{
    const char *order = NULL, *p = NULL;
    attr = (char)toupper((unsigned char)attr);
    if (band == 1) order = "CSLXPWYM";
    else if (band == 2) order = "CDLXPWYM";
    else return -1;
    p = strchr(order, attr);
    return p ? (int)(p - order) : -1;
}

static int select_l1_l2_indices(const obsd_t *obs, const nav_t *nav, int sat,
                                int *k1, int *k2, double *f1, double *f2)
{
    int i;
    int best1 = -1, best2 = -1;
    int best1_fb = -1, best2_fb = -1;
    int rank1 = -1, rank2 = -1;
    double f1_fb = 0.0, f2_fb = 0.0;
    double fi;

    if (!obs || !nav || !k1 || !k2 || !f1 || !f2) return 0;

    for (i = 0; i < NFREQ + NEXOBS; i++) {
        int r;
        if (!obs->code[i] || obs->L[i] == 0.0 || obs->P[i] == 0.0) continue;
        fi = sat2freq(sat, obs->code[i], nav);
        if (fi <= 0.0) continue;

        /* Keep fallback to avoid holes if code attr is unknown. */
        if (fabs(fi - FREQ1) <= 1e6) {
            if (best1_fb < 0) { best1_fb = i; f1_fb = fi; }
            r = gps_attr_rank(1, code_attr(obs->code[i]));
            if (r > rank1 || (r == rank1 && best1 >= 0 && i < best1)) {
                rank1 = r;
                best1 = i;
                *f1 = fi;
            }
        }
        else if (fabs(fi - FREQ2) <= 1e6) {
            if (best2_fb < 0) { best2_fb = i; f2_fb = fi; }
            r = gps_attr_rank(2, code_attr(obs->code[i]));
            if (r > rank2 || (r == rank2 && best2 >= 0 && i < best2)) {
                rank2 = r;
                best2 = i;
                *f2 = fi;
            }
        }
    }

    if (best1 < 0 && best1_fb >= 0) { best1 = best1_fb; *f1 = f1_fb; }
    if (best2 < 0 && best2_fb >= 0) { best2 = best2_fb; *f2 = f2_fb; }
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

enum {
    DCB_REQ_P1P2 = 1 << 0,
    DCB_REQ_P1C1 = 1 << 1,
    DCB_REQ_P2C2 = 1 << 2
};

static int map_code_to_corr_req(uint8_t code, const double cor[5], double *corr, int *req_mask)
{
    int band;
    char attr;
    if (!corr || !req_mask) return 0;
    *req_mask = 0;
    band = code_band(code);
    attr = code_attr(code);

    if (band == 1) {
        if (attr == 'C') {                                      /* C1C */
            *corr = cor[2];
            *req_mask = DCB_REQ_P1P2 | DCB_REQ_P1C1;
            return 1;
        }
        if (attr == 'P' || attr == 'Y' || attr == 'W') {       /* C1P/C1Y/C1W */
            *corr = cor[0];
            *req_mask = DCB_REQ_P1P2;
            return 1;
        }
        *corr = 0.0;
        return 0;
    }
    if (band == 2) {
        if (attr == 'C') {                                      /* C2C */
            *corr = cor[3];
            *req_mask = DCB_REQ_P1P2 | DCB_REQ_P2C2;
            return 1;
        }
        if (attr == 'P' || attr == 'Y' || attr == 'W') {       /* C2P/C2Y/C2W */
            *corr = cor[1];
            *req_mask = DCB_REQ_P1P2;
            return 1;
        }
        if (attr == 'X') {                                      /* C2X */
            *corr = cor[4];
            *req_mask = DCB_REQ_P1P2 | DCB_REQ_P1C1;
            return 1;
        }
        if (attr == 'D') { *corr = 0.0; return 1; }            /* C2D */
        *corr = 0.0;
        return 0;
    }
    *corr = 0.0;
    return 0;
}

static int mw_dcb_corr_cycles(const obsd_t *obs, const nav_t *nav, int sat,
                              int k1, int k2, double lam1, double lam2,
                              double *corr_cyc, double *corr_sat_cyc,
                              double *corr_rcv_cyc, dcb_mode_t *mode)
{
    double p1p2 = 0.0, p1c1 = 0.0, p2c2 = 0.0;
    double cor[5];
    double corr1 = 0.0, corr2 = 0.0;
    double fact, c1, c2;
    int has_p1p2, has_p1c1, has_p2c2;
    int ok1, ok2;
    int req1 = 0, req2 = 0;
    int req_mask = 0, have_mask = 0;

    if (mode) *mode = DCB_MODE_NONE;
    if (!obs || !nav || !corr_cyc || lam1 <= 0.0 || lam2 <= 0.0) return 0;
    *corr_cyc = 0.0;
    if (corr_sat_cyc) *corr_sat_cyc = 0.0;
    if (corr_rcv_cyc) *corr_rcv_cyc = 0.0;

    /* Strict GREAT sat-only DCB for MW correction. Receiver DSB is diagnostics-only. */
    has_p1p2 = get_dsb(nav, sat, CODE_L1W, CODE_L2W, &p1p2);
    has_p1c1 = get_dsb(nav, sat, CODE_L1W, CODE_L1C, &p1c1);
    has_p2c2 = get_dsb(nav, sat, CODE_L2W, CODE_L2C, &p2c2);

    fact = lam1 / lam2;
    if (fabs(1.0 - fact * fact) < 1E-12) return 0;
    c1 = 1.0 / (1.0 - fact * fact);
    c2 = fact * fact * c1;

    cor[0] = c2 * p1p2;
    cor[1] = c1 * p1p2;
    cor[2] = cor[0] + p1c1;
    cor[3] = cor[1] + p2c2;
    cor[4] = cor[1] + p1c1;

    ok1 = map_code_to_corr_req(obs->code[k1], cor, &corr1, &req1);
    ok2 = map_code_to_corr_req(obs->code[k2], cor, &corr2, &req2);
    g_dcb_diag.total++;
    if (!ok1 || !ok2) {
        g_dcb_diag.map_fail++;
        if (mode) *mode = DCB_MODE_MAP_FAIL;
        return 0;
    }

    req_mask = req1 | req2;
    if (has_p1p2) have_mask |= DCB_REQ_P1P2;
    if (has_p1c1) have_mask |= DCB_REQ_P1C1;
    if (has_p2c2) have_mask |= DCB_REQ_P2C2;

    if ((req_mask & have_mask) == req_mask) {
        g_dcb_diag.full++;
        if (mode) *mode = DCB_MODE_FULL;
    }
    else {
        g_dcb_diag.fallback++;
        if (mode) *mode = DCB_MODE_FALLBACK;
    }

    *corr_cyc = -(corr1 / lam1 + corr2 / lam2) * (1.0 - fact) / (1.0 + fact);
    if (corr_sat_cyc) {
        *corr_sat_cyc = *corr_cyc;
    }
    if (corr_rcv_cyc) {
        *corr_rcv_cyc = 0.0;
    }
    return 1;
}

/* ---------------------------------------------------------------------------
 * mw_raw() — compute Melbourne-Wübbena combination from RAW observations.
 *
 * Implements the same structure as GREAT-UPD:
 *   1) MW from raw code/phase
 *   2) additional satellite-only DCB correction term based on P1P2/P1C1/P2C2 model.
 *
 * Returns MW in cycles, or 0.0 on failure (sets *ok=0).
 * --------------------------------------------------------------------------*/
static double mw_raw(const obsd_t *obs, const nav_t *nav,
                     int sat, int k1, int k2,
                     double f1, double f2, int *ok,
                     double *mw_base, double *mw_dcb,
                     double *mw_dcb_sat, double *mw_dcb_rcv,
                     dcb_mode_t *dcb_mode)
{
    double lam1, lam2, fact;
    double C1, C2, L1, L2;
    double mw, dcb_corr = 0.0, dcb_sat = 0.0, dcb_rcv = 0.0;

    *ok = 0;
    if (mw_base) *mw_base = 0.0;
    if (mw_dcb) *mw_dcb = 0.0;
    if (mw_dcb_sat) *mw_dcb_sat = 0.0;
    if (mw_dcb_rcv) *mw_dcb_rcv = 0.0;
    if (dcb_mode) *dcb_mode = DCB_MODE_NONE;

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
    if (mw_base) *mw_base = mw;

    /* GREAT-like DCB term (if required pairs are available). */
    (void)mw_dcb_corr_cycles(obs, nav, sat, k1, k2, lam1, lam2,
                             &dcb_corr, &dcb_sat, &dcb_rcv, dcb_mode);
    if (mw_dcb) *mw_dcb = dcb_corr;
    if (mw_dcb_sat) *mw_dcb_sat = dcb_sat;
    if (mw_dcb_rcv) *mw_dcb_rcv = dcb_rcv;
    mw += dcb_corr;

    *ok = 1;
    return mw;
}

static gtime_t normalize_epoch_time(gtime_t t)
{
    double ep[6], sod;
    time2epoch(t, ep);
    sod = ep[3] * 3600.0 + ep[4] * 60.0 + ep[5];
    if (sod >= 86399.9995) {
        /* Prevent 24:00:00 style carry-over into previous file/day key. */
        t = timeadd(t, 0.001);
    }
    return t;
}

static void trace_dcb_diag(const char *tag)
{
    trace(3,
          "ambupd: DCB diag [%s] total=%lu full=%lu fallback=%lu map_fail=%lu\n",
          tag ? tag : "-", g_dcb_diag.total, g_dcb_diag.full,
          g_dcb_diag.fallback, g_dcb_diag.map_fail);
}

static const char *dcb_mode_str(dcb_mode_t mode)
{
    switch (mode) {
    case DCB_MODE_FULL:     return "full";
    case DCB_MODE_FALLBACK: return "fallback";
    case DCB_MODE_MAP_FAIL: return "map_fail";
    default:                return "none";
    }
}

static int arc_rows_push(int sat, const amb_epoch_row_t *row)
{
    int idx, nmax;
    amb_epoch_row_t *tmp;
    if (!row || sat <= 0 || sat > MAXSAT) return 0;
    idx = sat - 1;
    if (g_arc_rows_n[idx] >= g_arc_rows_nmax[idx]) {
        nmax = g_arc_rows_nmax[idx] + 128;
        tmp = (amb_epoch_row_t *)realloc(g_arc_rows[idx], sizeof(amb_epoch_row_t) * nmax);
        if (!tmp) return 0;
        g_arc_rows[idx] = tmp;
        g_arc_rows_nmax[idx] = nmax;
    }
    g_arc_rows[idx][g_arc_rows_n[idx]++] = *row;
    return 1;
}

static int pending_push(const amb_epoch_row_t *row, double wl_mean, double sigma_wl)
{
    amb_pending_row_t *tmp;
    int nmax;
    if (!row) return 0;
    if (g_pending_n >= g_pending_nmax) {
        nmax = g_pending_nmax + 4096;
        tmp = (amb_pending_row_t *)realloc(g_pending, sizeof(amb_pending_row_t) * nmax);
        if (!tmp) return 0;
        g_pending = tmp;
        g_pending_nmax = nmax;
    }
    g_pending[g_pending_n].row      = *row;
    g_pending[g_pending_n].wl_mean  = wl_mean;
    g_pending[g_pending_n].sigma_wl = sigma_wl;
    g_pending_n++;
    return 1;
}

static int pending_cmp(const void *a, const void *b)
{
    const amb_pending_row_t *pa = (const amb_pending_row_t *)a;
    const amb_pending_row_t *pb = (const amb_pending_row_t *)b;
    if (pa->row.mjd < pb->row.mjd) return -1;
    if (pa->row.mjd > pb->row.mjd) return  1;
    if (pa->row.sod < pb->row.sod) return -1;
    if (pa->row.sod > pb->row.sod) return  1;
    return strcmp(pa->row.satid, pb->row.satid);
}

static void pending_flush_to_file(void)
{
    int i;
    if (g_pending_n <= 0) return;
    if (!g_fp_amb) { g_pending_n = 0; return; }

    qsort(g_pending, g_pending_n, sizeof(amb_pending_row_t), pending_cmp);
    for (i = 0; i < g_pending_n; i++) {
        const amb_epoch_row_t *r = &g_pending[i].row;
        double wl_mean  = g_pending[i].wl_mean;
        double sigma_wl = g_pending[i].sigma_wl;
        fprintf(g_fp_amb, "%8.0f%10.1f%5s%4s%19.3f%19.3f%10.3f\n",
                floor(r->mjd), r->sod, r->station, r->satid,
                r->if_out, wl_mean, sigma_wl);
        if (g_fp_diag) {
            fprintf(g_fp_diag,
                    "%s,%.1f,%s,%s,%s,%d,%d,%.3f,%.3f,%.6f,%.6f,%.6f,%.6f,%s,%s,%d,%.6f,%.6f\n",
                    r->station, r->sod, r->satid, r->code1, r->code2,
                    r->k1, r->k2, r->f1, r->f2, r->mw_raw, r->mw_dcb,
                    r->mw_dcb_sat, r->mw_dcb_rcv, dcb_mode_str(r->dcb_mode),
                    "sat-only", r->bsx_rcv_rows,
                    wl_mean, r->if_out);
        }
    }
    g_pending_n = 0;
}

static void ambupd_flush_sat_arc(int sat, const char *reason)
{
    int idx, i;
    double wl_mean, sigma_wl, arc_sec;
    if (sat <= 0 || sat > MAXSAT) return;
    idx = sat - 1;
    if (g_arc_rows_n[idx] <= 0) return;

    if (!g_fp_amb) {
        g_arc_rows_n[idx] = 0;
        return;
    }

    wl_mean  = wl_arc[idx].mean;
    sigma_wl = arc_sigma_mean(&wl_arc[idx]);    /* GREAT writes SE-of-mean */
    arc_sec  = wl_arc[idx].n * g_sample_dt;

    /* Ge 2008 quality gates: drop short or noisy arcs entirely */
    if (arc_sec < g_min_arc_sec) {
        g_drop_arc_sh++;
        trace(2, "ambupd: drop arc sat=%d arc=%.0fs std=%.3f cy reason=%s "
                 "(min_arc=%.0fs)\n",
              sat, arc_sec, sigma_wl, reason ? reason : "-", g_min_arc_sec);
        g_arc_rows_n[idx] = 0;
        return;
    }
    if (sigma_wl > g_max_wl_std) {
        g_drop_arc_sd++;
        trace(2, "ambupd: drop arc sat=%d arc=%.0fs std=%.3f cy reason=%s "
                 "(max_std=%.3f)\n",
              sat, arc_sec, sigma_wl, reason ? reason : "-", g_max_wl_std);
        g_arc_rows_n[idx] = 0;
        return;
    }
    g_arc_kept++;

    for (i = 0; i < g_arc_rows_n[idx]; i++) {
        if (!pending_push(&g_arc_rows[idx][i], wl_mean, sigma_wl)) {
            trace(1, "ambupd: pending buffer alloc failed sat=%d\n", sat);
            break;
        }
    }
    trace(4, "ambupd: pend sat=%d rows=%d reason=%s wl=%.4f sig=%.4f\n",
          sat, g_arc_rows_n[idx], reason ? reason : "-", wl_mean, sigma_wl);
    g_arc_rows_n[idx] = 0;
}

static void ambupd_flush_all_arcs(const char *reason)
{
    int sat;
    for (sat = 1; sat <= MAXSAT; sat++) {
        ambupd_flush_sat_arc(sat, reason);
    }
    pending_flush_to_file();
    if (g_fp_amb) fflush(g_fp_amb);
    if (g_fp_diag) fflush(g_fp_diag);
}

static void ambupd_close_outputs(void)
{
    if (g_fp_diag) {
        fclose(g_fp_diag);
        g_fp_diag = NULL;
    }
    if (g_fp_amb) {
        fclose(g_fp_amb);
        g_fp_amb = NULL;
    }
    g_amb_file[0] = '\0';
    g_diag_file[0] = '\0';
}

static void ambupd_gate_diag(const char *tag)
{
    trace(2,
          "ambupd: gate diag [%s] kept=%lu drop_if=%lu drop_short=%lu "
          "drop_std=%lu drop_outlier=%lu  (max_if_sd=%.3f m, min_arc=%.0fs, "
          "max_wl_std=%.3f cy, dt=%.1fs)\n",
          tag ? tag : "-",
          g_arc_kept, g_drop_if_sd, g_drop_arc_sh, g_drop_arc_sd, g_drop_outlr,
          g_max_if_sd, g_min_arc_sec, g_max_wl_std, g_sample_dt);
}

static void ambupd_atexit_hook(void)
{
    int i;
    ambupd_flush_all_arcs("atexit");
    trace_dcb_diag("atexit");
    ambupd_gate_diag("atexit");
    ambupd_close_outputs();
    bsx_cache_reset(&g_bsx_rcv);
    for (i = 0; i < MAXSAT; i++) {
        if (g_arc_rows[i]) {
            free(g_arc_rows[i]);
            g_arc_rows[i] = NULL;
        }
        g_arc_rows_n[i] = 0;
        g_arc_rows_nmax[i] = 0;
    }
    if (g_pending) {
        free(g_pending);
        g_pending = NULL;
    }
    g_pending_n = 0;
    g_pending_nmax = 0;
}

/* ---------------------------------------------------------------- */
extern void ambupd_write_epoch(const rtk_t *rtk, const obsd_t *obs, int n,
                               const nav_t *nav)
{
    int week, year, doy, dow;
    double tow, mjd;
    double ep[6], time_of_day;
    gtime_t t_epoch;

    char station_name[32];
    const char *src_name;
    char *p;

    char temp_amb_file[1024];
    char temp_diag_file[1024];
    char satid[8];

    int iobs;
    int sat, sys;
    int k1, k2;
    double f1 = 0.0, f2 = 0.0;
    double gap_reset;
    double Nwl_mw;
    double Nwl_mw_base = 0.0, Nwl_mw_dcb = 0.0, Nwl_mw_dcb_sat = 0.0, Nwl_mw_dcb_rcv = 0.0;
    double N_if;
    int idx_if;
    int mw_ok;
    dcb_mode_t dcb_mode = DCB_MODE_NONE;
    amb_epoch_row_t row;
    const char *c1, *c2;

    if (!(rtk->opt.ionoopt == IONOOPT_IFLC && rtk->opt.nf >= 2)) return;
    if (n <= 0) return;
    t_epoch = normalize_epoch_time(obs[0].time);

    /* Track actual sampling interval for arc-length gate */
    {
        double dt = fabs(rtk->tt);
        if (dt > 0.001 && dt < 3600.0) g_sample_dt = dt;
    }
    apply_gate_overrides(&rtk->opt);

    if (!g_ambupd_atexit_registered) {
        if (atexit(ambupd_atexit_hook) == 0) {
            g_ambupd_atexit_registered = 1;
        }
    }

    src_name = rtk->opt.station_name[0] ? rtk->opt.station_name : "UNKN";
    strncpy(station_name, src_name, sizeof(station_name) - 1);
    station_name[sizeof(station_name) - 1] = '\0';
    for (p = station_name; *p; p++) {
        *p = (char)toupper((unsigned char)*p);
    }

    tow = time2gpst(t_epoch, &week);
    time2epoch(t_epoch, ep);
    year = (int)ep[0];
    doy = time2doy(t_epoch);
    sprintf(temp_amb_file, "%s_ambupd_%4d%03d", station_name, year, doy);
    sprintf(temp_diag_file, "%s_ambdiag_%4d%03d.csv", station_name, year, doy);

    bsx_prepare_from_opt(&rtk->opt, t_epoch, station_name);

    if (g_fp_amb == NULL || strcmp(temp_amb_file, g_amb_file) != 0) {
        if (g_fp_amb) {
            trace_dcb_diag("close-file");
            ambupd_reset();
            ambupd_close_outputs();
        }
        strncpy(g_amb_file, temp_amb_file, sizeof(g_amb_file) - 1);
        g_amb_file[sizeof(g_amb_file) - 1] = '\0';
        g_fp_amb = fopen(g_amb_file, "w");
        if (!g_fp_amb) {
            trace(1, "ambupd: cannot create file: %s\n", g_amb_file);
            g_amb_file[0] = '\0';
            return;
        }

        strncpy(g_diag_file, temp_diag_file, sizeof(g_diag_file) - 1);
        g_diag_file[sizeof(g_diag_file) - 1] = '\0';
        g_fp_diag = fopen(g_diag_file, "w");
        if (!g_fp_diag) {
            trace(1, "ambupd: cannot create debug file: %s\n", g_diag_file);
            g_diag_file[0] = '\0';
        }
        else {
            fprintf(g_fp_diag,
                    "station,sod,sat,code1,code2,k1,k2,f1,f2,mw_raw,mw_dcb,mw_dcb_sat,mw_dcb_rcv,dcb_mode,mw_dcb_source,bsx_rcv_rows,wl_mean,if_state\n");
        }
        trace(1, "ambupd: creating file: %s\n", g_amb_file);
    }

    dow = (int)floor(tow / 86400.0);
    time_of_day = tow - dow * 86400.0;
    if (time_of_day >= 86399.9995) {
        time_of_day = 0.0;
        dow += 1;
    }
    while (time_of_day < 0.0) {
        time_of_day += 86400.0;
        dow -= 1;
    }
    while (dow < 0) {
        dow += 7;
        week -= 1;
    }
    while (dow > 6) {
        dow -= 7;
        week += 1;
    }
    mjd = 44244.0 + (week * 7.0) + dow;

    gap_reset = MAX(3.0 * fabs(rtk->tt), MIN_GAP_RESET_S);

    for (iobs = 0; iobs < n; iobs++) {
        sat = obs[iobs].sat;
        if (sat <= 0 || sat > MAXSAT) continue;

        sys = satsys(sat, NULL);
        if (sys != SYS_GPS) continue;

        if (!select_l1_l2_indices(obs + iobs, nav, sat, &k1, &k2, &f1, &f2)) {
            trace(4, "ambupd: skip sat=%d no valid L1/L2 pair in obs record\n", sat);
            continue;
        }
        if (fabs(f1 - f2) < 1.0) continue;

        idx_if = IB(sat, 0, &rtk->opt);
        if (idx_if < 0 || idx_if >= rtk->nx) continue;
        N_if = rtk->x[idx_if];
        if (fabs(N_if) < 1e-12) continue;

        {
            int arc_reset_flag = 0;
            if (!ambflag_is_valid(sat, obs[iobs].time, &g_ambflag, &arc_reset_flag)) {
                trace(4, "ambupd: sat=%d excluded by ambflag DEL/BAD\n", sat);
                continue;
            }
            if (arc_reset_flag) {
                ambupd_flush_sat_arc(sat, "ambflag-reset");
                arc_reset(sat);
                wl_last_init[sat - 1] = 0;
                trace(3, "ambupd: WL arc reset by ambflag AMB sat=%d\n", sat);
            }
        }

        if (!wl_arc_init[sat - 1]) arc_reset(sat);
        if (wl_last_init[sat - 1] &&
            fabs(timediff(obs[iobs].time, wl_last_obs[sat - 1])) > gap_reset &&
            g_ambflag.n == 0) {
            ambupd_flush_sat_arc(sat, "gap-reset");
            arc_reset(sat);
            wl_last_init[sat - 1] = 0;
        }

        Nwl_mw = mw_raw(obs + iobs, nav, sat, k1, k2, f1, f2, &mw_ok,
                        &Nwl_mw_base, &Nwl_mw_dcb,
                        &Nwl_mw_dcb_sat, &Nwl_mw_dcb_rcv, &dcb_mode);
        if (!mw_ok) continue;

        /* --- IF maturity check (Song et al. 2026; avoids leaking
         *     unconverged float-PPP IF state into ambupd output). */
        {
            int idx_p = idx_if + idx_if * rtk->nx;
            double sd_if = (idx_p >= 0 && idx_p < rtk->nx * rtk->nx)
                           ? sqrt(MAX(rtk->P[idx_p], 0.0)) : 1e9;
            if (sd_if > g_max_if_sd) {
                g_drop_if_sd++;
                trace(4, "ambupd: sat=%d IF not mature sd=%.3f m (limit=%.3f)\n",
                      sat, sd_if, g_max_if_sd);
                continue;   /* don't write this epoch yet */
            }
        }

        /* --- Elevation-weighted MW update + robust outlier check ----- */
        {
            double w = elev_weight(rtk->ssat[sat-1].azel[1]);
            if (w <= 0.0) continue;
            if (arc_mw_outlier(sat, Nwl_mw)) {
                g_drop_outlr++;
                trace(4, "ambupd: sat=%d MW outlier rejected raw=%.3f mean=%.3f n=%.0f\n",
                      sat, Nwl_mw, wl_arc[sat - 1].mean, wl_arc[sat - 1].n);
                continue;
            }
            arc_update(sat, Nwl_mw, w);
        }
        wl_last_obs[sat - 1] = obs[iobs].time;
        wl_last_init[sat - 1] = 1;

        satno2id(sat, satid);
        trace(3, "ambupd: sat=%-4s idx_if=%3d IF=%10.4f m  "
              "WL_raw=%7.4f cy WL_dcb=%7.4f cy sat=%7.4f cy rcv=%7.4f cy dcb=%s n=%.0f\n",
              satid, idx_if, N_if, Nwl_mw_base, Nwl_mw_dcb,
              Nwl_mw_dcb_sat, Nwl_mw_dcb_rcv, dcb_mode_str(dcb_mode),
              wl_arc[sat - 1].n);

        memset(&row, 0, sizeof(row));
        row.mjd = floor(mjd);
        row.sod = time_of_day;
        row.if_out = N_if;
        row.mw_raw = Nwl_mw_base;
        row.mw_dcb = Nwl_mw_dcb;
        row.mw_dcb_sat = Nwl_mw_dcb_sat;
        row.mw_dcb_rcv = Nwl_mw_dcb_rcv;
        row.dcb_mode = dcb_mode;
        row.bsx_rcv_rows = g_bsx_rcv.n;
        row.k1 = k1;
        row.k2 = k2;
        row.f1 = f1;
        row.f2 = f2;
        strncpy(row.station, station_name, sizeof(row.station) - 1);
        row.station[sizeof(row.station) - 1] = '\0';
        strncpy(row.satid, satid, sizeof(row.satid) - 1);
        row.satid[sizeof(row.satid) - 1] = '\0';
        c1 = code2obs(obs[iobs].code[k1]);
        c2 = code2obs(obs[iobs].code[k2]);
        strncpy(row.code1, (c1 && *c1) ? c1 : "", sizeof(row.code1) - 1);
        row.code1[sizeof(row.code1) - 1] = '\0';
        strncpy(row.code2, (c2 && *c2) ? c2 : "", sizeof(row.code2) - 1);
        row.code2[sizeof(row.code2) - 1] = '\0';
        if (!arc_rows_push(sat, &row)) {
            trace(1, "ambupd: arc buffer allocation failed sat=%d\n", sat);
        }
    }

    if (g_dcb_diag.total > 0UL && (g_dcb_diag.total % 50000UL) == 0UL) {
        trace_dcb_diag("periodic");
    }
}
