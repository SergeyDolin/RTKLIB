/* ambupd.c : write ambupd file (GREAT-UPD style) from the PPP state
 * The file contains:
 *     1) Modified Julian Day
 *     2) Time of day, (unit: second)
 *     3) Station name
 *     4) PRN number of GNSS satellite
 *     5) IF ambiguity derived from PPP float solution, (unit: meter)
 *     6) WL ambiguity derived from Melbourne-Wübbena combination, (unit: cycle)
 *     7) Standard deviation of WL ambiguity.
 * The ambupd_write_epoch() function is intended to be called from pppos() once per epoch. */

#include "rtklib.h"
#include <ctype.h>
#include <string.h>
#include <math.h>

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

extern void ambupd_write_epoch(const rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
    /* Output IF and WL ambiguities */
    if (rtk->opt.ionoopt == IONOOPT_IFLC && rtk->opt.nf >= 2) {
        int week;
        double tow = time2gpst(rtk->sol.time, &week);
        char satid[8];
        int sat, i, idx1, idx2;
        double N1, N2, f1, f2, N_if, N_wl, sigma_if, sigma_wl, sigma;
        const obsd_t *pobs;

        static FILE *fp_amb = NULL;
        static char amb_file[1024] = "";
        static int last_week = -1;
        static double last_tow = -1.0;
        double ep[6];
        double mjd;

        /* Use station name from options if provided, otherwise default to "UNKN" */
        char station_name[32];
        const char *src_name = rtk->opt.station_name[0] ? rtk->opt.station_name : "UNKN";
        strncpy(station_name, src_name, sizeof(station_name) - 1);
        station_name[sizeof(station_name) - 1] = '\0';
        char *p;
        for (p = station_name; *p; p++) {
            *p = (char)toupper((unsigned char)*p);
        }

        /* Threshold for extreme/incorrect values */
        const double threshold = 1000000.0;

        /* Create filename with date if not already created */
        int year, doy;
        time2epoch(rtk->sol.time, ep);
        year = (int)ep[0];
        doy = time2doy(rtk->sol.time);
        char temp_amb_file[1024];
        sprintf(temp_amb_file, "%s_ambupd_%4d%03d", station_name, year, doy);

        if (fp_amb == NULL || strcmp(temp_amb_file, amb_file) != 0) {
            /* Close existing file if open */
            if (fp_amb != NULL) {
                fclose(fp_amb);
                fp_amb = NULL;
            }

            strcpy(amb_file, temp_amb_file);

            fp_amb = fopen(amb_file, "w");
            if (fp_amb == NULL) {
                trace(1, "Cannot create ambiguity file: %s\n", amb_file);
            } else {
                trace(1, "Creating ambiguity file: %s\n", amb_file);
            }
        }

        /* Convert GPS time to Modified Julian Day */
        mjd = 44244.0 + (week * 7.0) + (tow / 86400.0);

        /* Extract time of day from epoch */
        time2epoch(rtk->sol.time, ep);
        double time_of_day = ep[3] * 3600.0 + ep[4] * 60.0 + ep[5];

        /* Output all satellites for this epoch */
        for (sat = 1; sat <= MAXSAT; sat++) {
            if (!rtk->ssat[sat-1].vs) continue;

            idx1 = IB(sat, 0, &rtk->opt); /* 0 = L1 */
            idx2 = IB(sat, 2, &rtk->opt); /* 2 = L5 */

            N1 = rtk->x[idx1];
            N2 = rtk->x[idx2];

            /* if (N1 == 0.0 || N2 == 0.0) continue; */
            if (fabs(N1) < 1e-12 && fabs(N2) < 1e-12) continue;

            /* Find observation for this satellite */
            pobs = NULL;
            int i;
            for (i = 0; i < n; i++) {
                if (obs[i].sat == sat) {
                    pobs = &obs[i];
                    break;
                }
            }
            if (!pobs) continue;

            /* Get carrier frequencies */
            f1 = sat2freq(sat, pobs->code[0], nav);
            f2 = sat2freq(sat, pobs->code[2], nav);
            if (f1 <= 0.0 || f2 <= 0.0) continue;

            /* Cycle to meters */
            double lambda1 = CLIGHT / f1;
            double lambda2 = CLIGHT / f2;
            double phi1_m = N1 * lambda1;  /* L1  */
            double phi2_m = N2 * lambda2;  /* L5  */

            /* IF  */
            double denom_if = SQR(f1) - SQR(f2);
            N_if = (SQR(f1) * phi1_m - SQR(f2) * phi2_m) / denom_if;

            /* Calculate WL */
            double denom_wl = f1 - f2;
            N_wl = (f1 * N1 - f2 * N2) / denom_wl;

            /* Check if values are reasonable*/
            /*if (fabs(N_if) > threshold || fabs(N_wl) > threshold) {
                continue;
            } */

            /* Compute approximate standard deviation */
            double var1 = rtk->P[idx1 + idx1 * rtk->nx];
            double var2 = rtk->P[idx2 + idx2 * rtk->nx];
            double cov12 = rtk->P[idx1 + idx2 * rtk->nx];  /* off-diagonal */

            /* Standard deviation for WL ambiguity */
            double coef_wl_1 = f1 / denom_wl;
            double coef_wl_2 = -f2 / denom_wl;
            sigma = SQRT(MAX(0.0, SQR(coef_wl_1)*var1 + SQR(coef_wl_2)*var2 + 2.0*coef_wl_1*coef_wl_2*cov12));

            satno2id(sat, satid);

            /* Output to file */
            if (fp_amb != NULL) {
                fprintf(fp_amb, "%8.0f%10.1f%5s%4s%19.3f%19.3f%10.3f\n",
                    floor(mjd), time_of_day, station_name, satid, N_if, N_wl, sigma);
            }
        }

        /* Flush file to ensure data is written */
        if (fp_amb != NULL) {
            fflush(fp_amb);
        }
    }
}