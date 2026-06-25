/*------------------------------------------------------------------------------
* gpt3.c : Global Pressure and Temperature 3 (GPT3) model with VMF3 mapping
*
* references :
*   [1] Landskron D, Böhm J (2018) VMF3/GPT3: Refined discrete and empirical
*       troposphere mapping functions. J Geodesy 92(4):349-360.
*       doi:10.1007/s00190-017-1066-2
*   [2] Askne J, Nordius (1987) Estimation of tropospheric delay for
*       microwaves from surface weather data. Radio Sci 22(3):379-386.
*   [3] Niell AE (1996) Global mapping functions for the atmosphere delay
*       at radio wavelengths. J Geophys Res 101(B2):3227-3246.
*
* Grid file gpt3_5.grd: 36 lat x 72 lon = 2592 cells, 64 cols per row.
*   Cell centers: lat -87.5..87.5 (step 5), lon 2.5..357.5 (step 5).
*
*   Columns per row (after lat, lon):
*     p[5]  : pressure (Pa)           [/100 -> hPa]
*     T[5]  : temperature (K)
*     Q[5]  : specific humidity (g/kg) [/1000 -> kg/kg]
*     dT[5] : lapse rate (K/km)        [/1000 -> K/m]
*     u     : geoid undulation (m)
*     Hs    : orthometric height (m)
*     ah[5] : hydrostatic VMF3 coeff *1000  [/1000]
*     aw[5] : wet VMF3 coeff *1000          [/1000]
*     la[5] : water vapour decrease factor
*     Tm[5] : mean temperature of water vapour (K)
*     Gn_h[5], Ge_h[5], Gn_w[5], Ge_w[5]: gradients (ignored)
*
*   All harmonic fields have 5 terms: A0, A1cos, B1sin, A2cos, B2sin
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

/* 5-degree grid — cell centers */
#define GPT3_DLAT  5.0
#define GPT3_DLON  5.0
#define GPT3_NLAT  36      /* -87.5 to 87.5 */
#define GPT3_NLON  72      /* 2.5 to 357.5   */
#define GPT3_NGRID (GPT3_NLAT*GPT3_NLON)
#define GPT3_LAT0  (-87.5) /* first latitude  */
#define GPT3_LON0  (2.5)   /* first longitude */

/* physical constants */
#define GPT3_G0    9.80665    /* standard gravity (m/s^2) */
#define GPT3_MDRY  0.028965   /* molar mass of dry air (kg/mol) */
#define GPT3_RG    8.3143     /* universal gas constant (J/mol/K) */
#define GPT3_RW    461.52     /* specific gas const for water vapour (J/kg/K) */
/* k2' = k2 - k1*Mw/Md (Askne & Nordius 1987) */
#define GPT3_K2D   16.52      /* refractivity coeff k2' (K/hPa) */
#define GPT3_K3    377600.0   /* refractivity coeff k3 (K^2/hPa) */
#define GPT3_GMR   (GPT3_G0*GPT3_MDRY/GPT3_RG) /* g*M/R (K/m) */

/* VMF3 mapping function constants */
#define VMF3_BH    0.0029
#define VMF3_BW    0.00146
#define VMF3_CW    0.04391
#define VMF3_AHT   2.53e-5   /* height correction, Niell 1996 */
#define VMF3_BHT   5.49e-3
#define VMF3_CHT   1.14e-3

/* grid cell — all parameters have 5 harmonic terms */
typedef struct {
    float p[5];    /* pressure (Pa) */
    float T[5];    /* temperature (K) */
    float Q[5];    /* specific humidity (g/kg) */
    float dT[5];   /* lapse rate (K/km) */
    float u;       /* geoid undulation (m) */
    float Hs;      /* orthometric height (m) */
    float ah[5];   /* hydrostatic MF coeff *1000 */
    float aw[5];   /* wet MF coeff *1000 */
    float la[5];   /* water vapour decrease factor */
    float Tm[5];   /* mean temp of water vapour (K) */
} gpt3cell_t;

static gpt3cell_t gpt3grid[GPT3_NGRID];
static int gpt3_loaded = 0;

/* continued fraction mapping function */
static double mfcf(double el, double a, double b, double c)
{
    double se = sin(el);
    return (1.0 + a/(1.0 + b/(1.0 + c))) / (se + a/(se + b/(se + c)));
}

/* 5-term harmonic: A0 + A1*cos(t) + B1*sin(t) + A2*cos(2t) + B2*sin(2t) */
static double hval5(const float *c, double t)
{
    return (double)c[0] + (double)c[1]*cos(t)     + (double)c[2]*sin(t)
                        + (double)c[3]*cos(2.0*t)  + (double)c[4]*sin(2.0*t);
}

/* bilinear interpolation of n-element float array */
static void interpNf(float *out, int n, const float *a, const float *b,
                     const float *c, const float *d, const double w[4])
{
    int i;
    for (i = 0; i < n; i++)
        out[i] = (float)(w[0]*a[i] + w[1]*b[i] + w[2]*c[i] + w[3]*d[i]);
}

static float interp1f(float a, float b, float c, float d, const double w[4])
{
    return (float)(w[0]*a + w[1]*b + w[2]*c + w[3]*d);
}

/* latitude-dependent c_h coefficient for VMF3 */
static double vmf3_ch(double lat_rad)
{
    double phi = fabs(lat_rad) * R2D / 90.0;
    return 0.062766 + phi*(0.001202 + phi*0.005730);
}

/*------------------------------------------------------------------------------
* read GPT3 5-degree grid file (gpt3_5.grd)
* return : 1:ok, 0:error
*-----------------------------------------------------------------------------*/
extern int gpt3_read(const char *file)
{
    FILE *fp;
    char buf[4096];
    int nread, ilat, ilon, idx;
    double lat, lon;
    gpt3cell_t cell;
    float ign[20];

    trace(3, "gpt3_read: file=%s\n", file);

    if (!file || !*file) return 0;
    if (!(fp = fopen(file, "r"))) {
        trace(1, "gpt3_read: cannot open %s\n", file);
        return 0;
    }
    memset(gpt3grid, 0, sizeof(gpt3grid));
    nread = 0;

    while (fgets(buf, sizeof(buf), fp)) {
        int n;
        if (buf[0] == '%' || buf[0] == '#') continue;

        n = sscanf(buf,
            "%lf %lf "
            "%f %f %f %f %f "    /* p[5]  (Pa)    */
            "%f %f %f %f %f "    /* T[5]  (K)     */
            "%f %f %f %f %f "    /* Q[5]  (g/kg)  */
            "%f %f %f %f %f "    /* dT[5] (K/km)  */
            "%f %f "             /* u, Hs (m)     */
            "%f %f %f %f %f "    /* ah[5] (*1000) */
            "%f %f %f %f %f "    /* aw[5] (*1000) */
            "%f %f %f %f %f "    /* la[5]         */
            "%f %f %f %f %f "    /* Tm[5] (K)     */
            "%f %f %f %f %f "    /* Gn_h (ignored)*/
            "%f %f %f %f %f "    /* Ge_h (ignored)*/
            "%f %f %f %f %f "    /* Gn_w (ignored)*/
            "%f %f %f %f %f",    /* Ge_w (ignored)*/
            &lat, &lon,
            &cell.p[0],  &cell.p[1],  &cell.p[2],  &cell.p[3],  &cell.p[4],
            &cell.T[0],  &cell.T[1],  &cell.T[2],  &cell.T[3],  &cell.T[4],
            &cell.Q[0],  &cell.Q[1],  &cell.Q[2],  &cell.Q[3],  &cell.Q[4],
            &cell.dT[0], &cell.dT[1], &cell.dT[2], &cell.dT[3], &cell.dT[4],
            &cell.u, &cell.Hs,
            &cell.ah[0], &cell.ah[1], &cell.ah[2], &cell.ah[3], &cell.ah[4],
            &cell.aw[0], &cell.aw[1], &cell.aw[2], &cell.aw[3], &cell.aw[4],
            &cell.la[0], &cell.la[1], &cell.la[2], &cell.la[3], &cell.la[4],
            &cell.Tm[0], &cell.Tm[1], &cell.Tm[2], &cell.Tm[3], &cell.Tm[4],
            &ign[0],  &ign[1],  &ign[2],  &ign[3],  &ign[4],
            &ign[5],  &ign[6],  &ign[7],  &ign[8],  &ign[9],
            &ign[10], &ign[11], &ign[12], &ign[13], &ign[14],
            &ign[15], &ign[16], &ign[17], &ign[18], &ign[19]);

        if (n < 44) continue;   /* need at least lat+lon+p+T+Q+dT+u+Hs+ah+aw+la+Tm */

        /* cell-center indexing: lat in [-87.5, 87.5], lon in [2.5, 357.5] */
        ilat = (int)round((lat - GPT3_LAT0) / GPT3_DLAT);
        ilon = ((int)round((fmod(lon + 360.0, 360.0) - GPT3_LON0) / GPT3_DLON)
                + GPT3_NLON) % GPT3_NLON;

        if (ilat < 0 || ilat >= GPT3_NLAT) continue;

        idx = ilat * GPT3_NLON + ilon;
        gpt3grid[idx] = cell;
        nread++;
    }
    fclose(fp);
    gpt3_loaded = nread > 0 ? 1 : 0;
    trace(2, "gpt3_read: %d/%d grid points loaded\n", nread, GPT3_NGRID);
    return gpt3_loaded;
}

/*------------------------------------------------------------------------------
* GPT3 tropospheric model
* args   : gtime_t time   I   GPST
*          double *pos    I   {lat(rad), lon(rad), h(m ellipsoidal)}
*          double *pres   O   pressure (hPa)
*          double *temp   O   temperature (K)
*          double *e      O   water vapour pressure (hPa)
*          double *ah     O   hydrostatic VMF3 coeff
*          double *aw     O   wet VMF3 coeff
*          double *zhd    O   zenith hydrostatic delay (m)
*          double *zwd    O   zenith wet delay (m)
* return : 1:ok, 0:error
*-----------------------------------------------------------------------------*/
extern int gpt3(gtime_t time, const double *pos,
                double *pres, double *temp, double *e,
                double *ah, double *aw, double *zhd, double *zwd)
{
    gpt3cell_t c;
    double lat, lon_d, hell, doy, t;
    double fi, fj, dx, dy, w[4];
    int i0, i1, j0, j1, idx4[4];
    double p0, T0, Q0, dT0, u0, Hs0, ah0, aw0, la0, Tm0;
    double dh, T1, p1, e0, e1, gm;

    if (!gpt3_loaded) {
        trace(2, "gpt3: grid not loaded\n");
        return 0;
    }
    lat  = pos[0];
    hell = pos[2];

    if (hell!=hell || hell<-1000.0 || hell>9000.0) return 0;
    if (lat!=lat   || lat<-PI/2.0  || lat>PI/2.0)  return 0;

    lon_d = pos[1] * R2D;
    lon_d = fmod(lon_d + 360.0, 360.0); /* 0..360 */

    /* day-of-year harmonic argument */
    doy = time2doy(time);
    t   = 2.0 * PI * doy / 365.25;

    /* fractional cell-center indices */
    fi = (pos[0] * R2D - GPT3_LAT0) / GPT3_DLAT;  /* 0..35 */
    fj = (lon_d       - GPT3_LON0)  / GPT3_DLON;  /* 0..71 */

    i0 = (int)floor(fi);
    i1 = i0 + 1;
    if (i0 < 0)          { i0 = 0;             i1 = 1; }
    if (i1 >= GPT3_NLAT) { i1 = GPT3_NLAT - 1; i0 = i1 - 1; }

    j0 = ((int)floor(fj) % GPT3_NLON + GPT3_NLON) % GPT3_NLON;
    j1 = (j0 + 1) % GPT3_NLON;

    dy = fi - floor(fi);
    dx = fj - floor(fj);
    if (dx < 0.0) dx += 1.0;

    w[0] = (1.0-dy)*(1.0-dx);
    w[1] = (1.0-dy)*dx;
    w[2] = dy*(1.0-dx);
    w[3] = dy*dx;

    idx4[0] = i0*GPT3_NLON + j0;
    idx4[1] = i0*GPT3_NLON + j1;
    idx4[2] = i1*GPT3_NLON + j0;
    idx4[3] = i1*GPT3_NLON + j1;

    interpNf(c.p,  5, gpt3grid[idx4[0]].p,  gpt3grid[idx4[1]].p,
                      gpt3grid[idx4[2]].p,  gpt3grid[idx4[3]].p,  w);
    interpNf(c.T,  5, gpt3grid[idx4[0]].T,  gpt3grid[idx4[1]].T,
                      gpt3grid[idx4[2]].T,  gpt3grid[idx4[3]].T,  w);
    interpNf(c.Q,  5, gpt3grid[idx4[0]].Q,  gpt3grid[idx4[1]].Q,
                      gpt3grid[idx4[2]].Q,  gpt3grid[idx4[3]].Q,  w);
    interpNf(c.dT, 5, gpt3grid[idx4[0]].dT, gpt3grid[idx4[1]].dT,
                      gpt3grid[idx4[2]].dT, gpt3grid[idx4[3]].dT, w);
    interpNf(c.ah, 5, gpt3grid[idx4[0]].ah, gpt3grid[idx4[1]].ah,
                      gpt3grid[idx4[2]].ah, gpt3grid[idx4[3]].ah, w);
    interpNf(c.aw, 5, gpt3grid[idx4[0]].aw, gpt3grid[idx4[1]].aw,
                      gpt3grid[idx4[2]].aw, gpt3grid[idx4[3]].aw, w);
    interpNf(c.la, 5, gpt3grid[idx4[0]].la, gpt3grid[idx4[1]].la,
                      gpt3grid[idx4[2]].la, gpt3grid[idx4[3]].la, w);
    interpNf(c.Tm, 5, gpt3grid[idx4[0]].Tm, gpt3grid[idx4[1]].Tm,
                      gpt3grid[idx4[2]].Tm, gpt3grid[idx4[3]].Tm, w);
    c.u  = interp1f(gpt3grid[idx4[0]].u,  gpt3grid[idx4[1]].u,
                    gpt3grid[idx4[2]].u,   gpt3grid[idx4[3]].u,  w);
    c.Hs = interp1f(gpt3grid[idx4[0]].Hs, gpt3grid[idx4[1]].Hs,
                    gpt3grid[idx4[2]].Hs,  gpt3grid[idx4[3]].Hs, w);

    /* apply unit conversions matching gpt3_5.grd encoding */
    p0  = hval5(c.p,  t) / 100.0;   /* Pa -> hPa */
    T0  = hval5(c.T,  t);            /* K */
    Q0  = hval5(c.Q,  t) / 1000.0;  /* g/kg -> kg/kg */
    dT0 = hval5(c.dT, t) / 1000.0;  /* K/km -> K/m */
    u0  = (double)c.u;
    Hs0 = (double)c.Hs;
    ah0 = hval5(c.ah, t) / 1000.0;  /* stored *1000 */
    aw0 = hval5(c.aw, t) / 1000.0;  /* stored *1000 */
    la0 = hval5(c.la, t);
    Tm0 = hval5(c.Tm, t);            /* K */

    /* basic sanity on interpolated meteo */
    if (T0 < 150.0 || T0 > 350.0 || p0 < 500.0 || p0 > 1100.0 ||
        Tm0 < 150.0 || Tm0 > 350.0) {
        trace(2, "gpt3: bad grid values T0=%.1f p0=%.1f Tm0=%.1f lat=%.1f\n",
              T0, p0, Tm0, lat*R2D);
        return 0;
    }

    /* height difference: station ellipsoidal minus grid ellipsoidal (Hs+u) */
    dh = hell - (Hs0 + u0);

    /* temperature at station height */
    T1 = T0 + dT0 * dh;
    if (T1 < 100.0) T1 = 100.0;

    /* pressure at station height (hypsometric) */
    if (fabs(dT0) > 1.0e-6)
        p1 = p0 * pow(T1/T0, -GPT3_GMR/dT0);
    else
        p1 = p0 * exp(-GPT3_GMR * dh / T0);

    /* water vapour pressure */
    e0 = Q0 * p0 / (0.622 + 0.378*Q0);
    e1 = (p0 > 0.0) ? e0 * pow(p1/p0, la0 + 1.0) : 0.0;
    if (e1 < 0.0) e1 = 0.0;

    /* site gravity */
    gm = 9.784 * (1.0 - 0.00266*cos(2.0*lat) - 0.00028*hell/1000.0);

    /* zenith hydrostatic delay */
    *zhd = 2.2768e-3 * p1 / (1.0 - 0.00266*cos(2.0*lat) - 0.00028*hell/1000.0);

    /* zenith wet delay (Askne & Nordius 1987) */
    if (la0 + 1.0 > 0.1 && gm > 0.0)
        *zwd = 1.0e-6 * (GPT3_K2D + GPT3_K3/Tm0) * GPT3_RW * e1
               / (gm * (la0 + 1.0));
    else
        *zwd = 0.0;
    if (*zwd < 0.0) *zwd = 0.0;

    *pres = p1;
    *temp = T1;
    *e    = e1;
    *ah   = ah0;
    *aw   = aw0;

    if (*zhd!=*zhd || *zhd<0.5 || *zhd>3.5) {
        trace(2,"gpt3: invalid zhd=%.4f p1=%.2f lat=%.1f h=%.1f\n",
              *zhd, p1, lat*R2D, hell);
        return 0;
    }
    if (*zwd!=*zwd || *zwd<0.0 || *zwd>0.8) {
        trace(2,"gpt3: invalid zwd=%.4f e1=%.4f la=%.3f Tm=%.1f\n",
              *zwd, e1, la0, Tm0);
        *zwd = 0.05;
    }
    trace(4,"gpt3: lat=%.2f h=%.1f p=%.2f T=%.1f e=%.3f "
          "zhd=%.4f zwd=%.4f ah=%.6f aw=%.6f\n",
          lat*R2D, hell, p1, T1, e1, *zhd, *zwd, *ah, *aw);
    return 1;
}

/*------------------------------------------------------------------------------
* VMF3 mapping function
* return : hydrostatic mapping function; *mfw = wet MF (if mfw!=NULL)
*-----------------------------------------------------------------------------*/
extern double vmf3(double ah, double aw, double el, double lat, double hgt,
                   double *mfw)
{
    double ch, mfh, mfh_ht;

    if (el <= 0.0) {
        if (mfw) *mfw = 0.0;
        return 0.0;
    }
    ch  = vmf3_ch(lat);
    mfh = mfcf(el, ah, VMF3_BH, ch);

    /* height correction (Niell 1996) */
    mfh_ht = mfcf(el, VMF3_AHT, VMF3_BHT, VMF3_CHT);
    mfh   += (1.0/sin(el) - mfh_ht) * hgt / 1000.0;

    if (mfw) *mfw = mfcf(el, aw, VMF3_BW, VMF3_CW);
    return mfh;
}
