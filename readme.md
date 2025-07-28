# RTKLIB 2.4.3 d3

Added the ability to process multi-system GNSS measurements in PPP, taking into account differential code delays. L1/L5 combination is supported.

*For post-processing use BSX file by CAS (https://cddis.nasa.gov/Data_and_Derived_Products/GNSS/gnss_differential_code_bias_product.html)*

*For real-time you can use any SSR-stream with DCB or use SSRA00CNE1 from IGS*

The algorithm of variational Bayesian-based robust adaptive Kalman filter is implemented.
In config file you may chouse base algorithm or vbbra (pos1-kalman).
