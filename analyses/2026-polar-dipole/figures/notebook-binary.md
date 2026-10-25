Paired Campbell comparisons from `results/binary-xy-summary.csv`. Pressure deviations are %, composition deviations are mole fraction; validation means withheld from this study's binary interaction fit.

| model                 | variant                                        |                   k_ij | role        |   points |   failures |   mean_abs_dy |   max_abs_dy |   AAD_P_percent |
|:----------------------|:-----------------------------------------------|-----------------------:|:------------|---------:|-----------:|--------------:|-------------:|----------------:|
| polar/GV2006          | k_ij=0                                         |              0         | calibration |        8 |          0 |    0.0363908  |    0.0787182 |      12.6815    |
| polar/GV2006          | k_ij=0                                         |              0         | validation  |       19 |          0 |    0.0283152  |    0.0603453 |       9.83968   |
| polar/GV2006          | k_ij calibrated                                |              0.0333387 | calibration |        8 |          0 |    0.00976427 |    0.0288948 |       1.79945   |
| polar/GV2006          | k_ij calibrated                                |              0.0333387 | validation  |       19 |          0 |    0.0101825  |    0.0307724 |       1.65381   |
| nonpolar/refit        | k_ij=0                                         |              0         | calibration |        8 |          0 |    0.0692881  |    0.127348  |      20.9282    |
| nonpolar/refit        | k_ij=0                                         |              0         | validation  |       19 |          0 |    0.0552071  |    0.109086  |      16.9901    |
| nonpolar/refit        | k_ij calibrated                                |              0.0592    | calibration |        8 |          0 |    0.0170833  |    0.026743  |       3.34642   |
| nonpolar/refit        | k_ij calibrated                                |              0.0592    | validation  |       19 |          0 |    0.0159381  |    0.038246  |       2.85063   |
| Margules/Campbell1986 | four-suffix eq 8 (Barker P-x fit) per isotherm | not applicable         | calibration |        8 |          0 |    0.008375   |    0.029     |       0.0896472 |
| Margules/Campbell1986 | four-suffix eq 8 (Barker P-x fit) per isotherm | not applicable         | validation  |       19 |          0 |    0.00678947 |    0.016     |       0.192155  |
