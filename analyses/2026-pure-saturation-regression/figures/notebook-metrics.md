Average absolute relative deviation (AARD), computed from `results/propane-replay-2026-10-01.csv`; unavailable predictions would stop presentation, rather than disappear from the mean.

| Role     | Model   |   States |   Minimum T / K |   Maximum T / K |   Pressure AARD / % |   Liquid density AARD / % |
|:---------|:--------|---------:|----------------:|----------------:|--------------------:|--------------------------:|
| stress   | fitted  |        2 |             110 |             340 |          22.0196    |                  0.381817 |
| stress   | packet  |        2 |             110 |             340 |          22.1228    |                  0.442421 |
| held_out | fitted  |       19 |             120 |             320 |           1.00068   |                  0.179524 |
| held_out | packet  |       19 |             120 |             320 |           1.07295   |                  0.228264 |
| training | fitted  |        3 |             210 |             330 |           0.0138239 |                  0.134111 |
| training | packet  |        3 |             210 |             330 |           0.0669128 |                  0.260725 |
