# Matched paper experiments



Target coverage reference: Beta(183, 18); alpha=0.089552239.

hat: Romao final disk, ell=r/d. check: one-fit residual threshold.

bar_L6: batches (3,3,3,3,3,3); bar_L2: batches (9,9) at the default rank 18.

All schedules are fixed before observing data. Intersection alpha is the next float above r/(N+1).

The hat coverage CDF is bounded above by the Beta CDF; equality is not assumed.

The check/bar laws are exact subject to the assumptions in main.tex.

Disk problems can have fewer than d supports: bar uses residual ranking, hat uses label padding.

Fit times include standalone solves; area quadrature is timed separately.

The plotted reference is the theoretical true-coverage Beta CDF. Empirical coverage includes finite-validation noise.

See summary.csv for run-level 95% mean intervals and paired.csv for all pairwise differences.

Failures: 0; inspect runs.csv. Successful-only summaries are conditional.


## Results

| Horizon | Method | Coverage | Area | Fit ms | Query ms | Failures |
|---|---|---:|---:|---:|---:|---:|
| 2 | bar_L2 | 0.91072 | 3.57535e-05 | 19.533 | 0.375 | 0 |
| 2 | bar_L6 | 0.91118 | 2.88685e-05 | 66.167 | 0.813 | 0 |
| 2 | check | 0.91069 | 4.8598e-05 | 2.276 | 0.207 | 0 |
| 2 | hat | 0.92535 | 3.22461e-05 | 66.077 | 0.245 | 0 |
| 30 | bar_L2 | 0.91291 | 5.60698 | 19.628 | 0.375 | 0 |
| 30 | bar_L6 | 0.91107 | 5.35077 | 60.933 | 0.864 | 0 |
| 30 | check | 0.91283 | 5.94011 | 2.508 | 0.232 | 0 |
| 30 | hat | 0.92987 | 5.88262 | 58.674 | 0.272 | 0 |
| 50 | bar_L2 | 0.91191 | 54.7529 | 21.186 | 0.421 | 0 |
| 50 | bar_L6 | 0.90994 | 53.1566 | 63.675 | 0.835 | 0 |
| 50 | check | 0.91293 | 56.1931 | 2.460 | 0.219 | 0 |
| 50 | hat | 0.92638 | 56.6955 | 60.927 | 0.236 | 0 |

Mean intervals use independent runs. Horizons share trajectories; do not pool them as independent runs. Coverage quantiles include test-sample noise.

## Paired differences

Differences are method minus reference, using runs where both succeeded.

| Horizon | Method | Reference | Metric | Difference | 95% interval | Pairs |
|---|---|---|---|---:|---|---:|
| 2 | bar_L6 | bar_L2 | coverage | 0.000456 | [-0.00211073, 0.00302273] | 100 |
| 2 | check | bar_L2 | coverage | -3.4e-05 | [-0.00214455, 0.00207655] | 100 |
| 2 | hat | bar_L2 | coverage | 0.014632 | [0.0115231, 0.0177409] | 100 |
| 2 | check | bar_L6 | coverage | -0.00049 | [-0.00373528, 0.00275528] | 100 |
| 2 | hat | bar_L6 | coverage | 0.014176 | [0.0120788, 0.0162732] | 100 |
| 2 | hat | check | coverage | 0.014666 | [0.0110539, 0.0182781] | 100 |
| 2 | bar_L6 | bar_L2 | area | -6.88494e-06 | [-9.87241e-06, -3.89747e-06] | 100 |
| 2 | check | bar_L2 | area | 1.28446e-05 | [6.20522e-06, 1.94839e-05] | 100 |
| 2 | hat | bar_L2 | area | -3.50732e-06 | [-6.5291e-06, -4.85532e-07] | 100 |
| 2 | check | bar_L6 | area | 1.97295e-05 | [1.0588e-05, 2.8871e-05] | 100 |
| 2 | hat | bar_L6 | area | 3.37762e-06 | [2.78919e-06, 3.96605e-06] | 100 |
| 2 | hat | check | area | -1.63519e-05 | [-2.54539e-05, -7.24989e-06] | 100 |
| 2 | bar_L6 | bar_L2 | fit_seconds | 0.0466342 | [0.0443559, 0.0489125] | 100 |
| 2 | check | bar_L2 | fit_seconds | -0.0172568 | [-0.0181782, -0.0163353] | 100 |
| 2 | hat | bar_L2 | fit_seconds | 0.0465439 | [0.0441863, 0.0489015] | 100 |
| 2 | check | bar_L6 | fit_seconds | -0.063891 | [-0.0664784, -0.0613035] | 100 |
| 2 | hat | bar_L6 | fit_seconds | -9.0313e-05 | [-0.00221176, 0.00203114] | 100 |
| 2 | hat | check | fit_seconds | 0.0638006 | [0.0612139, 0.0663873] | 100 |
| 2 | bar_L6 | bar_L2 | query_seconds | 0.000438102 | [0.00039318, 0.000483024] | 100 |
| 2 | check | bar_L2 | query_seconds | -0.00016815 | [-0.000193668, -0.000142632] | 100 |
| 2 | hat | bar_L2 | query_seconds | -0.000130141 | [-0.000154064, -0.000106218] | 100 |
| 2 | check | bar_L6 | query_seconds | -0.000606252 | [-0.000648586, -0.000563918] | 100 |
| 2 | hat | bar_L6 | query_seconds | -0.000568243 | [-0.00060837, -0.000528116] | 100 |
| 2 | hat | check | query_seconds | 3.8009e-05 | [2.36742e-05, 5.23437e-05] | 100 |
| 30 | bar_L6 | bar_L2 | coverage | -0.001844 | [-0.00424657, 0.000558573] | 100 |
| 30 | check | bar_L2 | coverage | -8.4e-05 | [-0.00150232, 0.00133432] | 100 |
| 30 | hat | bar_L2 | coverage | 0.016956 | [0.0142747, 0.0196373] | 100 |
| 30 | check | bar_L6 | coverage | 0.00176 | [-0.00108446, 0.00460446] | 100 |
| 30 | hat | bar_L6 | coverage | 0.0188 | [0.0164924, 0.0211076] | 100 |
| 30 | hat | check | coverage | 0.01704 | [0.0139874, 0.0200926] | 100 |
| 30 | bar_L6 | bar_L2 | area | -0.256209 | [-0.316655, -0.195762] | 100 |
| 30 | check | bar_L2 | area | 0.333133 | [0.107024, 0.559241] | 100 |
| 30 | hat | bar_L2 | area | 0.275641 | [0.206567, 0.344715] | 100 |
| 30 | check | bar_L6 | area | 0.589341 | [0.351294, 0.827389] | 100 |
| 30 | hat | bar_L6 | area | 0.531849 | [0.475214, 0.588485] | 100 |
| 30 | hat | check | area | -0.0574919 | [-0.281871, 0.166887] | 100 |
| 30 | bar_L6 | bar_L2 | fit_seconds | 0.0413053 | [0.0389519, 0.0436586] | 100 |
| 30 | check | bar_L2 | fit_seconds | -0.0171195 | [-0.0182272, -0.0160118] | 100 |
| 30 | hat | bar_L2 | fit_seconds | 0.0390462 | [0.0367009, 0.0413915] | 100 |
| 30 | check | bar_L6 | fit_seconds | -0.0584248 | [-0.0613159, -0.0555336] | 100 |
| 30 | hat | bar_L6 | fit_seconds | -0.00225904 | [-0.00504501, 0.000526921] | 100 |
| 30 | hat | check | fit_seconds | 0.0561657 | [0.0535873, 0.0587442] | 100 |
| 30 | bar_L6 | bar_L2 | query_seconds | 0.000489027 | [0.000428098, 0.000549956] | 100 |
| 30 | check | bar_L2 | query_seconds | -0.00014212 | [-0.000170779, -0.000113461] | 100 |
| 30 | hat | bar_L2 | query_seconds | -0.000102421 | [-0.000127554, -7.72885e-05] | 100 |
| 30 | check | bar_L6 | query_seconds | -0.000631147 | [-0.000695205, -0.000567089] | 100 |
| 30 | hat | bar_L6 | query_seconds | -0.000591448 | [-0.000654176, -0.00052872] | 100 |
| 30 | hat | check | query_seconds | 3.9699e-05 | [1.57605e-05, 6.36375e-05] | 100 |
| 50 | bar_L6 | bar_L2 | coverage | -0.00197 | [-0.00448666, 0.000546665] | 100 |
| 50 | check | bar_L2 | coverage | 0.001024 | [-0.000498652, 0.00254665] | 100 |
| 50 | hat | bar_L2 | coverage | 0.01447 | [0.011827, 0.017113] | 100 |
| 50 | check | bar_L6 | coverage | 0.002994 | [0.000229975, 0.00575803] | 100 |
| 50 | hat | bar_L6 | coverage | 0.01644 | [0.0139867, 0.0188933] | 100 |
| 50 | hat | check | coverage | 0.013446 | [0.0105589, 0.0163331] | 100 |
| 50 | bar_L6 | bar_L2 | area | -1.5963 | [-2.01382, -1.17879] | 100 |
| 50 | check | bar_L2 | area | 1.44024 | [1.01481, 1.86567] | 100 |
| 50 | hat | bar_L2 | area | 1.94261 | [1.51641, 2.3688] | 100 |
| 50 | check | bar_L6 | area | 3.03654 | [2.44876, 3.62433] | 100 |
| 50 | hat | bar_L6 | area | 3.53891 | [3.09846, 3.97935] | 100 |
| 50 | hat | check | area | 0.502366 | [-0.0385387, 1.04327] | 100 |
| 50 | bar_L6 | bar_L2 | fit_seconds | 0.042489 | [0.0400945, 0.0448835] | 100 |
| 50 | check | bar_L2 | fit_seconds | -0.018726 | [-0.0201058, -0.0173462] | 100 |
| 50 | hat | bar_L2 | fit_seconds | 0.0397407 | [0.0375998, 0.0418816] | 100 |
| 50 | check | bar_L6 | fit_seconds | -0.0612149 | [-0.0638052, -0.0586247] | 100 |
| 50 | hat | bar_L6 | fit_seconds | -0.00274831 | [-0.00472202, -0.000774599] | 100 |
| 50 | hat | check | fit_seconds | 0.0584666 | [0.0561005, 0.0608328] | 100 |
| 50 | bar_L6 | bar_L2 | query_seconds | 0.0004133 | [0.000354118, 0.000472482] | 100 |
| 50 | check | bar_L2 | query_seconds | -0.000201974 | [-0.000255819, -0.000148129] | 100 |
| 50 | hat | bar_L2 | query_seconds | -0.000185475 | [-0.000239387, -0.000131562] | 100 |
| 50 | check | bar_L6 | query_seconds | -0.000615274 | [-0.000665818, -0.00056473] | 100 |
| 50 | hat | bar_L6 | query_seconds | -0.000598775 | [-0.000644558, -0.000552992] | 100 |
| 50 | hat | check | query_seconds | 1.6499e-05 | [-2.05183e-06, 3.50498e-05] | 100 |
