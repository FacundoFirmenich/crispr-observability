# MP27 source-allele endpoint by absolute error budget

All counts below use padding 5. Source-allele candidate sets are verified equal across paddings. For HG002 these are concordances with conditional genotype labels, not biological accuracy. Unknown genotype records remain in the raw denominator.

| Budget | Held-out allele correct | Held-out allele discordant | Held-out allele U | HG002 concordant | HG002 discordant | HG002 U (label available) | HG002 label unavailable |
|---|---:|---:|---:|---:|---:|---:|---:|
| 0 | 2880 | 1 | 7919 | 8462 | 21 | 6470 | 1863 |
| 1 | 4359 | 3 | 6438 | 9570 | 27 | 5356 | 1863 |
| 2 | 5330 | 8 | 5462 | 7826 | 29 | 7098 | 1863 |
| 3 | 5574 | 3 | 5223 | 6823 | 35 | 8095 | 1863 |
| 5 | 4580 | 0 | 6220 | 5014 | 25 | 9914 | 1863 |
| 8 | 2906 | 0 | 7894 | 3019 | 12 | 11922 | 1863 |
| 12 | 1451 | 0 | 9349 | 1616 | 8 | 13329 | 1863 |

Held-out rows include histories outside each budget. All observed discordances in the simulated panel occur outside the declared literal-history budget. This does not make those failures irrelevant: the literal error count is unknown in real data. No budget is selected as a calibrated operating point.
