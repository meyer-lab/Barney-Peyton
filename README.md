# Barney-Peyton

Code for modeling performed in Barney et al.

[![Build Status](https://travis-ci.org/meyer-lab/Barney-Peyton.svg?branch=master)](https://travis-ci.org/meyer-lab/Barney-Peyton)

You can find the output built through Travis CI [here](https://meyer-lab.github.io/Barney-Peyton/Figures.md.html).








# Methods text

## Multilinear regression modeling

Multilinear regression modeling was performed using ‘stats::lm’ within R, regressing each compendium of signaling measurements against the corresponding proliferation measured in the same conditions. All relationships were assumed to be additive, and no interaction or intercept terms were included in the models. Each measurement was z-score normalized before regression. The two-sided p-values presented were calculated based upon the significance of each parameter being non-zero. The overall performance of each model (as calculated by the F-statistic) corresponded well to the significance of individual terms.
