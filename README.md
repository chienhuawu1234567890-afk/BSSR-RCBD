#SSR-RBD: Simulation R Codes for Blinded Sample Size Re-Estimation in RCBD

This repository contains R simulation codes for the manuscript:

**"An Unbiased Blinded Approach to Sample Size Re-Estimation in Randomized Complete Block Designs"**

## Overview

The R codes implement simulation studies for blinded sample size re-estimation (SSR) in randomized complete block designs (RCBDs). The simulations compare three variance estimation approaches:

1. Naïve blinded variance estimator
2. Adjusted blinded variance estimator
3. Unblinded benchmark estimator

The goal is to evaluate how different variance estimators affect re-estimated sample sizes under balanced and unbalanced RCBD settings.

## Simulation Settings

The simulation studies cover the following RCBD frameworks:

- Balanced RCBD with fixed block effects
- Balanced RCBD with random block effects
- Unbalanced RCBD with fixed block effects
- Unbalanced RCBD with random block effects

The simulations consider different values of:

- Interim information fraction
- Residual variance
- Block variance
- Treatment effect structure
- Allocation balance or imbalance

## Main Methods

The codes implement:

- Exact sample size planning using noncentral F-distribution inversion
- Interim blinded variance estimation
- Adjusted blinded variance estimation using planned structured variation
- Unblinded residual variance estimation
- Re-estimated sample size calculation
- Monte Carlo simulation summaries

## Output

The simulation outputs include summary tables of re-estimated sample sizes, including:

- Mean
- Standard deviation
- Median
- 2.5% quantile
- 97.5% quantile

These summaries are used to compare the performance of naïve blinded, adjusted blinded, and unblinded SSR procedures.

## Repository Structure

```text
SSR-RBD/
├── README.md
├── balanced_fixed_rcbd.R
├── balanced_random_rcbd.R
├── unbalanced_fixed_rcbd.R
├── unbalanced_random_rcbd.R
├── functions/
│   ├── sample_size_planning.R
│   ├── variance_estimators.R
│   └── simulation_summary.R
└── results/
    ├── tables/
    └── figures/
