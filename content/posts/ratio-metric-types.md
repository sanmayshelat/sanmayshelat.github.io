---
title: "Comparing Ratio Metric Types"
date: 2023-12-04
description: "Comparing ratio-of-averages versus average-of-ratios"
tags: [experimentation]
---

# Ratio-of-averages v. Average-of-ratios
To calculate the variance of ratio-of-averages we have to use the [Delta method](/posts/cuped-derivation/#1-delta-method-formulation) while the variance of average-of-ratios can be estimated using the sample variance.

Let $i$ be the independent randomization unit. Let there be $K$ clusters.

Then in each unit, $Y_{j}$ are independent such that, $\mathbb{E}(Y_{j})=\mu_i$ and $var(Y_{j})=\sigma_{i}^2$. Thus, the between-cluster variance, $\tau^2=var(\mu_i)$. Further assume that $N_i=m, \sigma_{i}^2=\sigma^2; \forall{i}$.

Based on the [Variance property](/posts/cuped-derivation/#33-cuped-on-ratio-metrics) for the Delta method, 

$$
\begin{align*}
\text{var}\left(\sum_{i}{a_{i}X_{i}}\right) &= \sum_{i}{a_{i}^{2}\text{var}(x_{i})} + 2\sum_{i\lt{j}\lt{N}}{a_{i}a_{j}\text{cov}(x_{i},x_{j})}\\\
\text{var}\left(\frac{\sum_{j=1}^{N}{Y_{ij}}}{N}\right) &= \frac{1}{N^2} \text{var}\left(\sum_{j=1}^{N}{Y_{ij}}\right) \\\
&= \frac{1}{N^2} \left[ \sum_{j=1}^{N}{\text{var}(Y_{ij})} + 2\sum_{p\lt{q}\lt{\frac{N}{2}}}{\text{cov}(Y_{ip},Y_{iq})} \right]\\\
&= \frac{1}{N^2} \left[ \sum_{j=1}^{N}{(\sigma_b + \sigma_w)} + \sum_{p=1}^{N}\sum_{q=1,q\neq p}^{N}{\sigma_b} \right]\\\
&= \frac{1}{N} \left[ (\sigma_b + \sigma_w) + (N-1)\sigma_b \right]\\\
&= \frac{(\sigma_b + \sigma_w)}{N} \left[ 1 + (N-1)\frac{\sigma_b}{(\sigma_b + \sigma_w)} \right]\\\
&= \frac{(\sigma_b + \sigma_w)}{N} \left[ 1 + (N-1)\rho \right]\\\
\text{var}\left(\frac{\sum_{i=1}^{K}\sum_{j=1}^{N}{Y_{ij}}}{KN}\right) &= \frac{\text{var}\left(\frac{\sum_{j=1}^{N}{Y_{ij}}}{N}\right)}{K}\\\
&= \frac{(\sigma_b + \sigma_w)}{KN} \left[ 1 + (N-1)\rho \right]\\\
\end{align*}
$$