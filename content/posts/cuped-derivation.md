---
title: "CUPED Derived for Ratio Metrics"
date: 2023-12-20
description: "Derivation of CUPED for ratio-type metrics with each step explained—from Delta method to final formulation."
tags: [experimentation, causal]
---

# Table of Contents <!-- omit in toc -->
- [1. Introduction](#1-introduction)
- [2. Delta method formulation](#2-delta-method-formulation)
  - [2.1. Univariate](#21-univariate)
  - [2.2. Multivariate](#22-multivariate)
    - [2.2.1. Sidebar: Variance property](#221-sidebar-variance-property)
- [3. Delta method for ratio metrics](#3-delta-method-for-ratio-metrics)
- [4. CUPED](#4-cuped)
  - [4.1. What is CUPED?](#41-what-is-cuped)
  - [4.2. Why does it work?](#42-why-does-it-work)
  - [4.3. CUPED on ratio metrics](#43-cuped-on-ratio-metrics)
- [5. Sources](#5-sources)
  - [5.1. Delta method](#51-delta-method)
  - [5.2. CUPED](#52-cuped)


# 1. Introduction
CUPED—controlled-experiment using pre-experiment data—is a technique to reduce variance in A/B tests by using information about the test metric from a pre-experimental period. While the method itself is well-documented in the [original paper](https://www.exp-platform.com/Documents/2013-02-CUPED-ImprovingSensitivityOfControlledExperiments.pdf) and blog articles based on it (e.g., [1](https://www.statsig.com/blog/cuped), [2](https://matteocourthoud.github.io/post/cuped/)), I missed a complete derivation for ratio metrics. In this article, I derive this formulation. For relevant Python code please see `ate.py` in my [CUPED simulation project](https://github.com/sanmayshelat/cuped).



# 2. Delta method formulation
We start with the Delta method which is critical for obtaining the variance of ratio metrics. To avoid confusion, note that, here, I refer to metrics that are ratios-of-averages as ratio metrics; for example, (in a website) click-through rate defined as the average number of clicks divided by the average number of page-views or (in a ride-hailing marketplace) the driver acceptance rate defined as the average number of offers accepted divided by the average number of offers received.

The Delta method allows us to derive the asymptotic (i.e., for a large number of observations) distribution of a function on a random variable given that the underlying random variable is asymptotically normal. The later is typically true for many metrics of interest as by the central limit theorem, sample averages _are_ asymptotically normal. First we'll give the proof in the univariate case and then do the same for multivariate.

## 2.1. Univariate
If $X_{n}$ is a sequence of averages where $n\to\infty$:

$$
\begin{aligned}
X_{n} &\sim \mathcal{N}(\mu, \frac{\sigma^{2}}{n}) \\\
X_{n} - \mu &\sim \mathcal{N}(0, \frac{\sigma^{2}}{n}) \\\
X_{n} - \mu &\sim \frac{\sigma}{\sqrt{n}}\mathcal{N}(0, 1) \\\
\frac{\sqrt{n}(X_{n} - \mu)}{\sigma} &\sim \mathcal{N}(0, 1) \\\
\end{aligned}
$$


The Delta method then gives that for any continuous transformation $g(x)$:

$$
\begin{align*}
\frac{\sqrt{n}(g(X_{n})-g(\mu))}{|g'(\mu)|\sigma} &\sim \mathcal{N}(0, 1) \\\
g(X_{n}) &\sim \mathcal{N} \left( g(\mu), \frac{g'(\mu)^{2}\sigma^{2}}{n} \right) \\\
\end{align*}
$$


The above applies because from the Taylor series:

$$
\begin{align*}
g(X_{n}) &= g(\mu) + \frac{g'(\mu)}{1!}(X_{n}-\mu) + \frac{g''(\mu)}{2!}(X_{n}-\mu)^{2} + \dots\\\
\frac{g(X_{n})-g(\mu)}{g'(\mu)} &\approx (X_{n}-\mu)\\\
\frac{\sqrt{n}(g(X_{n})-g(\mu))}{g'(\mu)\sigma} &\approx \frac{\sqrt{n}}{\sigma}(X_{n}-\mu) \sim {\sf N}(0, 1)\\\
\end{align*}
$$



## 2.2. Multivariate
Similarly, for the multivariate case:

$$
\frac{\sqrt{n}(B-\beta)}{\Sigma} \sim \mathcal{N}(0,1)
$$
where:
* $B$: vector of $k$ random variables
* $\beta$: vector of $k$ constants
* $\Sigma$: covariance matrix of $B$

For any continuous function, $h(B)$, first order Taylor series implies:

$$
h(B) \approx h(\beta) + \nabla h(\beta)^{T}(B-\beta)
$$

where:
* $\nabla h(B)$: gradient; i.e., vector of partial differentiation of $h(B)$ by each element of $B$

The variance of both sides should also be equal:

$$
\text{var}(h(B)) = \text{var}(\nabla h(\beta)^{T}B)
$$

$\because$ all other terms on the _rhs_ are constants. Then by variance property:

$$
\begin{align*}
\text{var}(h(B)) &= \nabla h(\beta)^{T}\text{var}(B)\nabla h(\beta) \\\
&= \nabla h(\beta)^{T} \Sigma \nabla h(\beta) \\\
\end{align*}
$$


### 2.2.1. Sidebar: Variance property
Variance when multiplying scalar with random variable.

$$
\begin{align*}
\text{var}\left(\sum_{i}a_{i}X_{i}\right)
&= \sum_{i,j}a_{i}a_{j}\text{cov}(x_{i},x_{j})\\\
&= \sum_{i}a_{i}^{2}\text{var}(x_{i}) + \sum_{i,j, i \neq j}a_{i}a_{j}\text{cov}(x_{i},x_{j})\\\
&= \sum_{i}a_{i}^{2}\text{var}(x_{i}) + 2\sum_{i<j<N}a_{i}a_{j}\text{cov}(x_{i},x_{j})\\\
\end{align*}
$$

Equivalently, in matrix notation:

$$
\text{var}(a^{T}X) = a^{T} \Sigma a
$$




# 3. Delta method for ratio metrics
Now that we can use the Delta method for the multivariate case, let's calculate the treatment effect for a ratio (ratio-of-averages) metric with the Delta method. In the derivation below, $\bar{Y}$ and $\bar{N}$ are, for instance, average accepted offers and received offers for drivers in a ride-hailing marketplace.

* Random variables: $B = [\bar{Y}, \bar{N}]$ (thus, $k=2$)
* Let mean vector: $\beta = [\mu_{Y}, \mu_{N}]$
* Covariance: $\Sigma$

Averages of randomization unit level quantities (which are $\therefore$ i.i.d):
* $\bar{Y} = \frac{\sum_{i}{Y_{i+}}}{n}$
* $\bar{N} = \frac{\sum_{i}{N_{i+}}}{n}$

Let

$$
\begin{align*}
h(B) = h(\bar{Y}, \bar{N}) &= \frac{\bar{Y}}{\bar{N}} \\\\\\
\therefore \text{var}(h(B)) &= \nabla h(\beta)^{T} \Sigma \nabla h(\beta) \\\
\end{align*}
$$

where

$$
\begin{align*}
\nabla h(\beta) &= \left[ 
    \frac{\partial{\frac{\bar{Y}}{\bar{N}}}}{\partial{\bar{Y}}},
    \frac{\partial{\frac{\bar{Y}}{\bar{N}}}}{\partial{\bar{N}}} 
    \right]\\\
    &= \left[ 
    \frac{1}{\mu_{N}},
    \frac{-\mu_{Y}}{\mu_{N}^{2}} 
    \right]\\\
\end{align*}
$$

This way calculate $h(B)$ and $var(h(B))$ for the treatment and control sets. The difference of the former gives the mean effect (_delta_) and the latter can be combined to get the pooled variance. The delta and pooled variance can be combined to calculate  the z-statistic.

# 4. CUPED
In this section, we first discuss what CUPED is and why it works, deriving its formulation for count type metrics. Then we apply the Delta method to obtain the matrix formulation that can be used for count metrics.

## 4.1. What is CUPED?
Control variates is a variance reduction method commonly applied in Monte Carlo procedures to reduce the error of estimates. Let's say we are interested in estimating $\mu_{Y}=\mathbb{E}(\bar{Y})$. A typical estimator would be the average of sample observations, $\bar{Y} = \frac{\sum_{i}{Y_{i+}}}{n}$.

Let $\bar{X}$ be some other metric where $\mathbb{E}(\bar{X})$ is known. Then $Y_{cv}=\bar{Y} + \theta(\bar{X} - \mathbb{E}(\bar{X}))$, where $\theta$ is a fixed value, is also an unbiased estimator of $\mu_{Y}$. To prove this, take the expectation on both sides:

$$
\begin{align*}
\mathbb{E}(Y_{cv}) &= \mathbb{E}(\bar{Y} + \theta(\bar{X} - \mathbb{E}(\bar{X}))) \\\
&= \mathbb{E}(\bar{Y}) + \theta\mathbb{E}(\bar{X}) - \theta\mathbb{E}(\bar{X}) \\\
&= \mathbb{E}(\bar{Y})
\end{align*}
$$

The variance of this estimator is:

$$
\begin{align*}
\text{var}(Y_{cv}) &= \text{var}(\bar{Y} + \theta(\bar{X} - \mathbb{E}(X))) \\\
&= \text{var}(\bar{Y} + \theta\bar{X}) \\\
&= \text{var}(\bar{Y}) + \text{var}(\theta\bar{X}) + \text{cov}(\bar{Y},\theta\bar{X}) \\\
&= \text{var}(\bar{Y}) + \theta^{2} \text{var}(\bar{X}) + 2\theta \text{cov}(\bar{Y},\bar{X})
\end{align*}
$$


We would like to minimize the variance of any estimator. As $\theta$ is a free parameter we can optimize to minimize the variance:

$$
\begin{align*}
\frac{\mathrm{d}\left[\text{var}(Y_{cv})\right]}{\mathrm{d}\theta} &= 0 \\\
\frac{\mathrm{d}\left[\text{var}(\bar{Y}) + \theta^{2} \text{var}(\bar{X}) + 2\theta \text{cov}(\bar{Y},\bar{X})\right]}{\mathrm{d}\theta} &= 0 \\\
2\theta \text{var}(\bar{X}) + 2\text{cov}(\bar{Y},\bar{X}) &= 0 \\\
\theta &= - \frac{\text{cov}(\bar{Y},\bar{X})}{\text{var}(\bar{X})}
\end{align*}
$$

Replacing theta in the variance formula with this optimal value:

$$
\begin{align*}
\text{var}(Y_{cv}) &= \text{var}(\bar{Y}) + 
\frac{\text{cov}(\bar{Y},\bar{X})^{2}}{\text{var}(\bar{X})}\\\
&- 2\frac{\text{cov}(\bar{Y},\bar{X})^{2}}{\text{var}(\bar{X})}\\\
&= \text{var}(\bar{Y}) - \frac{\text{cov}(\bar{Y},\bar{X})^{2}}{\text{var}(\bar{X})}\\\
&= \text{var}(\bar{Y}) \left[1-\frac{\text{cov}(\bar{Y},\bar{X})^{2}}{\text{var}(\bar{X})\text{var}(\bar{Y})} \right]\\\
&= \text{var}(\bar{Y}) (1-\rho_{\bar{Y},\bar{X}}^{2})
\end{align*}
$$

This reduced variance can be used to calculate the pooled variance of the treatment and control groups and used with the average treatment effect to obtain the test statistic.

Points of note:
* We don't know the terms required to calculate $\theta$ exactly: the sample estimates for the terms can be used.
* Which X to use: usually the same metric but pre-experiment will give the most variance reduction because it will have the highest correlation with the experiment metric.


## 4.2. Why does it work?
Typically, we don't have an $X$ for which we know $\mathbb{E}(X)$. If we just use the sample mean to estimate it then we would just get the old metric back (i.e., sample mean of $\bar{Y}$).

We could use some other data than $X_1, X_2, \dots$ to estimate $\mathbb{E}(\bar{X})$ with some variance. But remember that for A/B tests this is not a problem if $\bar{X}$ is some pre-experimental data because, we are typically interested in estimating $\mathbb{E}(\bar{Y_t}) - \mathbb{E}(\bar{Y_c})$ (referring to the metric in treatment and control).

With the new estimator that becomes as shown below; and $\mathbb{E}(\bar{X_c}) - \mathbb{E}(\bar{X_t}) = 0$ because using pre-experimental data, there should be no treatment effect and random assignment would assure that there are no differences.

$$
\begin{align*}
\mathbb{E}(\bar{Y}\_{cv,t}) - \mathbb{E}(\bar{Y}\_{cv,c}) &= 
(\mathbb{E}(\bar{Y}\_{t}) + \theta(\bar{X}\_{t} - \mathbb{E}(\bar{X}\_{t})))\\\ 
&-(\mathbb{E}(\bar{Y}\_{c}) + \theta(\bar{X}\_{c} - \mathbb{E}(\bar{X}\_{c})))\\\
&= (\mathbb{E}(\bar{Y}\_{t}) + \theta\bar{X}\_{t}) - 
(\mathbb{E}(\bar{Y}\_{c}) + \theta\bar{X}\_{c})\\\
&+
(\mathbb{E}(\bar{X}\_{c}) - \mathbb{E}(\bar{X}\_{t}))
\end{align*}
$$


## 4.3. CUPED on ratio metrics
We now combine CUPED and the Delta method. Two additional random variables are going to be included from the pre-experimental data.

* Random variables: $B = [\bar{Y}, \bar{N}, \bar{X}, \bar{M}]$ (thus, $k=4$)
    * Here $\bar{Y}$ and $\bar{X}$ refer to the numerator of the ratio metric in the experimental and pre-experimental data, respectively. Correspondingly, $\bar{N}$ and $\bar{M}$ refer to the denominator.
* Let mean vector: $\beta = [\mu_{Y}, \mu_{N}, \mu_{X}, \mu_{M}]$
* Covariance: $\Sigma$

To apply CUPED:
* Let $s(B)=\frac{\bar{Y}}{\bar{N}}$ and $r(B)=\frac{\bar{X}}{\bar{M}}$ then using the Delta method:
    * $var(s(B)) = \nabla{s(\beta)}^{T} \Sigma \nabla{s(\beta)}$
    * $var(r(B)) = \nabla{r(\beta)}^{T} \Sigma \nabla{r(\beta)}$
    * $covar(s(B), r(B)) = \nabla{s(\beta)}^{T} \Sigma \nabla{r(\beta)}$

* The gradients calculations are thus:

$$
\begin{align*}
\nabla s(\beta) &= \left[ 
    \frac{\partial{\frac{\bar{Y}}{\bar{N}}}}{\partial{\bar{Y}}},
    \frac{\partial{\frac{\bar{Y}}{\bar{N}}}}{\partial{\bar{N}}},
    \frac{\partial{\frac{\bar{Y}}{\bar{N}}}}{\partial{\bar{X}}},
    \frac{\partial{\frac{\bar{Y}}{\bar{N}}}}{\partial{\bar{M}}}
    \right]\\\
    &= \left[ 
    \frac{1}{\mu_{N}},
    \frac{-\mu_{Y}}{\mu_{N}^{2}},
    0,
    0
    \right]\\\
\nabla r(\beta) &= \left[ 
    \frac{\partial{\frac{\bar{X}}{\bar{M}}}}{\partial{\bar{Y}}},
    \frac{\partial{\frac{\bar{X}}{\bar{M}}}}{\partial{\bar{N}}},
    \frac{\partial{\frac{\bar{X}}{\bar{M}}}}{\partial{\bar{X}}},
    \frac{\partial{\frac{\bar{X}}{\bar{M}}}}{\partial{\bar{M}}}
    \right]\\\
    &= \left[
    0,
    0,
    \frac{1}{\mu_{M}},
    \frac{-\mu_{X}}{\mu_{M}^{2}}
    \right]
\end{align*}
$$


* $\theta = -\frac{covar(s(B), r(B))}{var(r(B))}$
* Replacing theta in the variance formula:

$$
\begin{align*}
\text{var} \left( \left[\frac{\bar{Y}}{\bar{N}}\right]\_{cv} \right) &= \text{var} \left( \frac{\bar{Y}}{\bar{N}} \right) + 
\frac{\text{cov}\left( \frac{\bar{Y}}{\bar{N}}, \frac{\bar{X}}{\bar{M}} \right)^{2}}{\text{var} \left( \frac{\bar{X}}{\bar{M}} \right)} - 
2\frac{\text{cov}\left( \frac{\bar{Y}}{\bar{N}}, \frac{\bar{X}}{\bar{M}} \right)^{2}}{\text{var} \left( \frac{\bar{X}}{\bar{M}} \right)}\\\
&= \text{var} \left( \frac{\bar{Y}}{\bar{N}} \right) \left(1-\rho_{\frac{\bar{Y}}{\bar{N}},\frac{\bar{X}}{\bar{M}}}^{2} \right)
\end{align*}
$$

Just as with count metrics, this variance can be calculated for both experiment groups and used to obtain the test statistic.


# 5. Sources
Here are some sources I found useful:

## 5.1. Delta method
* [Applying the Delta Method in Metric Analytics: A Practical Guide with Novel Ideas](https://arxiv.org/pdf/1803.06336)
* [10 Fundamental Theorems for Econometrics](https://bookdown.org/ts_robinson1994/10_fundamental_theorems_for_econometrics/dm.html)

## 5.2. CUPED
* [Improving the Sensitivity of Online Controlled Experiments by Utilizing Pre-Experiment Data](https://doi.org/10.1145/2433396.2433413)
* [Improving the Sensitivity of Online Controlled Experiments: Case Studies at Netflix](https://doi.org/10.1145/2939672.2939733)