JM demonstration
================
Nathan Constantine-Cooke

  - [Introduction](#introduction)
  - [Running JM](#running-jm)
  - [Interpretting the output](#interpretting-the-output)
  - [Table of `method` arguments](#table-of-method-arguments)
  - [Competing risks](#competing-risks)

## Introduction

``` r
library(lattice)
library(JM, quietly = TRUE, warn.conflicts = FALSE)
```

`JM` by Dimitris Rizopoulos is one of the most popular packages for
joint modeling of longitudinal and time-to-event data. `JM` includes a
dataset, `aids` which will be used for this demonstration.

``` r
data(aids)
knitr::kable(head(aids, 8))
```

| patient |  Time | death |       CD4 | obstime | drug | gender | prevOI | AZT         | start |  stop | event |
| :------ | ----: | ----: | --------: | ------: | :--- | :----- | :----- | :---------- | ----: | ----: | ----: |
| 1       | 16.97 |     0 | 10.677078 |       0 | ddC  | male   | AIDS   | intolerance |     0 |  6.00 |     0 |
| 1       | 16.97 |     0 |  8.426150 |       6 | ddC  | male   | AIDS   | intolerance |     6 | 12.00 |     0 |
| 1       | 16.97 |     0 |  9.433981 |      12 | ddC  | male   | AIDS   | intolerance |    12 | 16.97 |     0 |
| 2       | 19.00 |     0 |  6.324555 |       0 | ddI  | male   | noAIDS | intolerance |     0 |  6.00 |     0 |
| 2       | 19.00 |     0 |  8.124038 |       6 | ddI  | male   | noAIDS | intolerance |     6 | 12.00 |     0 |
| 2       | 19.00 |     0 |  4.582576 |      12 | ddI  | male   | noAIDS | intolerance |    12 | 18.00 |     0 |
| 2       | 19.00 |     0 |  5.000000 |      18 | ddI  | male   | noAIDS | intolerance |    18 | 19.00 |     0 |
| 3       | 18.53 |     1 |  3.464102 |       0 | ddI  | female | AIDS   | intolerance |     0 |  2.00 |     0 |

## Running JM

JM is easy to use. We simply fit a LME using the `nlme` package (which
is loaded by JM).

``` r
fitLME <- lme(sqrt(CD4) ~ obstime + obstime:drug,
              random = ~ obstime | patient,
              data = aids)
summary(fitLME)
```

    ## Linear mixed-effects model fit by REML
    ##  Data: aids 
    ##        AIC      BIC    logLik
    ##   2699.069 2735.789 -1342.535
    ## 
    ## Random effects:
    ##  Formula: ~obstime | patient
    ##  Structure: General positive-definite, Log-Cholesky parametrization
    ##             StdDev     Corr  
    ## (Intercept) 0.87143264 (Intr)
    ## obstime     0.03617033 -0.015
    ## Residual    0.36844785       
    ## 
    ## Fixed effects: sqrt(CD4) ~ obstime + obstime:drug 
    ##                      Value  Std.Error  DF  t-value p-value
    ## (Intercept)      2.5118005 0.04258901 936 58.97766  0.0000
    ## obstime         -0.0375070 0.00440225 936 -8.51997  0.0000
    ## obstime:drugddI  0.0082141 0.00632277 936  1.29912  0.1942
    ##  Correlation: 
    ##                 (Intr) obstim
    ## obstime         -0.118       
    ## obstime:drugddI  0.000 -0.687
    ## 
    ## Standardized Within-Group Residuals:
    ##           Min            Q1           Med            Q3           Max 
    ## -4.2480426446 -0.4082420034 -0.0002391741  0.4336550881  3.7150583354 
    ## 
    ## Number of Observations: 1405
    ## Number of Groups: 467

And fit a Cox proportional hazards model using the standard `survival`
approach. Note the use of the `aids.id` dataset which uses only baseline
measurements.

``` r
fitSURV <- coxph(Surv(Time, death) ~ drug, data = aids.id, x = TRUE)
summary(fitSURV)
```

    ## Call:
    ## coxph(formula = Surv(Time, death) ~ drug, data = aids.id, x = TRUE)
    ## 
    ##   n= 467, number of events= 188 
    ## 
    ##           coef exp(coef) se(coef)     z Pr(>|z|)
    ## drugddI 0.2102    1.2339   0.1462 1.437    0.151
    ## 
    ##         exp(coef) exp(-coef) lower .95 upper .95
    ## drugddI     1.234     0.8104    0.9264     1.643
    ## 
    ## Concordance= 0.531  (se = 0.019 )
    ## Likelihood ratio test= 2.07  on 1 df,   p=0.2
    ## Wald test            = 2.07  on 1 df,   p=0.2
    ## Score (logrank) test = 2.07  on 1 df,   p=0.1

We use the `JM::jointModel()` function to build our joint model. We pass
the lme and coxph objects we have created, and provide the time variable
from the LME. We provide a method argument to specify the type of joint
model to fit. A full table of method arguments is included below. We can
apply a time lag, and there is also support for competing risks when
using `method = "spline-PH-GH"`. A full demonstration for competing
risks can be found in the appropriate section.

``` r
fit.JM <- jointModel(fitLME,
                     fitSURV,
                     timeVar = "obstime",
                     method = "piecewise-PH-GH",
                     lag = 0, 
                     CompRisk = FALSE)
summary(fit.JM)
```

    ## 
    ## Call:
    ## jointModel(lmeObject = fitLME, survObject = fitSURV, timeVar = "obstime", 
    ##     method = "piecewise-PH-GH", lag = 0, CompRisk = FALSE)
    ## 
    ## Data Descriptives:
    ## Longitudinal Process     Event Process
    ## Number of Observations: 1405 Number of Events: 188 (40.3%)
    ## Number of Groups: 467
    ## 
    ## Joint Model Summary:
    ## Longitudinal Process: Linear mixed-effects model
    ## Event Process: Relative risk model with piecewise-constant
    ##      baseline risk function
    ## Parameterization: Time-dependent 
    ## 
    ##    log.Lik      AIC      BIC
    ##  -2107.647 4247.295 4313.636
    ## 
    ## Variance Components:
    ##              StdDev    Corr
    ## (Intercept)  0.8660  (Intr)
    ## obstime      0.0388  0.0680
    ## Residual     0.3754        
    ## 
    ## Coefficients:
    ## Longitudinal Process
    ##                   Value Std.Err z-value p-value
    ## (Intercept)      2.5558  0.0372 68.7961 <0.0001
    ## obstime         -0.0423  0.0046 -9.1931 <0.0001
    ## obstime:drugddI  0.0051  0.0065  0.7821  0.4342
    ## 
    ## Event Process
    ##             Value Std.Err z-value p-value
    ## drugddI    0.3511  0.1537  2.2839  0.0224
    ## Assoct    -1.1016  0.1180 -9.3388 <0.0001
    ## log(xi.1) -1.6489  0.2498 -6.6000        
    ## log(xi.2) -1.3393  0.2394 -5.5940        
    ## log(xi.3) -1.0231  0.2861 -3.5758        
    ## log(xi.4) -1.5802  0.3736 -4.2299        
    ## log(xi.5) -1.4722  0.3500 -4.2069        
    ## log(xi.6) -1.4383  0.4283 -3.3584        
    ## log(xi.7) -1.4780  0.5455 -2.7094        
    ## 
    ## Integration:
    ## method: Gauss-Hermite
    ## quadrature points: 15 
    ## 
    ## Optimization:
    ## Convergence: 0

## Interpretting the output

We are immediately presented with the call, typical for `summary()`. We
are given some descriptive data: primarily the number of observations
and the total number of events.

We are then given statistics for the sub models. For the LME, we are
given log-likelihood, AIC, and BIC as well as variance components. We
are also given betas for the LME, alongside their respective standard
errors, their z-values (beta / standard error) and corresponding
p-values.

We are also provided with similar statistics for the survival model.

Convergence indicates if convergence was successful. 0 indicates
successful convergence, whilst 1 indicates an issue.

## Table of `method` arguments

| method            | description                                                                                                          |
| ----------------- | -------------------------------------------------------------------------------------------------------------------- |
| “weibull-AFT-GH”  | time-dependent Weibull model under the accelerated failure time formulation                                          |
| “weibull-PH-GH”   | time-dependent relative risk model postulated with a Weibull baseline risk function                                  |
| “piecewise-PH-GH” | time-dependent relative risk model postulated with a piecewise constant baseline risk function                       |
| “spline-PH-GH”    | time-dependent relative risk model in which the log baseline risk function is approximated using B-splines           |
| “ch-Laplace”      | additive model on the log cumulative hazard scale (see Rizopoulos et al., 2009)                                      |
| “Cox-PH-GH”       | time-dependent relative risk model where the baseline risk function is left unspecified (Wulfsohn and Tsiatis, 1997) |

-----

## Competing risks

It is possible to fit a Joint model with competing risks using the `JM`.
For this example, we will use the `pbc` data-set included with `JM`.

``` r
data(pbc)
knitr::kable(head(pbc, 8))
```

| id | time | status | trt |      age | sex | ascites | hepato | spiders | edema | bili | chol | albumin | copper | alk.phos |    ast | trig | platelet | protime | stage |
| -: | ---: | -----: | --: | -------: | :-- | ------: | -----: | ------: | ----: | ---: | ---: | ------: | -----: | -------: | -----: | ---: | -------: | ------: | ----: |
|  1 |  400 |      2 |   1 | 58.76523 | f   |       1 |      1 |       1 |   1.0 | 14.5 |  261 |    2.60 |    156 |   1718.0 | 137.95 |  172 |      190 |    12.2 |     4 |
|  2 | 4500 |      0 |   1 | 56.44627 | f   |       0 |      1 |       1 |   0.0 |  1.1 |  302 |    4.14 |     54 |   7394.8 | 113.52 |   88 |      221 |    10.6 |     3 |
|  3 | 1012 |      2 |   1 | 70.07255 | m   |       0 |      0 |       0 |   0.5 |  1.4 |  176 |    3.48 |    210 |    516.0 |  96.10 |   55 |      151 |    12.0 |     4 |
|  4 | 1925 |      2 |   1 | 54.74059 | f   |       0 |      1 |       1 |   0.5 |  1.8 |  244 |    2.54 |     64 |   6121.8 |  60.63 |   92 |      183 |    10.3 |     4 |
|  5 | 1504 |      1 |   2 | 38.10541 | f   |       0 |      1 |       1 |   0.0 |  3.4 |  279 |    3.53 |    143 |    671.0 | 113.15 |   72 |      136 |    10.9 |     3 |
|  6 | 2503 |      2 |   2 | 66.25873 | f   |       0 |      1 |       0 |   0.0 |  0.8 |  248 |    3.98 |     50 |    944.0 |  93.00 |   63 |       NA |    11.0 |     3 |
|  7 | 1832 |      0 |   2 | 55.53457 | f   |       0 |      1 |       0 |   0.0 |  1.0 |  322 |    4.09 |     52 |    824.0 |  60.45 |  213 |      204 |     9.7 |     3 |
|  8 | 2466 |      2 |   2 | 53.05681 | f   |       0 |      0 |       0 |   0.0 |  0.3 |  280 |    4.00 |     52 |   4651.2 |  28.38 |  189 |      373 |    11.0 |     3 |

We fit a linear mixed model - similar to before.

``` r
lmeFit.pbc <- lme(log(serBilir) ~ drug * ns(year, 3),
                  random = list(id = pdDiag(form = ~ ns(year, 3))),
                  data = pbc2)
```

However for the survival model, we will use `JM::crLong` to convert the
data to a long format. We also use the competing risks indicator as a
stratification factor.

``` r
pbc2.idCR <- crLong(pbc2.id, "status", "alive")
coxCRFit.pbc <- coxph(Surv(years, status2) ~ (drug + sex) * strata + strata(strata),
                      data = pbc2.idCR,
                      x = TRUE)
```

Then, we use `jointModel()`, setting `method = "spline-PH-aGH"` and
`CompRisk = TRUE`. We include strata as an interaction factor to allow
our longitudinal outcome to have a different effect for each of the two
competing risks.

``` r
jmCRFit.pbc <- jointModel(lmeFit.pbc,
                          coxCRFit.pbc,
                          timeVar = "year",
                          method = "spline-PH-aGH",
                          interFact = list(value = ~ strata, data = pbc2.idCR),
                          CompRisk = TRUE)
summary(jmCRFit.pbc)
```

    ## 
    ## Call:
    ## jointModel(lmeObject = lmeFit.pbc, survObject = coxCRFit.pbc, 
    ##     timeVar = "year", method = "spline-PH-aGH", interFact = list(value = ~strata, 
    ##         data = pbc2.idCR), CompRisk = TRUE)
    ## 
    ## Data Descriptives:
    ## Longitudinal Process     Event Process
    ## Number of Observations: 1945 Number of Events: 169 (54.2%)
    ## Number of Groups: 312
    ## 
    ## Joint Model Summary:
    ## Longitudinal Process: Linear mixed-effects model
    ## Event Process: Competing risks relative risk model with spline-approximated
    ##      baseline risk function
    ## Parameterization: Time-dependent 
    ## 
    ##    log.Lik      AIC      BIC
    ##  -1994.844 4063.688 4202.179
    ## 
    ## Variance Components:
    ##              StdDev
    ## (Intercept)  1.0017
    ## ns(year, 3)1 1.4856
    ## ns(year, 3)2 1.0765
    ## ns(year, 3)3 1.2694
    ## Residual     0.2895
    ## 
    ## Coefficients:
    ## Longitudinal Process
    ##                              Value Std.Err z-value p-value
    ## (Intercept)                 0.5813  0.0718  8.0944 <0.0001
    ## drugD-penicil              -0.1095  0.1016 -1.0781  0.2810
    ## ns(year, 3)1                0.8663  0.1690  5.1275 <0.0001
    ## ns(year, 3)2                1.5466  0.1440 10.7386 <0.0001
    ## ns(year, 3)3                1.4059  0.2243  6.2681 <0.0001
    ## drugD-penicil:ns(year, 3)1  0.2772  0.2342  1.1833  0.2367
    ## drugD-penicil:ns(year, 3)2 -0.4096  0.2041 -2.0065  0.0448
    ## drugD-penicil:ns(year, 3)3 -0.6555  0.3209 -2.0426  0.0411
    ## 
    ## Event Process
    ##                             Value Std.Err z-value p-value
    ## drugD-penicil             -0.4417  0.3858 -1.1449  0.2523
    ## sexfemale                  0.2551  0.6038  0.4225  0.6727
    ## drugD-penicil:stratadead   0.4398  0.4194  1.0486  0.2944
    ## sexfemale:stratadead      -0.5791  0.6400 -0.9048  0.3655
    ## Assoct                     1.1842  0.2041  5.8008 <0.0001
    ## Assoct:stratadead          0.1185  0.2247  0.5272  0.5980
    ## bs1(transplanted)        -10.7415  6.0500 -1.7755  0.0758
    ## bs2(transplanted)         -8.0150  2.4681 -3.2474  0.0012
    ## bs3(transplanted)         -4.5890  1.5932 -2.8804  0.0040
    ## bs4(transplanted)         -6.1349  1.1876 -5.1659 <0.0001
    ## bs5(transplanted)         -4.1988  1.1395 -3.6848  0.0002
    ## bs6(transplanted)         -6.1544  1.6546 -3.7196  0.0002
    ## bs7(transplanted)         -8.5434  5.4371 -1.5713  0.1161
    ## bs8(transplanted)         -8.6613 11.9534 -0.7246  0.4687
    ## bs9(transplanted)         -7.2987 13.1393 -0.5555  0.5786
    ## bs1(dead)                 -3.8781  0.5396 -7.1876 <0.0001
    ## bs2(dead)                 -4.8013  0.6567 -7.3112 <0.0001
    ## bs3(dead)                 -3.7555  0.6523 -5.7577 <0.0001
    ## bs4(dead)                 -4.5481  0.5673 -8.0176 <0.0001
    ## bs5(dead)                 -4.1473  0.5936 -6.9869 <0.0001
    ## bs6(dead)                 -4.4641  0.6588 -6.7760 <0.0001
    ## bs7(dead)                 -1.5794  1.1395 -1.3861  0.1657
    ## bs8(dead)                 -8.2645  2.2495 -3.6739  0.0002
    ## bs9(dead)                 -2.0341  1.9005 -1.0703  0.2845
    ## 
    ## Integration:
    ## method: (pseudo) adaptive Gauss-Hermite
    ## quadrature points: 3 
    ## 
    ## Optimization:
    ## Convergence: 0
