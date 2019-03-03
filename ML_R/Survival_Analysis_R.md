---
title: "Survival Analysis With R"
author: "Meichen Lu"
date: "March 3, 2019"
output: 
  html_document:
    keep_md: true
    toc: true
    toc_depth: 3
---
## Setup

```r
library(survival)
library(survminer)
```

```
## Loading required package: ggplot2
```

```
## Loading required package: ggpubr
```

```
## Loading required package: magrittr
```

```r
library(splines)
library(lattice)
data("lung", package = "survival")
lung$sex <- factor(lung$sex, levels = 1:2, labels = c("male", "female"))
```
We use the lung dataset from the survival model, consisting of data from 228 patients. 

## Non-parametric model
The `Surv()` function from the `survival` package create a survival object, which is used in many other functions. In most cases, the first argument the observed survival times, and as second the event indicator. The type of censoring is also specified in this function.

To calculate the Kaplan-Meier (KM) estimate, we use `survfit()` function and its first argument is an R formula. The left-hand side of this formula is a survival object created by `Surv()`, and the right-hand side specifies the grouping variable. Here we would like to estimate a single survival curve of the entire population, so we put $1$ as below. 

We can call the plot() method on the survfit object to display the estimated curve, where the 95% confidence interval is included by default:


```r
KM_fit <- survfit(Surv(time, status) ~ 1, data = lung)
plot(KM_fit, xlab = "Time", ylab = "Survival Probability", 
    main = "Kaplan-Meier Estimate of S(t) for the lung Data")
```

![](Survival_Analysis_R_files/figure-html/unnamed-chunk-1-1.png)<!-- -->
The ggsurvplot() function from `survminer` creates ggplot2 plots from survfit objects.


```r
ggsurvplot(KM_fit, data = lung)
```

![](Survival_Analysis_R_files/figure-html/unnamed-chunk-2-1.png)<!-- -->

We can also plot the cumulative events or cumulative hazard by specifying  `fun = "event"` or `fun = "cumhaz"` in the ggsurvplot function.

The Nelson-Aalen estimator is also known as the Breslow and can be estimated by specifying the `type` argument to be `"fleming-harrington"`. You can see below that the two estimators give nearly identical survival curves.

```r
NA_fit <- survfit(Surv(time, status) ~ 1, data = lung, 
                  type = "fleming-harrington")
plot(KM_fit, conf.int = FALSE, xlab = "Time", ylab = "Survival Probability")
lines(NA_fit, conf.int = FALSE, col = 'red')
```

![](Survival_Analysis_R_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

## Cox proportional hazard model

### Model fitting and significance test
`coxph()` fits a Cox proportional hazard model to the data and the syntax is similar to `survfit()`. Here, we fit a model using only the age predictor and called `summary()` to examine the details of the coxph fit. From the output, we can see that the coefficient for age is greater than $0$ and $\exp(\text{coef}) > 1$, meaning that the age contributes to increasing hazard.

```r
cox_fit <- coxph(Surv(time, status) ~ age, data = lung)
summary(cox_fit)
```

```
## Call:
## coxph(formula = Surv(time, status) ~ age, data = lung)
## 
##   n= 228, number of events= 165 
## 
##         coef exp(coef) se(coef)     z Pr(>|z|)  
## age 0.018720  1.018897 0.009199 2.035   0.0419 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
##     exp(coef) exp(-coef) lower .95 upper .95
## age     1.019     0.9815     1.001     1.037
## 
## Concordance= 0.55  (se = 0.025 )
## Rsquare= 0.018   (max possible= 0.999 )
## Likelihood ratio test= 4.24  on 1 df,   p=0.04
## Wald test            = 4.14  on 1 df,   p=0.04
## Score (logrank) test = 4.15  on 1 df,   p=0.04
```
To test whether gender is  an important factor, we can conduct a log rank test by calling `survdiff()` function. In this example, the test statistics $\chi^2 = 10.3$ is equivalent to $p = 0.001$, indicating that sex is a significant factor.


```r
survdiff(Surv(time, status) ~ sex, data = lung)
```

```
## Call:
## survdiff(formula = Surv(time, status) ~ sex, data = lung)
## 
##              N Observed Expected (O-E)^2/E (O-E)^2/V
## sex=male   138      112     91.6      4.55      10.3
## sex=female  90       53     73.4      5.68      10.3
## 
##  Chisq= 10.3  on 1 degrees of freedom, p= 0.001
```


We can also use anova to compare two coxph model fitted to the same data. Here, the test results evaluate how significant sex is, given that we already have the age predictor.

```r
cox_fit2 <- coxph(Surv(time, status) ~ sex + age, data = lung)
anova(cox_fit, cox_fit2)
```

```
## Analysis of Deviance Table
##  Cox model: response is  Surv(time, status)
##  Model 1: ~ age
##  Model 2: ~ sex + age
##    loglik  Chisq Df P(>|Chi|)   
## 1 -747.79                       
## 2 -742.85 9.8822  1  0.001669 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

### Prediction
To use coxph to make prediction, we can derive the following formula from the definition. In coxph model
$$S(t|X) = \exp[-\Lambda(t)\exp(X\beta)] = \exp[-\Lambda(t)]^{\exp(X\beta)} = S_{KM}(t)^{\exp(X\beta)}$$
where $S_{KM}$ is the KM-estimator. In the following computation, we computed $X\beta$ for the first patient in the `lung` dataset. The mean (`cox_fit2$means`) is subtracted from $X$ and the linear combination of mean-adjusted $X$ given $\beta$ is also called linear predictors (lp). As we can see, lp computed manually using the definition is exactly the same in `cox_fit2$linear.predictors`. 


```r
data <- lung[1, c('sex', 'age')]
data$sex = (data$sex == 'female')
lp <- sum(t(data - cox_fit2$means) * cox_fit2$coefficients)
data.frame(lp_manual = lp, 
           lp_model = cox_fit2$linear.predictors[1])
```

```
##   lp_manual  lp_model
## 1 0.3995047 0.3995047
```

To obtain the individual survival curve, we can extract $S_{KM}$ from a KM estimate and the survival curve is calculated using the previous equation.

```r
pred <- KM_fit$surv ^ exp(lp)
plot(KM_fit, conf.int = FALSE, xlab = "Time", ylab = "Survival Probability")
lines(KM_fit$time, pred, type = 's', col = "red")
legend("topright",
       legend = c('baseline', 'first patient'),
       lty = c(1,1),
       col = c('black', 'red'))
```

![](Survival_Analysis_R_files/figure-html/unnamed-chunk-8-1.png)<!-- -->



The alternative is to use cumulative hazard function, which can be obtained from calling `basehaz()` function. Here, we use the cumulative hazard function to verify the calculation of martingale residual, which is defined as
$$r_{Mi} = \delta_i - \Lambda_0(T_i)\exp(X_j\beta)$$
where $\delta_i$ is the state ($0$ is alive/censored and $1$ is dead) and $\Lambda_0(T_i)\exp(\hat{\beta}Z_j)$ is the cumulative hazard function at time $T_i$. The $r_{Mi}$ is dependent on time but the one we are most interested in is the maximum survival time. In the following example, we use `basehaz()` to obtain $\Lambda_0$ and `cox_fit2$linear.predictors` for $X_j\beta$. Martingale residual is also accessible at `cox_fit2$residuals`.


```r
# baseline hazard is the cumulative hazard function
cum_haz <- basehaz(cox_fit2)
# For the first patient in lung, he died at day = 306, the index is 109
# Thus the prediction is
pred1 <- cum_haz[109,1]*exp(cox_fit2$linear.predictors[1])
data.frame(residual_manual = 1 - pred1, 
           residual_model = cox_fit2$residuals[1])
```

```
##   residual_manual residual_model
## 1     0.004389994    0.004389994
```


### Residual and effect plots
As martingal residual evaluates the goodness of fit, it can be used to examine non-linearity by plotting residuals against one of the predictors. In this example, we suspect that the survival may have a non-linear relationship with age.

```r
lung$martingale <- residuals(cox_fit2, type = "martingale")
# this is equivalent to lung$martingale <- cox_fit2$residuals
xyplot(martingale ~ age | sex, data = lung, main='martingale residual of coxph')
```

![](Survival_Analysis_R_files/figure-html/residual_linear-1.png)<!-- -->

Here, we also show the effect plots, further confirming the linear relationship.


```r
new_data <- with(lung, expand.grid(
    age = seq(45, 75, 5),
    sex = levels(sex)
))
predictions <- predict(cox_fit2, newdata = new_data, type = "lp", se.fit = TRUE)
new_data$predict <- predictions$fit
new_data$se <- predictions$se.fit
new_data$lower <- predictions$fit - 1.96 * predictions$se.fit
new_data$upper <- predictions$fit + 1.96 * predictions$se.fit
head(new_data,3)
```

```
##   age  sex      predict         se      lower     upper
## 1  45 male -0.094809928 0.17730643 -0.4423305 0.2527107
## 2  50 male -0.009583269 0.13559528 -0.2753500 0.2561835
## 3  55 male  0.075643391 0.09791848 -0.1162768 0.2675636
```


```r
xyplot(predict + lower + upper ~ age | sex, data = new_data, 
       type = "l", col = "black", lwd = 2, lty = c(1, 2, 2),
       abline = list(h = 0, lty = 2, lwd = 2, col = "red"),
       xlab = "Age (years)", ylab = "log Hazard Ratio")
```

![](Survival_Analysis_R_files/figure-html/effect_plot_linear-1.png)<!-- -->

Then we introduce non-linearity using a natural cubic spline on age variable. 

```r
cox_fit2_ns <- coxph(Surv(time, status) ~ sex + ns(age, 3), data = lung)
cox_fit2_ns
```

```
## Call:
## coxph(formula = Surv(time, status) ~ sex + ns(age, 3), data = lung)
## 
##                coef exp(coef) se(coef)      z       p
## sexfemale   -0.5037    0.6043   0.1679 -3.001 0.00269
## ns(age, 3)1  0.1285    1.1371   0.3595  0.357 0.72074
## ns(age, 3)2  1.7182    5.5744   1.1618  1.479 0.13916
## ns(age, 3)3  1.1237    3.0762   0.4796  2.343 0.01914
## 
## Likelihood ratio test=16.05  on 4 df, p=0.002956
## n= 228, number of events= 165
```

From the model summary of `cox_fit2_ns`, the non-linear terms of the age predictor are not all significant. We can see that there is not significant improvement on the residual plot.

```r
lung$martingale <- residuals(cox_fit2_ns, type = "martingale")
xyplot(martingale ~ age | sex, data = lung, main='martingale residual of coxph')
```

![](Survival_Analysis_R_files/figure-html/residual_nonlinear-1.png)<!-- -->

Regardless, it is helpful to visualise the nonlinear relationships.

```r
predictions <- predict(cox_fit2_ns, newdata = new_data, type = "lp", se.fit = TRUE)
new_data$predict <- predictions$fit
new_data$se <- predictions$se.fit
new_data$lower <- predictions$fit - 1.96 * predictions$se.fit
new_data$upper <- predictions$fit + 1.96 * predictions$se.fit

xyplot(predict + lower + upper ~ age | sex, data = new_data, 
       type = "l", col = "black", lwd = 2, lty = c(1, 2, 2),
       abline = list(h = 0, lty = 2, lwd = 2, col = "red"),
       xlab = "Age (years)", ylab = "log Hazard Ratio")
```

![](Survival_Analysis_R_files/figure-html/effect_plot_nonlinear-1.png)<!-- -->

### Testing proportionality assumption
Schonfeld residual is often used to test the proportionality assumption and we can use `ggcoxzph()` from `survminer` package to produce the diagnostic plot. As there is slight non-linear trend in the diagnostic plot with respect to `sex`, we can hypothesise that the proportionality may not hold for `sex`.


```r
ftest <- cox.zph(cox_fit2)
ggcoxzph(ftest)
```

![](Survival_Analysis_R_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

To address this potential violation of proportionality assumption, one way is to use stratification. If we startify using `sex`, we are essentially fitting two models regarding different sex.

```r
cox_fit_strat <- coxph(Surv(time, status) ~ strata(sex)*age, data = lung)
predictions <- predict(cox_fit_strat, newdata = new_data, type = "lp", se.fit = TRUE)
new_data$predict <- predictions$fit
new_data$se <- predictions$se.fit
new_data$lower <- predictions$fit - 1.96 * predictions$se.fit
new_data$upper <- predictions$fit + 1.96 * predictions$se.fit
ggplot(new_data, aes(x = age, y = predict, col = sex)) +
    geom_line() + xlab("Age (years)") + ylab("log Hazard Ratio")
```

![](Survival_Analysis_R_files/figure-html/stratified_model-1.png)<!-- -->




### ROC and AUC
So far, we have looked at several ways to diagnose the fitted coxph model and tricks to improve the model by adding complexity, such as non-linearity and stratification. One missing piece is to have a single number that evaluates the goodness-of-fit. As I pointed out last time, survival analysis is equivalent to a classification problem at a fixed time point. Therefore, a useful metric is AUC for different models. Here, I used `survivalROC` package as a demonstration and many advanced AUC measures are implemented in `survAUC`.


```r
library(survivalROC)
```


```r
par(mfrow=c(2,2))
for (i in seq(1,4)){
    t = 200 * i
    roc <- survivalROC(Stime = lung$time,
                       status = lung$status,
                       marker = cox_fit2$linear.predictors,
                       predict.time = t,
                       method = "NNE",
                       span = 0.25 * nrow(lung)^(-0.20))
    fig_title <- sprintf("AUC at time %d: %.3f",t, roc$AUC)
    plot(roc$FP, roc$TP, xlab = "FP", ylab = "TP", main =  fig_title)
}
```

![](Survival_Analysis_R_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

We then plot the ROC for `cox_fit2_ns` model to examine the difference.

```r
par(mfrow=c(2,2))
for (i in seq(1,4)){
    t = 200 * i
    roc <- survivalROC(Stime = lung$time,
                       status = lung$status,
                       marker = cox_fit2_ns$linear.predictors,
                       predict.time = t,
                       method = "NNE",
                       span = 0.25 * nrow(lung)^(-0.20))
    fig_title <- sprintf("AUC at time %d: %.3f",t, roc$AUC)
    plot(roc$FP, roc$TP, xlab = "FP", ylab = "TP", main =  fig_title)
}
```

![](Survival_Analysis_R_files/figure-html/unnamed-chunk-13-1.png)<!-- -->

Comparing these two models, we can see that adding non-linearity improves the AUC for long-term prediction (see t = 600, 800) while the performance is worse at short-term (t = 200, 400).

