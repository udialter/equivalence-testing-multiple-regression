# Equivalence Testing for Multiple Regression

Equivalence testing can be applied to evaluate whether an observed
effect from an individual predictor in a multiple regression model is
small enough to be considered statistically and practically negligible
([Alter & Counsell, 2021](https://psyarxiv.com/5bv2p/)). For more
information, please refer to the [OSF
page](https://osf.io/w96xe/ "Equivalence Testing for Multiple Regression OSF Page")
and/or a freely available [preprint on
PsyArXiv](https://psyarxiv.com/5bv2p/ "PsyArXiv preprint").

The following functions offer appropriate equivalence-based alternatives
for concluding negligible effect between a predictor and outcome in
multiple regression

These R functions are designed to accommodate multiple research contexts
effortlessly, with or without access to the full dataset. The two
functions, `reg.equiv.fd()` and `reg.equiv()`, provide similar output
but differ on the type of input information required by the user.

Specifically, the first function, `reg.equiv.fd()`, requires the full
dataset and model in R (`lm` object), whereas the second does not. The
`reg.equiv()` is intended for researchers who do not have access to the
complete dataset but still wish to evaluate a certain predictor’s lack
of association with the outcome variable in multiple regression, for
example, using information typically presented in a results section or
table reported in a published article.

## `reg.equiv.fd()`: full dataset required

### *Required* input:

-   `datfra=` a data frame (e.g., mtcars)
-   `model=` the model, an lm object (e.g., `mod1`, where
    `mod1<- mpg~hp+cyl`)
-   `delta=` the smallest effect size of interest (SESOI), minimally
    meaningful effect size (MMES), or upper bound of the equivalence
    interval (𝛅) (e.g., .15)
-   `predictor=` the name of the predictor to be tested (e.g., `"cyl"`)

### Default settings:

-   `test=` test type is set automatically to Two One-Sided Test (TOST;
    Schuirmann, 1987), the other option is the Anderson-Hauck (AH;
    Anderson & Hauck, 1983)
-   `std=` the delta (or, SESOI) is the set as as standardized by
    default. Indicate `std=FALSE` to assume unstandardized units
-   `alpha=` the nominal Type I error rate is set to .05 by default. To
    change, simply indicate the alpha level. E.g., `alpha=.10`

`reg.equiv.fd()` example: ![alt
text](Figure 4 reg.equiv.fd .png)

## `reg.equiv()`: full dataset _not_ required

### *Required* input:

-   `b=` the estimated effect size associated with the predictor of
    interest, this could be either standardized or unstandardized (e.g.,
    .02)
-   `se=` the standard error associated with the effect size of the
    predictor of interest (if the effect size is standardized, make sure
    the `se` value is tied to the standardized and not the raw effect)
-   `p=` the number of total predictors in the regression model
    (excluding intercept)
-   `n=` sample size
-   `delta=` the smallest effect size of interest (SESOI), minimally
    meaningful effect size (MMES), or upper bound of the equivalence
    interval (𝛅) (e.g., .15)
-   `predictor=` the name of the predictor to be tested (e.g., `"cyl"`)

### Default settings:

-   `test=` test type is set automatically to Two One-Sided Test (TOST;
    Schuirmann, 1987), the other option is the Anderson-Hauck (AH;
    Anderson & Hauck, 1983)
-   `std=` the delta (or, SESOI) and the indicated effect size are set
    as as standardized by default. Indicate `std=FALSE` to assume
    unstandardized units
-   `alpha=` the nominal Type I error rate is set to .05 by default. To
    change, simply indicate the alpha level. E.g., `alpha=.10`

`reg.equiv()` example: ![alt
text](~/Documents/Master's%20thesis/Figure%205%20reg.equiv.png)

### Recommendations for Test Result Interpretation

Equivalence testing is a method designed within the null-hypothesis
significance testing (NHST) framework. NHST has been heavily criticized
for its overreliance on the dichotomous results of *p* values with
little, or no consideration of the effect’s magnitude or its
implications in practice (e.g., Cumming, 2012; Fidler & Loftus, 2009;
Harlow, 1997; Kirk, 2003; Lee, 2016 2014). Researchers must be mindful
of the limitations of NHST, and disentangle the practical and
statistical aspects of the test results.

To minimize the limitations of *p* values, it is more informative to
**interpret the observed effect’s magnitude and precision** beyond the
conclusion of “negligible effects” or “insufficient evidence for
negligible effects.” **Observed effects should be construed in relation
to the equivalence bounds, the extent of their uncertainty, and their
practical implications (or lack thereof)**. For this reason, **the two R
functions offered here also include a graphical representation of the
observed effect and its associated uncertainty in relation to the
equivalence interval. The resulting plot aids in illustrating how close
or far and wide or narrow the observed effect and its margin of error
are from the equivalence bounds; inferring about the proportion and
position of the confidence band in relation to the equivalence interval
can help interpret the results over and above *p* values**.

#### References

1.  Alter, U., & Counsell, A. (2021, June 17). Equivalence Testing for
    Multiple Regression. <https://doi.org/10.17605/OSF.IO/W96XE>
2.  Anderson, S., & Hauck, W. W. (1983). A new procedure for testing
    equivalence in comparative bioavailability and other clinical
    trials. Statistics and Communications-Theory and Methods, 12
    ,2663-2692. <https://doi.org/10.1080/03610928308828634>
3.  Cumming, G. (2012). Understanding The New Statistics: Effect Sizes,
    Confidence Intervals, and Meta-Analysis. New York, NY: Taylor &
    Francis Group, LLC
4.  Cumming, G. (2014). The New Statistics: Why and How. Psychological
    Science, 25(1), 7–29. <https://doi.org/10.1177/0956797613504966>
5.  Fidler, F., & Loftus, G. (2009). Why figures with error bars should
    replace p values: Some conceptual arguments and empirical
    demonstrations. Zeitschrift Für Psychologie/Journal of Psychology,
    217(1), 27-37. <https://doi.org/10.1027/0044-3409.217.1.27>
6.  Harlow, L.L. (1997). Significance Testing in Introduction and
    Overview. In L.L. Harlow, S.A. Muliak & J.H. Steiger (Eds.). What If
    There Were No Significance Tests? (pp.1-17). Mahwah, NJ, USA:
    Lawrence Erlbaum.
7.  Hauck, W. W., & Anderson, S. (1984). A new statistical procedure for
    testing equivalence in two-group comparative bioavailability trials.
    Journal of Pharmacokinetics and Biopharmaceutics, 12(1), 83-91.
    <https://doi.org/10.1007/BF01063612>
8.  Kirk, R. E. (2003). The importance of effect magnitude. In S. F.
    Davis (Ed.), Handbook of research methods in experimental psychology
    (pp. 83–105). Malden, MA: Blackwell.
9.  Lee, D. K. (2016). Alternatives to P value: Confidence interval and
    effect size. Korean Journal of Anesthesiology, 69(6), 555-562.
    <https://doi.org/10.4097/kjae.2016.69.6.555>
10. Schuirmann, D. J. (1987). A comparison of the two one-sided tests
    procedure and the power approach for assessing the equivalence of
    average bioavailability. Journal of Pharmacokinetics and
    Biopharmaceutics, 15 ,657-680. <https://doi.org/10.1007/BF01068419>
