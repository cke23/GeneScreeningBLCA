Gene Screening for Prognosis of Non-Muscle-Invasive Bladder Carcinoma
under Competing Risks Endpoints
================

Bladder cancer is a common type of cancer associated with high morbidity
and mortality rates, if not treated optimally. Early diagnosis with
personalized treatment and follow-up is the key to a successful outcome.
Developments in microarray and sequencing technologies have allowed
collection of massive genomic information that substantially advances
the understanding of molecular mechanisms, biomarker discovery, and
personalized medicine. High-throughput data produced by those techniques
are characterized by a large number of features that far exceeds the
sample size
(![p\>\>n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p%3E%3En "p>>n")).
In clinical studies, a vital research task is to find predictive
features for survival outcomes and build prognostic models for cancer
patients, which often requires techniques that were developed for
specific time-to-event responses. In bladder cancer, one primary
endpoint is time-to-progression, but competing events such as death from
non-cancer causes can also be observed. The proportional subdistribution
hazard (PSH) model proposed by Fine and Gray (Fine and Gray 1999) has
become a popular semi-parametric model for competing risks data, which
can be regarded as an adaption of the Cox PH model for the
subdistribution hazard related to the cumulative incidence function. In
this document, we demonstrate the implementation of a PSH-based gene
screening procedure for high-throughput competing risks data.

## The Gene Screening Procedure

Let
![T](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;T "T")
and
![C](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;C "C")
be the failure and censoring times, and
![\\epsilon\\in\\{1,...,K\\}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cepsilon%5Cin%5C%7B1%2C...%2CK%5C%7D "\epsilon\in\{1,...,K\}")
be cause of failure. Let
![\\mathbf{X}\\in \\mathbb{R}^p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BX%7D%5Cin%20%5Cmathbb%7BR%7D%5Ep "\mathbf{X}\in \mathbb{R}^p")
denote the vector of
![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p")
covariates subject to selection and
![\\mathbf{Z}\\in \\mathbb{R}^{p_0}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BZ%7D%5Cin%20%5Cmathbb%7BR%7D%5E%7Bp_0%7D "\mathbf{Z}\in \mathbb{R}^{p_0}")
denote the vector of
![p_0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p_0 "p_0")
covariates to be controlled for in the analysis. For typical
right-censored data, we observe
![Y=\\min(T,C)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;Y%3D%5Cmin%28T%2CC%29 "Y=\min(T,C)")
and
![\\delta=I(T\\leq C)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cdelta%3DI%28T%5Cleq%20C%29 "\delta=I(T\leq C)"),
where
![I(\\cdot)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;I%28%5Ccdot%29 "I(\cdot)")
is the indicator function. Our goal is to model the cumulative incidence
function (CIF) for failure from the cause of interest
(![\\epsilon=1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cepsilon%3D1 "\epsilon=1"))
conditional on the covariates:

![
F_1(t;\\mathbf{X},\\mathbf{Z}) = Pr(T\\leq t, \\epsilon=1\|\\mathbf{X},\\mathbf{Z}),
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0AF_1%28t%3B%5Cmathbf%7BX%7D%2C%5Cmathbf%7BZ%7D%29%20%3D%20Pr%28T%5Cleq%20t%2C%20%5Cepsilon%3D1%7C%5Cmathbf%7BX%7D%2C%5Cmathbf%7BZ%7D%29%2C%0A "
F_1(t;\mathbf{X},\mathbf{Z}) = Pr(T\leq t, \epsilon=1|\mathbf{X},\mathbf{Z}),
")

i.e., the probability of experiencing event 1 before time
![t](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;t "t")
and before the occurrence of any other types of event. The
subdistribution hazard (Gray 1988) associated with event 1 is defined as
which measures the instantaneous risk of failure from event 1 for
patients who have not yet experienced the event. Note that this risk set
includes those who are currently event free as well as who have
previously experienced a competing event. The subdistribution hazard for
event 1 is assumed to follow a proportional hazard model (Fine and Gray
1999)

![
\\lambda_1(t;\\mathbf{X},\\mathbf{Z}) = \\lambda\_{1,0}(t)\\exp(\\boldsymbol{\\beta}^T\\mathbf{X}+\\boldsymbol{\\gamma}^T\\mathbf{Z}),
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Clambda_1%28t%3B%5Cmathbf%7BX%7D%2C%5Cmathbf%7BZ%7D%29%20%3D%20%5Clambda_%7B1%2C0%7D%28t%29%5Cexp%28%5Cboldsymbol%7B%5Cbeta%7D%5ET%5Cmathbf%7BX%7D%2B%5Cboldsymbol%7B%5Cgamma%7D%5ET%5Cmathbf%7BZ%7D%29%2C%0A "
\lambda_1(t;\mathbf{X},\mathbf{Z}) = \lambda_{1,0}(t)\exp(\boldsymbol{\beta}^T\mathbf{X}+\boldsymbol{\gamma}^T\mathbf{Z}),
")

where
![\\lambda\_{1,0}(t)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Clambda_%7B1%2C0%7D%28t%29 "\lambda_{1,0}(t)")
is an unspecified baseline hazard, and
![\\boldsymbol{\\beta}\\in \\mathbb{R}^p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cboldsymbol%7B%5Cbeta%7D%5Cin%20%5Cmathbb%7BR%7D%5Ep "\boldsymbol{\beta}\in \mathbb{R}^p")
and
![\\boldsymbol{\\gamma} \\in \\mathbb{R}^{p_0}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cboldsymbol%7B%5Cgamma%7D%20%5Cin%20%5Cmathbb%7BR%7D%5E%7Bp_0%7D "\boldsymbol{\gamma} \in \mathbb{R}^{p_0}")
are regression coefficients. Given a finite sample
![\\{\\mathbf{X}\_i,\\mathbf{Z}\_i,Y_i,\\delta_i,\\epsilon_i\\}\_{i=1}^n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5C%7B%5Cmathbf%7BX%7D_i%2C%5Cmathbf%7BZ%7D_i%2CY_i%2C%5Cdelta_i%2C%5Cepsilon_i%5C%7D_%7Bi%3D1%7D%5En "\{\mathbf{X}_i,\mathbf{Z}_i,Y_i,\delta_i,\epsilon_i\}_{i=1}^n"),
the coefficients can be estimated by maximizing the log partial
likelihood function

![
l_n(\\boldsymbol{\\beta},\\boldsymbol{\\gamma}) = \\sum\_{i=1}^n\\int_0^\\infty
\\{\\boldsymbol{\\beta}^T\\mathbf{x}\_i+\\boldsymbol{\\gamma}^T\\mathbf{z}\_i - 
\\log\\sum_j w_j(t)R_j(t)\\exp(\\boldsymbol{\\beta}^T\\mathbf{x}\_j+\\boldsymbol{\\gamma}^T\\mathbf{z}\_j)\\} w_i(t)dN_i(t),
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0Al_n%28%5Cboldsymbol%7B%5Cbeta%7D%2C%5Cboldsymbol%7B%5Cgamma%7D%29%20%3D%20%5Csum_%7Bi%3D1%7D%5En%5Cint_0%5E%5Cinfty%0A%5C%7B%5Cboldsymbol%7B%5Cbeta%7D%5ET%5Cmathbf%7Bx%7D_i%2B%5Cboldsymbol%7B%5Cgamma%7D%5ET%5Cmathbf%7Bz%7D_i%20-%20%0A%5Clog%5Csum_j%20w_j%28t%29R_j%28t%29%5Cexp%28%5Cboldsymbol%7B%5Cbeta%7D%5ET%5Cmathbf%7Bx%7D_j%2B%5Cboldsymbol%7B%5Cgamma%7D%5ET%5Cmathbf%7Bz%7D_j%29%5C%7D%20w_i%28t%29dN_i%28t%29%2C%0A "
l_n(\boldsymbol{\beta},\boldsymbol{\gamma}) = \sum_{i=1}^n\int_0^\infty
\{\boldsymbol{\beta}^T\mathbf{x}_i+\boldsymbol{\gamma}^T\mathbf{z}_i - 
\log\sum_j w_j(t)R_j(t)\exp(\boldsymbol{\beta}^T\mathbf{x}_j+\boldsymbol{\gamma}^T\mathbf{z}_j)\} w_i(t)dN_i(t),
")

where
![N_i(t)=I(T_i\\leq t, \\epsilon_i=1)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;N_i%28t%29%3DI%28T_i%5Cleq%20t%2C%20%5Cepsilon_i%3D1%29 "N_i(t)=I(T_i\leq t, \epsilon_i=1)"),
![R_i(t)=1-N_i(t-)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;R_i%28t%29%3D1-N_i%28t-%29 "R_i(t)=1-N_i(t-)"),
and
![w_i(t)=I(C_i\\geq T_i\\wedge t)\\hat{G}(t)/\\hat{G}(T_i\\wedge t)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;w_i%28t%29%3DI%28C_i%5Cgeq%20T_i%5Cwedge%20t%29%5Chat%7BG%7D%28t%29%2F%5Chat%7BG%7D%28T_i%5Cwedge%20t%29 "w_i(t)=I(C_i\geq T_i\wedge t)\hat{G}(t)/\hat{G}(T_i\wedge t)")
is the inverse probability of censoring weight with
![\\hat{G}(t)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7BG%7D%28t%29 "\hat{G}(t)")
being the Kaplan-Meier estimate for
![G(t)=Pr(C\\geq t)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;G%28t%29%3DPr%28C%5Cgeq%20t%29 "G(t)=Pr(C\geq t)").
Denote the maximizer by
![(\\hat{\\boldsymbol{\\beta}},\\hat{\\boldsymbol{\\gamma}}) = \\arg\\max\_{\\boldsymbol{\\beta},\\boldsymbol{\\gamma}} l_n(\\boldsymbol{\\beta},\\boldsymbol{\\gamma})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%28%5Chat%7B%5Cboldsymbol%7B%5Cbeta%7D%7D%2C%5Chat%7B%5Cboldsymbol%7B%5Cgamma%7D%7D%29%20%3D%20%5Carg%5Cmax_%7B%5Cboldsymbol%7B%5Cbeta%7D%2C%5Cboldsymbol%7B%5Cgamma%7D%7D%20l_n%28%5Cboldsymbol%7B%5Cbeta%7D%2C%5Cboldsymbol%7B%5Cgamma%7D%29 "(\hat{\boldsymbol{\beta}},\hat{\boldsymbol{\gamma}}) = \arg\max_{\boldsymbol{\beta},\boldsymbol{\gamma}} l_n(\boldsymbol{\beta},\boldsymbol{\gamma})").
Having obtained the estimated regression coefficients, the estimated CIF
is obtained by

![
\\hat{F}\_1(t) = 1-\\exp(-\\hat{H}\_{1,0}(t)\\exp(\\hat{\\boldsymbol{\\beta}}^T\\mathbf{x}\_j+\\hat{\\boldsymbol{\\gamma}}^T\\mathbf{z}\_j)),
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Chat%7BF%7D_1%28t%29%20%3D%201-%5Cexp%28-%5Chat%7BH%7D_%7B1%2C0%7D%28t%29%5Cexp%28%5Chat%7B%5Cboldsymbol%7B%5Cbeta%7D%7D%5ET%5Cmathbf%7Bx%7D_j%2B%5Chat%7B%5Cboldsymbol%7B%5Cgamma%7D%7D%5ET%5Cmathbf%7Bz%7D_j%29%29%2C%0A "
\hat{F}_1(t) = 1-\exp(-\hat{H}_{1,0}(t)\exp(\hat{\boldsymbol{\beta}}^T\mathbf{x}_j+\hat{\boldsymbol{\gamma}}^T\mathbf{z}_j)),
")

where
![\\hat{H}\_{1,0}(t)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7BH%7D_%7B1%2C0%7D%28t%29 "\hat{H}_{1,0}(t)")
is the Breslow estimator of the cumulative baseline subdistribution
hazard. The prediction error of the estimated CIF can then be calculated
as

![
Err(t) = \\frac{1}{n}\\sum\_{i=1}^n(N(t)-\\hat{F}\_1(t))^2w_i(t)
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0AErr%28t%29%20%3D%20%5Cfrac%7B1%7D%7Bn%7D%5Csum_%7Bi%3D1%7D%5En%28N%28t%29-%5Chat%7BF%7D_1%28t%29%29%5E2w_i%28t%29%0A "
Err(t) = \frac{1}{n}\sum_{i=1}^n(N(t)-\hat{F}_1(t))^2w_i(t)
")

after accounting for censoring, which is often used to evaluate the
performance of the fitted PSH model. The estimation, prediction and
evaluation of the PSH model can be achieved via the R packages and .
However, the partial likelihood estimation is no longer applicable when
![p+p_0\>n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p%2Bp_0%3En "p+p_0>n"),
demanding new techniques to be developed for high dimensional settings.

In ultrahigh dimensional sparse PSH model with
![p+p_0\>\>n](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p%2Bp_0%3E%3En "p+p_0>>n"),
we realistically assume that the true parameter
![\\boldsymbol{\\beta}=(\\beta_1,...,\\beta_p)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cboldsymbol%7B%5Cbeta%7D%3D%28%5Cbeta_1%2C...%2C%5Cbeta_p%29 "\boldsymbol{\beta}=(\beta_1,...,\beta_p)")
is sparse. In other words, the subset

![
\\mathbf{X}\_\\mathcal{A}=\\{X_j:\\beta_j\\neq 0, j=1,...,p\\}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cmathbf%7BX%7D_%5Cmathcal%7BA%7D%3D%5C%7BX_j%3A%5Cbeta_j%5Cneq%200%2C%20j%3D1%2C...%2Cp%5C%7D%0A "
\mathbf{X}_\mathcal{A}=\{X_j:\beta_j\neq 0, j=1,...,p\}
")

is small. Our aim is therefore to identify the active subset
![\\mathbf{X}\_\\mathcal{A}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BX%7D_%5Cmathcal%7BA%7D "\mathbf{X}_\mathcal{A}")
and estimate
![\\boldsymbol{\\beta}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cboldsymbol%7B%5Cbeta%7D "\boldsymbol{\beta}"),
as well as predicting the survival outcome. For each
![j=1,...,p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;j%3D1%2C...%2Cp "j=1,...,p"),
consider the PSH model containing an individual predictor
![X_j](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;X_j "X_j")
in addition to the control variables

![
\\lambda_1(t;X_j,\\mathbf{Z}) = \\lambda\_{1,0}(t)\\exp(\\beta_jX_j+\\boldsymbol{\\gamma}^T\\mathbf{Z}),
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Clambda_1%28t%3BX_j%2C%5Cmathbf%7BZ%7D%29%20%3D%20%5Clambda_%7B1%2C0%7D%28t%29%5Cexp%28%5Cbeta_jX_j%2B%5Cboldsymbol%7B%5Cgamma%7D%5ET%5Cmathbf%7BZ%7D%29%2C%0A "
\lambda_1(t;X_j,\mathbf{Z}) = \lambda_{1,0}(t)\exp(\beta_jX_j+\boldsymbol{\gamma}^T\mathbf{Z}),
")

and let
![\\hat{u}\_j=\\max\_{\\beta_j,\\boldsymbol{\\gamma}}l_n(\\beta_j,\\boldsymbol{\\gamma})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7Bu%7D_j%3D%5Cmax_%7B%5Cbeta_j%2C%5Cboldsymbol%7B%5Cgamma%7D%7Dl_n%28%5Cbeta_j%2C%5Cboldsymbol%7B%5Cgamma%7D%29 "\hat{u}_j=\max_{\beta_j,\boldsymbol{\gamma}}l_n(\beta_j,\boldsymbol{\gamma})")
denote its estimated maximum log partial likelihood. Then
![\\hat{u}\_j](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7Bu%7D_j "\hat{u}_j")
can be regarded as the change in likelihood associated with
![X_j](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;X_j "X_j")
given
![\\mathbf{Z}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BZ%7D "\mathbf{Z}"),
after ignoring the common constant
![\\max\_{\\boldsymbol{\\gamma}}l_n(\\boldsymbol{\\gamma})](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmax_%7B%5Cboldsymbol%7B%5Cgamma%7D%7Dl_n%28%5Cboldsymbol%7B%5Cgamma%7D%29 "\max_{\boldsymbol{\gamma}}l_n(\boldsymbol{\gamma})").
That is,
![\\hat{u}\_j](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%7Bu%7D_j "\hat{u}_j")
measures the marginal contribution of
![X_j](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;X_j "X_j")
to the survival outcome after adjusting for the effect of
![\\mathbf{Z}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cmathbf%7BZ%7D "\mathbf{Z}"),
and a large value suggests
![\\beta_j\\neq 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Cbeta_j%5Cneq%200 "\beta_j\neq 0").
We therefore propose to recruit variables

![
\\widehat{\\mathbf{X}}\_\\mathcal{A}=\\{X_j:\\hat{u}\_j \~\\text{is among the first}\~d\~\\text{largest of all}\\}
](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%0A%5Cwidehat%7B%5Cmathbf%7BX%7D%7D_%5Cmathcal%7BA%7D%3D%5C%7BX_j%3A%5Chat%7Bu%7D_j%20~%5Ctext%7Bis%20among%20the%20first%7D~d~%5Ctext%7Blargest%20of%20all%7D%5C%7D%0A "
\widehat{\mathbf{X}}_\mathcal{A}=\{X_j:\hat{u}_j ~\text{is among the first}~d~\text{largest of all}\}
")

for a pre-specified model size
![d](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;d "d").
We henceforth refer to the above procedure as conditional sure
independence screening for PSH model, or PSH-CSIS for short. According
to the sure screening property (Fan and Lv 2008), the choice of
![d](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;d "d")
can be relatively generous to ensure that all the important predictors
are preserved with high probability. Conventional choices of
![d](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;d "d")
are
![\[n/\\log(n)\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Bn%2F%5Clog%28n%29%5D "[n/\log(n)]"),
![2\[n/\\log(n)\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;2%5Bn%2F%5Clog%28n%29%5D "2[n/\log(n)]"),
![3\[n/\\log(n)\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;3%5Bn%2F%5Clog%28n%29%5D "3[n/\log(n)]"),
and
![n-p_0-1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;n-p_0-1 "n-p_0-1")
(Fan and Lv 2008; Li, Zhong, and Zhu 2012). Once the dataset is
sufficiently downsized by PSH-CSIS, existing lower dimensional methods
can be used afterwards for more precise variable selection and
statistical inference.

## Non-muscle-invasive Bladder Carcinoma Data

Gene expression data and clinical data for patients with
non-muscle-invasive bladder carcinoma were acquired from GEO database
GSE5479 (accession number GSE5479, Dyrskjøt et al. 2007). In total, 300
patients with complete information on 1,381 microarray features and 5
important clinical covariates (age, sex, reevaluated WHO grade,
reevaluated pathological disease stage and BCG/MMC treatment) were
included for analysis. The primary endpoint is the time to progression.
Progression or death from bladder cancer, which is the event of
interest, was observed for 83 patients. Besides, 33 patients died from
other or unknown causes and 184 patients were censored during follow-up.
The progression-free survival time ranges from 0 to 185 months, with a
median of 47 months. Gene expressions were standardized prior to
analysis.

``` r
### load the cleaned dataset and ESIS functions
blca <- read.table(file = "BLCA.csv", header = T)
```

## Statistical Analysis

We first divide the patients into a training subgroup and a testing
subgroup in a ratio of 4:1 such that the two cohorts share similar
clinical characteristics.

``` r
### model matrix
blca$SEX <- as.factor(blca$SEX)
blca$treatment <- as.factor(blca$treatment)
blca$grade <- as.factor(blca$grade)
blca$pT <- as.factor(blca$pT)
blca1 <- data.frame(model.matrix(~., blca))[,-1]

### training - testing split
set.seed(123)
library(SPlit)
ttind <- SPlit(cbind(blca1[,1:6], as.factor(blca1[,7])), splitRatio = 0.2)
```

    ## Optimizing <1 thread> [                    ] 1%Optimizing <1 thread> [                    ] 2%Optimizing <1 thread> [                    ] 3%Optimizing <1 thread> [                    ] 4%Optimizing <1 thread> [+                   ] 5%Optimizing <1 thread> [+                   ] 6%Optimizing <1 thread> [+                   ] 7%Optimizing <1 thread> [+                   ] 8%Optimizing <1 thread> [+                   ] 9%Optimizing <1 thread> [++                  ] 10%Optimizing <1 thread> [++                  ] 11%Optimizing <1 thread> [++                  ] 12%Optimizing <1 thread> [++                  ] 13%Optimizing <1 thread> [++                  ] 14%Optimizing <1 thread> [+++                 ] 15%Optimizing <1 thread> [+++                 ] 16%Optimizing <1 thread> [+++                 ] 17%Optimizing <1 thread> [+++                 ] 18%Optimizing <1 thread> [+++                 ] 19%Optimizing <1 thread> [++++                ] 20%Optimizing <1 thread> [++++                ] 21%Optimizing <1 thread> [++++                ] 22%Optimizing <1 thread> [++++                ] 23%Optimizing <1 thread> [++++                ] 24%Optimizing <1 thread> [+++++               ] 25%Optimizing <1 thread> [+++++               ] 26%Optimizing <1 thread> [+++++               ] 27%Optimizing <1 thread> [+++++               ] 28%Optimizing <1 thread> [+++++               ] 29%Optimizing <1 thread> [++++++              ] 30%Optimizing <1 thread> [++++++              ] 31%Optimizing <1 thread> [++++++              ] 32%Optimizing <1 thread> [++++++              ] 33%Optimizing <1 thread> [++++++              ] 34%Optimizing <1 thread> [+++++++             ] 35%Optimizing <1 thread> [+++++++             ] 36%Optimizing <1 thread> [+++++++             ] 37%Optimizing <1 thread> [+++++++             ] 38%Optimizing <1 thread> [+++++++             ] 39%Optimizing <1 thread> [++++++++            ] 40%Optimizing <1 thread> [++++++++            ] 41%Optimizing <1 thread> [++++++++            ] 42%Optimizing <1 thread> [++++++++            ] 43%Optimizing <1 thread> [++++++++            ] 44%Optimizing <1 thread> [+++++++++           ] 45%Optimizing <1 thread> [+++++++++           ] 46%Optimizing <1 thread> [+++++++++           ] 47%Optimizing <1 thread> [+++++++++           ] 48%Optimizing <1 thread> [+++++++++           ] 49%Optimizing <1 thread> [++++++++++          ] 50%Optimizing <1 thread> [++++++++++          ] 51%Optimizing <1 thread> [++++++++++          ] 52%Optimizing <1 thread> [++++++++++          ] 53%Optimizing <1 thread> [++++++++++          ] 54%Optimizing <1 thread> [+++++++++++         ] 55%Optimizing <1 thread> [+++++++++++         ] 56%Optimizing <1 thread> [+++++++++++         ] 57%Optimizing <1 thread> [+++++++++++         ] 58%Optimizing <1 thread> [+++++++++++         ] 59%Optimizing <1 thread> [++++++++++++        ] 60%Optimizing <1 thread> [++++++++++++        ] 61%Optimizing <1 thread> [++++++++++++        ] 62%Optimizing <1 thread> [++++++++++++        ] 63%Optimizing <1 thread> [++++++++++++        ] 64%Optimizing <1 thread> [+++++++++++++       ] 65%Optimizing <1 thread> [+++++++++++++       ] 66%Optimizing <1 thread> [+++++++++++++       ] 67%Optimizing <1 thread> [+++++++++++++       ] 68%Optimizing <1 thread> [+++++++++++++       ] 69%Optimizing <1 thread> [++++++++++++++      ] 70%Optimizing <1 thread> [++++++++++++++      ] 71%Optimizing <1 thread> [++++++++++++++      ] 72%Optimizing <1 thread> [++++++++++++++      ] 73%Optimizing <1 thread> [++++++++++++++      ] 74%Optimizing <1 thread> [+++++++++++++++     ] 75%Optimizing <1 thread> [+++++++++++++++     ] 76%Optimizing <1 thread> [+++++++++++++++     ] 77%Optimizing <1 thread> [+++++++++++++++     ] 78%Optimizing <1 thread> [+++++++++++++++     ] 79%Optimizing <1 thread> [++++++++++++++++    ] 80%Optimizing <1 thread> [++++++++++++++++    ] 81%Optimizing <1 thread> [++++++++++++++++    ] 82%Optimizing <1 thread> [++++++++++++++++    ] 83%Optimizing <1 thread> [++++++++++++++++    ] 84%Optimizing <1 thread> [+++++++++++++++++   ] 85%Optimizing <1 thread> [+++++++++++++++++   ] 86%Optimizing <1 thread> [+++++++++++++++++   ] 87%Optimizing <1 thread> [+++++++++++++++++   ] 88%Optimizing <1 thread> [+++++++++++++++++   ] 89%Optimizing <1 thread> [++++++++++++++++++  ] 90%Optimizing <1 thread> [++++++++++++++++++  ] 91%Optimizing <1 thread> [++++++++++++++++++  ] 92%Optimizing <1 thread> [++++++++++++++++++  ] 93%Optimizing <1 thread> [++++++++++++++++++  ] 94%Optimizing <1 thread> [+++++++++++++++++++ ] 95%Optimizing <1 thread> [+++++++++++++++++++ ] 96%Optimizing <1 thread> [+++++++++++++++++++ ] 97%Optimizing <1 thread> [+++++++++++++++++++ ] 98%Optimizing <1 thread> [+++++++++++++++++++ ] 99%Optimizing <1 thread> [++++++++++++++++++++] 100%

``` r
Xtr <- blca1[-ttind,]
Xte <- blca1[ttind,]
```

In the next we perform PSH-CSIS on the training subset and pre-selected
![240-1-5=234](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;240-1-5%3D234 "240-1-5=234")
genes after adjusting for the effect of the clinical covariates.

``` r
library(pec)
library(riskRegression)

p <- ncol(blca1[,-(1:7)]) # no. of genes
w <- numeric(p)
for(j in 1:p){
  fm <- as.formula(paste("Hist(time,cause)~",
                         paste(colnames(blca1)[c(1:5,7+j)], collapse= "+")))
  mod <- FGR(fm, cause=1, data=Xtr) # individual PSH model
  w[j] <- mod$crrFit$loglik
}
ow <- order(-w)

d <- 234
ix <- c(1:5,7+ow[1:d]) # variables for further selection
```

The likelihood-based boosting approach (CoxBoost, Binder et al. 2009)
can then be applied to the reduced training data for further gene
selection and prognostic modeling simultaneously through the R package .
The clinical covariates remain unpenalized in the boosting procedure and
the optimal tuning parameters can be determined through cross
validation.

``` r
library(CoxBoost)
set.seed(9753)
optim.res <- optimCoxBoostPenalty(Xtr$time, Xtr$cause, as.matrix(Xtr[,ix]),
                                  unpen.index=1:5, standardize=F,
                                  minstepno=50, maxstepno=200,
                                  stepsize.factor=1, start.penalty=sum(Xtr$cause==1))
fm.csis <- as.formula(paste("Hist(time,cause)~",
                            paste(colnames(blca1)[ix], collapse= "+")))
fit.cb <- coxboost(fm.csis, data=Xtr, cause=1, cv=F,
                   CoxBoost.unpen.index=1:5, CoxBoost.standardize=F,
                   CoxBoost.stepsize.factor=1, CoxBoost.penalty=optim.res$penalty,
                   CoxBoost.stepno=optim.res$cv.res$optimal.step)
```

A patient’s risk score is defined as the linear combination of the
selected genes where the coefficients are extracted from the above
fitted model (i.e. the gene signature values). The risk score is used to
classify the patient as having high or low risk, with the median score
of the training group being the cutoff. The same cutoff value is also
applied when assigning the test samples. The two risk groups are
contrasted by cumulative incidence analysis.

``` r
coef.cb = fit.cb$coxboost$coefficients[fit.cb$stepno+1,]
fitted.lp <- drop(as.matrix(Xtr[,ix][,-(1:5)])%*%matrix(coef.cb[-(1:5)])) # training risk scores 
pred.lp <- drop(as.matrix(Xte[,ix][,-(1:5)])%*%matrix(coef.cb[-(1:5)])) # testing risk scores
risk.test <- pred.lp>median(fitted.lp)
library(cmprsk)
ciftest <- cuminc(Xte$time, Xte$cause, risk.test)
ciftest
```

    ## Tests:
    ##         stat         pv df
    ## 1 4.53265651 0.03325396  1
    ## 2 0.02744928 0.86841018  1
    ## Estimates and Variances:
    ## $est
    ##                 20         40        60        80       100       120
    ## FALSE 1 0.06451613 0.06451613 0.1396051 0.1396051 0.1396051        NA
    ## TRUE 1  0.33142857 0.33142857 0.4290612 0.4290612 0.4290612 0.4290612
    ## FALSE 2 0.00000000 0.03598015 0.1095048 0.1095048 0.1095048        NA
    ## TRUE 2  0.07428571 0.07428571 0.1363673 0.1363673 0.1363673 0.1363673
    ## 
    ## $var
    ##                  20          40          60          80         100         120
    ## FALSE 1 0.002012949 0.002012949 0.004412615 0.004412615 0.004412615          NA
    ## TRUE 1  0.008593248 0.008593248 0.010839726 0.010839726 0.010839726 0.010839726
    ## FALSE 2 0.000000000 0.001297549 0.003719307 0.003719307 0.003719307          NA
    ## TRUE 2  0.002678058 0.002678058 0.006343877 0.006343877 0.006343877 0.006343877

The cumulative incidence analysis revealed that the gene signature
identified by the PSH-CSIS+CoxBoost model provides effective risk
stratification (p-value = 0.033 for the testing cohort).

Lastly, we fit a PSH model to the entire dataset to make inference about
independent prognostic factors associated with progression, where the
gene signature identified by the PSH-CSIS+CoxBoost model and the five
clinical covariates are used.

``` r
blca.final <- data.frame(blca1[,1:7],
                         gene = drop(as.matrix(blca1[,ix][,-(1:5)])%*%matrix(coef.cb[-(1:5)])))
fm.final <- as.formula(paste("Hist(time,cause)~",
                             paste(colnames(blca.final)[c(1:5,8)], collapse= "+")))
fit.final <- FGR(fm.final,cause=1,data=blca.final)
fit.final
```

    ## 
    ## Right-censored response of a competing.risks model
    ## 
    ## No.Observations: 300 
    ## 
    ## Pattern:
    ##          
    ## Cause     event right.censored
    ##   1          83              0
    ##   2          33              0
    ##   unknown     0            184
    ## 
    ## 
    ## Fine-Gray model: analysis of cause 1 
    ## 
    ## Competing Risks Regression
    ## 
    ## Call:
    ## FGR(formula = fm.final, data = blca.final, cause = 1)
    ## 
    ##                  coef exp(coef) se(coef)      z p-value
    ## AGE            0.0285     1.029   0.0123  2.317 2.1e-02
    ## SEXM          -0.1313     0.877   0.3007 -0.437 6.6e-01
    ## treatmentnone  0.9712     2.641   0.3037  3.198 1.4e-03
    ## gradeLOW      -0.8486     0.428   0.2960 -2.867 4.1e-03
    ## pTTa           0.5331     1.704   0.2726  1.955 5.1e-02
    ## gene           2.5087    12.289   0.3604  6.961 3.4e-12
    ## 
    ##               exp(coef) exp(-coef)  2.5%  97.5%
    ## AGE               1.029     0.9719 1.004  1.054
    ## SEXM              0.877     1.1403 0.486  1.581
    ## treatmentnone     2.641     0.3786 1.457  4.789
    ## gradeLOW          0.428     2.3363 0.240  0.765
    ## pTTa              1.704     0.5868 0.999  2.908
    ## gene             12.289     0.0814 6.064 24.905
    ## 
    ## Num. cases = 300
    ## Pseudo Log-likelihood = -404 
    ## Pseudo likelihood ratio test = 78.7  on 6 df,
    ## 
    ## Convergence: TRUE

The 6-gene signature selected by PSH-CSIS+CoxBoost was a strong
predictor with an hazard ratio of 10.74 (p-value
![\<](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%3C "<")
0.001), adjusted for other clinical covariates. There was a 3% increase
in the expected hazard relative to a one year increase in age
(p-value=0.021). Patients with low grade tumor experienced reduction of
hazard by 57% (p-value=0.004), compared to those with high grade tumor.
Receiving BCG/MMC treatment reduced the hazard of progression by 62%
(p-value=0.001).

## Reference

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-binder2009boosting" class="csl-entry">

Binder, Harald, Arthur Allignol, Martin Schumacher, and Jan Beyersmann.
2009. “Boosting for High-Dimensional Time-to-Event Data with Competing
Risks.” *Bioinformatics* 25 (7): 890–96.

</div>

<div id="ref-dyrskjot2007gene" class="csl-entry">

Dyrskjøt, Lars, Karsten Zieger, Francisco X Real, Núria Malats, Alfredo
Carrato, Carolyn Hurst, Sanjeev Kotwal, et al. 2007. “Gene Expression
Signatures Predict Outcome in Non–Muscle-Invasive Bladder Carcinoma: A
Multicenter Validation Study.” *Clinical Cancer Research* 13 (12):
3545–51.

</div>

<div id="ref-fan2008sure" class="csl-entry">

Fan, Jianqing, and Jinchi Lv. 2008. “Sure Independence Screening for
Ultrahigh Dimensional Feature Space.” *Journal of the Royal Statistical
Society: Series B (Statistical Methodology)* 70 (5): 849–911.

</div>

<div id="ref-fine1999proportional" class="csl-entry">

Fine, Jason P, and Robert J Gray. 1999. “A Proportional Hazards Model
for the Subdistribution of a Competing Risk.” *Journal of the American
Statistical Association* 94 (446): 496–509.

</div>

<div id="ref-gray1988class" class="csl-entry">

Gray, Robert J. 1988. “A Class of k-Sample Tests for Comparing the
Cumulative Incidence of a Competing Risk.” *The Annals of Statistics*,
1141–54.

</div>

<div id="ref-li2012feature" class="csl-entry">

Li, Runze, Wei Zhong, and Liping Zhu. 2012. “Feature Screening via
Distance Correlation Learning.” *Journal of the American Statistical
Association* 107 (499): 1129–39.

</div>

</div>
