# ivxPredictive 

## Unified Inference with General Autoregressive Roots

### Description 

The R package 'ivxPredictive' (under development) implements robust econometric inference methodologies for predictive regression models with autoregressive roots across the spectrum of stationarity and nonstationarity, with persistence types as defined by Magdalinos and Phillips (2020). In particular, the R implementation of the 'ivxPredictive' R package extends the original IVX instrumentation of [Phillips and Magdalinos (2009)](https://ideas.repec.org/p/skb/wpaper/cofie-06-2009.html) and [Kostakis, Magdalinos and Stamatogiannis (2015)](https://academic.oup.com/rfs/article/28/5/1506/1867633?login=true), to general autoregressive roots for parametric functional forms. Specifically, here we focus on linear specifications for the predictive regression model and implement the proposed procedure for the case of a conditional mean loss function (see, Magdalinos and Petrova (2022)) as well as a conditional quantile loss function (see, [Katsouris, C. (2022)](https://arxiv.org/abs/2204.02073)).  

<p align="center">
  
<img src="https://github.com/christiskatsouris/ivxPredictive/blob/main/data/persistence.jpg" width="460"/>

</p>  

We consider the following persistence classes:
   
$\textbf{(P1). nearly stable}$ processes: 

$$\text{if} \ \ \left( \theta_n \right)_{ n \in \mathbb{N} } \ \ \text{is such that} \ \zeta = - \infty \ \ \text{and it holds that} \ \ \theta_n \to | \theta | < 1.$$
    
$\textbf{(P2). nearly unstable}$ processes:   

$$\text{if} \ \ \left( \theta_n \right)_{ n \in \mathbb{N} } \ \ \text{is such that} \ \ \zeta \equiv c \in \mathbb{R} \ \ \text{and it holds that} \ \ \theta_n \to \theta = 1.$$

    
$\textbf{(P3). nearly explosive}$ processes:   

$$\text{if} \ \ \left( \theta_n \right)_{ n \in \mathbb{N} } \ \ \text{is such that} \ \ \zeta = + \infty \ \ \text{and it holds that} \ \ \theta_n \to | \theta | > 1.$$

  
### Methodology  
  
This R package implements a novel endogenous instrumentation approach based on the IVX estimator examined by [Phillips and Magdalinos (2009)](https://ideas.repec.org/p/skb/wpaper/cofie-06-2009.html) and [Kostakis, Magdalinos and Stamatogiannis (2015)](https://academic.oup.com/rfs/article/28/5/1506/1867633?login=true). The current procedure has a similar construction to the IV instrumentation proposed in the recent working paper of Magdalinos and Petrova (2022), with the aim to provide uniform and robust to the nuisance parameter of persistence inference across the spectrum of stationary and nonstationary roots, specifically for quantile autoregressive processes. We call this variant of the original IVX estimator, IVX-P, which can be employed to both conditional mean and conditional quantile functional forms when the model includes either univariate or multivariate regressors. The novelty of the IVX-P estimator is that is a 'hybrid estimator' which in contrast to the classical least squares estimator has desirable asymptotic theory properties and is constructed based on the underline nonstationary stochastic processes using information both within the admissible parameter space as well as outside the usual parameter space.  
  
## Installation (under development) 

The R package ‘ivxPredictive’ will be able to be installed from Github.

## Usage 

```R

# After development the package will be able to be installed using
install.packages("ivxPredictive")
library("ivxPredictive")

```
## Notes

1. The implementation of the original (current) IVX instrumentation method proposed by [Phillips and Magdalinos (2009)](https://ideas.repec.org/p/skb/wpaper/cofie-06-2009.html) along with related statistical inference and empirical applications has been prepared as an R package [ivx](https://github.com/kvasilopoulos/ivx) developed by [Kostas Vasilopoulos](https://github.com/kvasilopoulos).

2. The Matlab code of the original IVX instrumentation method can be found on the website of [Michalis Stamatogiannis](https://sites.google.com/site/mpstamatogiannis/home).

3. The author would like to thank [Ji Hyung Lee](https://economics.illinois.edu/profile/jihyung) for kindly sharing the Matlab code for the implementation of the original IVX instrumentation method for the quantile predictive regression. 

## Key References:

- Katsouris, C. (2022b). "Asymptotic Theory for Moderate Deviations from the Unit Boundary in Quantile Autoregression and Predictive Regression Models". University of Southampton, Working paper.  

- Katsouris, C. (2022a). "Asymptotic Theory for Moderate Deviations from the Unit Boundary in Quantile Autoregressive Time Series". arXiv preprint [arXiv:2204.02073](https://arxiv.org/abs/2204.02073).

- Katsouris, C. (2021d). "Testing for Structural Breaks in Predictive Regression Models". University of Southampton, Working paper.  

- Katsouris, C. (2021e). "Bootstrapping Nonstationary Autoregressive Processes in Predictive Regression". University of Southampton, Working paper.   

- Kostakis, A., Magdalinos, T., & Stamatogiannis, M. P. (2015). "Robust econometric inference for stock return predictability". The Review of Financial Studies, 28(5), 1506-1553.

- Lee, J. H. (2016). "Predictive quantile regression with persistent covariates: IVX-QR approach". Journal of Econometrics, 192(1), 105-118.

- Fan, R., & Lee, J. H. (2019). Predictive quantile regressions under persistence and conditional heteroskedasticity. Journal of Econometrics, 213(1), 261-280.

- Magdalinos, T. (2021). Least squares and IVX limit theory in systems of predictive regressions with GARCH innovations. Econometric Theory, 1-38.

- Magdalinos, T., & Petrova, K. (2022). "Uniform and distribution-free inference with general autoregressive processes". University of Southampton, Working paper. 

- Magdalinos, T., & Phillips, P. C. B. (2020). "Econometric inference in matrix vicinities of unity and stationarity". University of Southampton, Working paper.  

- Magdalinos, T., & Phillips, P. C. B. (2009). "Limit theory for cointegrated systems with moderately integrated and moderately explosive regressors". Econometric Theory, 25(2), 482-526.

- Phillips, P. C. B., & Magdalinos, T. (2009). "Econometric inference in the vicinity of unity". Singapore Management University, CoFie Working paper, 7.

- Phillips, P. C. B., & Magdalinos, T. (2007). Limit theory for moderate deviations from a unit root. Journal of Econometrics, 136(1), 115-130.

## References on IVX Variants:

- Cai, Z., Chen, H., & Liao, X. (2022). A new robust inference for predictive quantile regression. Journal of Econometrics.
- Demetrescu, M., Georgiev, I., Rodrigues, P. M., & Taylor, A. R. (2022). Extensions to IVX methods of inference for return predictability. Journal of Econometrics.
- Demetrescu, M., & Rodrigues, P. M. (2020). Residual-augmented IVX predictive regression. Journal of Econometrics.
- Hosseinkouchack, M., & Demetrescu, M. (2021). Finite-sample size control of IVX-based tests in predictive regressions. Econometric Theory, 37(4), 769-793.
- Liu, Y., & Phillips, P. C. (2022). Robust inference with stochastic local unit root regressors in predictive regressions. Journal of Econometrics.
- Ren, Y., Tu, Y., & Yi, Y. (2019). Balanced predictive regressions. Journal of Empirical Finance, 54, 118-142.
- Yang, B., Long, W., & Yang, Z. (2022). Testing predictability of stock returns under possible bubbles (IVX-BUB). Journal of Empirical Finance.
- Yang, B., Long, W., Peng, L., & Cai, Z. (2020). Testing the predictability of US housing price index returns based on an IVX-AR model (IVX-AR). Journal of the American Statistical Association, 115(532), 1598-1619.

## Seminal Studies:

### 1950's

- Anderson, T. W. (1959). On asymptotic distributions of estimates of parameters of stochastic difference equations. The Annals of Mathematical Statistics, 676-687.
- Cramér, H. (1951). A contribution to the theory of stochastic processes. In Proceedings of the Second Berkeley Symposium on Mathematical Statistics and Probability (Vol. 2, pp. 329-340). University of California Press.
- Mann, H. B., & Wald, A. (1943). On the statistical treatment of linear stochastic difference equations. Econometrica, Journal of the Econometric Society, 173-220.
- White, J. S. (1959). The limiting distribution of the serial correlation coefficient in the explosive case II. The Annals of Mathematical Statistics, 831-834.
- White, J. S. (1958). The limiting distribution of the serial correlation coefficient in the explosive case. The Annals of Mathematical Statistics, 1188-1197.

### 1980's

- Chan, N. H., & Wei, C. Z. (1987). Asymptotic inference for nearly nonstationary AR (1) processes. The Annals of Statistics, 1050-1063.
- Chan, N. H. (1988). The parameter inference for nearly nonstationary time series. Journal of the American Statistical Association, 83(403), 857-862.
- Dickey, D. A., & Fuller, W. A. (1979). Distribution of the estimators for autoregressive time series with a unit root. Journal of the American statistical association, 74(366a), 427-431.
- Dickey, D. A., & Fuller, W. A. (1981). Likelihood ratio statistics for autoregressive time series with a unit root. Econometrica: journal of the Econometric Society, 1057-1072.
- Phillips, P. C. B. (1987). Towards a unified asymptotic theory for autoregression. Biometrika, 74(3), 535-547.
- Phillips, P. C. B. & Perron, P. (1988). Testing for a unit root in time series regression. Biometrika, 75(2), 335-346.
- Phillips, P. C. B. (1988). Regression theory for near-integrated time series. Econometrica: Journal of the Econometric Society, 1021-1043.

### 1990's

- Cavanagh, C. L., Elliott, G., & Stock, J. H. (1995). Inference in models with nearly integrated regressors. Econometric theory, 11(5), 1131-1147.

### 2000's
- Aue, A., & Horváth, L. (2007). A limit theorem for mildly explosive autoregression with stable errors. Econometric Theory, 23(2), 201-220.
- Jansson, M., & Moreira, M. J. (2006). Optimal inference in regression models with nearly integrated regressors. Econometrica, 74(3), 681-714.


## Code of Coduct

Please note that the ‘ivxPredictive’ project will be released with a Contributor Code of Coduct (under construction). By contributing to this project, you agree to abide by its terms.

## Declarations

The author declares no conflicts of interest. 

In particular, the author declares that has no affiliations with or involvement in any organization or entity with any financial interest (such as honoraria; educational grants; participation in speakers’ bureaus; membership, employment, consultancies, stock ownership, or other equity interest; and expert testimony or patent-licensing arrangements), or non-financial interest (such as personal or professional relationships, affiliations, knowledge or beliefs) in the subject matter or materials discussed in the manuscript and implemented in the R package.

# Acknowledgments

The author greatfully acknowledges financial support from Graduate Teaching Assistantships at the School of Economic, Social and Political Sciences of the University of Southampton as well as from the Vice-Chancellor's PhD Scholarship of the University of Southampton, for the duration of the academic years 2018 to 2021. The author also gratefully acknowledges previously obtained funding from Research Grants of interdisciplinary Centers of research excellence based at the University of Cyprus (UCY) as well as at the University College London (UCL).

Part of the aspects implemented in the R package 'ivxPredictive', are discussed in the PhD thesis of this R package creator (Christis G. Katsouris) titled: "Aspects of Estimation and Inference for Predictive Regression Models", completed under the supervision of  [Professor Jose Olmo](https://www.southampton.ac.uk/people/5xb29b/professor-jose-olmo) and [Professor Anastasios Magdalinos](https://www.southampton.ac.uk/people/5x7sz7/professor-anastasios-magdalinos#supervision).

$\textbf{News:}$ The author will be joining the [University of Exeter Business School](http://business-school.exeter.ac.uk/) as a Visiting Lecturer in Economics (Education and Scholarship) at the Department of Economics in September 2022.

# Historical Background

$\textbf{Charles Max Stein}$  (March 22, 1920 – November 24, 2016) was an American mathematical statistician and professor of statistics at Stanford University. 
He received his Ph.D in 1947 at Columbia University with advisor Abraham Wald. He held faculty positions at Berkeley and the University of Chicago before moving permanently to Stanford in 1953. He is known for Stein's paradox in decision theory, which shows that ordinary least squares estimates can be uniformly improved when many parameters are estimated; for Stein's lemma, giving a formula for the covariance of one random variable with the value of a function of another when the two random variables are jointly normally distributed; and for Stein's method, a way of proving theorems such as the Central Limit Theorem that does not require the variables to be independent and identically distributed (Source: Wikepedia).

## Remark: 

In contrast to the ordinary least squares estimates, the opposite result seems to hold for the original IVX estimator, that is, it seems to not uniformly be improving when many parameters are estimated. The particular aspect, is currently an open research question in the literature which worths further investigation. 
