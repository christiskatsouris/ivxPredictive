# ivxPredictive 

## Unified Inference with General Autoregressive Roots

### Description 

The R package 'ivxPredictive' (under development) implements robust econometric inference methodologies for predictive regression models with autoregressive roots across the spectrum of stationarity and nonstationarity, with persistence types as defined by Magdalinos and Phillips (2020). In particular, the 'ivxPredictive' package extends the original IVX instrumentation of [Phillips and Magdalinos (2009)](https://ideas.repec.org/p/skb/wpaper/cofie-06-2009.html) and [Kostakis, Magdalinos and Stamatogiannis (2015)](https://academic.oup.com/rfs/article/28/5/1506/1867633?login=true), to general autoregressive roots for parametric functional forms. Specifically, here we focus on linear specifications for the predictive regression model and implement the proposed procedure for the case of a conditional mean loss function (see, Magdalinos and Petrova (2022)) as well as a conditional quantile loss function (see, [Katsouris, C. (2022)](https://arxiv.org/abs/2204.02073)).  

<p align="center">
  
<img src="https://github.com/christiskatsouris/ivxPredictive/blob/main/data/persistence.jpg" width="460"/>

</p>  

We consider the following persistence classes:
   
#### (P1). nearly stable processes: 

$$\text{if} \ \ \left( \theta_n \right)_{ n \in \mathbb{N} } \ \ \text{is such that} \ \zeta = - \infty \ \ \text{and it holds that} \ \ \theta_n \to | \theta | < 1.$$
    
#### (P2). nearly unstable processes:   

$$\text{if} \ \ \left( \theta_n \right)_{ n \in \mathbb{N} } \ \ \text{is such that} \ \ \zeta \equiv c \in \mathbb{R} \ \ \text{and it holds that} \ \ \theta_n \to \theta = 1.$$

    
#### (P3). nearly explosive processes:   

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
## Notes:

1. The implementation of the original (current) IVX instrumentation method proposed by [Phillips and Magdalinos (2009)](https://ideas.repec.org/p/skb/wpaper/cofie-06-2009.html) along with statistical inference procedures and empirical applications can be found in the R package ['ivx'](https://github.com/kvasilopoulos/ivx) developed by [Kostas Vasilopoulos](https://github.com/kvasilopoulos).

2. The Matlab code of the original IVX instrumentation method can be found on the website of [Michalis Stamatogiannis](https://sites.google.com/site/mpstamatogiannis/home).

3. The author would like to thank [Ji Hyung Lee](https://economics.illinois.edu/profile/jihyung) for kindly sharing the Matlab code for the implementation of the original IVX instrumentation method for the quantile predictive regression. 

## Empirical Illustrations

### Example 1: (Autocorrelation Matrix Estimation)

Consider the predictive regression model with multiple predictors. Then, the Aols vector of parameters (not a matrix in this case, since the response is a scalar vector) and the Rols matrix of coefficients can be obtained as below

```R

lm1  <- lm(yt ~ Xlag)
Aols <- coefficients(lm1)  

Rn   <- matrix(0, p, p)
for (i in 1:p) 
{
   Rn[i, i] <- lm( Xt[, i] ~ 0 + Xlag[, i] )$coefficients
}
  
> Rn
          [,1]      [,2]       [,3]      [,4]     [,5]        [,6]     
[1,] 0.9891264 0.0000000 0.00000000 0.0000000 0.000000  0.00000000  
[2,] 0.0000000 0.8461566 0.00000000 0.0000000 0.000000  0.00000000  
[3,] 0.0000000 0.0000000 0.07760762 0.0000000 0.000000  0.00000000  
[4,] 0.0000000 0.0000000 0.00000000 0.9985889 0.000000  0.00000000  
[5,] 0.0000000 0.0000000 0.00000000 0.0000000 1.000366  0.00000000 
[6,] 0.0000000 0.0000000 0.00000000 0.0000000 0.000000  1.04066777  

```

### Remarks:

(i) The $R_n$ which is a $(p \times p)$ matrix corresponds to the autocorrelation coefficients of the nonstationary regressors of the predictive regression model. In the theoretical speciciation of the predictive regression model, these are modeled using an AR(1) model with a local-unit-root expression to capture persistence and near the unity characteristics. From the empirical example above, we see that some of these coefficients are close to the unit boundary from below and some reach the unit boundary from the explosive side. 

(ii) The autocorrelation matrix that corresponds to the autoregressive specification of regressors is estimated based on least-squares optimization. The OLS estimator of the autocorrelation matrix similar to the case of univariate autoregressions is not a consistent estimator when the local-to-unity specification is imposed. However, the particular estimator is required when estimating statistics from fitted predictive regression models.   

(iii) Another important aspect when fitting predictive regression models and especially when constructing test statistics such as testing for parameter instability in predictive regressions is to apply correctly the demeaning of random variables. Specifically, when the true model is assumed to have a non-zero model intercept then demeaning the random variables of the model can alter the asymptotic theory of estimators and test statistics depending on the underline persistence properties of regressors.  

### Example 2: (Covariance Matrix Estimation)

The estimation procedure for the covariance matrix of estimators and test statistics in predictive regression models does not follow conventional formulations as in the case of stationary time series regression models. Due to the local-to-unity specification of the autocorrelation coefficient(s), which captures the unknown 'degree of persistence' in time series regressors, the predictive regression specification permits to include both 'nonstationary predictors' and/or 'cointegrated regressors'; therefore the estimation of the covariance matrix is obtained using nonparametric kernel methodologies in order to introduce an endogeneity bias correction.   

```R

# Step 1: Estimation of the Rn autocorrelation coefficient 
for i=1:l
  rn(i,i) = regress( xt(:,i), xlag(:,i) );
end

# autoregressive residual estimation 
u = xt-xlag*rn;

# residuals' correlation matrix
corrmat = corrcoef([epshat u]);

# covariance matrix estimation (predictive regression)
covepshat = epshat'*epshat/nn;
covu = zeros(l,l);

for t=1:nn
    covu = covu+u(t,:)'*u(t,:);
end

# covariance matrix estimation (autoregression)
covu=covu/nn;
covuhat=zeros(1,l);

for i=1:l
    covuhat(1,i)=sum(epshat.*u(:,i));
end

```

### Remarks:

(i) 

(ii) 

### Example 3: (Instrumental Variables Estimation)

Firstly, we consider estimation procedure for the original IVX instrumentation of [Phillips and Magdalinos (2009)](https://ideas.repec.org/p/skb/wpaper/cofie-06-2009.html) and [Kostakis, Magdalinos and Stamatogiannis (2015)](https://academic.oup.com/rfs/article/28/5/1506/1867633?login=true) and discuss the main features. Secondly, we present the estimation procedure of the IVX-P estimator (i.e., the IVX instrumentation that covers the whole spectrum of persistence classes in a unified way) and discuss the main features.  

``` Matlab

# Estimation of the original IVX instruments 





```

### Remarks:

Notice that various Software packages now exist which can be employed for monitoring time series data for explosive bubbles, although the main aim of developing the R package 'ivxPredictive' is for estimation and inference purposes, we briefly mention some of them here. Specifically, the PWY (2011) and PSY (2015) statistical tests for monitoring explosive behaviour are implemented both in MATLAB as well as in Eviews. Furthermore, the R package ['exuber'](https://github.com/kvasilopoulos/exuber)   


## Key References:

- Katsouris, C. (2022b). "Asymptotic Theory for Moderate Deviations from the Unit Boundary in Quantile Autoregression and Predictive Regression Models". University of Southampton, Working paper.  
- Katsouris, C. (2022a). "Asymptotic Theory for Moderate Deviations from the Unit Boundary in Quantile Autoregressive Time Series". arXiv preprint [arXiv:2204.02073](https://arxiv.org/abs/2204.02073).
- Katsouris, C. (2021d). "Testing for Structural Breaks in Predictive Regression Models". University of Southampton, Working paper.  
- Katsouris, C. (2021e). "Bootstrapping Nonstationary Autoregressive Processes in Predictive Regression". University of Southampton, Working paper.   
- Kostakis, A., Magdalinos, T., & Stamatogiannis, M. P. (2015). "Robust econometric inference for stock return predictability". The Review of Financial Studies, 28(5), 1506-1553.
- Lee, J. H. (2016). "Predictive quantile regression with persistent covariates: IVX-QR approach". Journal of Econometrics, 192(1), 105-118.
- Fan, R., & Lee, J. H. (2019). Predictive quantile regressions under persistence and conditional heteroskedasticity. Journal of Econometrics, 213(1), 261-280.
- Magdalinos, T. (2016). Least squares and IVX limit theory in systems of predictive regressions with GARCH innovations. Econometric Theory, 1-38.
- Magdalinos, T., & Petrova, K. (2022). "Uniform and distribution-free inference with general autoregressive processes". University of Southampton, Working paper. 
- Magdalinos, T., & Phillips, P. C. B. (2020). "Econometric inference in matrix vicinities of unity and stationarity". University of Southampton, Working paper.  
- Magdalinos, T., & Phillips, P. C. B. (2009). "Limit theory for cointegrated systems with moderately integrated and moderately explosive regressors". Econometric Theory, 25(2), 482-526.
- Phillips, P. C. B., & Magdalinos, T. (2009). "Econometric inference in the vicinity of unity". Singapore Management University, CoFie Working paper, 7.
- Phillips, P. C. B., & Magdalinos, T. (2008). Limit theory for explosively cointegrated systems. Econometric Theory, 24(4), 865-887.
- Phillips, P. C. B., & Magdalinos, T. (2007). Limit theory for moderate deviations from a unit root. Journal of Econometrics, 136(1), 115-130.

## References on IVX Variants:

- Christou, C., Gupta, R., Hassapis, C., & Suleman, T. (2018). The role of economic uncertainty in forecasting exchange rate returns and realized volatility: Evidence from quantile predictive regressions. Journal of Forecasting, 37(7), 705-719.
- Cai, Z., Chen, H., & Liao, X. (2022). A new robust inference for predictive quantile regression. Journal of Econometrics.
- Cai, Z., & Chang, S. Y. (2020). A new test on asset return predictability with structural breaks. University of Kansas. Working paper. 
- Demetrescu, M., Georgiev, I., Rodrigues, P. M., & Taylor, A. R. (2022). Extensions to IVX methods of inference for return predictability. Journal of Econometrics.
- Demetrescu, M., & Rodrigues, P. M. (2020). Residual-augmented IVX predictive regression. Journal of Econometrics.
- Hosseinkouchack, M., & Demetrescu, M. (2021). Finite-sample size control of IVX-based tests in predictive regressions. Econometric Theory, 37(4), 769-793.
- Ren, Y., Tu, Y., & Yi, Y. (2019). Balanced predictive regressions. Journal of Empirical Finance, 54, 118-142.
- Xu, K. L. (2020). Testing for multiple-horizon predictability: Direct regression based versus implication based. The Review of Financial Studies, 33(9), 4403-4443.
- Yang, B., Long, W., & Yang, Z. (2022). Testing predictability of stock returns under possible bubbles (IVX-BUB). Journal of Empirical Finance.
- Yang, B., Long, W., Peng, L., & Cai, Z. (2020). Testing the predictability of US housing price index returns based on an IVX-AR model (IVX-AR). Journal of the American Statistical Association, 115(532), 1598-1619.

## Seminal Studies:

### 1930's

- Slutzky, E. (1937). The summation of random causes as the source of cyclic processes. Econometrica: Journal of the Econometric Society, 105-146.

### 1940's

- Anderson, R. L. (1942). Distribution of the serial correlation coefficient. The Annals of Mathematical Statistics, 13(1), 1-13.
- Haavelmo, T. (1943). The statistical implications of a system of simultaneous equations. Econometrica, Journal of the Econometric Society, 1-12.
- Mann, H. B., & Wald, A. (1943). On the statistical treatment of linear stochastic difference equations. Econometrica, Journal of the Econometric Society, 173-220.
- Kac, M., & Erdos, P. (1946). On certain limit theorems of the theory of probability. Bull. Amer. Math. Soc, 52, 292-302.
- Wald, A. (1949). Note on the consistency of the maximum likelihood estimate. The Annals of Mathematical Statistics, 20(4), 595-601.

### 1950's

- Anderson, T. W. (1959). On asymptotic distributions of estimates of parameters of stochastic difference equations. The Annals of Mathematical Statistics, 676-687.
- Cramér, H. (1951). A contribution to the theory of stochastic processes. In Proceedings of the Second Berkeley Symposium on Mathematical Statistics and Probability (Vol. 2, pp. 329-340). University of California Press.
- Donsker, M. D. (1951). An invariance principle for certain probability limit theorems. Mem. Amer. Math. Soc. 6 1-12.
- Grenander, U. (1954). On the estimation of regression coefficients in the case of an autocorrelated disturbance. The Annals of Mathematical Statistics, 25(2), 252-272.
- Rubin, H. (1950). Consistency of maximum likelihood estimates in the explosive case. Statistical Inference in Dynamic Economic Models, 356-364.
- White, J. S. (1958). The limiting distribution of the serial correlation coefficient in the explosive case. The Annals of Mathematical Statistics, 1188-1197.
- White, J. S. (1959). The limiting distribution of the serial correlation coefficient in the explosive case II. The Annals of Mathematical Statistics, 831-834.


### 1960's

- Mandelbrot, B. (1963). The Variation of Certain Speculative Prices. The Journal of Business, 36(4), 394-419.
- Rao, M. M. (1961). Consistency and limit distributions of estimators of parameters in explosive stochastic difference equations. The Annals of Mathematical Statistics, 32(1), 195-218.
- Rao, M. M. (1966). Inference in stochastic processes. II. Zeitschrift für Wahrscheinlichkeitstheorie und Verwandte Gebiete, 5(4), 317-335.
- Hee, O. (1966). Tests for Predictability of Statistical Models. Journal of Farm Economics, 48(5), 1479-1484.

### 1970's

- Dickey, D. A., & Fuller, W. A. (1979). Distribution of the estimators for autoregressive time series with a unit root. Journal of the American statistical association, 74(366a), 427-431.
- Rao, M. M. (1978). Asymptotic distribution of an estimator of the boundary parameter of an unstable process. The Annals of Statistics, 185-190.

### 1980's

- Cavanagh, C. (1985). Roots local to unity. Manuscript, Harvard University, Department of Economics.
- Cumberland, W. G., & Sykes, Z. M. (1982). Weak convergence of an autoregressive process used in modeling population growth. Journal of Applied Probability, 19(2), 450-455.
- Chan, N. H., & Wei, C. Z. (1987). Asymptotic inference for nearly nonstationary AR (1) processes. The Annals of Statistics, 1050-1063.
- Chan, N. H. (1988). The parameter inference for nearly nonstationary time series. Journal of the American Statistical Association, 83(403), 857-862.
- Dickey, D. A., & Fuller, W. A. (1981). Likelihood ratio statistics for autoregressive time series with a unit root. Econometrica: journal of the Econometric Society, 1057-1072.
- Lai, T. L., & Wei, C. Z. (1982). Least squares estimates in stochastic regression models with applications to identification and control of dynamic systems. The Annals of Statistics, 10(1), 154-166.
- Phillips, P. C. B. (1987). Towards a unified asymptotic theory for autoregression. Biometrika, 74(3), 535-547.
- Phillips, P. C. B. (1988). Regression theory for near-integrated time series. Econometrica: Journal of the Econometric Society, 1021-1043.
- Phillips, P. C. B. & Perron, P. (1988). Testing for a unit root in time series regression. Biometrika, 75(2), 335-346.

### 1990's

- Chan, N. H. (1990). Inference for near-integrated time series with infinite variance. Journal of the American Statistical Association, 85(412), 1069-1074.
- Cavanagh, C. L., Elliott, G., & Stock, J. H. (1995). Inference in models with nearly integrated regressors. Econometric theory, 11(5), 1131-1147.
- Elliott, G., Rothenberg, T. J., & Stock, J. H. (1996). Efficient tests for an autoregressive unit root. Econometrica, 64 (4), 813-836.
- Jeganathan, P. (1991). On the asymptotic behavior of least-squares estimators in AR time series with roots near the unit circle. Econometric Theory, 7(3), 269-306.
- Phillips, P. C. B. (1991). Optimal inference in cointegrated systems. Econometrica: Journal of the Econometric Society, 283-306.
- Phillips, P. C. B. (1995). Fully modified least squares and vector autoregression. Econometrica: Journal of the Econometric Society, 1023-1078.

### 2000's
- Aue, A., & Horváth, L. (2007). A limit theorem for mildly explosive autoregression with stable errors. Econometric Theory, 23(2), 201-220.
- Buchmann, B., & Chan, N. H. (2007). Asymptotic theory of least squares estimators for nearly unstable processes under strong dependence. The Annals of Statistics, 35(5), 2001-2017.
- Campbell, J. Y., & Yogo, M. (2006). Efficient tests of stock return predictability. Journal of financial economics, 81(1), 27-60.
- Jansson, M., & Moreira, M. J. (2006). Optimal inference in regression models with nearly integrated regressors. Econometrica, 74(3), 681-714.
- Mikusheva, A. (2007). Uniform inference in autoregressive models. Econometrica, 75(5), 1411-1452.
- Müller, U. K., & Elliott, G. (2003). Tests for unit roots and the initial condition. Econometrica, 71(4), 1269-1286.
- Phillips, P. C. B., Moon, H. R., & Xiao, Z. (2001). How to estimate autoregressive roots near unity. Econometric Theory, 17(1), 29-69.

### Remark: 

In contrast to the ordinary least squares estimates, the opposite result seems to hold for the original IVX estimator, that is, it seems to not uniformly be improving when many parameters are estimated. The particular aspect, is currently an open question in the literature which worths further investigation. Some studies in this direction include: 

- Lee, J. H., Shi, Z., & Gao, Z. (2022). On LASSO for predictive regression. Journal of Econometrics, 229(2), 322-349.
- Ke-Li Xu & Junjie Guo (2022). A New Test for Multiple Predictive Regression. Journal of Financial Econometrics.

# Code of Coduct

Please note that the ‘ivxPredictive’ project will be released with a Contributor Code of Coduct (under construction). By contributing to this project, you agree to abide by its terms.

# Declarations

The author (Christis G. Katsouris) declares no conflicts of interest. 

In particular, the author declares that has no affiliations with or involvement in any organization or entity with any financial interest (such as honoraria; educational grants; participation in speakers’ bureaus; membership, employment, consultancies, stock ownership, or other equity interest; and expert testimony or patent-licensing arrangements), or non-financial interest (such as personal relationships) in the subject matter or materials discussed in the manuscript and implemented in the R package.

# Acknowledgments

The author greatfully acknowledges financial support from Graduate Teaching Assistantships at the School of Economic, Social and Political Sciences of the University of Southampton as well as from the Vice-Chancellor's PhD Scholarship of the University of Southampton, for the duration of the academic years 2018 to 2021. The author also gratefully acknowledges previously obtained funding from Research Grants of interdisciplinary Centers of research excellence based at the University of Cyprus (UCY) as well as at the University College London (UCL).

Part of the aspects implemented in the R package 'ivxPredictive', are discussed in the PhD thesis of the author (Christis G. Katsouris) titled: "Aspects of Estimation and Inference for Predictive Regression Models", completed under the supervision of  [Professor Jose Olmo](https://www.southampton.ac.uk/people/5xb29b/professor-jose-olmo) and [Professor Anastasios Magdalinos](https://www.southampton.ac.uk/people/5x7sz7/professor-anastasios-magdalinos#supervision). In his PhD thesis, Katsouris C. generalises the work of Phillips P.C.B (under the guidance of his advisor Magdalinos T.) on robust econometric inference for nonstationary time series models, under more general persistence conditions than his predecessors. Furthermore, he proposes testing methodologies for structural break detection in predictive regressions.

####  An Econometrician amongst Statisticians or A Statistician amongst Econometricians? 

Christis G. Katsouris has participated to various international academic conferences such as: (i) the International Symposium  on Nonparametric Statistics (ISNS), June 2022 in Paphos, Cyprus, (ii) the Hausdorff Center for Mathematics Summer School in High-Dimensional Statistics, July 2021 in Bonn, Germany and (iii) the International Association for Applied Econometrics, Annual Conference, June 2019 in Nicosia, Cyprus. 


# Historical Background

> Standing on the shoulders of giants.
> 
> ''_If I have been able to see further, it was only because I stood on the shoulders of giants._"
> Isaac Newton, 1676 

#### Aleksandr Lyapunov 

Aleksandr Lyapunov  (6 June 1857 – 3 November 1918) was a Russian mathematician, mechanician and physicist. Lyapunov is known for his development of the stability theory of a dynamical system, as well as for his many contributions to mathematical physics and probability theory. A major theme in Lyapunov's research was the stability of a rotating fluid mass with possible astronomical application. His main preoccupations were the stability of equilibria and the motion of mechanical systems, especially rotating fluid masses, and the study of particles under the influence of gravity. Lyapunov's work in the field of mathematical physics regarded the boundary value problem of the equation of Laplace. In the theory of potential, his work from 1897: "On some questions connected with Dirichlet's problem" clarified several important aspects of the theory. His work in this field is in close connection with the work of Steklov. Lyapunov developed many important approximation methods. His methods, which he developed in 1899, make it possible to define the stability of sets of ordinary differential equations. In the theory of probability, he generalised the works of Chebyshev and Markov, and proved the Central Limit Theorem under more general conditions than his predecessors. The method of characteristic functions he used for the proof later found widespread use in probability theory (Source: Wikepedia).

#### Harald Cramér 

Harald Cramér  (25 September 1893 – 5 October 1985) was a Swedish mathematician, actuary, and statistician, specializing in mathematical statistics and probabilistic number theory. John Kingman described him as "one of the giants of statistical theory". A large portion of Cramér's work concerned the field of actuarial science and insurance mathematics. In 1929, Cramér was appointed to a newly created chair in Stockholm University, becoming the first Swedish professor of Actuarial Mathematics and Mathematical Statistics. Cramér retained this position up until 1958. During his tenure at Stockholm University, Cramér was a PhD advisor for 10 students, most notably Herman Wold and Kai Lai Chung. In 1950 he was elected as a Fellow of the American Statistical Association. Starting in 1950, Cramér took on the additional responsibility of becoming the President of Stockholm University. In 1958, he was also appointed to be Chancellor of the entire Swedish university system (Source: Wikepedia). 

#### Charles Stein

Charles Stein  (22 March 1920 – 24 November 2016) was an American mathematical statistician and professor of statistics at Stanford University. 
He received his Ph.D in 1947 at Columbia University with advisor Abraham Wald. He held faculty positions at Berkeley and the University of Chicago before moving permanently to Stanford in 1953. He is known for Stein's paradox in decision theory, which shows that ordinary least squares estimates can be uniformly improved when many parameters are estimated; for Stein's lemma, giving a formula for the covariance of one random variable with the value of a function of another when the two random variables are jointly normally distributed; and for Stein's method, a way of proving theorems such as the Central Limit Theorem that does not require the variables to be independent and identically distributed (Source: Wikepedia).


