# ivxPredictive 

## Unified Inference with General Autoregressive Roots

### Description 

The R package 'ivxPredictive' (under development) implements robust econometric inference methodologies for predictive regression models with autoregressive roots across the spectrum of stationarity and nonstationarity, with persistence classes as defined by Magdalinos and Phillips (2020). In particular, our R implementation extends the original IVX instrumentation of [Phillips and Magdalinos (2009)](https://ideas.repec.org/p/skb/wpaper/cofie-06-2009.html) and [Kostakis, Magdalinos and Stamatogiannis (2015)](https://academic.oup.com/rfs/article/28/5/1506/1867633?login=true) to general autoregressive roots for parametric functional forms. In particular, we focus on linear specifications for the predictive regression model and implement the proposed procedure for the case of a conditional mean loss function (see, Magdalinos and Petrova (2022)) as well as a conditional quantile loss function (see, [Katsouris, C. (2022)](https://arxiv.org/abs/2204.02073)).  

<p align="center">
  
<img src="https://github.com/christiskatsouris/ivxPredictive/blob/main/data/persistence.jpg" width="460"/>

</p>  

We consider the following persistence classes:
   
$\textbf{(P1): nearly stable}$ processes: 

$$\text{if} \ \ \left( \theta_n \right)_{ n \in \mathbb{N} } \ \ \text{is such that} \ \zeta = - \infty \ \ \text{and it holds that} \ \ \theta_n \to | \theta | < 1.$$
    
$\textbf{(P2): nearly unstable}$ processes:   

$$\text{if} \ \ \left( \theta_n \right)_{ n \in \mathbb{N} } \ \ \text{is such that} \ \ \zeta \equiv c \in \mathbb{R} \ \ \text{and it holds that} \ \ \theta_n \to \theta = 1.$$

    
$\textbf{(P3): nearly unstable}$ processes:   

$$\text{if} \ \ \left( \theta_n \right)_{ n \in \mathbb{N} } \ \ \text{is such that} \ \ \zeta = + \infty \ \ \text{and it holds that} \ \ \theta_n \to | \theta | > 1.$$

  
### Methodology  
  
This R package implements a novel endogenous instrumentation approach based on the IVX estimator examined by [Phillips and Magdalinos (2009)](https://ideas.repec.org/p/skb/wpaper/cofie-06-2009.html) and [Kostakis, Magdalinos and Stamatogiannis (2015)](https://academic.oup.com/rfs/article/28/5/1506/1867633?login=true). The current procedure has a similar construction to the IV instrumentation proposed in the recent working paper of Magdalinos and Petrova (2022), with the aim to provide uniform and robust to the nuisance parameter of persistence inference across the spectrum of stationary and nonstationary roots, specifically for quantile autoregressive processes. We call our novel estimator IVX-P, which can be employed to both conditional mean and conditional quantile functional forms when the model includes either univariate or multivariate regressors. The novelty of the IVX-P estimator is that is a 'hybrid estimator' which in contrast to the classical least squares estimator has desirable asymptotic theory properties and is constructed based on the underline nonstationary stochastic processes using information both within the admissible parameter space as well as outside the usual parameter space.  
  
## Installation (under development) 

The R package ‘ivxPredictive’ will be able to be installed from Github.

## Usage 

```R

# After development the package will be able to be installed using
install.packages("ivxPredictive")
library("ivxPredictive")

```

## Main References

- Katsouris, C. (2022b). "Asymptotic Theory for Moderate Deviations from the Unit Boundary in Quantile Autoregression and Predictive Regression Models". University of Southampton, Working paper.  

- Katsouris, C. (2022a). "Asymptotic Theory for Moderate Deviations from the Unit Boundary in Quantile Autoregressive Time Series". arXiv preprint [arXiv:2204.02073](https://arxiv.org/abs/2204.02073).

- Kostakis, A., Magdalinos, T., & Stamatogiannis, M. P. (2015). "Robust econometric inference for stock return predictability". The Review of Financial Studies, 28(5), 1506-1553.

- Lee, J. H. (2016). "Predictive quantile regression with persistent covariates: IVX-QR approach". Journal of Econometrics, 192(1), 105-118.

- Magdalinos, T., & Petrova, K. (2022). "Uniform and distribution-free inference with general autoregressive processes". University of Southampton, Working paper. 

- Magdalinos, T., & Phillips, P. C. B. (2020). "Econometric inference in matrix vicinities of unity and stationarity". University of Southampton, Working paper.  

- Magdalinos, T., & Phillips, P. C. B. (2009). "Limit theory for cointegrated systems with moderately integrated and moderately explosive regressors". Econometric Theory, 25(2), 482-526.

- Phillips, P. C. B., & Magdalinos, T. (2009). "Econometric inference in the vicinity of unity". Singapore Management University, CoFie Working paper, 7.

## Notes

1. The implementation of the original (current) IVX instrumentation method proposed by [Phillips and Magdalinos (2009)](https://ideas.repec.org/p/skb/wpaper/cofie-06-2009.html) along with related statistical inference and empirical applications has been prepared as an R package [ivx](https://github.com/kvasilopoulos/ivx) developed by [Kostas Vasilopoulos](https://github.com/kvasilopoulos).

2. The Matlab code of the original IVX instrumentation method can be found on the website of [Michalis Stamatogiannis](https://sites.google.com/site/mpstamatogiannis/home).

3. The author would like to thank [Ji Hyung Lee](https://economics.illinois.edu/profile/jihyung) for kindly sharing the Matlab code for the implementation of the original IVX instrumentation method for the quantile predictive regression. 

## Code of Coduct

Please note that the ‘ivxPredictive’ project will be released with a Contributor Code of Coduct (under construction). By contributing to this project, you agree to abide by its terms.

## Declarations

The author declares no conflicts of interest. 

In particular, the author declares that has no affiliations with or involvement in any organization or entity with any financial interest (such as honoraria; educational grants; participation in speakers’ bureaus; membership, employment, consultancies, stock ownership, or other equity interest; and expert testimony or patent-licensing arrangements), or non-financial interest (such as personal or professional relationships, affiliations, knowledge or beliefs) in the subject matter or materials discussed in the manuscript and implemented in the R package.

# Acknowledgments

The author greatfully acknowledges financial support from Graduate Teaching Assistantships at the School of Economic, Social and Political Sciences of the University of Southampton as well as from the Vice-Chancellor's PhD Scholarship of the University of Southampton, for the duration of the academic years 2018 to 2021. Furthermore, the author gratefully acknowledges funding from Research Grants of interdisciplinary Centers of research excellence based at the University of Cyprus (UCY) as well as at the University College London (UCL).

Part of the aspects implemented in the R package 'ivxPredictive', were discussed in the PhD thesis of the R package creator (Christis G. Katsouris) titled: "Aspects of Estimation and Inference for Predictive Regression Models".
