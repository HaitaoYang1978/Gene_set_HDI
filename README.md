# Simulation-_real_SNP
## Description
Our proposed procedure aims to detect the set-based genetic variation (e.g., SNP-set, gene-set) associated with complex disease trait under a high-dimensional inference frame. Furthermore, We propose an omnibus testing procedure that employs a robust and powerful p-value combination method to enhance the power of set-based genetic variation association.
## Dependency
Our proposed procedure relies on the `glmnet`  `future.apply` `MASS` `mvtnorm` `devtools` `ACAT` `screening`. All required package can be installed automatically in our procedure. In addition, all corresponding R functions of de-sparsified lasso method are in the *'./de-sparsified lasso functions'* folder.
## Usage
A simple working example. We use a simulation data in our manuscript’s Simulation study’s CaseⅡ. Please run the R file, **main_HOLP_HDI_Omnibus.R** in which the analysis procedure of a random subset was illustrated. 
## Author
Haitao Yang and Fuzhao Chen. We appreciate your feedback under issues of this repository.
## References
1. Fan, Jianqing, and Jinchi Lv. "Sure independence screening for ultrahigh dimensional feature space." Journal of the Royal Statistical Society: Series B (Statistical Methodology) 70.5 (2008): 849-911.
2. Wang X, Leng C. High dimensional ordinary least squares projection for screening variables, Journal of the Royal Statistical Society: Series B (Statistical Methodology) 2016;78:589-611.3. Li, Gaorong, et al. "Robust rank correlation based screening." The Annals of Statistics 40.3 (2012): 1846-1877.
3. Liu Y, Xie J. Cauchy combination test: a powerful test with analytic p-value calculation under arbitrary dependency structures, Journal of the American Statistical Association 2019:1-18.
4. Vsevolozhskaya OA, Hu F, Zaykin DV. Detecting weak signals by combining small P-values in genetic association studies, Frontiers in Genetics 2019;10:1051.
5. Zhang CH, Zhang SS. Confidence intervals for low dimensional parameters in high dimensional linear models, Journal of the Royal Statistical Society: Series B (Statistical Methodology) 2014;76:217-242.
6. Van de Geer S, Bühlmann P, Ritov Ya et al. On asymptotically optimal confidence regions and tests for high-dimensional models, The Annals of Statistics 2014;42:1166-1202.

