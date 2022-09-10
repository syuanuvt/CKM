# Cardinality K-Means (CKM): A Brief Tutorial
The package implements the Cardinality K-Means algorithm that performs simultaneous variable selection and clustering. Compared to existing approaches, CKM is tested to be faster and yield at least comparable performance. The algorithm is especially suitable to tackle the clustering problem with a large numer of variables; in such cases picking up signaling variables is vital to recovering the clusters

## Corresponding address
Issues on the packages could be filed as pull requests (recommended) or sent to s.yuan@uva.nl

## New version
A new version of CKM is currently under-development with updated functions for model selection. Stay tuned!

## Installation
To install the Package, make sure that the R package "devtool" and "remotes" have been installed. Then, use the following code to install the Package CKM. 

```{r eval = FALSE}
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS="true")
remotes::install_github("syuanuvt/CKM")
```

## Data requirement
A dataset (stored in R with the class "matrix" or "data frame") should be processed to have all variables mean-centered and standardized before putting into CKM. A simulation function could be also used to create artificial datasets.

## Data simulation

### Function
**DataGenCKM**: Data simulation that accords to the CKM models, with the option to select some most important features, including (1) the number of observations, (2) the number of clusters, (3) the number of masking variables, (4) the number of signaling variables, (5) the mean difference in centroids between two neighboring clusters, and (6) the standard normal distribution that the observations are drawn from.

The output of the function is a list containing two elements: (1) the dataset as a matrix, (2) the cluster assignment as a vector

### Example
```{r eval = FALSE}
nobs <- 60
ncluster <- 3
nnoisevar <- 100
nsigvar <- 20
mu <- 1 # the mean differences of the centroids between two neighboring clusters
var <- 1 # therefore, each response is drawn from N(mean, 1)
sim.data <- DataGenCKM(nobs, ncluster, nvalidvar, nnoisevar, mu, var)
dataset <- sim.data[[1]]
true.cluster.assign <- sim.data[[2]]
```

## Execute the CKM algorithm
Depending on available knowledge and assumptions, one has three options to carry out the CKM analysis: (1)assume both parameters are known and therefore do not perform any model selection; (2) assume the number of clusters is known, therefore only select the number of signaling variables; (3) assume both parameters are unknown and rely on the algorithm to select the number of signaling variables as well as the number of clusters. The package CKM offers three separate functions for these three options. Although the inputs of three functions may differ, the output is always an object of class "ckm", which could be called subsequently.

### No model selection
**CKMSelNo**: the function could be used when the number of masking variables and the number of clusters is known. Both of them are treated as input when using the function. In addition, the user could specify the number of starts for the inherent conventional K-Means algorithm (a higher number of starts takes more computational time but likely yield more accurate estimation).

```{r eval = FALSE}
ncluster <- 3 #assumed to be the true value
nnoisevar <- 100 #assumed to be the true value
ckm.sel.no <- CKMSelNo(dataset, ncluster, nnoisevar) #the defaults are used for other input arguments
```

### Selection of the number of masking variables
**CKMSelVar**: the function could be used when the number of clusters has been determined (from existing knowledge or other model selection approaches) and the number of masking variables has to be selected. In this case, only the number of clusters will be used as the input of the function. In addition, the user could select the strategy applied to select over the grid of the parameter ("all" = the full grid, "sub" = always using the "grid search with a zoom" strategy, "dep" = the strategy employed should depend on the number of variables), and the number of re-samples in model selection (a higher number of re-samples takes more computational time but likely yield more accurate estimation)

```{r eval = FALSE}
ncluster <- 3
ckm.sel.var <- CKMSelVar(dataset, ncluster) #the defaults are used for other input arguments
```

### Selection of both the number of masking variables and the number of clusters
**CKMSelAll**: the function should be used when both parameters have to be determined by the algorithm. In this case, the user only needs to supply the expected maximal number of clusters. In addition, since the model selection approach with the gap statistic, as implemented here, could use various criteria, the user could specify the preferred criterion, as done in the package \code{cluster}. Other arguments of the function are similar to those of the function *CKMSelVar*. Note that the results of the selection of the number of clusters are reported in the console, which should be taken by the user manually and is useful for the subsequent reporting and plotting.

```{r eval = FALSE}
maxcluster <- 10
ckm.sel.all <- CKMSelAll(dataset, mincluster, maxcluster) #the defaults are used for other input arguments
```

### Summary and plot
## summary
The "ckm" object returned from any of the above three functions will be used for the reports of the plots. In addition, the original dataset, as well as the assumed or selected number of clusters are used as input arguments.
A total of four elements, compiled in a list, will be returned. The first is an integer that indicates the number of noise variables. The second is a vector that contains the indices of signaling variables. The third demonstrates the cluster assignment. The last one, a matrix, represent the centroids of the signaling variables (columns) for all clusters (rows).

```{r eval = FALSE}
ncluster <- 3
## the output of all of the above three functions is of a S3 class "ckm", which could be directly used as an argument in the generic function "summary"
estimates.no <- summary(dataset, ckm.sel.no)
estimates.var <- summary(dataset, ckm.sel.var)
estimates.all <- summary(dataset, ckm.sel.all)
```

## Plot
the function plots a graph with x-axis representing various clusters and y-axis referring to the centroids of signaling variables.

```{r eval = FALSE}
## the output of all of the above three functions is of a S3 class "ckm", which could be directly used as an argument in the generic function "plot"
plot(dataset, ckm.sel.no)
plot(dataset, ckm.sel.var)
plot(dataset, ckm.sel.all)
```

### Analysis Pipeline
Data preprocess (variable standardization and centering) -> Select the number of clusters (with function CKMSelAll) -> 
Select the number of irrelevant variables (with function CKMSelVar) -> Final computation (with function CKMSelNo) ->
Summary and Plot

## Depends
RSpectra, cluster, GGally

## License
GPL-3
