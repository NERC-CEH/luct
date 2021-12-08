--- 
title: "Improving Land Use Change Tracking in the UK Greenhouse Gas Inventory. Report on final outputs November 2021"
author: "Peter Levy, Sam Tomlinson, and Beth Raine"
date: "2021-11-10"
site: bookdown::bookdown_site
documentclass: book
bibliography: [luc.bib]
biblio-style: apalike
link-citations: yes
github-repo: NERC-CEH/luc_track
description: "Code and documentation for the BEIS project Improving Land Use Change Tracking in the UK Greenhouse Gas Inventory"
---

# Summary {-}
This report describes work on the project "Improving Land Use Change Tracking in the UK Greenhouse Gas Inventory" for the Department for Business, Energy & Industrial Strategy (reference TRN 2384/05/2020). The aim of the project was to make improved estimates of land-use change in the UK, using multiple sources of data. We applied a method for estimating land-use change using a Bayesian data assimilation approach. This allows us to constrain estimates of gross land-use change with national-scale census data, whilst retaining the detailed information available from several other sources. Previous reports covered work with existing data sets (WP-A) and on developing new data sets (WP-B). This report describes subsequent work, focussed on quantifying uncertainties and improving the representation of crop-grass rotation.

The random uncertainty term $\sigma$ for each data source has been estimated, based on agreement with the net change observed in a reference data set, comprising those data sets we believe to be most reliable for net change. The random uncertainty is now represented with both absolute and relative components, as a linear function with intercept and slope terms.	

Systematic uncertainties represent biases in the data, and can be charaterised as the false positive and false negative rates (i.e. detecting apparent change where none really occurs, and failing to detect actual change, respectively). We compared three different methods for estimating these rates, which gave similar results. These showed very high false positive rates, which mean that large corrections have to be applied to the data. The data assimilation algorithm now includes terms to account for: data-source-specific random uncertainty (sigma); data-source-specific systematic uncertainty (false positive and false negative rates); sampling error due to low survey frequency (e.g. decadal cf. annual); and uncertainty in extrapolation outwith (pre/post) survey time period.

Additional data sets, which were not previously available have also been included in the analysis. All the holdings-level agricultural census data available in the UK from 1990 to 2019 have now been included. The CROME data set, covering crops in England, has been added.
 
We have improved the accuracy of how we represent crop-grass rotations, using the idea of "life tables" or age-specific transition probabilities, from population modelling. The life table probabilities are based on an analysis of the data available data from IACS, CROME, LCM, and LCC, and these all show consistent patterns. Of these, IACS has much the longest time span of data, and is the dominant source of information. Using this method, the observed frequencies of transitions from crop to grassland as a function of cropland age (and v.v. for grassland) are now reproduced in the land-use vectors, and thereby approximate the observed frequency of crop-grass rotations.

Since land-use change is now represented spatially, we can separately identify that occurring on mineral soils and organic soils. In summary form, land-use vectors and the matrices of land-use change are now provided separately for mineral and organic soils for each DA.

The land-use change output data are provided as space-time data cubes at 100-m and 1000-m resolution for the period 1950-2020. These are also provided in summary (non-spatial) form as the set of unique land-use vectors, and the matrices of land-use change.  In principle, many thousand samples of the posterior distribution are available (currently 20000 samples of the matrices are stored per DA). Many samples of the mapped and vector data can be provided from these, given practical computing and storage constraints, but the optimal way to produce these depends on how they are used in the LULUCF inventory work, and this is open to discussion.  Storing land use as a 71-digit character string at 91 million locations (for the time period 1950-2020 at 100-m resolution) requires ~50 MB, or ~5MB in vector summary form, and much less at 1-km resolution.

The workflow for the project uses a "Make"-like paradigm to maintain a reproducible workflow, implemented in the R "targets" package (https://books.ropensci.org/targets/). This means that the dependencies among functions and input data are analysed dynamically and stored in a hash table. Any changes to source code functions or input data are detected automatically, and only the invalidated components re-run as necessary. This has several advantages: forcing the workflow to be declared at a higher level of abstraction; only running the necessary computation, so saving run-time for tasks that are already up to date; and most importantly, providing tangible evidence that the results match the underlying code and data, and confirm the computation is reproducible. All the code is under version control on GitHub.

MCMC algorithms typically have to run for many iterations to produce an acceptable sample of the posterior, so that parameter estimates (the Beta matrices) are robust. For each year and each DA, the algorithm was run with nine MCMC chains for 200,000 iterations.  Convergence was assessed by a standard method - comparing within-chain variance with among-chain variance via the Gelman-Rubin statistic. As chains move towards convergence, the value of this tends towards 1, and the standard test requires that this value is less than 1.1. This criterion is well met with this number of iterations. It is standard to discard "burn-in" samples prior to convergence, and to thin chains by a factor of ten to remove any autocorrelation.  This yields around 16000 posterior Beta samples, which gives potentially the same number of posterior samples of the land-use vectors. However, it may be more practical to combine the calculation of fluxes into this code, without storing all the underlying data.

We produced a time series of maps describing our best estimate of land-use change given the available data, as well as the full posterior distribution of this space-time data cube. This quantifies the joint probability distribution of the parameters, and properly propagates the uncertainty from input data to final output. The output data has been summarised in the form of land-use vectors. The results show that we can provide improved estimates of past land-use change using this method. The main advantage of the approach is that it provides a coherent, generalised framework for combining multiple disparate sources of data, and adding additional sources of data is straightforward.


<!--- { rendering -->

<!--- } -->

