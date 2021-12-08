




# Quantifying Random Uncertainty in Land-Use Data Sources

## Introduction

Several different data sources provide observations of the transition matrix $B$ and the net and gross changes in area of each land use ($D, G, L$). The method used in WP-A treated all these data sources as equally uncertain, and assumed the same relative error for all observations. However, in reality, we know that these data sets have different levels of uncertainty: some data sets are closer to direct observations, are more plausible, and we have greater faith in these. We want to reflect this in the methodology by quantitatively associating different uncertainties with each data set. This is straightforward in principle, but there are several considerations when doing this in practice:

1. We can consider increasing levels of detail:

    -   variable-specific uncertainties (i.e. different for $B, G, L$ & $D$)

    -   data source type-specific uncertainties (i.e. different for ground-based vs EO data)

    -   data set-specific uncertainties (i.e. different for CS, IACS, LCM etc.)

    -   land-use type-specific uncertainties (i.e. different for woods, crops, grass etc.)

    -   time-specific uncertainties (i.e. different for 1990, 2000 ... 2019)

2. Rather than continuous data with a simple $\sigma$ error term, the $B$ observations are count data in a 6 x 6-way classification. When considering land-use *change*, we can compare the 36 elements of this classification from one data source with another (or the truth), so we have a 36 x 36 error matrix (or "confusion" matrix). Various metrics can be calculated which summarise the agreement measured by this matrix. 
3. We can specify uncertainty with greater or lesser rigour: there are several possibilities for how we represent "uncertainty" in the mathematical model. 
4. We can estimate uncertainty subjectively, or base it more closely on data. There are also several possibilities for how we translate measures of uncertainty in the data in to the mathematical model.

A limitation is that none of the data sources represents absolute truth, and we have no definitive data set against which to calibrate.

### Representation of random uncertainty
The data sources are assimilated in the Bayesian method via the likelihood function, which includes a term $\sigma^{obs}$, representing the standard deviation in the probability density function for the observation. The observation is thus not assumed to be the true value, but subject to errors which make it deviate from this. Random uncertainty is represented by the magnitude of $\sigma^{obs}$ - large values of $\sigma^{obs}$ represent high uncertainty (systematic uncertainty is considered in later Sections).
For each observation, a likelihood is calculated, assuming that measurement errors show a Gaussian distribution and are independent of each other.
In mathematical notation, the likelihood of observing the area changing from land use $i$ to land use $j$, $\beta_{ij}^{\mathrm{obs}}$, is

\begin{equation} \label{eq:likBeta}
 \mathcal{L} = \frac{1}{\sigma_{ij}^\mathrm{obs} \sqrt{2\pi}} \mathrm{exp}(-(\beta_{ij}^{\mathrm{obs}} - \beta_{ij}^{\mathrm{pred}})^2/2 {\sigma_{ij}^\mathrm{obs}}^2)
\end{equation}

where $\beta_{ij}^{\mathrm{pred}}$ is the corresponding prediction, and $\sigma_{ij}^\mathrm{obs}$ is the uncertainty in the observation. There are analagous terms for $G, L$ and $D$ which can all be multiplied. For example, the term for the likelihood of observing the net change in land use $u$, $D_{u}^{\mathrm{obs}}$, is
 
\begin{equation} \label{eq:likNet}
 \mathcal{L} = \frac{1}{\sigma_{u}^\mathrm{obs} \sqrt{2\pi}} \mathrm{exp}(-(D_{u}^{\mathrm{obs}} - D_{u}^{\mathrm{pred}})^2/2 {\sigma_{u}^\mathrm{obs}}^2)
\end{equation}
 

We previously assumed that *relative* measurement uncertainty was the same for all observations, i.e. a constant proportion of the observed value. Thus, observations of large areas come with larger absolute uncertainty. Here, we estimate more specific uncertainties $\sigma^\mathrm{obs}$ for the different data sources, and potentially this can be extended to be specific for each individual observation.

Two additional issues concern the specification of the random uncertainty when accounting for the effects of the frequency of surveys, and avoiding step changes when data sets begin and end. For example, Countryside Survey data come from approximately decadal surveys, but are interpolated and used as if they were constant annual values within each decade. We incorporated some simple methods to include these effects appropriately.

## Methods
### Estimating $\sigma^\mathrm{obs}$ from reference data
A pre-cursor step is to define a reference data set with which we compare each data source. In the absence of ideal ground-truth data, some subjectivity is inevitable here, and we use the data sets which we believe to be the most plausible or closest to the truth, based on judgement and prior knowledge. For agricultural land, we defined the reference data set as the June Agricultural Census data for crops, grass and rough grazing, as this is has a very long record and is submitted annually as part of reporting to FAO.
For forests, the pre-existing record of afforestation and deforestation based on FC statistics was used, as this is also long-running and has been submitted as official national data.

We applied a simple method, basing all the $\sigma^\mathrm{obs}$ terms on the "lowest common denominator" data set, the time series of net area change $D$. Although they vary in the higher levels of detail, all data sources produce estimates of $D$, so we can calculate a comparable metric of agreement across all data sets. Suitable metrics are the root-mean-square error (RMSE), mean-absolute error (MAE) and the correlation coefficient.
We used the RMSE and the correlation coefficient, as measures of absolute and relative agreement, respectively.
We multiplied these (rescaling the correlation coefficient as $1/(r^2+1)$) to give a single scaling factor for $\sigma^\mathrm{obs}$ ("CV" in the table). 
This metric is used as a scaling factor in estimating $\sigma^\mathrm{obs}$, such that data sources with poor agreement receive high $\sigma^\mathrm{obs}$ (high uncertainty). Some subjectivity comes into which metrics to use, and the absolute values of $\sigma^\mathrm{obs}$, but the relative uncertainties (and therefore weighting) are based quantitatively on data.

### Effects of survey frequency and survey start/stop dates
Each pair of surveys gives an estimate of land-use change over some time interval.
When survey observations are not available annually, we are effectively trying to estimate the change in a population (all years) from a sample (the years when surveys were carried out).
Using the analogy with conventional sampling theory, the standard error in our estimate should be lower in more frequent surveys, because more samples are included.
We can apply the same logic to infer the appropriate correction to apply to $\sigma^\mathrm{obs}$, such that it reflects the uncertainty about the rate of land-use change in any given year within the inter-survey interval:

\begin{equation} \label{eq:sigma_ann}
  \sigma_\mathrm{ann}^\mathrm{obs} = \sigma^\mathrm{obs} \sqrt{n}
\end{equation}

where $n$ is the survey interval length in years. This scaling means that the correction evaluates to 1 for annual data (no effect).

All data sources have a limited period over which they were collected, and this would introduce an artefact at the boundaries.
When a data set which tends to provide lower estimates begins, this would tend to pull the mean estimate down at this point, and introduce a step change purely as an effect of data availability.
To counter this effect, we apply a similar logic to that above - the available data is an imperfect sample of the surrounding time interval, and the uncertainty increases with fewer samples and extended distance in time. 
We therefore extend each data set to the limits of the time period considered here (1950-2020), assuming a constant rate, but we force the uncertainty to increase with the square of distance in time ($\Delta_t^2$) beyond the boundaries where data was actually available. The random uncertainty in any given year is

\begin{equation} \label{eq:sigma_t}
  \sigma_{t}^\mathrm{obs} = \sigma^\mathrm{obs} \Delta_t^2
\end{equation}

again scaled so that there is no effect within the time bounds of the observed data.


## Results
The table below shows the RMSE, the correlation coefficient and the scaling factor for $\sigma^\mathrm{obs}$ (expressed as a coefficient of variation (CV), $\sigma^\mathrm{obs}$ as a fraction of the observed value. 
The ranking shows that CS has the lowest uncertainty and the CROME has the highest (bearing in mind all the imperfections in the reference data). This produces a quantitative means of accounting for the different relative uncertainties in these data sources. 












|       |      RMSE|         r2|        CV|
|:------|---------:|----------:|---------:|
|CS     |  654.9496|  0.2850220| 0.2000000|
|CORINE |  830.2048|  0.3533137| 0.2293023|
|LCM    |  892.9640|  0.1650214| 0.3184481|
|IACS   | 1170.0198|  0.1882285| 0.4056547|
|LCC    | 1340.3314|  0.1350960| 0.4951190|
|CROME  | 7899.3774| -0.0819340| 3.6502540|

## Discussion

Representing data source-specific random uncertainty is relatively straightforward in principle. We need to estimate appropriate $\sigma$ values for each data source and use these in the likelihood function.
The most fundamental problem is accurately estimating $\sigma$ in the absence of any data which we regard as "true", particularly for the $B$ matrices which are key. There is no immediate solution to this, and a pragmatic approach is to define a reference data set, with more or less subjectivity/expert judgement, and potentially with some cross-validation.
Here we implemented a simple method, whereby $\sigma$ for each data set is scaled according to metrics measuring its correspondence with reference data. Currently, this is based only on the net change data, as this makes cross-comparison simplest, but this could be extended to include the gross changes $G, L$ and $B$.

An alternative approach would be to estimate the uncertainties as part of the data assimilation. This avoids the sticking point of subjectively defining a suitable reference data-set, when all of the available data sources, including ground-truthed data, are flawed in some way. The downside of this approach is that it is more complicated, involves estimating more parameters, and will have greater computation time, but merits some exploration.



```
## Computation time (excl. render): 0.2 sec elapsed
```
