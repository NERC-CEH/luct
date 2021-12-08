



# Estimating false-positive rates in detection of land-use change based on classification accuracy

## Introduction
In the "Tracking Land-Use Change" project, several data sources provide a time series of maps of land use. An obvious approach is to estimate land-use change as the difference between these maps over time. However, any error in land classification will also be included in the estimate of land-use change, so it is important to quantify these errors properly and note how they propagate.
Here, we show the calculations which propagate the error in classification through to its resulting effect on the estimate of land-use change.


## Methods

The accuracy of land-use maps is usually estimated by comparison with some reference data set. We can then calculate confusion matrices and metrics of overall agreement, of which there are several. Common choices are overall accuracy ($\alpha$, the fraction of locations where estimated land use agrees with the reference data set) and the $\kappa$ statistic, which corrects for the  probability of chance agreement. $\kappa$ therefore gives a more robust measure, typically 5-20 % lower than simple percentage agreement. The probability of misclassification can be estimated simply as 1 - $\alpha$ or more stringently as 1 - $\kappa$. For maps at times $t_1$ and $t_2$, the probabilities of misclassification are denoted $p_1$ and $p_2$. Estimating land-use change involves calculating the difference between maps, and the errors are additive in the result. The probability of estimating erroneous land-use change because of misclassification in a pair of maps can be written as

\begin{equation} \label{eq:probmis}
  p_{1 \cup 2} = p_1 + p_2 - p_{1 \cap 2}.
\end{equation}

That is, the probability of error is the union of two events (misclassification occurring at time 1 or at time 2, minus their intersection $p_{1 \cap 2}$, which is the probability of misclassification at both time 1 and time 2, which would otherwise be double-counted.  $p_{1 \cap 2}$ can be estimated as $p_1 p_2$ assuming that the errors leading to misclassification at times 1 and 2 are independent of each other.
In practice, our estimates of the probabilities of misclassification at times 1 and 2 are usually the same ($p_1 \simeq p_2$), so this simplifies to:

\begin{equation} \label{eq:probmis2}
  p_{1 \cup 2} = 2 p_1 - p_1^2.
\end{equation}


## Results
Estimates of $\alpha$ and $\kappa$ from some of the data sets used in the LUC Tracking project are shown in the table below. 


\begin{tabular}{l|r|r}
\hline
Data source & $\alpha$ & $\kappa$\\
\hline
Corine & 0.80 & 0.64\\
\hline
LCC & 0.91 & 0.82\\
\hline
LCM & 0.88 & NA\\
\hline
\end{tabular}

For the purposes of the examples below, we use the value of 0.88, the overall accuracy of the LCM, as a relatively optimistic metric.
The value of $p_1 \simeq p_2$ is 1 - $\alpha$, and therefore = 0.12. 



Using this value in Equation \ref{eq:probmis2} yields a probability of estimating erroneous land-use change because of misclassification of 0.226.
Because this probability applies at every location on the map, multiplying by the total area yields the expected area of erroneous land-use change.
So, when comparing two UK maps which each have a classification accuracy of 88 %, 22.6 % of the area, around 55000 km$^2$, will show land-use change where none actually occurs. This provides a huge amount of measurement noise when we are attempting to detect a very small signal: the expected magnitude of actual land use change in the UK is of the order of a few hundred km$^2$, and at most a few thousand km$^2$, based on Forestry Commission planting rates, Agricultural Census, and urban expansion data. The area of land changing use is thus less than 1 % of the total area, and we would therefore need the probability of misclassification error to be less than this in order to accurately detect true change (meaning the accuracy needs to be > 99 %).

We can extend this to calculate the false positive rates for terms in the $\beta$ matrix and gross gains and losses, given the appropriate denominators and estimates of the true extent of land-use change.  For example, the area of cropland in England is approximately 45000 km$^2$, and the area of gross gains and losses are estimated to be in the range 300-800 km$^2$ y$^{-1}$ based on CS.  Based on the June Agricultural Census, we might estimate these rates to be higher, perhaps reaching 1000-3000 km$^2$ y$^{-1}$.
 
 


Expressing the estimated true land-use change $A_{true}$ as a fraction of the total recorded land-use change (i.e. true + erroneous), we can calculate the relevant false positive rate, $F_P$.  

\begin{equation} \label{eq:falsepos}
  F_P = 1 - \frac{A_{true}}{A_{true} + A_{false}}
\end{equation}

If CS rates of land-use change are correct (300-800 km$^2$ y$^{-1}$), the false positive rate is in the range 92.7 to 97.1 %.



If the Agricultural Census rates of land-use change are correct (1000-3000 km$^2$ y$^{-1}$), the false positive rate is in the range 77.2 to 91 %.

Given this relationship between the classification accuracy of individual maps and the resulting false positive rates in detecting land-use change, we can examine the improvement needed to obtain false positive rates below a given level.  The figure below shows the change in false positive rate with classification error (1- $\alpha$), using the example of cropland gains in England as above, assuming a true rate of change of 1000 km$^2$ y$^{-1}$.

![](m_fp_from_accuracy_files/figure-latex/unnamed-chunk-5-1.pdf)<!-- --> 

The figure shows asymptotic behaviour because of the form of Equation \ref{eq:falsepos}, with a constant true area expressed with an increasing $A_{false}$ term in the denominator.  The result is that it takes a substatial decrease in misclassification (or increase in accuracy) from present values to achieve a marked increase in the false positive rate.  For example, to reduce the false positive rate to 0.5 requires an accuracy of 0.988. The basic problem is that the true areas of change are very small compared to current error rates, and it would require an order of magnitude improvement in accuracy to reduce the measurement noise to a similar level to the signal we want to detect.


## Conclusions

* Accuracy of land-use classification is in the range 0.8-0.9 in the available data sets. The corresponding probability of misclassification is 10-20 %.
* The errors in individual maps can be propagated to calculate the error in their differences i.e. estimates of land-use change.
* Error rates for land-use change propagating from misclassification will typically be greater than 20 %. This is 1-2 orders of magnitude larger than the expected land-use change.
* The probability of misclassification can be used to calculate the false positive rate of land-use change detection, and values are typically around 90 %.
* Whilst these errors are predictable and can be accounted for, directly detecting the expected land-use change rates <1 % is not currently practicable.



