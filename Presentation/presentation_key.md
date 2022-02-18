# Slide 1:
Hi, my name is Guillaume, I am a Phd. Candidate from the University of Geneva.
Today, we'll talk about how to deal with multivariate data with missing values.

# Slide 2:
Because this presentation is short, we'll take a very simple setting:

* on the left-hand sidey, we have a dataset with 3 variables, Y1, Y2, and Y3. The first variable, Y1, has some missing values marked as "NA".
* on the right-hand side, we have what is called the missingness "Mask", where a 1 encodes a missing data and a 0 an observed data.
* My **Goal** is to estimate the mean of Y1: to do that, I propose to impute the missing values by the average of the observed values, and then take the average of the whole as estimate for the mean.
* This estimator will be unbiased if the missingness mechanism is what we call "missing completely at random", or MCAR: that is, the probability that a data is missing does not depend on its value.
* We'll take a second example where the missingness mechanism is NOT MCAR: say, for instance, the data is truncated and only the 50%  smallest values are observed. In this case, since we do not observed the larger values, imputing by the mean will result in a biased estimator.
* We can use Monte Carlo simulations to get a rough idea of the bias of our estimates.

# Slide 3:
* Here, I simulated 1000 times a dataset, applied the two missingness mechanisms and computed the bias of the estimate. 
* We can clearly see that imputation by the mean is unbiased for MCAR data, and has a large bias for truncated data.
  
* If we KNOW the missigness mechanism, we **estimate** and correct for the bias. In our example, if we know how the data are truncated, then we can generate some data, truncate accordingly, and observe the effect of the truncation on the bias of our estimate, and then correct for it. **This is precisely our strategy.**

For that, we need a **Generative model** for the data.

# Slide 4: 
For that, we not to model not only the mean, but also the variance of the data, and since we're in multivariate data, it is useful to model the covariance, too. 

We'll use factor analysis or probabilistic PCA to model the multivariate data, but you could similarly use a **Variational Autoencoder** or a **GAN**, but you'll have to be extra careful not to overfit the data. Our proposed models are parameterized by a loading matrix lambda, the mean of the responses mu and their variance psi. 

It is assumed that the covariance of the data can be accounted for by a low-dimensional, unobserved, factor $Z$.

There are some more conditional independence assumption on the right-hand side but don't worry about that, suffices to say that this is a useful model large datasets, and perform dimension reduction.

Here, imputation by the mean will not only affect mu, but also the other parameters: regardless, they can all be corrected using our strategy and the **knowledge of the missingness mechanism**.

# Slide 5:
Here we have the same experiment as before with the de-biased estimator using this strategy.

* Big caveat: in practice we hardly ever know the Missingness mechanism. If we don't know it, we can estimate it. That's our final model and my contribution to this problem.

# Slide 6:
* The idea is to model not only the covariance between the data, but also between the data AND the mask that encodes the missingness. We now modify the model to account for this. 
* In factor analysis, we used the variable Z to model the dependence within the data: now, it also encodes the dependence between the mask and the variable Y.
* Here, the ones and zeroes of the mask are modeled as bernoulli random variables.
* This model is far from easy to estimate: it has been a big part of my research: the goal being to handle large datasets with not three, but tens of thausands of variables.

# Slide 7
* Let's go back to the experiment with truncated data. The rightmost panel shows the result when the mechanism is estimated from the incomplete data. We can see that it can be used without any harm on MCAR, and that, for the truncated data,
* it is able to correct for most of the bias of the estimate. Of course, we can never hope to do as well as when the mechanism is known.

# That's it.
