# Package index

## Preparing Data

Prepare your data for ploidal estimation with nQuack.

- [`prepare_data()`](http://mlgaynor.com/nQuack/reference/prepare_data.md)
  : Prepare Data - Step 1
- [`process_data()`](http://mlgaynor.com/nQuack/reference/process_data.md)
  : Process Data - Step 2
- [`process_nquire()`](http://mlgaynor.com/nQuack/reference/process_nquire.md)
  : Use nQuire's Data
- [`denoise_data()`](http://mlgaynor.com/nQuack/reference/denoise_data.md)
  : Denoise Data
- [`Bclean()`](http://mlgaynor.com/nQuack/reference/Bclean.md) : Remove
  Noise with the Beta Distribution

## Simulating data

Simulate data for a given ploidal level.

- [`sim.ind.BB()`](http://mlgaynor.com/nQuack/reference/sim.ind.BB.md) :
  Simulate Allele Counts for Single Individual - Beta-Binomial
  Distribution
- [`sim.ind.BB.tau()`](http://mlgaynor.com/nQuack/reference/sim.ind.BB.tau.md)
  : Simulate Allele Counts for Single Individual - Beta-Binomial
  Distribution with Overdispersion and Error
- [`sim.ind.simple()`](http://mlgaynor.com/nQuack/reference/sim.ind.simple.md)
  : Simulate Allele Counts for Single Individual - Simple Approach

## Basic Mixture Models

Wrapper functions to run model selection for diploids, triploids,
tetraploids, pentaploids, and hexaploids.

- [`quackNormal()`](http://mlgaynor.com/nQuack/reference/quackNormal.md)
  : Model Selection - Expectation Maximization - Normal Mixture
- [`quackNormalNQ()`](http://mlgaynor.com/nQuack/reference/quackNormalNQ.md)
  : Model Selection - Expectation Maximization - Normal Mixture (nQuire)
- [`quackBeta()`](http://mlgaynor.com/nQuack/reference/quackBeta.md) :
  Model Selection - Expectation Maximization - Beta Mixture
- [`quackBetaBinom()`](http://mlgaynor.com/nQuack/reference/quackBetaBinom.md)
  : Model Selection - Expectation Maximization - Beta-Binomial Mixture

## More on Mixture Models

Wrapper functions for model selection, running specific sets of mixture
models, and bootstrap replication.

- [`quackit()`](http://mlgaynor.com/nQuack/reference/quackit.md) : Model
  Selection - Based on BIC or Log-Likelihood
- [`bestquack()`](http://mlgaynor.com/nQuack/reference/bestquack.md) :
  Model Selection - Expectation Maximization - Optimal Distribution and
  Type
- [`quackNboots()`](http://mlgaynor.com/nQuack/reference/quackNboots.md)
  : Bootstrapping - Expectation Maximization - Optimal Distribution and
  Type

## Expectation Maximization Base Functions

Functions used within our mixture models

- [`emstepN()`](http://mlgaynor.com/nQuack/reference/emstepN.md) :
  Expectation maximization - Normal Distribution
- [`emstepNU()`](http://mlgaynor.com/nQuack/reference/emstepNU.md) :
  Expectation maximization - Normal and Uniform Distribution
- [`emstepNA()`](http://mlgaynor.com/nQuack/reference/emstepNA.md) :
  Expectation maximization - Normal Distribution
- [`emstepNUA()`](http://mlgaynor.com/nQuack/reference/emstepNUA.md) :
  Expectation maximization - Normal Distribution
- [`emstepB()`](http://mlgaynor.com/nQuack/reference/emstepB.md) :
  Expectation maximization - Beta Distribution
- [`emstepBU()`](http://mlgaynor.com/nQuack/reference/emstepBU.md) :
  Expectation maximization - Beta and Uniform Distributions
- [`emstepBB()`](http://mlgaynor.com/nQuack/reference/emstepBB.md) :
  Expectation maximization - Beta-Binomial Distribution
- [`emstepBBU()`](http://mlgaynor.com/nQuack/reference/emstepBBU.md) :
  Expectation maximization - Beta-Binomial and Uniform Distributions
- [`estepB3()`](http://mlgaynor.com/nQuack/reference/estepB3.md) :
  E-Step for Expectation Maximization - Beta + Beta + Beta Distribution
- [`emstepB3()`](http://mlgaynor.com/nQuack/reference/emstepB3.md) :
  Expectation maximization - Beta + Beta + Beta Distribution

## Helper functions

Everything else.

- [`SetupBasicExample()`](http://mlgaynor.com/nQuack/reference/SetupBasicExample.md)
  : Setup Basic Example
- [`alphabetacalc()`](http://mlgaynor.com/nQuack/reference/alphabetacalc.md)
  : Calculate Alpha and Beta from Mean and Variance
- [`alphabetacalctau()`](http://mlgaynor.com/nQuack/reference/alphabetacalctau.md)
  : Calculate Alpha and Beta from Mean, Tau, and Error rate.
- [`alphabetacalcvec()`](http://mlgaynor.com/nQuack/reference/alphabetacalcvec.md)
  : Vector-based - Calculate Alpha and Beta from Mean and Variance
- [`alphabetacalctauvec()`](http://mlgaynor.com/nQuack/reference/alphabetacalctauvec.md)
  : Vector-based - Calculate Alpha and Beta from Mean, Tau, and Error
  rate.
- [`muvarcalcvec()`](http://mlgaynor.com/nQuack/reference/muvarcalcvec.md)
  : Variance calculation from Mean, Tau, and Sequencing Error
- [`process_rcpp()`](http://mlgaynor.com/nQuack/reference/process_rcpp.md)
  : Data Preparation - Matrix Filtering
- [`nQuire_reformat()`](http://mlgaynor.com/nQuack/reference/nQuire_reformat.md)
  : Data Preparation - Use nQuire's Data
- [`setconvert()`](http://mlgaynor.com/nQuack/reference/setconvert.md) :
  Calculate Variance from Mean, Tau, and Sequencing Error
- [`resample_xm()`](http://mlgaynor.com/nQuack/reference/resample_xm.md)
  : Calculate Alpha and Beta from Mean and Variance
