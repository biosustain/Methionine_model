# Methionine_model

This repository contains code for case studies in the paper “Bayesian
regression facilitates quantitative modelling of cell metabolism”.

To reproduce our results, run the following shell commands from the root of
this repository:

1. `make build`
2. `make sample_methionine`
3. `make sample_ode`
4. `make validate_example_ode`
5. `make figures`

This may take some time, especially the `sample` commands!

You may hit an error at Step 3, as at this point the code attempts to perform
Laplace approximation for a complex model, which is unstable and generally not
recommended. In fact, the reason we include this method is to make the case
that you shouldn't do it! If this happens we recommend simply repeating the
command until it runs without an error.

