Bayesian Optimization Fundamentals
==================================

Bayesian optimization is a powerful technique for optimizing expensive-to-evaluate functions, making it ideal for crystal field parameter estimation where each theoretical spectrum calculation can be computationally costly. This guide introduces the key concepts behind Bayesian optimization and explains why it's particularly well-suited for spectroscopic parameter fitting.

What is Bayesian Optimization?
------------------------------

Bayesian optimization is a sequential model-based approach for finding the global optimum of a black-box function. Unlike traditional optimization methods that may require many function evaluations, Bayesian optimization is designed to find good solutions with as few evaluations as possible.

The core idea is to:

1. **Build a probabilistic model** of the objective function based on observed data points
2. **Use this model to predict** where the next evaluation should occur
3. **Update the model** with each new observation
4. **Repeat** until convergence or budget exhaustion

Why Use Bayesian Optimization for Spectroscopy?
-----------------------------------------------

Crystal field parameter estimation presents several challenges that make Bayesian optimization particularly attractive:

**Expensive Function Evaluations**
    Computing theoretical X-ray emission spectra requires solving complex quantum mechanical calculations. Each parameter set evaluation can take seconds to minutes, making exhaustive search impractical.

**Black-box Nature**
    The relationship between crystal field parameters and spectral features is highly nonlinear and difficult to express analytically. We can evaluate the function (compute a spectrum) but cannot easily compute gradients.

**Noisy Observations**
    Experimental spectra contain noise, and numerical calculations may have small uncertainties, requiring robust optimization approaches.

**Multi-modal Landscapes**
    The parameter space often contains multiple local minima, requiring global optimization techniques.

**Limited Budget**
    In practice, we want to find good parameters with as few spectrum calculations as possible.

Key Components
--------------

Gaussian Process Surrogate Model
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

At the heart of Bayesian optimization is a **Gaussian Process (GP)**, which serves as a probabilistic surrogate model for the true objective function.

**What it provides:**
    - **Mean prediction**: Best estimate of the function value at any point
    - **Uncertainty quantification**: Confidence intervals around predictions
    - **Probabilistic framework**: Full posterior distribution over possible functions

**Why it's useful:**
    The GP allows us to reason about uncertainty. In regions where we have few observations, the uncertainty is high, suggesting these areas might be worth exploring. In well-sampled regions, uncertainty is low, and we can confidently use the model predictions.

**Mathematical foundation:**
    A Gaussian process is defined by its mean function μ(x) and covariance function k(x, x'). For any finite set of points, the function values follow a multivariate Gaussian distribution:

    .. math::

        f(\\mathbf{x}) \\sim \\mathcal{GP}(\\mu(\\mathbf{x}), k(\\mathbf{x}, \\mathbf{x}'))

Acquisition Functions
~~~~~~~~~~~~~~~~~~~~

The acquisition function determines where to sample next by balancing **exploration** (sampling where uncertainty is high) and **exploitation** (sampling where the predicted function value is good).

**Expected Improvement (EI)**
    Measures the expected improvement over the current best observation:

    .. math::

        \\text{EI}(\\mathbf{x}) = \\mathbb{E}[\\max(f(\\mathbf{x}) - f^+, 0)]

    where f^+ is the best observed value so far.

**Upper Confidence Bound (UCB)**
    Balances mean prediction with uncertainty:

    .. math::

        \\text{UCB}(\\mathbf{x}) = \\mu(\\mathbf{x}) + \\kappa \\sigma(\\mathbf{x})

    where κ controls the exploration-exploitation trade-off.

**Probability of Improvement (PI)**
    Probability that a point will improve upon the current best:

    .. math::

        \\text{PI}(\\mathbf{x}) = P(f(\\mathbf{x}) > f^+ + \\xi)

Application to Crystal Field Parameters
--------------------------------------

Problem Formulation
~~~~~~~~~~~~~~~~~~

In crystal field parameter estimation, we want to minimize the difference between experimental and theoretical spectra:

.. math::

    \\mathbf{x}^* = \\arg\\min_{\\mathbf{x}} \\| I_{\\text{exp}}(E) - I_{\\text{theo}}(E; \\mathbf{x}) \\|^2

where:
    - **x** represents the crystal field parameters (e.g., 10Dq, B, C)
    - I_exp(E) is the experimental spectrum
    - I_theo(E; **x**) is the theoretical spectrum for parameters **x**

Parameter Space Considerations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Dimensionality**
    Crystal field problems typically involve 3-10 parameters, which is well-suited for Bayesian optimization.

**Bounds and Constraints**
    Physical constraints (e.g., positive values, Racah parameter relationships) can be incorporated as bounds or constraint functions.

**Scale Sensitivity**
    Parameters may have different scales and sensitivities. Proper normalization or kernel design helps the GP model these relationships effectively.

The Optimization Process
-----------------------

Initial Design
~~~~~~~~~~~~~

Bayesian optimization begins with an initial set of parameter evaluations to train the initial GP model:

**Latin Hypercube Sampling**
    Provides good coverage of the parameter space

**Random Sampling**
    Simple but effective for initial exploration

**Physics-Informed Initialization**
    Use prior knowledge about reasonable parameter ranges

Iterative Refinement
~~~~~~~~~~~~~~~~~~~

1. **Fit GP model** to all observed data points
2. **Optimize acquisition function** to find the most promising next point
3. **Evaluate objective** at the selected point (compute theoretical spectrum)
4. **Update dataset** with new observation
5. **Repeat** until convergence criteria are met

Convergence Criteria
~~~~~~~~~~~~~~~~~~~

**Maximum Iterations**
    Stop after a fixed number of evaluations

**Improvement Threshold**
    Stop when improvement falls below a threshold

**Acquisition Value**
    Stop when the acquisition function value becomes very small

Advantages and Limitations
-------------------------

Advantages
~~~~~~~~~

- **Sample Efficiency**: Finds good solutions with fewer function evaluations than grid search or random search
- **Global Optimization**: Designed to avoid local minima through exploration
- **Uncertainty Quantification**: Provides confidence estimates for the optimal parameters
- **Flexible**: Can incorporate prior knowledge and constraints
- **Robust to Noise**: Gaussian processes naturally handle noisy observations

Limitations
~~~~~~~~~~

- **Computational Overhead**: GP fitting scales as O(n³) with number of observations
- **High-Dimensional Challenges**: Performance may degrade in very high-dimensional spaces (>20 parameters)
- **Hyperparameter Sensitivity**: GP kernel hyperparameters need proper tuning
- **Local Minima in Acquisition**: Acquisition function optimization can itself get stuck

Practical Considerations
-----------------------

Kernel Selection
~~~~~~~~~~~~~~~

The choice of kernel function affects how the GP interpolates between points:

**RBF (Radial Basis Function)**
    Assumes smooth functions; good default choice

**Matérn Kernels**
    More flexible than RBF; can model less smooth functions

**Additive Kernels**
    Useful when parameters have different characteristics

Hyperparameter Optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~

GP hyperparameters (kernel parameters, noise level) significantly impact performance:

**Maximum Likelihood Estimation**
    Standard approach for learning hyperparameters

**Cross-Validation**
    More robust but computationally expensive

**Bayesian Treatment**
    Integrates over hyperparameter uncertainty

Multi-Objective Considerations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In some cases, you might want to optimize multiple objectives simultaneously:

- **Spectral fit quality** vs. **parameter physical reasonableness**
- **Multiple spectroscopic techniques** (XES, XAS, etc.)
- **Computational cost** vs. **accuracy**

Multi-objective Bayesian optimization extends these concepts to Pareto frontier exploration.

Conclusion
----------

Bayesian optimization provides an elegant framework for crystal field parameter estimation by efficiently navigating the complex relationship between parameters and spectral properties. By maintaining a probabilistic model of the objective function, it can make informed decisions about where to sample next, leading to faster convergence than traditional optimization approaches.

The method's ability to balance exploration and exploitation makes it particularly valuable when computational resources are limited, as is often the case in theoretical spectroscopy calculations.

Further Reading
--------------

- Shahriari, B., et al. "Taking the human out of the loop: A review of Bayesian optimization." *Proceedings of the IEEE*, 104.1 (2016): 148-175.
- Frazier, P. I. "A tutorial on Bayesian optimization." *arXiv preprint arXiv:1807.02811* (2018).
- Garnett, R. "Bayesian Optimization." Cambridge University Press (2023).
