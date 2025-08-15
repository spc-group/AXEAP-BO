Sequential vs. Simultaneous Optimization
======================================

When fitting crystal field parameters to experimental X-ray emission spectra, you have a fundamental choice in optimization strategy: should parameters be optimized one at a time (sequentially) or all together (simultaneously)? This choice significantly impacts both the quality of your results and the computational efficiency of the fitting process.

This guide explores both approaches, their trade-offs, and provides guidance on when to use each strategy for crystal field parameter estimation.

Overview of Optimization Strategies
-----------------------------------

Sequential Optimization
~~~~~~~~~~~~~~~~~~~~~~~

Sequential optimization involves optimizing parameters one at a time or in small groups, typically holding other parameters fixed during each optimization step.

**Basic approach:**
    1. Fix all parameters except one (or a small subset)
    2. Optimize the unfixed parameter(s) using your chosen method
    3. Fix the newly optimized parameter(s) and move to the next
    4. Repeat until all parameters have been optimized
    5. Optionally, perform multiple passes through all parameters

Simultaneous Optimization
~~~~~~~~~~~~~~~~~~~~~~~~

Simultaneous optimization treats all parameters as variables in a single, multi-dimensional optimization problem.

**Basic approach:**
    1. Define the full parameter vector **x** = [10Dq, B, C, ζ, ...]
    2. Optimize all parameters together in one unified optimization procedure
    3. The optimizer explores the full parameter space simultaneously

Advantages and Disadvantages
---------------------------

Sequential Optimization
~~~~~~~~~~~~~~~~~~~~~~~

**Advantages:**

*Computational Simplicity*
    Each optimization step is lower-dimensional, making it easier to visualize and understand the optimization landscape.

*Reduced Memory Requirements*
    Gaussian processes scale as O(n³) with the number of observations. Sequential optimization in 1D requires much less memory than simultaneous optimization in high dimensions.

*Interpretability*
    You can observe how each parameter individually affects the spectral fit, providing physical insight.

*Debugging and Validation*
    Easier to identify problematic parameters or understand why optimization is failing.

*Parameter Constraints*
    Simpler to enforce physical constraints when dealing with one parameter at a time.

**Disadvantages:**

*Parameter Coupling Ignored*
    Crystal field parameters are often strongly correlated. Sequential optimization cannot capture these correlations effectively.

*Suboptimal Solutions*
    May converge to local minima that would be avoided in simultaneous optimization.

*Order Dependence*
    Results may depend on the order in which parameters are optimized.

*Inefficient Sampling*
    May require more total function evaluations to reach convergence.

*False Convergence*
    May appear to converge when parameters are individually optimal but the combination is suboptimal.

Simultaneous Optimization
~~~~~~~~~~~~~~~~~~~~~~~~

**Advantages:**

*Global Optimality*
    Better chance of finding the global minimum by exploring the full parameter space.

*Parameter Correlations*
    Naturally accounts for correlations and trade-offs between parameters.

*Theoretical Rigor*
    Mathematically sound approach that treats the problem as it truly is: multi-dimensional.

*Efficient Exploration*
    Bayesian optimization can efficiently explore the multi-dimensional space using acquisition functions.

*No Order Dependence*
    Results are independent of any arbitrary parameter ordering.

**Disadvantages:**

*Computational Complexity*
    Higher-dimensional optimization is computationally more expensive, especially for Gaussian processes.

*Curse of Dimensionality*
    As the number of parameters increases, the volume of the parameter space grows exponentially.

*Less Interpretable*
    Harder to understand individual parameter contributions to the fit.

*Hyperparameter Sensitivity*
    GP kernel hyperparameters become more critical and harder to tune in high dimensions.

*Visualization Challenges*
    Cannot easily visualize or understand the optimization landscape.

Parameter Coupling in Crystal Field Theory
------------------------------------------

Understanding Parameter Relationships
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Crystal field parameters exhibit complex interdependencies that sequential optimization may miss:

**10Dq and Racah Parameters (B, C)**
    The crystal field splitting (10Dq) and electron-electron repulsion parameters are often correlated. Changes in 10Dq can be partially compensated by adjustments in B and C.

**Spin-Orbit Coupling (ζ)**
    Spin-orbit coupling affects peak positions and intensities in ways that can trade off with crystal field effects.

**Covalency Parameters**
    In ligand field theory, parameters describing covalent bonding can correlate with both crystal field and Racah parameters.

Mathematical Perspective
~~~~~~~~~~~~~~~~~~~~~~~

Consider a simplified 2-parameter case where the objective function is:

.. math::

    f(x_1, x_2) = (I_{\\text{exp}} - I_{\\text{theo}}(x_1, x_2))^2

If parameters are correlated, the contours of constant f are elliptical rather than circular:

.. math::

    f(x_1, x_2) \\approx f_0 + \\mathbf{x}^T \\mathbf{H} \\mathbf{x}

where **H** is the Hessian matrix. Off-diagonal terms in **H** represent parameter coupling.

**Sequential optimization** follows axis-aligned steps, which can be highly inefficient for correlated parameters.

**Simultaneous optimization** can follow the natural gradient direction, accounting for parameter correlations.

Hybrid Approaches
-----------------

Block Sequential Optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Optimize parameters in physically meaningful groups:

**Crystal Field Block**
    Optimize [10Dq, Dt, Ds] together while holding Racah parameters fixed

**Electronic Repulsion Block**
    Optimize [B, C] together while holding crystal field parameters fixed

**Complete Parameter Block**
    After block optimization, perform simultaneous optimization on all parameters

Hierarchical Optimization
~~~~~~~~~~~~~~~~~~~~~~~~

Use different strategies at different stages:

1. **Coarse Sequential Phase**: Get rough parameter estimates quickly
2. **Fine Simultaneous Phase**: Polish the solution with full optimization
3. **Validation Phase**: Check results with additional sequential refinement

Adaptive Strategy Selection
~~~~~~~~~~~~~~~~~~~~~~~~~~

Choose strategy based on problem characteristics:

.. code-block:: python

    def choose_optimization_strategy(n_params, budget, correlations):
        if n_params <= 3:
            return "simultaneous"
        elif budget < 50:
            return "sequential"  
        elif correlations > 0.7:
            return "simultaneous"
        else:
            return "hybrid"

Practical Recommendations
-------------------------

When to Use Sequential Optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Exploration and Understanding**
    When you're learning about your system and want to understand individual parameter effects.

**Limited Computational Budget**
    When you have very few function evaluations available (<50 total).

**High-Dimensional Problems**
    When you have more than 10 parameters and simultaneous optimization becomes intractable.

**Debugging**
    When your simultaneous optimization is failing and you need to identify problematic parameters.

**Well-Separated Parameters**
    When you have physical or chemical reasons to believe parameters are weakly coupled.

When to Use Simultaneous Optimization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Production Fitting**
    When you need the best possible fit quality and have sufficient computational resources.

**Known Parameter Coupling**
    When you know or suspect strong correlations between parameters.

**Moderate Dimensionality**
    For problems with 3-8 parameters, where simultaneous optimization is computationally feasible.

**Multiple Spectra**
    When fitting multiple related spectra where parameter correlations are expected.

**High-Quality Initial Guess**
    When you have good starting parameter estimates from previous work or physical intuition.

Implementation Considerations
----------------------------

Sequential Implementation Details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Parameter Ordering Strategies:**

*Physical Significance*
    Start with the most important parameters (e.g., 10Dq first)

*Sensitivity Analysis*
    Optimize parameters in order of their spectral sensitivity

*Iterative Refinement*
    Make multiple passes through all parameters

**Convergence Criteria:**

.. math::

    \\max_i |x_i^{(k+1)} - x_i^{(k)}| < \\epsilon

where k is the iteration number and ε is the tolerance.

Simultaneous Implementation Details
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Initialization Strategies:**

*Latin Hypercube Sampling*
    Generate diverse initial points across the parameter space

*Sequential Bootstrap*
    Use sequential optimization results as starting points

*Multi-Start Approach*
    Run simultaneous optimization from multiple initial points

**Dimensionality Reduction:**

*Principal Component Analysis*
    Identify the most important parameter combinations

*Active Subspace Methods*
    Find low-dimensional parameter subspaces that capture most spectral variation

Case Study: Transition Metal Complex
------------------------------------

Consider fitting parameters for a Cr³⁺ complex with the following parameters:
- 10Dq (crystal field splitting)
- B, C (Racah parameters)  
- ζ (spin-orbit coupling)

Sequential Approach Results
~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    Iteration 1: 10Dq = 18,500 cm⁻¹, B = 1,030 cm⁻¹, C = 3,850 cm⁻¹, ζ = 230 cm⁻¹
    Final RMSE: 0.15
    Function evaluations: 45
    Compute time: 12 minutes

Simultaneous Approach Results
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

::

    Final: 10Dq = 18,200 cm⁻¹, B = 1,050 cm⁻¹, C = 3,900 cm⁻¹, ζ = 225 cm⁻¹  
    Final RMSE: 0.09
    Function evaluations: 67
    Compute time: 25 minutes

The simultaneous approach achieved better fit quality but required more evaluations and time.

Best Practices
--------------

General Guidelines
~~~~~~~~~~~~~~~~~

**Start Simple**
    Begin with sequential optimization to understand your system, then move to simultaneous if needed.

**Monitor Correlations**
    Calculate parameter correlation matrices to understand coupling strength.

**Use Physical Intuition**
    Let crystal field theory guide your optimization strategy choices.

**Validate Results**
    Cross-check results between different optimization approaches.

**Document Strategy**
    Keep detailed records of which approach works best for different types of systems.

Quality Metrics
~~~~~~~~~~~~~~

**Spectral Fit Quality**
    RMSE, R², χ² between experimental and theoretical spectra

**Parameter Uncertainty**
    Bootstrap or Bayesian confidence intervals

**Physical Reasonableness**
    Do the optimized parameters make chemical sense?

**Robustness**
    How sensitive are results to different initial conditions?

Conclusion
----------

The choice between sequential and simultaneous optimization represents a fundamental trade-off between computational efficiency and fit quality. Sequential optimization offers simplicity and interpretability, making it ideal for exploration and understanding. Simultaneous optimization provides superior results for coupled parameter systems but at greater computational cost.

For crystal field parameter estimation, the optimal strategy often depends on:
- The number of parameters (dimensionality)
- Available computational resources
- Degree of parameter coupling
- Required fit quality

A pragmatic approach combines both strategies: use sequential optimization for initial exploration and parameter understanding, then apply simultaneous optimization for final, high-quality results. This hybrid strategy leverages the strengths of both approaches while mitigating their individual limitations.

The key is to remain flexible and choose the approach that best matches your specific problem requirements, computational constraints, and accuracy goals.

Further Reading
--------------

- Nocedal, J. & Wright, S. J. "Numerical Optimization." Springer (2006).
- Powell, M. J. D. "A Direct Search Optimization Method That Models the Objective and Constraint Functions by Linear Interpolation." Advances in Optimization and Numerical Analysis (1994).
- Conn, A. R., Scheinberg, K. & Vicente, L. N. "Introduction to Derivative-Free Optimization." SIAM (2009).
