Fitting Theoretical Spectra with Pseudo-Voigt Profiles
=====================================================

A critical step in crystal field parameter estimation is converting stick spectra from theoretical calculations into realistic lineshapes that can be compared with experimental X-ray emission spectra. The pseudo-Voigt profile provides an excellent approximation for the natural lineshapes observed in X-ray spectroscopy, combining the physical realism of Lorentzian profiles with the computational advantages of Gaussian functions.

This guide explains how to implement pseudo-Voigt profile fitting for spectral comparison and optimization.

Understanding Peak Profiles in X-ray Spectroscopy
-------------------------------------------------

Physical Origins of Peak Shapes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Real X-ray emission peaks are never infinitely sharp due to several broadening mechanisms:

**Natural Linewidth (Lorentzian)**
    Core-hole lifetime broadening follows the uncertainty principle. Shorter core-hole lifetimes lead to broader peaks with Lorentzian profiles.

**Instrumental Broadening (Gaussian)**
    Spectrometer resolution, detector effects, and sample-related factors contribute Gaussian broadening.

**Inhomogeneous Broadening**
    Sample heterogeneity, crystal imperfections, and environmental variations add additional Gaussian character.

**Thermal Effects**
    Phonon interactions and thermal motion contribute to both Gaussian and Lorentzian components.

The combination of these effects typically results in peak profiles that are neither purely Gaussian nor purely Lorentzian, but rather a convolution of both.

Mathematical Foundation
-----------------------

Pure Profile Functions
~~~~~~~~~~~~~~~~~~~~~~

**Gaussian Profile**

.. math::

    G(E; E_0, \\sigma) = \\frac{1}{\\sigma\\sqrt{2\\pi}} \\exp\\left(-\\frac{(E - E_0)^2}{2\\sigma^2}\\right)

where E₀ is the peak center and σ is the standard deviation.

**Lorentzian Profile**

.. math::

    L(E; E_0, \\gamma) = \\frac{1}{\\pi} \\frac{\\gamma}{(E - E_0)^2 + \\gamma^2}

where E₀ is the peak center and γ is the half-width at half-maximum (HWHM).

**Voigt Profile**

The true Voigt profile is the convolution of Gaussian and Lorentzian functions:

.. math::

    V(E) = \\int_{-\\infty}^{\\infty} G(E'; \\sigma) L(E - E'; \\gamma) dE'

This integral has no closed-form solution and requires numerical evaluation, making it computationally expensive for optimization.

Pseudo-Voigt Approximation
~~~~~~~~~~~~~~~~~~~~~~~~~~

The pseudo-Voigt function approximates the Voigt profile as a weighted sum:

.. math::

    \\text{pV}(E; E_0, \\text{FWHM}, \\eta) = \\eta L(E; E_0, \\gamma_L) + (1-\\eta) G(E; E_0, \\sigma_G)

where:
    - **η** is the mixing parameter (0 ≤ η ≤ 1)
    - **η = 0**: pure Gaussian
    - **η = 1**: pure Lorentzian
    - **FWHM**: full-width at half-maximum of the pseudo-Voigt profile

**Parameter Relationships**

To maintain consistent FWHM across the profile:

.. math::

    \\gamma_L = \\frac{\\text{FWHM}}{2}

.. math::

    \\sigma_G = \\frac{\\text{FWHM}}{2\\sqrt{2\\ln(2)}} \\approx \\frac{\\text{FWHM}}{2.355}

**Normalized Pseudo-Voigt**

For spectral fitting, we typically use the normalized form:

.. math::

    \\text{pV}_{\\text{norm}}(E) = \\frac{A}{\\text{FWHM}} \\left[ \\frac{\\eta}{\\pi} \\frac{\\gamma_L}{(E - E_0)^2 + \\gamma_L^2} + \\frac{1-\\eta}{\\sigma_G\\sqrt{2\\pi}} \\exp\\left(-\\frac{(E - E_0)^2}{2\\sigma_G^2}\\right) \\right]

where A is the integrated area under the peak.

Implementation for Spectral Comparison
--------------------------------------

From Stick Spectrum to Realistic Spectrum
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Theoretical calculations typically produce "stick spectra" - discrete transitions with specific energies and intensities:

::

    Transition    Energy (eV)    Intensity
    ²E → ⁴A₂        7110.5         0.85
    ²T₁ → ⁴A₂       7112.3         0.12
    ²T₂ → ⁴A₂       7115.8         0.03

**Conversion Process:**

1. **Peak Assignment**: Each calculated transition becomes a pseudo-Voigt peak
2. **Parameter Setting**: Assign FWHM and η values to each peak
3. **Profile Generation**: Compute the pseudo-Voigt profile for each transition
4. **Spectrum Assembly**: Sum all individual profiles to create the complete theoretical spectrum

Mathematical Implementation
~~~~~~~~~~~~~~~~~~~~~~~~~~

For a theoretical spectrum with n transitions:

.. math::

    I_{\\text{theo}}(E) = \\sum_{i=1}^{n} A_i \\cdot \\text{pV}(E; E_i, \\text{FWHM}_i, \\eta_i) + I_{\\text{bg}}(E)

where:
    - Aᵢ is the intensity of transition i
    - Eᵢ is the energy of transition i
    - FWHMᵢ is the linewidth of transition i
    - ηᵢ is the mixing parameter for transition i
    - I_bg(E) is the background function

Profile Parameter Strategies
---------------------------

Fixed Parameters Approach
~~~~~~~~~~~~~~~~~~~~~~~~~

**Uniform Parameters**
    Use the same FWHM and η for all peaks:
    
    - Advantages: Simple, fewer fitting parameters
    - Disadvantages: May not capture energy-dependent broadening

**Energy-Dependent FWHM**
    Model natural linewidth variation with energy:

.. math::

    \\text{FWHM}(E) = \\text{FWHM}_0 + \\alpha E

where FWHM₀ and α are constants determined from tabulated core-hole lifetimes.

**Transition-Specific Parameters**
    Different η values for different types of transitions:
    
    - Main lines: η = 0.3 (more Gaussian due to instrumental broadening)
    - Satellite lines: η = 0.7 (more Lorentzian due to shorter lifetimes)

Variable Parameters Approach
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Fitting Profile Parameters**
    Include FWHM and/or η as additional optimization parameters:

**Advantages:**
    - More flexible fitting
    - Can account for unknown broadening mechanisms
    - May improve spectral agreement

**Disadvantages:**
    - More parameters to optimize
    - Risk of overfitting
    - May compensate for poor crystal field parameters

**Constraint Considerations:**
    - Physical bounds: FWHM > 0, 0 ≤ η ≤ 1
    - Energy scaling: FWHM should increase monotonically with energy
    - Transition type consistency: Similar transitions should have similar parameters

Objective Function Formulation
------------------------------

Spectral Comparison Metrics
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Least Squares (L2 norm)**

.. math::

    \\chi^2 = \\sum_{j=1}^{N} \\left(\\frac{I_{\\text{exp}}(E_j) - I_{\\text{theo}}(E_j)}{\\sigma_j}\\right)^2

where σⱼ is the experimental uncertainty at energy Eⱼ.

**Robust Alternatives**

*L1 norm (less sensitive to outliers):*

.. math::

    \\text{MAE} = \\sum_{j=1}^{N} |I_{\\text{exp}}(E_j) - I_{\\text{theo}}(E_j)|

*Huber loss (combines L1 and L2):*

.. math::

    L_{\\text{Huber}} = \\begin{cases}
    \\frac{1}{2}(I_{\\text{exp}} - I_{\\text{theo}})^2 & \\text{if } |I_{\\text{exp}} - I_{\\text{theo}}| \\leq \\delta \\\\
    \\delta |I_{\\text{exp}} - I_{\\text{theo}}| - \\frac{1}{2}\\delta^2 & \\text{otherwise}
    \\end{cases}

Data Preprocessing Considerations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Energy Alignment**
    Ensure experimental and theoretical energy scales are properly aligned:

.. math::

    E_{\\text{aligned}} = E_{\\text{raw}} + \\Delta E_{\\text{shift}}

**Intensity Scaling**
    Account for experimental factors affecting absolute intensities:

.. math::

    I_{\\text{scaled}} = s \\cdot I_{\\text{raw}} + I_{\\text{offset}}

**Background Subtraction**
    Remove experimental background before comparison:

.. math::

    I_{\\text{corrected}} = I_{\\text{raw}} - I_{\\text{background}}

Advanced Considerations
----------------------

Peak Shape Analysis
~~~~~~~~~~~~~~~~~~

**Asymmetric Profiles**
    Real peaks may show asymmetry due to:
    - Post-collision interaction effects
    - Crystal field splitting beyond the model
    - Sample charging or energy referencing issues

*Asymmetric pseudo-Voigt:*

.. math::

    \\text{FWHM}(E) = \\begin{cases}
    \\text{FWHM}_L & \\text{if } E < E_0 \\\\
    \\text{FWHM}_R & \\text{if } E > E_0
    \\end{cases}

**Peak Correlation Effects**
    Nearby peaks may show interference patterns or correlation effects that require careful modeling.

Multi-Component Fitting
~~~~~~~~~~~~~~~~~~~~~~

**Multiplet Structure**
    Crystal field calculations often predict closely spaced transitions that appear as a single broadened peak:

.. math::

    I_{\\text{multiplet}}(E) = \\sum_{k} A_k \\cdot \\text{pV}(E; E_0 + \\Delta E_k, \\text{FWHM}, \\eta)

**Shake-up/Shake-off Features**
    Additional satellite structure may require extra pseudo-Voigt components.

Computational Implementation
---------------------------

Efficient Calculation
~~~~~~~~~~~~~~~~~~~~

**Vectorized Implementation**
    For computational efficiency, implement pseudo-Voigt calculation using vectorized operations:

.. code-block:: python

    def pseudo_voigt(E, E0, FWHM, eta, A):
        """
        Compute pseudo-Voigt profile
        E: energy array
        E0: peak center
        FWHM: full width at half maximum
        eta: mixing parameter (0=Gaussian, 1=Lorentzian)
        A: peak area
        """
        gamma = FWHM / 2
        sigma = FWHM / (2 * np.sqrt(2 * np.log(2)))
        
        # Lorentzian component
        lorentzian = (1/np.pi) * gamma / ((E - E0)**2 + gamma**2)
        
        # Gaussian component  
        gaussian = (1/(sigma * np.sqrt(2*np.pi))) * np.exp(-0.5 * ((E - E0)/sigma)**2)
        
        # Pseudo-Voigt combination
        profile = eta * lorentzian + (1 - eta) * gaussian
        
        return A * profile / FWHM

**Energy Grid Considerations**
    - Use fine energy grids (0.1-0.5 eV spacing) for accurate profile representation
    - Extend grid beyond experimental range to avoid edge effects
    - Consider adaptive grids with higher density near peaks

Optimization Integration
~~~~~~~~~~~~~~~~~~~~~~~

**Parameter Vector Structure**
    When including profile parameters in optimization:

.. code-block:: python

    # Parameter vector organization
    params = [
        # Crystal field parameters
        10*Dq, B, C, zeta,
        # Profile parameters (optional)
        FWHM_global, eta_global,
        # Alignment parameters (optional)
        energy_shift, intensity_scale
    ]

**Gradient Information**
    For gradient-based optimizers, analytical derivatives of pseudo-Voigt profiles can significantly improve convergence speed.

Quality Assessment and Validation
---------------------------------

Goodness-of-Fit Metrics
~~~~~~~~~~~~~~~~~~~~~~

**Statistical Measures**

*R-squared:*

.. math::

    R^2 = 1 - \\frac{\\sum_i (I_{\\text{exp},i} - I_{\\text{theo},i})^2}{\\sum_i (I_{\\text{exp},i} - \\bar{I}_{\\text{exp}})^2}

*Reduced chi-squared:*

.. math::

    \\chi_{\\text{red}}^2 = \\frac{\\chi^2}{N - p}

where N is the number of data points and p is the number of fitted parameters.

**Physical Reasonableness**
    - Do the fitted FWHM values match expected core-hole lifetimes?
    - Are η values consistent with the expected Gaussian/Lorentzian balance?
    - Do peak positions match transition energies from theory?

Residual Analysis
~~~~~~~~~~~~~~~~

**Systematic Patterns**
    Look for systematic deviations in residuals (experimental - theoretical):
    - Energy-dependent trends may indicate incorrect profile parameters
    - Peak-specific patterns may suggest missing transitions or incorrect assignments

**Statistical Tests**
    - Durbin-Watson test for autocorrelation in residuals
    - Runs test for randomness
    - Normality tests (Shapiro-Wilk) for residual distribution

Best Practices and Recommendations
----------------------------------

Profile Parameter Selection
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Starting Values**
    - FWHM: Use tabulated core-hole lifetimes as starting estimates
    - η: Start with 0.5 (equal Gaussian/Lorentzian mix) and refine based on experimental comparison

**Physical Constraints**
    - Enforce monotonic increase of FWHM with energy
    - Limit η to physically reasonable ranges (typically 0.2-0.8 for X-ray emission)
    - Use prior knowledge about instrumental resolution

**Sensitivity Analysis**
    Test how sensitive your crystal field parameters are to profile parameter choices:

.. code-block:: python

    # Pseudo-code for sensitivity analysis
    for eta in [0.2, 0.4, 0.6, 0.8]:
        for FWHM in [1.0, 1.5, 2.0, 2.5]:  # eV
            optimized_params = optimize_crystal_field(eta, FWHM)
            sensitivity_matrix[eta, FWHM] = optimized_params

Common Pitfalls and Solutions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Over-parameterization**
    *Problem:* Including too many profile parameters can lead to overfitting
    *Solution:* Use physical constraints and cross-validation to prevent overfitting

**Profile Parameter Compensation**
    *Problem:* Profile parameters may compensate for incorrect crystal field parameters
    *Solution:* Use independent constraints on profile parameters from known physics

**Energy Calibration Issues**
    *Problem:* Incorrect energy alignment can lead to poor fits regardless of profile quality
    *Solution:* Implement robust energy alignment procedures using reference peaks

**Background Treatment**
    *Problem:* Inadequate background modeling affects both profile fitting and parameter estimation
    *Solution:* Use appropriate background functions (polynomial, exponential, etc.) and fit simultaneously

Conclusion
----------

Pseudo-Voigt profiles provide an excellent balance between physical realism and computational efficiency for fitting theoretical X-ray emission spectra to experimental data. The key to successful implementation lies in:

1. **Understanding the physical origins** of peak broadening in your system
2. **Choosing appropriate profile parameters** based on known spectroscopic properties
3. **Implementing efficient computational algorithms** for optimization
4. **Validating results** through residual analysis and physical reasonableness checks
5. **Avoiding over-parameterization** that can mask underlying crystal field parameter issues

When properly implemented, pseudo-Voigt profile fitting enables accurate comparison between theory and experiment, forming the foundation for reliable crystal field parameter estimation through Bayesian optimization or other advanced fitting techniques.

The choice of profile parameters should always be guided by physical understanding rather than purely empirical fitting considerations. This ensures that the resulting crystal field parameters reflect genuine electronic structure properties rather than artifacts of the fitting procedure.
