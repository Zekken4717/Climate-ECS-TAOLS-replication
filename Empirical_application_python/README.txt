Updated: 2026-06-09
Directory: Empirical_application_python\

==========================================================================
RECOMMENDED RUN ORDER
==========================================================================

  1. Plot_Data_Unit_Root_Normality_Tests.py   (pre-estimation diagnostics)
  2. Empirical_Analysis_Climate_trendType_1_in_mc.py   (run 3 times: model I, II, III)
  3. extract_ranges_and_combine.py   (combined figures + numerical ranges)

All three scripts are Python. They read the same data file (OHC Data.xls).

Python dependencies (all scripts):
  numpy, pandas, scipy, statsmodels, matplotlib, xlrd (>= 2.0.1)
  Install with: pip install numpy pandas scipy statsmodels matplotlib xlrd

==========================================================================
SCRIPT 1 — DIAGNOSTIC TESTS (run first)
==========================================================================

Plot_Data_Unit_Root_Normality_Tests.py

  Pre-estimation diagnostic tests and data plots for the climate
  multicointegration analysis. Tests all four input series:
    TotalRF   (full-efficacy RF)         [Model I]
    MarvelRF  (partial-efficacy RF)      [Models II & III]
    Berkeley  (Berkeley Earth temp)      [Models I & II]
    HADCRUT   (HadCRUT temperature)      [Model III]

  Figures produced:
    Figure 1 -- Time series of all four data series (1 x 2 panel)
                Left:  both RF series (full-efficacy, partial-efficacy)
                Right: both surface temperature series (Berkeley, HadCRUT)
    Figure 2 -- Histograms of AR(3) residuals with fitted normal density;
                JB statistic and p-value annotated on each panel
    Figure 3 -- Normal Q-Q plots of AR(3) residuals

  Tests performed:
    Table I:  ADF unit root on levels (ARD model, BIC lag selection, max 3 lags)
              KPSS stationarity on levels (trend specification)
    Table II: Jarque-Bera normality on AR(3) residuals of first differences
              Model: Dy_t = c + phi_1*Dy_{t-1} + phi_2*Dy_{t-2} + phi_3*Dy_{t-3}
              T_ar = 161 residuals (1854-2014)

  Key findings:
    - All four level series fail to reject ADF and reject KPSS => I(1).
    - RF first differences reject normality even after AR(3) filtering
      (volcanic outliers), motivating TAOLS over Gaussian MLE.
    - Temperature first differences are consistent with normality.

  Outputs saved to Results_Diagnostics\ :
    unit_root_normality_tests.out        -- full printed tables
    rf_temperature_combined.eps/.png     -- time series plot (Figure 1)
    fd_histograms.eps/.png               -- histograms (Figure 2)
    fd_qqplots.eps/.png                  -- Q-Q plots (Figure 3)

  Also copies rf_temperature_combined.eps/.png to ..\Paper_June_2026\ for LaTeX.

  Implementation notes:
    - ADF lag selected by BIC using a custom adf_bic_select() function
      (ARD model: constant + lagged level + lagged differences).
    - KPSS bandwidth: Newey-West (1994) automatic selection (statsmodels default).

==========================================================================
SCRIPT 2 — MAIN ESTIMATION (run once per model)
==========================================================================

Empirical_Analysis_Climate_trendType_1_in_mc.py

  Estimates lambda (climate feedback parameter) and phi (multicointegration
  parameter) using TAOLS, for two nested models:
    (1) Multicointegration model: F on [S, s, ds, trend]
    (2) Regular cointegration model: f on [s, ds, trend]

  Key user settings at the top of the script:
    model            = 'I', 'II', or 'III'  (see table below)
    trend_type_coint = 0 (constant only) or 1 (constant + linear trend)

  Data source: OHC Data.xls, sheet "Data" (annual, 1850-2014)

  Model definitions (following Bruns, Csereklyei & Stern, JOE 2020):
    Model I   -- full-efficacy RF (TotalRF)    + Berkeley Earth temperature
    Model II  -- partial-efficacy RF (MarvelRF) + Berkeley Earth temperature
    Model III -- partial-efficacy RF (MarvelRF) + HadCRUT temperature

  ECS = (5.35 * log(2)) / lambda_hat   [Myhre et al. 1998]
  %   = 100 * 0.31 / phi_hat           [atmospheric heat capacity proxy]

  Figures produced:
    Figure 1    -- TAOLS estimates of lambda and ECS (multicointegration model)
    Figure 2    -- TAOLS estimates of phi and % heat content (multicointegration)
    Figure 3    -- TAOLS estimates of trend coefficient(s) (multicointegration)
    Figure 4    -- TAOLS estimates of lambda and ECS (cointegration model)
    Figure 5    -- TAOLS estimates of trend coefficient (cointegration model,
                   only produced when trend_type_coint = 1)
    fig_ohc     -- Predicted OHC (TAOLS) vs observed OHC series (1940-2014).
                   Compares model-implied cumulative heat content Q_t against
                   Cheng et al. (2017) OHC 0-700m and 0-2000m, and the
                   Marvel et al. (2016) simulated series. Q_t scaled to
                   10^22 Joules and to the top-2000m ocean (x 0.81).

  Run this script three times (setting model = 'I', 'II', 'III' in turn)
  to produce results for all model specifications before running Script 3.

==========================================================================
SCRIPT 3 — COMBINED FIGURES AND NUMERICAL RANGES
==========================================================================

extract_ranges_and_combine.py

  Loads the .npz files produced by Script 2 for all three models and:
    (1) Prints TAOLS numerical ranges for lambda, ECS, phi, and % for each
        model and a cross-model summary.
    (2) Saves combined lambda/ECS figure (all three models, 2 panels) as
        combined_mc_lambda_ECS.png to ..\Paper_June_2026\
    (3) Saves combined phi/percent figure (all three models, 2 panels) as
        combined_mc_phi_percent.png to ..\Paper_June_2026\

  Prerequisite: Script 2 must have been run for models I, II, and III
  (i.e., all six .npz files in Results_Model_*\ must exist).

  Uses RGB colour coding:
    Model I   -- red   (0.8, 0.2, 0.2)
    Model II  -- blue  (0.2, 0.4, 0.8)
    Model III -- green (0.2, 0.7, 0.2)

==========================================================================
OUTPUT STRUCTURE
==========================================================================

Results are saved to a subfolder named:
  Results_Model_<model>_cointTrendType_<trend_type_coint>\

Each run of Script 2 produces the following files inside that subfolder:

  Model_<model>_cointTrendType_<N>.out
      Log of printed summary statistics (lambda, SE, ECS, phi, %).

  <prefix>_trendType_<M>_mc_lambda_ECS.eps / .pdf / .npz
      TAOLS estimates of lambda and ECS from the multicointegration model,
      with 95% CI bands. npz file contains: Ks, TA_OLS_mc_lambda, TA_OLS_mc_ECS.

  <prefix>_trendType_<M>_mc_phi_percent.eps / .pdf / .npz
      TAOLS estimates of phi and % heat content toward surface warming,
      with 95% CI bands. npz file contains: Ks, TA_OLS_mc_phi, TA_OLS_mc_percent.

  <prefix>_trendType_<M>_mc_trend.pdf
      TAOLS estimates of trend coefficient(s) from the multicointegration model.

  <prefix>_trendType_<M>_mc_OHC_comparison.eps / .pdf
      Predicted OHC (TAOLS, all K values, orange dotted) overlaid on
      observed OHC from Cheng et al. (2017) and Marvel et al. (2016).
      1940-2014, units: 10^22 Joules.

  <prefix>_trendType_<N>_coint_lambda_ECS.eps / .pdf / .npz
      TAOLS estimates of lambda and ECS from the regular cointegration model,
      with 95% CI bands. npz file contains: Ks, TA_OLS_lambda, TA_OLS_ECS.

  <prefix>_trendType_<N>_coint_trend.pdf
      TAOLS estimates of the linear trend coefficient from the cointegration
      model (only produced when trend_type_coint = 1).

  Where:
    <prefix> = Full_RF_BerkeleyT    (Model I)
               Partial_RF_BerkeleyT (Model II)
               Partial_RF_HADCRUT   (Model III)
    <M> = trend_type_mcoint = trend_type_coint + 1
    <N> = trend_type_coint

  Default: trend_type_coint = 0  -->  trend_type_mcoint = 1
           (linear trend in multicointegration model;
            constant only in cointegration model)

==========================================================================
CURRENT RESULTS (default settings: trend_type_coint = 0)
==========================================================================

Results_Model_I_cointTrendType_0\
  Full_RF_BerkeleyT_trendType_1_mc_lambda_ECS.eps/.pdf/.npz
  Full_RF_BerkeleyT_trendType_1_mc_phi_percent.eps/.pdf/.npz
  Full_RF_BerkeleyT_trendType_1_mc_OHC_comparison.eps/.pdf
  Full_RF_BerkeleyT_trendType_1_mc_trend.pdf
  Full_RF_BerkeleyT_trendType_0_coint_lambda_ECS.eps/.pdf/.npz
  Model_I_cointTrendType_0.out

Results_Model_II_cointTrendType_0\
  Partial_RF_BerkeleyT_trendType_1_mc_lambda_ECS.eps/.pdf/.npz
  Partial_RF_BerkeleyT_trendType_1_mc_phi_percent.eps/.pdf/.npz
  Partial_RF_BerkeleyT_trendType_1_mc_OHC_comparison.eps/.pdf
  Partial_RF_BerkeleyT_trendType_1_mc_trend.pdf
  Partial_RF_BerkeleyT_trendType_0_coint_lambda_ECS.eps/.pdf/.npz
  Model_II_cointTrendType_0.out

Results_Model_III_cointTrendType_0\
  Partial_RF_HADCRUT_trendType_1_mc_lambda_ECS.eps/.pdf/.npz
  Partial_RF_HADCRUT_trendType_1_mc_phi_percent.eps/.pdf/.npz
  Partial_RF_HADCRUT_trendType_1_mc_OHC_comparison.eps/.pdf
  Partial_RF_HADCRUT_trendType_1_mc_trend.pdf
  Partial_RF_HADCRUT_trendType_0_coint_lambda_ECS.eps/.pdf/.npz
  Model_III_cointTrendType_0.out

Results_Diagnostics\
  unit_root_normality_tests.out
  rf_temperature_combined.eps/.png
  fd_histograms.eps/.png
  fd_qqplots.eps/.png

==========================================================================
OTHER FILES
==========================================================================

OHC Data.xls     -- input data (Ocean Heat Content, RF, temperature series)
                    Sheet "Data": annual observations 1850-2014.
                    Variables: TotalRF, CTotalRF, MarvelRF, CMarvelRF,
                               Berkeley, CBerkeley, HADCRUT, CHadcrut,
                               OHCCheng700  (Cheng et al. 2017, 0-700m, 1940-2014),
                               OHCCheng2000 (Cheng et al. 2017, 0-2000m, 1940-2014),
                               OHCMarvel    (Marvel et al. 2016, 1940-2014).

OHC Figures.xlsx -- supplementary figures workbook


==========================================================================
REFERENCE
==========================================================================

Bruns, S. B., Csereklyei, Z., & Stern, D. I. (2020). A multicointegration
model of global climate change. Journal of Econometrics, 214(1), 175-197.
