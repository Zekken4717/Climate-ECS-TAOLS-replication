% =========================================================================
% Empirical_Analysis_Climate_trendType_1_in_mc.m
%
% Purpose:
%   Estimates the Equilibrium Climate Sensitivity (ECS) and the
%   multicointegration parameter phi using the Transformed and Augmented OLS (TAOLS)
%   applied to the multicointegration framework of Bruns, Csereklyei, and Stern (JOE 2020).
%
%   Two models are estimated sequentially:
%     (1) Multicointegration model  -- regresses cumulated RF (F) on
%         cumulated temperature (S), temperature (s), and its first
%         difference (ds), and a trend, recovering lambda and phi.
%     (2) Regular cointegration model -- regresses RF (f) on temperature
%         (s), its first difference (ds), and a (lower order) trend, recovering lambda.
%
%   ECS is computed from the TAOLS estimates as:
%     ECS = (5.35 * log(2)) / lambda_hat
%
%   phi is converted to a percentage of total heat content directed toward
%   surface warming as:
%     % = 100 * 0.31 / phi_hat
%
% Data:
%   'OHC Data.xls', sheet 'Data' -- annual data from 1850 to 2014.
%   Variables: TotalRF, CTotalRF, MarvelRF, CMarvelRF,
%              Berkeley, CBerkeley, HADCRUT, CHadcrut.
%
% Key user settings (lines below):
%   model            -- 'I', 'II', or 'III' (data/variable selection)
%   trend_type_coint -- 0 (constant only) or 1 (linear trend)
%   trend_type_mcoint -- automatically set to coint trend + 1 (linear or quadratic trend)
%
% Outputs (saved to results_dir subfolder):
%   <filename>_lambda_ECS.fig/.eps/.mat  -- TAOLS of lambda and ECS
%   <filename>_phi_percent.fig/.eps/.mat -- TAOLS of phi and percentage towards surface warming
%   Model_*_cointTrendType_*.out   -- diary of printed summary statistics
%
% Reference:
%   Bruns, S. B., Csereklyei, Z., & Stern, D. I. (2020). A multicointegration
%   model of global climate change. Journal of Econometrics, 214(1), 175-197.
%
% =========================================================================

clear;
close all;

cd(fileparts(mfilename('fullpath')));  % set working directory to the folder containing this script

% -------------------------------------------------------------------------
% USER SETTINGS
% -------------------------------------------------------------------------

% Trend specification for the regular cointegration model:
%   0 = constant only (default)
%   1 = constant + linear trend
trend_type_coint = 0;

% Trend specification for the multicointegration model is derived
% automatically: mcoint trend = coint trend + 1
%   1 = linear trend (default, paired with trend_type_coint = 0)
%   2 = quadratic trend (paired with trend_type_coint = 1)
trend_type_mcoint = trend_type_coint + 1;

% Data / variable selection (following Table 3 of Bruns et al. JOE 2020):
%   'I'   -- full-efficacy RF (TotalRF) and Berkeley Earth temperature
%   'II'  -- partial-efficacy RF (MarvelRF) and Berkeley Earth temperature
%   'III' -- partial-efficacy RF (MarvelRF) and HadCRUT temperature
model = 'I';

% -------------------------------------------------------------------------
% PHYSICAL CONSTANTS (do not modify)
% -------------------------------------------------------------------------
% Radiative forcing from a doubling of CO2 [W/m^2], from Myhre et al. (1998)
% eq. 3: DeltaF = 5.35 * ln(C/C0); doubling sets C/C0 = 2.
DELTA_F2X = 5.35 * log(2);

% Atmospheric heat capacity per unit Earth surface area [W-yr/m^2/degC].
% Derived from: mass of atmosphere (~5.15e18 kg) x cp of air (~1005 J/kg/K)
% divided by Earth surface area (5.1e14 m^2) and seconds per year (3.156e7).
% Used to convert phi estimates to the % of total heat content that warms
% the atmosphere: % = 100 * HC_ATMOS / phi_hat.
% See Bruns, Csereklyei & Stern (2020, p. 185).
HC_ATMOS = 0.31;

results_dir = ['Results_Model_' model '_cointTrendType_' num2str(trend_type_coint)];
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

diary_filename =  fullfile(results_dir, ['Model_' model '_cointTrendType_' num2str(trend_type_coint) '.out']);

diary(diary_filename)
diary on
% cleanup_diary = onCleanup(@() diary('off'));  % ensures diary is closed even if script errors

% Output the values to the diary file
disp('coint Trend Type:');
disp(trend_type_coint);

disp('Mcoint Trend Type:');
disp(trend_type_mcoint);

disp('Model:');
disp(model);

% Data file must be in the same directory as this script.
% (The cd() call above already sets the working directory to that location.)
data_fileName = 'OHC Data.xls';

% Load data from the 'Data' sheet
data = readtable(data_fileName, 'Sheet', 'Data');



%% Load variables from the data table
% Lowercase letters denote levels; uppercase letters denote cumulated sums.
% Suffix _m = partial-efficacy (Marvel et al.) version; _h = HadCRUT version.

f   = data.TotalRF;    % full-efficacy radiative forcing (Appendix A, Bruns et al. 2020)
F   = data.CTotalRF;   % cumulated full-efficacy RF

f_m = data.MarvelRF;   % partial-efficacy RF (Marvel et al. efficacy correction)
F_m = data.CMarvelRF;  % cumulated partial-efficacy RF

s   = data.Berkeley;   % global surface temperature -- Berkeley Earth
S   = data.CBerkeley;  % cumulated Berkeley Earth temperature

s_h = data.HADCRUT;    % global surface temperature -- HadCRUT (Hadley Centre / CRU)
S_h = data.CHadcrut;   % cumulated HadCRUT temperature

OHC700 = data.OHCCheng700;  % ocean heat content (Cheng et al. 2017)
OHC700 = OHC700(1940-1850+1:end);

OHC2000 = data.OHCCheng2000;  % ocean heat content (Cheng et al. 2017)
OHC2000 = OHC2000(1940-1850+1:end);

OHCMarvel = data.OHCMarvel;  % ocean heat content (Marvel et al. 2016)

OHCMarvel = OHCMarvel(1940-1850+1:end);  
OHCMarvel =  OHCMarvel - OHCMarvel(1);


year  = (1940:1:2014)'; % for plotting OHC

figure(1000)
  handle = plot(year, OHC2000,'k-',...
        year, OHC700,'ro-', ....
        year, OHCMarvel, 'b+-');

legend({'Cheng et al. (2017) OHC 0-2000m', ...
'Cheng et al. (2017) OHC 0-700m', 'Marvel et al. (2016)'}, 'Location', 'northwest');
set(handle, 'LineWidth', 1.5, 'markerSize', 3);
ylabel('10^{22} Joules');
xlabel('Year');
xticks(1940:10:2010);
grid on;
box on;


%% Pre-process the main data

T = length(f);  % full sample size (1850–2014 = 165 observations)

% First differences of temperature (needed as a regressor in both models)
ds   = s(2:T,1)   - s(1:T-1,1);
ds_h = s_h(2:T,1) - s_h(1:T-1,1);

% Drop the first observation from all level/cumulated series so that
% all series are aligned with the first-differenced temperature (ds, ds_h)
F   = F(2:T,1);    f   = f(2:T,1);
S   = S(2:T,1);    s   = s(2:T,1);
F_m = F_m(2:T,1);  f_m = f_m(2:T,1);
S_h = S_h(2:T,1);  s_h = s_h(2:T,1);

T = T - 1;       % effective sample after first-differencing (164 obs)
t = (1:1:T)';    % integer time index 1, 2, ..., T

%% Construct the cosine basis matrix (Phi)
% Phi is a T x K_end matrix whose k-th column is the k-th cosine basis
% function evaluated at t = 1, ..., T
% Dividing by sqrt(T) normalises so that Phi'*Phi / T --> I as T --> inf.
%
% TAOLS uses only the first K columns of Phi.  K is varied from K_start
% to K_end in steps of K_step to assess sensitivity to bandwidth choice.

K_end   = T - 10;                        % maximum K 
K_start = 8 + 2*trend_type_mcoint;       % minimum K (increases with trend order)
K_step  = 2;                             % increment (only even K values)
nK      = (K_end - K_start)/K_step + 1; % total number of K values evaluated

Phi = zeros(T, K_end);
for k = 1:K_end
    Phi(:,k) = sqrt(2) * sin(pi * t/T * (k - 1/2)) / sqrt(T);
end

%% Project all series onto the cosine basis
% V_X = Phi' * X  is the K_end-vector of cosine coefficients of series X.
% TAOLS then uses only the first K of these coefficients for each K value.

if strcmp(model, 'I')      % Model I (Table 3, Bruns et al. 2020): full-efficacy RF + Berkeley Earth temp
    V_F  = Phi'*F;        %#ok<UNRCH>
    V_S  = Phi'*S;
    V_f  = Phi'*f;
    V_s  = Phi'*s;
    V_ds = Phi'*ds;
    lambda_MLE = 1.709;    % MLE benchmark from Table 3 of Bruns et al. (2020)
    phi_MLE    = 41.482;   % MLE benchmark from Table 3 of Bruns et al. (2020)
    filename = ['Full_RF_BerkeleyT'  '_trendType_' num2str(trend_type_mcoint) '_mc'];
    F_use = F;  S_use = S;  s_use = s;

elseif strcmp(model, 'II') % Model II (Table 3, Bruns et al. 2020): partial-efficacy RF + Berkeley Earth temp
    V_F  = Phi'*F_m;      %#ok<UNRCH>
    V_S  = Phi'*S;
    V_f  = Phi'*f_m;
    V_s  = Phi'*s;
    V_ds = Phi'*ds;
    lambda_MLE = 1.326;    % MLE benchmark from Table 3 of Bruns et al. (2020)
    phi_MLE    = 33.355;   % MLE benchmark from Table 3 of Bruns et al. (2020)
    filename = ['Partial_RF_BerkeleyT' '_trendType_' num2str(trend_type_mcoint) '_mc'];
    F_use = F_m;  S_use = S;  s_use = s;

else                       % Model III (Table 3, Bruns et al. 2020): partial-efficacy RF + HadCRUT temp
    V_F  = Phi'*F_m;       %#ok<UNRCH>
    V_S  = Phi'*S_h;
    V_f  = Phi'*f_m;
    V_s  = Phi'*s_h;
    V_ds = Phi'*ds_h;
    lambda_MLE = 1.567;    % MLE benchmark from Table 3 of Bruns et al. (2020)
    phi_MLE    = 30.695;   % MLE benchmark from Table 3 of Bruns et al. (2020)
    filename = ['Partial_RF_HADCRUT'  '_trendType_' num2str(trend_type_mcoint) '_mc'];
    F_use = F_m;  S_use = S_h;  s_use = s_h;
end

% Trend regressors for the multicointegration model
if trend_type_mcoint == 1
    trend = [ones(T,1), t];          % constant + linear trend
else
    trend = [ones(T,1), t, t.^2];    % constant + linear + quadratic trend
end

V_ell = Phi'*trend;   % projected trend regressors

%% Estimation: multicointegration 
d_x = 1;                                          % the number of x's. 
TA_OLS_mc = zeros(nK, 3*d_x + trend_type_mcoint  + 1);    % for storing point estimates
TA_OLS_se_mc = zeros(nK, 3*d_x + trend_type_mcoint  + 1); % for storing standard errors
TA_OLS_cv_mc = zeros(nK, 1);                              % for storing critical values

for k = K_start:K_step:K_end
   V_F_k = V_F(1:k,:);
   V_S_k = V_S(1:k,:);
   V_s_k = V_s(1:k,:);
   V_ds_k = V_ds(1:k,:);
   V_ell_k = V_ell(1:k,:);

   reg = [V_S_k, V_s_k, V_ds_k, V_ell_k];  % regressors
    
   theta_hat = (reg\V_F_k);                % column vector of coefficient estimates
   
   residual = V_F_k - reg*theta_hat;
    
   omega_ee =  residual'*residual/(k-3*d_x-1-trend_type_mcoint);  
                                      % estimate of the variance of the regression residuals, 
                                      % with degree-of-freedom adjustment for finite-sample correction.
   se = sqrt(diag(omega_ee*inv(reg'*reg))); %#ok<MINV>

   TA_OLS_mc((k-K_start)/K_step+1,:) =   theta_hat' ;

   TA_OLS_se_mc((k-K_start)/K_step+1,:) =   se' ;

   TA_OLS_cv_mc((k-K_start)/K_step+1,1) = tinv(0.975, k-3*d_x-1-trend_type_mcoint );


    Q1 = F_use - theta_hat(1)*S_use;
    Q2 = theta_hat(2)*s_use;
    Q = (Q1 + Q2)/2;
    
    Q = Q(1940-1851+1:end);  % align Q with the ocean heat content series by starting in 1940
    Q = Q - Q(1);  % set the first value of Q to zero for better visual comparison with OHC series
    Q = Q *1.609*0.81;  % 1 Watt for one year/m^2 =  3.1536 x 10^7 Joules/m^2 
                        %  Account for Earth's Total Surface Area: 
                        %  3.1536 x 10^7  * 5.1007 * 10^14 = 1.609 *10^22 Joules
                        % Scaled by 0.81 to get the ocean heat content. 
   figure(1000)
   hold on
   if k == K_start
       plot(year, Q,':', 'Color', [1 0.5 0], 'LineWidth', 1.5, 'DisplayName', 'TAOLS Predictions for K=10:2:154');
   else
        plot(year, Q,':', 'Color', [1 0.5 0], 'LineWidth', 1.5, 'HandleVisibility', 'off');
   end


end


figure(1000)
set(gcf, 'PaperUnits', 'inches');


papersize = get(gcf, 'PaperSize');

width = 5;          % Initialize a variable for width.
height = 3.5;         % Initialize a variable for height.

left = (papersize(1) - width)/2;
bottom = (papersize(2)/2 - height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf, 'PaperPosition', myfiguresize);

saveas(figure(1000), fullfile(results_dir, [filename '_OHC_comparison.fig']));
print(figure(1000), fullfile(results_dir, [filename '_OHC_comparison.eps']), '-depsc', '-vector', '-r300');



%% print some summary statistic for \lambda_hat and ECS_hat

fprintf('\n==============================================\n');
fprintf(' Results for the multicointegration model \n');
fprintf('Summary statistics for lambda\n');
fprintf(' Min: %.4f,\n Max: %.4f, \n Mean: %.4f,\n Median: %.4f\n', ...
    min(TA_OLS_mc(:,1)), max(TA_OLS_mc(:,1)), mean(TA_OLS_mc(:,1)), median(TA_OLS_mc(:,1)));
fprintf('\n');

fprintf('Summary statistics for standard errors of lambda\n');
fprintf(' Min: %.4f,\n Max: %.4f, \n Mean: %.4f,\n Median: %.4f\n', ...
    min(TA_OLS_se_mc(:,1)), max(TA_OLS_se_mc(:,1)), mean(TA_OLS_se_mc(:,1)), median(TA_OLS_se_mc(:,1)));
fprintf('\n');

ECS = DELTA_F2X./TA_OLS_mc(:,1);  % compute the ECS based on the established formula
fprintf('Summary statistics for ECS\n');
fprintf(' Min: %.4f,\n Max: %.4f, \n Mean: %.4f,\n Median: %.4f\n', ...
    min(ECS), max(ECS), mean(ECS), median(ECS));
fprintf('\n');

fprintf('Summary statistics for phi\n');
fprintf(' Min: %.4f,\n Max: %.4f, \n Mean: %.4f,\n Median: %.4f\n', ...
    min(TA_OLS_mc(:,2)), max(TA_OLS_mc(:,2)), mean(TA_OLS_mc(:,2)), median(TA_OLS_mc(:,2)));
fprintf('\n');

fprintf('Summary statistics for percent (100 * HC_ATMOS / phi)\n');
fprintf(' Min: %.4f,\n Max: %.4f, \n Mean: %.4f,\n Median: %.4f\n', ...
    min(100*HC_ATMOS./TA_OLS_mc(:,2)), max(100*HC_ATMOS./TA_OLS_mc(:,2)), ...
    mean(100*HC_ATMOS./TA_OLS_mc(:,2)), median(100*HC_ATMOS./TA_OLS_mc(:,2)));


Ks = (K_start:K_step:K_end)';
TA_OLS_mc_lambda = [TA_OLS_mc(:,1) TA_OLS_mc(:,1) + TA_OLS_cv_mc .* TA_OLS_se_mc(:,1), TA_OLS_mc(:,1) - TA_OLS_cv_mc .* TA_OLS_se_mc(:,1)];

se_adj = DELTA_F2X ./ (TA_OLS_mc(:,1).^2);
% delta-method SE of ECS: |d(ECS)/d(lambda)| * se(lambda) = (DELTA_F2X/lambda^2) * se(lambda)

TA_OLS_mc_ECS = [DELTA_F2X./TA_OLS_mc(:,1), DELTA_F2X./TA_OLS_mc(:,1)+TA_OLS_cv_mc.*TA_OLS_se_mc(:,1).*se_adj, DELTA_F2X./TA_OLS_mc(:,1)-TA_OLS_cv_mc.*TA_OLS_se_mc(:,1).*se_adj];


%% plot the TAOLS estimates and Confidence Intervals of \lambda and ECS
figure(1);
subplot(1,2,1)
hold on;

% Fill the confidence bands and exclude them from the legend
h1 = fill([Ks', fliplr(Ks')], [TA_OLS_mc_lambda(:,2)' , fliplr(TA_OLS_mc_lambda(:,3)')], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
set(get(get(h1, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off'); 
                                                                                 % Columns 2 and 3 consist of the upper and lower limits of 95% CI
                                                                                 % no legend

plot(Ks, TA_OLS_mc_lambda(:,1), 'k-o', 'MarkerIndices', 1:5:length(Ks), 'LineWidth', 1.5);  % point estimate
plot(Ks, TA_OLS_mc_lambda(:,2), 'r:o', 'MarkerIndices', 1:5:length(Ks), 'LineWidth', 1.5);  % Upper CI limit

plot(Ks, ones(size(Ks))*lambda_MLE,'b:', 'LineWidth', 1.5);                                 % MLE from Table 3 of  Bruns et. al (JOE 2020)
legend({'$\hat{\lambda}_{\rm{TAOLS}}$', '95\% CI', '$\hat{\lambda}_{\rm{MLE}}$'}, 'Interpreter', 'latex','FontSize',9, 'Location','southeast');

h2 = plot(Ks, TA_OLS_mc_lambda(:,3), 'r:o', 'MarkerIndices', 1:5:length(Ks), 'LineWidth', 1.5);  % Lower CI limit
set(get(get(h2, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');      

hold off;

xlabel('$K$', 'Interpreter','latex','FontSize',9);
title('(a) TAOLS of $\lambda$ and 95\% CI', 'Interpreter', 'latex', 'FontSize',9);
xticks(10:20:150);
yticks(1:0.1:2.4);
axis([9 154 1.0 2.4]);
grid on;
box on;
set(gca, 'FontSize', 9)


subplot(1,2,2)

hold on 
h1 = fill([Ks', fliplr(Ks')], [TA_OLS_mc_ECS(:,2)' , fliplr(TA_OLS_mc_ECS(:,3)')], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
set(get(get(h1, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off'); 
                                                                                 % Columns 2 and 3 consist of the upper and lower limits of 95% CI
                                                                                 % no legend


plot(Ks, TA_OLS_mc_ECS(:,1), 'k-o', 'MarkerIndices', 1:5:length(Ks), 'LineWidth', 1.5); % point estimate
plot(Ks, TA_OLS_mc_ECS(:,2), 'r:o', 'MarkerIndices', 1:5:length(Ks), 'LineWidth', 1.5); % Upper CI limit

plot(Ks, ones(size(Ks))*DELTA_F2X/lambda_MLE,'b:', 'LineWidth', 1.5);   % Implied ECS from Table 3 of  Bruns et. al (JOE 2020)
legend({'TAOLS of ECS', '95\% CI', 'MLE of ECS'}, 'Interpreter', 'latex', 'Location', 'southeast', 'FontSize', 9);


h2 = plot(Ks, TA_OLS_mc_ECS(:,3), 'r:o', 'MarkerIndices', 1:5:length(Ks), 'LineWidth', 1.5);
set(get(get(h2, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');    % Lower CI limit with no legend

hold off

xlabel('$K$', 'Interpreter','latex','FontSize',9);
title('(b) TAOLS of ECS and 95\% CI', 'Interpreter', 'latex', 'FontSize',9);
xticks(10:20:150);
yticks(0.5:0.25:3.5);
axis([9 154 0.5 3.5]);
grid on;
box on;
set(gca, 'FontSize', 9)

%% Save figure(1)
figure(1)
set(gcf, 'PaperUnits', 'inches');

papersize = get(gcf, 'PaperSize');

width = 5;          % Initialize a variable for width.
height = 4;         % Initialize a variable for height.

left = (papersize(1)- width)/2;
bottom = (papersize(2)/2- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf, 'PaperPosition', myfiguresize);

saveas(gcf, fullfile(results_dir, [filename '_lambda_ECS.fig'])  );
print(gcf, fullfile(results_dir, [filename '_lambda_ECS.eps']), '-depsc', '-vector', '-r300'); % High-quality EPS

save(fullfile(results_dir, [filename '_lambda_ECS.mat']), 'Ks', 'TA_OLS_mc_lambda', 'TA_OLS_mc_ECS');


%% plot the estimates and confidence intervals of the multicointegration parameter \phi

TA_OLS_mc_phi = [TA_OLS_mc(:,2), TA_OLS_mc(:,2)+TA_OLS_cv_mc.*TA_OLS_se_mc(:,2), TA_OLS_mc(:,2)-TA_OLS_cv_mc.* TA_OLS_se_mc(:,2)];
                 % the second column of TA_OLS_mc and TA_OLS_se_mc are the
                 % point estimate and se for \phi

se_adj = HC_ATMOS ./ (TA_OLS_mc(:,2).^2);
% delta-method SE of (HC_ATMOS/phi): |d/d(phi)| * se(phi) = (HC_ATMOS/phi^2) * se(phi)
TA_OLS_mc_percent = 100*[HC_ATMOS./TA_OLS_mc(:,2), ...
                     HC_ATMOS./TA_OLS_mc(:,2)+TA_OLS_cv_mc.*TA_OLS_se_mc(:,2).*se_adj,....
                     HC_ATMOS./TA_OLS_mc(:,2)-TA_OLS_cv_mc.*TA_OLS_se_mc(:,2).*se_adj];


figure(2);
subplot(1,2,1)
hold on;

% Fill the confidence bands and exclude them from the legend
h1 = fill([Ks', fliplr(Ks')], [TA_OLS_mc_phi(:,2)' , fliplr(TA_OLS_mc_phi(:,3)')], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
set(get(get(h1, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off'); 
                                                                                 % Columns 2 and 3 of TA_OLS_mc_phi consist of the upper and lower limits of 95% CI
                                                                                 % no legend

plot(Ks, TA_OLS_mc_phi(:,1), 'k-o', 'MarkerIndices', 1:5:length(Ks), 'LineWidth', 1.5);  % point estimate
plot(Ks, TA_OLS_mc_phi(:,2), 'r:o', 'MarkerIndices', 1:5:length(Ks), 'LineWidth', 1.5);  % Upper CI limit


plot(Ks, ones(size(Ks))*phi_MLE,'b:', 'LineWidth', 1.5);                                 % MLE from Table 3 of  Bruns et. al (JOE 2020)
legend({'$\hat{\phi}_{\rm{TAOLS}}$', '95\% CI',  '$\hat{\phi}_{\rm{MLE}}$'}, 'Interpreter', 'latex','FontSize',9, 'Location','northeast');

h2 = plot(Ks, TA_OLS_mc_phi(:,3), 'r:o', 'MarkerIndices', 1:5:length(Ks), 'LineWidth', 1.5);  % Lower CI limit
set(get(get(h2, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');      

hold off;

xlabel('$K$', 'Interpreter','latex','FontSize',9);
title('(a) TAOLS of $\phi$ and 95\% CI', 'Interpreter', 'latex', 'FontSize',9);
xticks(10:20:150);
axis([9 154 -5 52]);
grid on;
box on;
set(gca, 'FontSize', 9)


subplot(1,2,2)
hold on;

% Fill the 95% CI band (excluded from legend)
h1 = fill([Ks', fliplr(Ks')], [TA_OLS_mc_percent(:,2)', fliplr(TA_OLS_mc_percent(:,3)')], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
set(get(get(h1, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');
                                                                                 % Columns 2 and 3 are the upper and lower 95% CI limits

plot(Ks, TA_OLS_mc_percent(:,1), 'k-o', 'MarkerIndices', 1:5:length(Ks), 'LineWidth', 1.5); % point estimate
plot(Ks, TA_OLS_mc_percent(:,2), 'r:o', 'MarkerIndices', 1:5:length(Ks), 'LineWidth', 1.5); % upper CI limit

plot(Ks, ones(size(Ks)) * 100*HC_ATMOS/phi_MLE, 'b:', 'LineWidth', 1.5);  % MLE implied percent = 100 * HC_ATMOS / phi_MLE

% '\%' escapes the percent sign for the LaTeX interpreter
legend({'\% Heat Content Directed to Warming the Atmosphere', '95\% CI', 'MLE of \%'}, 'Interpreter', 'latex', 'Location', 'southeast', 'FontSize', 9);

h2 = plot(Ks, TA_OLS_mc_percent(:,3), 'r:o', 'MarkerIndices', 1:5:length(Ks), 'LineWidth', 1.5); % lower CI limit
set(get(get(h2, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

hold off;

xlabel('$K$', 'Interpreter', 'latex', 'FontSize', 9);
title('(b) \% of Total Heat Content towards Surface Warming', 'Interpreter', 'latex', 'FontSize', 9);
xticks(10:20:150);
grid on;
box on;
set(gca, 'FontSize', 9)


%% Save figure(2)
figure(2)
set(gcf, 'PaperUnits', 'inches');

papersize = get(gcf, 'PaperSize');

width = 7.5;          % Initialize a variable for width.
height = 3;         % Initialize a variable for height.

left = (papersize(1)- width)/2;
bottom = (papersize(2)- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf, 'PaperPosition', myfiguresize);

saveas(gcf, fullfile(results_dir, [filename '_phi_percent.fig'])  );
print(gcf, fullfile(results_dir, [filename '_phi_percent.eps']), '-depsc', '-vector', '-r300'); % High-quality EPS

save(fullfile(results_dir, [filename '_phi_percent.mat']), 'Ks', 'TA_OLS_mc_phi', 'TA_OLS_mc_percent');


%% Figure 3: TAOLS of trend coefficient(s) with filled 95% CI band
% Mirrors the style of figure 1(a).
% trend_type_mcoint = 1 -> single panel (linear trend coefficient)
% trend_type_mcoint = 2 -> two panels (linear and quadratic trend coefficients)

if trend_type_mcoint == 1

    ind = 3*d_x + 2;   % column index of the linear trend coefficient in TA_OLS_mc
    pt  = TA_OLS_mc(:, ind);
    se  = TA_OLS_se_mc(:, ind);
    cv  = TA_OLS_cv_mc;
    upper = pt + cv .* se;
    lower = pt - cv .* se;

    figure(3)
    hold on;

    % Filled 95% CI band (excluded from legend)
    h1 = fill([Ks', fliplr(Ks')], [upper', fliplr(lower')], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    set(get(get(h1, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

    plot(Ks, pt,    'k-o', 'MarkerIndices', 1:5:length(Ks), 'LineWidth', 1.5);  % point estimate
    plot(Ks, upper, 'r:o', 'MarkerIndices', 1:5:length(Ks), 'LineWidth', 1.5);  % upper CI limit

    h2 = plot(Ks, lower, 'r:o', 'MarkerIndices', 1:5:length(Ks), 'LineWidth', 1.5);  % lower CI limit
    set(get(get(h2, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

    hold off;

    legend({'TAOLS of linear trend coeff', '95\% CI'}, 'Interpreter', 'latex', 'FontSize', 9, 'Location', 'best');
    xlabel('$K$', 'Interpreter', 'latex', 'FontSize', 9);
    ylabel('TAOLS', 'FontSize', 9);
    title('TAOLS of linear trend coefficient and 95\% CI', 'Interpreter', 'latex', 'FontSize', 9);
    grid on;
    box on;
    set(gca, 'FontSize', 9);

end

if trend_type_mcoint == 2

    figure(3)

    % --- Subplot (a): linear trend coefficient ---
    ind   = 3*d_x + 2;
    pt    = TA_OLS_mc(:, ind);
    se    = TA_OLS_se_mc(:, ind);
    cv    = TA_OLS_cv_mc;
    upper = pt + cv .* se;
    lower = pt - cv .* se;

    subplot(1,2,1)
    hold on;

    h1 = fill([Ks', fliplr(Ks')], [upper', fliplr(lower')], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    set(get(get(h1, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

    plot(Ks, pt,    'k-o', 'MarkerIndices', 1:5:length(Ks), 'LineWidth', 1.5);
    plot(Ks, upper, 'r:o', 'MarkerIndices', 1:5:length(Ks), 'LineWidth', 1.5);

    h2 = plot(Ks, lower, 'r:o', 'MarkerIndices', 1:5:length(Ks), 'LineWidth', 1.5);
    set(get(get(h2, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

    hold off;

    legend({'TAOLS of linear trend coeff', '95\% CI'}, 'Interpreter', 'latex', 'FontSize', 9, 'Location', 'best');
    xlabel('$K$', 'Interpreter', 'latex', 'FontSize', 9);
    ylabel('TAOLS', 'FontSize', 9);
    title('(a) Linear trend coefficient', 'Interpreter', 'latex', 'FontSize', 9);
    grid on;
    box on;
    set(gca, 'FontSize', 9);

    % --- Subplot (b): quadratic trend coefficient ---
    ind   = 3*d_x + 3;
    pt    = TA_OLS_mc(:, ind);
    se    = TA_OLS_se_mc(:, ind);
    upper = pt + cv .* se;
    lower = pt - cv .* se;

    subplot(1,2,2)
    hold on;

    h1 = fill([Ks', fliplr(Ks')], [upper', fliplr(lower')], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    set(get(get(h1, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

    plot(Ks, pt,    'k-o', 'MarkerIndices', 1:5:length(Ks), 'LineWidth', 1.5);
    plot(Ks, upper, 'r:o', 'MarkerIndices', 1:5:length(Ks), 'LineWidth', 1.5);

    h2 = plot(Ks, lower, 'r:o', 'MarkerIndices', 1:5:length(Ks), 'LineWidth', 1.5);
    set(get(get(h2, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

    hold off;

    legend({'TAOLS of quadratic trend coeff', '95\% CI'}, 'Interpreter', 'latex', 'FontSize', 9, 'Location', 'best');
    xlabel('$K$', 'Interpreter', 'latex', 'FontSize', 9);
    ylabel('TAOLS', 'FontSize', 9);
    title('(b) Quadratic trend coefficient', 'Interpreter', 'latex', 'FontSize', 9);
    grid on;
    box on;
    set(gca, 'FontSize', 9);

end




%% Estimation: cointegration 

if trend_type_coint == 0
    trend = ones(T,1);
else
    trend = [ones(T,1), t];
end

V_ell = Phi'*trend;


if strcmp(model, 'I')      % Model I: full-efficacy RF (f) and Berkeley Earth temperature (s)
    filename = ['Full_RF_BerkeleyT', '_trendType_' num2str(trend_type_coint) '_coint']; %#ok<UNRCH>

elseif strcmp(model, 'II') % Model II: partial-efficacy RF (f_m) and Berkeley Earth temperature (s)
    filename = ['Partial_RF_BerkeleyT' '_trendType_' num2str(trend_type_coint) '_coint']; %#ok<UNRCH>

else                       % Model III: partial-efficacy RF (f_m) and HadCRUT temperature (s_h)
    filename = ['Partial_RF_HADCRUT' '_trendType_' num2str(trend_type_coint) '_coint'];
end


d_x = 1;  % the number of x's. 
TA_OLS = zeros(nK, 2*d_x + trend_type_coint + 1);    % point estimates
TA_OLS_se = zeros(nK, 2*d_x + trend_type_coint + 1); % standard errors
TA_OLS_cv = zeros(nK, 1);   % critical values

for k = K_start:K_step:K_end
   
   V_f_k = V_f(1:k,:); 
   V_s_k = V_s(1:k,:);
   V_ds_k = V_ds(1:k,:);
   V_ell_k = V_ell(1:k,:);

   reg = [V_s_k, V_ds_k, V_ell_k]; % regressors
    
   theta_hat = (reg\V_f_k); % column vector of coefficient estimates
   
   residual = V_f_k - reg*theta_hat;
    
   omega_ee =  residual'*residual/(k-2*d_x-1-trend_type_coint);  
                                      % estimate of the variance of the regression residuals, 
                                      % with degree-of-freedom adjustment for finite-sample correction.
   se = sqrt(diag(omega_ee*inv(reg'*reg))); %#ok<MINV>

   TA_OLS((k-K_start)/K_step+1,:) =   theta_hat' ;

   TA_OLS_se((k-K_start)/K_step+1,:) =   se' ;

   TA_OLS_cv((k-K_start)/K_step+1,1) = tinv(0.975, k-2*d_x-1-trend_type_coint);

end


%% print some summary statistic for \lambda-hat and ECS_hat
fprintf('\n==============================================\n');
fprintf('Results for the regular cointegration model \n');
fprintf('Summary statistics for lambda\n');
fprintf(' Min: %.4f,\n Max: %.4f, \n Mean: %.4f,\n Median: %.4f\n', ...
    min(TA_OLS(:,1)), max(TA_OLS(:,1)), mean(TA_OLS(:,1)), median(TA_OLS(:,1)));
fprintf('\n');

fprintf('Summary statistics for standard errors of lambda\n');
fprintf(' Min: %.4f,\n Max: %.4f, \n Mean: %.4f,\n Median: %.4f\n', ...
    min(TA_OLS_se(:,1)), max(TA_OLS_se(:,1)), mean(TA_OLS_se(:,1)), median(TA_OLS_se(:,1)));
fprintf('\n');

ECS = DELTA_F2X./TA_OLS(:,1);  % compute the ECS based on the established formula
fprintf('Summary statistics for ECS\n');
fprintf(' Min: %.4f,\n Max: %.4f, \n Mean: %.4f,\n Median: %.4f\n', ...
    min(ECS), max(ECS), mean(ECS), median(ECS));



TA_OLS_lambda = [TA_OLS(:,1) TA_OLS(:,1)+TA_OLS_cv.* TA_OLS_se(:,1), TA_OLS(:,1)-TA_OLS_cv.* TA_OLS_se(:,1)];

se_adj = DELTA_F2X ./ (TA_OLS_lambda(:,1).^2);
% delta-method SE of ECS: (DELTA_F2X/lambda^2) * se(lambda)
TA_OLS_ECS = [DELTA_F2X./TA_OLS(:,1), DELTA_F2X./TA_OLS(:,1)+TA_OLS_cv.*TA_OLS_se(:,1).*se_adj, DELTA_F2X./TA_OLS(:,1)-TA_OLS_cv.*TA_OLS_se(:,1).*se_adj];



figure(4)
subplot(1,2,1)
hold on;

% Fill the confidence bands and exclude them from the legend
h1 = fill([Ks', fliplr(Ks')], [TA_OLS_lambda(:,2)' , fliplr(TA_OLS_lambda(:,3)')], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
set(get(get(h1, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off'); 

% Columns 2 and 3 consist of the upper and lower limits of 95% CI
% no legend

plot(Ks, TA_OLS_lambda(:,1), 'k-o', 'MarkerIndices', 1:5:length(Ks), 'LineWidth', 1.5);  % point estimate
plot(Ks, TA_OLS_lambda(:,2), 'r:o', 'MarkerIndices', 1:5:length(Ks), 'LineWidth', 1.5);  % Upper CI limit

legend({'$\hat{\lambda}_{\rm{TAOLS}}$', '95\% CI'}, 'Interpreter', 'latex','FontSize',9, 'Location','southeast');

h2 = plot(Ks, TA_OLS_lambda(:,3), 'r:o', 'MarkerIndices', 1:5:length(Ks), 'LineWidth', 1.5);  % Lower CI limit
set(get(get(h2, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');      

hold off;

xlabel('$K$', 'Interpreter','latex','FontSize',9);
title('(a) TAOLS of $\lambda$ and 95\% CI', 'Interpreter', 'latex', 'FontSize',9);
xticks(10:20:150);
yticks(1:0.2:3.3);
axis([9 154 1.5 3.3]);
grid on;
box on;
set(gca, 'FontSize', 9)


subplot(1,2,2)

hold on 
h1 = fill([Ks', fliplr(Ks')], [TA_OLS_ECS(:,2)' , fliplr(TA_OLS_ECS(:,3)')], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
set(get(get(h1, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off'); 
% Columns 2 and 3 consist of the upper and lower limits of 95% CI
% no legend


plot(Ks, TA_OLS_ECS(:,1), 'k-o', 'MarkerIndices', 1:5:length(Ks), 'LineWidth', 1.5); % point estimate
plot(Ks, TA_OLS_ECS(:,2), 'r:o', 'MarkerIndices', 1:5:length(Ks), 'LineWidth', 1.5); % Upper CI limit


legend({'TAOLS of ECS', '95\% CI'}, 'Interpreter', 'latex', 'Location', 'southeast', 'FontSize', 9);


h2 = plot(Ks, TA_OLS_ECS(:,3), 'r:o', 'MarkerIndices', 1:5:length(Ks), 'LineWidth', 1.5);
set(get(get(h2, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');    % Lower CI limit with no legend

hold off

xlabel('$K$', 'Interpreter','latex','FontSize',9);
title('(b) TAOLS of ECS and 95\% CI', 'Interpreter', 'latex', 'FontSize',9);
xticks(10:20:150);
yticks(0.25:0.25:3);
axis([9 154 0.25 3]);
grid on;
box on;
set(gca, 'FontSize', 9)

%% Save figure(4)
figure(4)
set(gcf, 'PaperUnits', 'inches');

papersize = get(gcf, 'PaperSize');

width = 7.5;          % Initialize a variable for width.
height = 3;         % Initialize a variable for height.

left = (papersize(1)- width)/2;
bottom = (papersize(2)/2- height)/2;
myfiguresize = [left, bottom, width, height];
set(gcf, 'PaperPosition', myfiguresize);

saveas(gcf, fullfile(results_dir, [filename '_lambda_ECS.fig'])  );
print(gcf, fullfile(results_dir, [filename '_lambda_ECS.eps']), '-depsc', '-vector', '-r300'); % High-quality EPS

save(fullfile(results_dir, [filename '_lambda_ECS.mat']), 'Ks', 'TA_OLS_lambda', 'TA_OLS_ECS');

if trend_type_coint == 1

    ind   = 2*d_x + 2;   % column index of the linear trend coefficient in TA_OLS
    pt    = TA_OLS(:, ind);
    se    = TA_OLS_se(:, ind);
    cv    = TA_OLS_cv;
    upper = pt + cv .* se;
    lower = pt - cv .* se;

    figure(5)
    hold on;

    % Filled 95% CI band (excluded from legend)
    h1 = fill([Ks', fliplr(Ks')], [upper', fliplr(lower')], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    set(get(get(h1, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

    plot(Ks, pt,    'k-o', 'MarkerIndices', 1:5:length(Ks), 'LineWidth', 1.5);  % point estimate
    plot(Ks, upper, 'r:o', 'MarkerIndices', 1:5:length(Ks), 'LineWidth', 1.5);  % upper CI limit

    h2 = plot(Ks, lower, 'r:o', 'MarkerIndices', 1:5:length(Ks), 'LineWidth', 1.5);  % lower CI limit
    set(get(get(h2, 'Annotation'), 'LegendInformation'), 'IconDisplayStyle', 'off');

    hold off;

    legend({'TAOLS of linear trend coeff', '95\% CI'}, 'Interpreter', 'latex', 'FontSize', 9, 'Location', 'best');
    xlabel('$K$', 'Interpreter', 'latex', 'FontSize', 9);
    ylabel('TAOLS', 'FontSize', 9);
    title('TAOLS of linear trend coefficient and 95\% CI', 'Interpreter', 'latex', 'FontSize', 9);
    grid on;
    box on;
    set(gca, 'FontSize', 9);

end


diary off