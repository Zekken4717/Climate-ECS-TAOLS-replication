%% extract_ranges_and_combine.m
% Run from Empirical_application_matlab\ to:
%   1. Extract TAOLS numerical ranges for Models I, II, III
%   2. Generate combined lambda/ECS and phi/percent figures
%   3. Export combined figures as PNG to Paper_June_2026\
%
% Mat files are in Results_Model_X_cointTrendType_0\ subfolders.
% Prerequisite: Empirical_Analysis_Climate_trendType_1_in_mc.m must have
%   been run with model = 'I', 'II', and 'III' (producing all six .mat files).

close all; clear; clc;

cd(fileparts(mfilename('fullpath')));  % set working directory to the folder containing this script

script_dir = fileparts(mfilename('fullpath'));           % .../Empirical_application_matlab
repo_root  = fileparts(script_dir);                     % .../A_Multicointegration_Global_Warming
matI   = fullfile(script_dir, 'Results_Model_I_cointTrendType_0');
matII  = fullfile(script_dir, 'Results_Model_II_cointTrendType_0');
matIII = fullfile(script_dir, 'Results_Model_III_cointTrendType_0');
out    = fullfile(repo_root, 'Paper_June_2026');

%% -----------------------------------------------------------------------
%  1.  Load lambda/ECS mat files
% -----------------------------------------------------------------------
d = load(fullfile(matI,   'Full_RF_BerkeleyT_trendType_1_mc_lambda_ECS.mat'));
Ks      = d.Ks;
lambda1 = d.TA_OLS_mc_lambda;   % N x 3: [estimate, lower CI, upper CI]
ECS1    = d.TA_OLS_mc_ECS;

d = load(fullfile(matII,  'Partial_RF_BerkeleyT_trendType_1_mc_lambda_ECS.mat'));
lambda2 = d.TA_OLS_mc_lambda;
ECS2    = d.TA_OLS_mc_ECS;

d = load(fullfile(matIII, 'Partial_RF_HADCRUT_trendType_1_mc_lambda_ECS.mat'));
lambda3 = d.TA_OLS_mc_lambda;
ECS3    = d.TA_OLS_mc_ECS;

%% -----------------------------------------------------------------------
%  2.  Load phi/percent mat files
% -----------------------------------------------------------------------
d = load(fullfile(matI,   'Full_RF_BerkeleyT_trendType_1_mc_phi_percent.mat'));
phi1 = d.TA_OLS_mc_phi;      % N x 3
pct1 = d.TA_OLS_mc_percent;

d = load(fullfile(matII,  'Partial_RF_BerkeleyT_trendType_1_mc_phi_percent.mat'));
phi2 = d.TA_OLS_mc_phi;
pct2 = d.TA_OLS_mc_percent;

d = load(fullfile(matIII, 'Partial_RF_HADCRUT_trendType_1_mc_phi_percent.mat'));
phi3 = d.TA_OLS_mc_phi;
pct3 = d.TA_OLS_mc_percent;

%% -----------------------------------------------------------------------
%  3.  Print numerical ranges
% -----------------------------------------------------------------------
fprintf('\n=== TAOLS Numerical Ranges ===\n');

fprintf('\n--- Model I (Full RF, Berkeley Earth) ---\n');
fprintf('lambda: [%.3f, %.3f],  mean = %.3f\n', min(lambda1(:,1)), max(lambda1(:,1)), mean(lambda1(:,1)));
fprintf('ECS:    [%.2f, %.2f] C,  mean = %.2f C\n', min(ECS1(:,1)), max(ECS1(:,1)), mean(ECS1(:,1)));
fprintf('phi:    [%.2f, %.2f] W-yr/m2\n', min(phi1(:,1)), max(phi1(:,1)));
fprintf('pct:    [%.2f, %.2f] %%\n', min(pct1(:,1)), max(pct1(:,1)));

fprintf('\n--- Model II (Partial RF, Berkeley Earth) ---\n');
fprintf('lambda: [%.3f, %.3f],  mean = %.3f\n', min(lambda2(:,1)), max(lambda2(:,1)), mean(lambda2(:,1)));
fprintf('ECS:    [%.2f, %.2f] C,  mean = %.2f C\n', min(ECS2(:,1)), max(ECS2(:,1)), mean(ECS2(:,1)));
fprintf('phi:    [%.2f, %.2f] W-yr/m2\n', min(phi2(:,1)), max(phi2(:,1)));
fprintf('pct:    [%.2f, %.2f] %%\n', min(pct2(:,1)), max(pct2(:,1)));

fprintf('\n--- Model III (Partial RF, HadCRUT) ---\n');
fprintf('lambda: [%.3f, %.3f],  mean = %.3f\n', min(lambda3(:,1)), max(lambda3(:,1)), mean(lambda3(:,1)));
fprintf('ECS:    [%.2f, %.2f] C,  mean = %.2f C\n', min(ECS3(:,1)), max(ECS3(:,1)), mean(ECS3(:,1)));
fprintf('phi:    [%.2f, %.2f] W-yr/m2\n', min(phi3(:,1)), max(phi3(:,1)));
fprintf('pct:    [%.2f, %.2f] %%\n', min(pct3(:,1)), max(pct3(:,1)));

fprintf('\n--- Cross-model summary ---\n');
all_lambda = [lambda1(:,1); lambda2(:,1); lambda3(:,1)];
all_ECS    = [ECS1(:,1);    ECS2(:,1);    ECS3(:,1)];
all_phi    = [phi1(:,1);    phi2(:,1);    phi3(:,1)];
all_pct    = [pct1(:,1);    pct2(:,1);    pct3(:,1)];
fprintf('lambda (all models): [%.3f, %.3f]\n', min(all_lambda), max(all_lambda));
fprintf('ECS    (all models): [%.2f, %.2f] C\n', min(all_ECS), max(all_ECS));
fprintf('phi    (all models): [%.2f, %.2f] W-yr/m2\n', min(all_phi), max(all_phi));
fprintf('pct    (all models): [%.2f, %.2f] %%\n', min(all_pct), max(all_pct));

%% -----------------------------------------------------------------------
%  4.  Combined lambda/ECS figure
% -----------------------------------------------------------------------
fig1 = figure(1);
set(fig1, 'Units','inches','Position',[1 1 9 3.5]);

% --- Panel (a): lambda ---
subplot(1,2,1); hold on;

h1a = fill([Ks', fliplr(Ks')], [lambda1(:,2)', fliplr(lambda1(:,3)')], ...
    [0.8 0.2 0.2], 'FaceAlpha',0.18,'EdgeColor','none');
set(get(get(h1a,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

h2a = fill([Ks', fliplr(Ks')], [lambda2(:,2)', fliplr(lambda2(:,3)')], ...
    [0.2 0.4 0.8], 'FaceAlpha',0.18,'EdgeColor','none');
set(get(get(h2a,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

h3a = fill([Ks', fliplr(Ks')], [lambda3(:,2)', fliplr(lambda3(:,3)')], ...
    [0.2 0.7 0.2], 'FaceAlpha',0.18,'EdgeColor','none');
set(get(get(h3a,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

p1 = plot(Ks, lambda1(:,1), 'r-',  'LineWidth',1.5);
p2 = plot(Ks, lambda2(:,1), 'b-s', 'LineWidth',1.5, 'MarkerIndices',1:5:length(Ks), 'MarkerSize',4);
p3 = plot(Ks, lambda3(:,1), 'g-o', 'LineWidth',1.5, 'MarkerIndices',1:5:length(Ks), 'MarkerSize',4);

hold off;
legend([p1 p2 p3], {'Model I','Model II','Model III'}, 'Location','southeast','FontSize',8);
xlabel('$K$','Interpreter','latex');
ylabel('$\hat{\lambda}$','Interpreter','latex');
title('(a) TAOLS of $\lambda$ and 95\% CI','Interpreter','latex','FontSize',10);
xticks(10:20:150); grid on; box on;

% --- Panel (b): ECS ---
subplot(1,2,2); hold on;

h1b = fill([Ks', fliplr(Ks')], [ECS1(:,2)', fliplr(ECS1(:,3)')], ...
    [0.8 0.2 0.2], 'FaceAlpha',0.18,'EdgeColor','none');
set(get(get(h1b,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

h2b = fill([Ks', fliplr(Ks')], [ECS2(:,2)', fliplr(ECS2(:,3)')], ...
    [0.2 0.4 0.8], 'FaceAlpha',0.18,'EdgeColor','none');
set(get(get(h2b,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

h3b = fill([Ks', fliplr(Ks')], [ECS3(:,2)', fliplr(ECS3(:,3)')], ...
    [0.2 0.7 0.2], 'FaceAlpha',0.18,'EdgeColor','none');
set(get(get(h3b,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

q1 = plot(Ks, ECS1(:,1), 'r-',  'LineWidth',1.5);
q2 = plot(Ks, ECS2(:,1), 'b-s', 'LineWidth',1.5, 'MarkerIndices',1:5:length(Ks), 'MarkerSize',4);
q3 = plot(Ks, ECS3(:,1), 'g-o', 'LineWidth',1.5, 'MarkerIndices',1:5:length(Ks), 'MarkerSize',4);

% MLE reference line from BCS (2020)
yline(2.80, 'k:', 'LineWidth',1.2, 'HandleVisibility','off');

hold off;
legend([q1 q2 q3], {'Model I','Model II','Model III'}, 'Location','northeast','FontSize',8);
xlabel('$K$','Interpreter','latex');
ylabel('ECS ($^\circ$C)','Interpreter','latex');
title('(b) TAOLS of ECS and 95\% CI','Interpreter','latex','FontSize',10);
xticks(10:20:150); grid on; box on;

% Export
set(fig1,'PaperUnits','inches','PaperPosition',[0 0 9 3.5]);
outfile1 = fullfile(out, 'combined_mc_lambda_ECS.png');
print(fig1, outfile1, '-dpng', '-r300');
fprintf('\nSaved: %s\n', outfile1);

%% -----------------------------------------------------------------------
%  5.  Combined phi/percent figure
% -----------------------------------------------------------------------
fig2 = figure(2);
set(fig2,'Units','inches','Position',[1 1 9 3.5]);

% --- Panel (a): phi ---
subplot(1,2,1); hold on;

h1c = fill([Ks', fliplr(Ks')], [phi1(:,2)', fliplr(phi1(:,3)')], ...
    [0.8 0.2 0.2], 'FaceAlpha',0.18,'EdgeColor','none');
set(get(get(h1c,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

h2c = fill([Ks', fliplr(Ks')], [phi2(:,2)', fliplr(phi2(:,3)')], ...
    [0.2 0.4 0.8], 'FaceAlpha',0.18,'EdgeColor','none');
set(get(get(h2c,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

h3c = fill([Ks', fliplr(Ks')], [phi3(:,2)', fliplr(phi3(:,3)')], ...
    [0.2 0.7 0.2], 'FaceAlpha',0.18,'EdgeColor','none');
set(get(get(h3c,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

r1 = plot(Ks, phi1(:,1), 'r-',  'LineWidth',1.5);
r2 = plot(Ks, phi2(:,1), 'b-s', 'LineWidth',1.5, 'MarkerIndices',1:5:length(Ks), 'MarkerSize',4);
r3 = plot(Ks, phi3(:,1), 'g-o', 'LineWidth',1.5, 'MarkerIndices',1:5:length(Ks), 'MarkerSize',4);

hold off;
legend([r1 r2 r3], {'Model I','Model II','Model III'}, 'Location','northeast','FontSize',8);
xlabel('$K$','Interpreter','latex');
ylabel('$\hat{\phi}$ (W-yr/m$^2$)','Interpreter','latex');
title('(a) TAOLS of $\phi$ and 95\% CI','Interpreter','latex','FontSize',10);
xticks(10:20:150); grid on; box on;

% --- Panel (b): percent ---
subplot(1,2,2); hold on;

h1d = fill([Ks', fliplr(Ks')], [pct1(:,2)', fliplr(pct1(:,3)')], ...
    [0.8 0.2 0.2], 'FaceAlpha',0.18,'EdgeColor','none');
set(get(get(h1d,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

h2d = fill([Ks', fliplr(Ks')], [pct2(:,2)', fliplr(pct2(:,3)')], ...
    [0.2 0.4 0.8], 'FaceAlpha',0.18,'EdgeColor','none');
set(get(get(h2d,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

h3d = fill([Ks', fliplr(Ks')], [pct3(:,2)', fliplr(pct3(:,3)')], ...
    [0.2 0.7 0.2], 'FaceAlpha',0.18,'EdgeColor','none');
set(get(get(h3d,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

s1 = plot(Ks, pct1(:,1), 'r-',  'LineWidth',1.5);
s2 = plot(Ks, pct2(:,1), 'b-s', 'LineWidth',1.5, 'MarkerIndices',1:5:length(Ks), 'MarkerSize',4);
s3 = plot(Ks, pct3(:,1), 'g-o', 'LineWidth',1.5, 'MarkerIndices',1:5:length(Ks), 'MarkerSize',4);

hold off;
legend([s1 s2 s3], {'Model I','Model II','Model III'}, 'Location','northeast','FontSize',8);
xlabel('$K$','Interpreter','latex');
ylabel('\% of Total Heat Content','Interpreter','latex');
title('(b) \% of Total Heat Content toward Surface Warming','Interpreter','latex','FontSize',10);
xticks(10:20:150); grid on; box on;

% Export
set(fig2,'PaperUnits','inches','PaperPosition',[0 0 9 3.5]);
outfile2 = fullfile(out, 'combined_mc_phi_percent.png');
print(fig2, outfile2, '-dpng', '-r300');
fprintf('Saved: %s\n', outfile2);

fprintf('\nDone.\n');
