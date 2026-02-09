clc; clear; close all;




Lambda_h = 908800;
Lambda_v = 7008000000;
mu_h = 0.0142;
mu_v = 0.5 ;
mu1 = 0.001;
mu2 = 0.073;
b = 0.6;
c = 0.3;
d = 0.45;
e = 0.8;

sigma = 0.34074;
phi = 0.81425;
rho = 0.52869;




%% India baseline parameter (median)
kappa_med  = 0.30449;
epsilon_med = 0.50368;
alpha_med   = 0.37665; 
eta_med     = 0.88633;
tau_med     = 0.50168;
delta_med   = 0.02077;

%% 95% CI bounds (from Table 2)
kappa_CI   = [0.282 0.327];
epsilon_CI = [0.466 0.541];
alpha_CI   = [0.348 0.405];
eta_CI     = [0.851 0.922];
tau_CI     = [0.464 0.540];
delta_CI   = [0.019 0.023];

%% Parameter grids
gamma = linspace(0,0.1773,150);   % recovery rate (rows)
beta  = linspace(0,0.1827,150);   % transmission rate (cols)
[Beta,Gamma] = meshgrid(beta,gamma);

%% Monte-Carlo simulation settings
Nsample = 200;      % for quick test; set to 500 or higher when satisfied
R0_all  = NaN([size(Beta), Nsample]);

rng(0); % for reproducibility

%% Monte Carlo loop
for n = 1:Nsample
    % sample parameters uniformly within 95% CI
    kappa   = kappa_CI(1)   + rand*(kappa_CI(2)-kappa_CI(1));
    epsilon = epsilon_CI(1) + rand*(epsilon_CI(2)-epsilon_CI(1));
    alpha   = alpha_CI(1)   + rand*(alpha_CI(2)-alpha_CI(1));
    eta     = eta_CI(1)     + rand*(eta_CI(2)-eta_CI(1));
    tau     = tau_CI(1)     + rand*(tau_CI(2)-tau_CI(1));
    delta   = delta_CI(1)   + rand*(delta_CI(2)-delta_CI(1));
    
    % compute R for this random parameter set
    % break the expression into parts for clarity & numerical safety
    A1 = (rho.*(b.*Beta - c.*Beta.*(eta-1))) ./ (mu_v.*(mu_v+rho));
    part1 = (Lambda_v.*d.*Beta.*epsilon.*kappa.*mu_h) ./ (Lambda_h.*mu_v.*(epsilon+mu_h).*(Gamma+mu_h+mu1));
    part2 = (Lambda_v.*e.*Beta.*epsilon.*mu_h.*(kappa-1)) ./ (Lambda_h.*mu_v.*(epsilon+mu_h).*(alpha+mu_h+mu2));
    R = A1 .* (part1 - part2);
    
    % Replace Inf with NaN (guard against divide by zero)
    R(~isfinite(R)) = NaN;
    
    % Clamp negative R to 0 before sqrt (R should be nonnegative for R0)
    R_clamped = R;
    R_clamped(R_clamped < 0) = 0;
    
    % If R_clamped is NaN everywhere for this n, skip storing to avoid polluting stats
    if all(isnan(R_clamped(:)))
        warning('Sample %d produced all NaNs - skipping', n);
        continue;
    end
    
    R0_all(:,:,n) = sqrt(R_clamped);
    
    % optional progress indicator every 50 samples
    if mod(n,50) == 0
        fprintf('Completed sample %d / %d\n', n, Nsample);
    end
end

%% Remove fully-NaN slices (if some samples were skipped)
valid_slices = squeeze(~all(all(isnan(R0_all),1),2));
R0_all = R0_all(:,:,valid_slices);

if isempty(R0_all)
    error('No valid Monte Carlo samples remained. Check your parameter ranges and expressions.');
end

%% Compute statistics across all samples
% Use median ignoring NaNs. For older MATLAB versions use nanmedian
try
    R0_med = median(R0_all,3,'omitnan');
    R0_p05 = prctile(R0_all,5,3);
    R0_p95 = prctile(R0_all,95,3);
catch
    % fallback for older MATLAB:
    R0_med = nanmedian(R0_all,3);
    R0_p05 = prctile(R0_all,5,3); % prctile handles NaN slices as well
    R0_p95 = prctile(R0_all,95,3);
end

%% Binary map (median)
regionMap = zeros(size(R0_med));
regionMap(R0_med<1)=1;
regionMap(R0_med>=1)=2;

%% Plot median contour with uncertainty shading
% figure('Color','w','Position',[100 100 800 600]);


figure;




% plot median safe/epidemic regions
hMap = pcolor(Beta,Gamma,regionMap);
shading flat;
colormap([0 1 0; 1 0 0]); % light green/red
caxis([1 2]);
hold on;

% contour for median threshold
[c1,h1] = contour(Beta,Gamma,R0_med,[1 1],'k','LineWidth',2);
clabel(c1,h1,'FontSize',11,'FontWeight','bold','Color','k','LabelSpacing',300);

% contour for uncertainty band (5th and 95th percentile thresholds)
[c2,h2] = contour(Beta,Gamma,R0_p05,[1 1],'--','Color',[0.3 0.3 0.3],'LineWidth',2.5);
[c3,h3] = contour(Beta,Gamma,R0_p95,[1 1],'--','Color',[0.3 0.3 0.3],'LineWidth',2.5);

xlabel('\bf \beta','FontSize',13);
ylabel('\bf \gamma','FontSize',13);

title('\bf Gujarat: Safe (Green) vs Epidemic (Red) Regions','FontSize',14,'FontWeight','bold');

% ---- Add invisible handles for legend ----
hSafe = patch(NaN, NaN, 'g', 'FaceAlpha',0.6, 'EdgeColor','none');
hEpi  = patch(NaN, NaN, 'r', 'FaceAlpha',0.6, 'EdgeColor','none');
hMed  = plot(NaN, NaN, 'k-', 'LineWidth', 2);
hCI   = plot(NaN, NaN, '--', 'Color',[0.3 0.3 0.3], 'LineWidth', 2);

legend([hSafe hEpi hMed hCI], ...
       {'Safe region (R_0 < 1)', ...
        'Epidemic region (R_0 > 1)', ...
        'Threshold (R_0 = 1)', ...
        '95% CI'}, ...
       'Location','northeast', ...
       'FontSize',12, ...
       'FontWeight','bold', ...
       'Box','on');

set(gca,'fontsize',14,'fontweight','bold','linewidth',1)

% set(gca,'FontSize',12,'FontWeight','bold','LineWidth',1.2);
grid on; box on;

hold off;

