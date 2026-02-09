clear all;
clc;

syms Lambda_v Lambda_h mu_v mu_h mu1 mu2 b c d e beta rho eta alpha kappa epsilon gamma


R = (rho*(b*beta-c*beta*(eta-1)))/(mu_v*(mu_v+rho)) * ((Lambda_v*d*beta*epsilon*kappa*mu_h)/(Lambda_h*mu_v*(epsilon+mu_h)*(gamma+mu_h+mu1)) - (Lambda_v*e*beta*epsilon*mu_h*(kappa-1))/(Lambda_h*mu_v*(epsilon+mu_h)*(alpha+mu_h+mu2)));
R0 = sqrt(R);



disp('Basic reproduction number R0:');
disp(R0);

% Compute Sensitivity Indices
a1 = diff(R0, Lambda_v) * (Lambda_v / R0);
a2 = diff(R0, Lambda_h) * (Lambda_h / R0);
a3 = diff(R0, mu_v) * (mu_v / R0);  % Corrected line
a4 = diff(R0, mu_h) * (mu_h / R0);
a5 = diff(R0, mu1) * (mu1 / R0);
a6 = diff(R0, mu2) * (mu2 / R0);
a7 = diff(R0, b) * (b / R0);
a8 = diff(R0, c) * (c / R0);
a9 = diff(R0, d) * (d / R0);
a10 = diff(R0, e) * (e / R0);
a11 = diff(R0, beta) * (beta / R0);
a12 = diff(R0, rho) * (rho / R0);
a13 = diff(R0, eta) * (eta / R0);
a14 = diff(R0, alpha) * (alpha / R0);
a15 = diff(R0, kappa) * (kappa / R0);
a16 = diff(R0, epsilon) * (epsilon / R0);
a17 = diff(R0, gamma) * (gamma / R0)


% Substitute numerical values
params = [Lambda_v Lambda_h mu_v mu_h mu1 mu2 b c d e beta rho eta alpha kappa epsilon gamma];


values = [7008000000, 908800, 0.5, 0.0142, 0.1, 0.073, 0.6, 0.3, 0.45, 0.8, 0.18719, 0.52869, 0.82674, 0.60235, 0.02521, 0.54630, 0.091015];

b_values = double(subs([a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17], params, values));

% Plot Sensitivity Indices
figure;
bar(b_values, 'g');

% Define parameter labels
parameter_labels = {'\Lambda_v', '\Lambda_h', '\mu_v', '\mu_h', '\mu_1', '\mu_2', 'b', 'c', 'd', 'e', '\beta', '\rho', ...
                    '\eta', '\alpha', '\kappa', ...
                    '\epsilon', '\gamma'};

xticks(1:length(parameter_labels));
xticklabels(parameter_labels);
xlabel('\bf Parameters', 'FontSize', 12);
ylabel('\bf Sensitivity Index', 'FontSize', 12);
title('\bf Karnataka', 'FontSize', 16);
ylim([-1.5, 1]);

set(gca,'fontsize',20,'fontweight','bold','linewidth',1)
