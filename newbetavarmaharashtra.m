clc;
clear all;

% Time span and initial conditions
tspan = [0 10];
xnot = [60439692, 1000000, 2500, 469, 207, 3000, 37510300, 49638138, 7802, 1300, 207];

% Solve the three models with different sigma & delta
[T1, X1] = ode15s(@(t,x) kar(t,x,0.12067), tspan, xnot); 
[T2, X2] = ode15s(@(t,x) kar(t,x,0.07067), tspan, xnot); 
[T3, X3] = ode15s(@(t,x) kar(t,x,0.02067), tspan, xnot); 

% Plot results
figure;
plot(T1, X1(:,5), 'r', 'LineWidth', 2); hold on;
plot(T2, X2(:,5), 'b', 'LineWidth', 2);
plot(T3, X3(:,5), 'g', 'LineWidth', 2);
xlabel('Time (days)');
ylabel('Infected Symptomatic (I_s)');
legend('\beta = 0.12067', '\beta = 0.12067', '\beta = 0.2067');
set(gca, 'FontWeight', 'bold'); % Set font weight of axes
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dx = kar(t, x, beta)
 Lambda_h = 908800; Lambda_v = 7008000000;
    mu_h = 0.0142; mu_v = 0.5; mu1 = 0.001; mu2 = 0.073;
    b = 0.6; c = 0.3; d = 0.45; e = 0.8;
    sigma = 0.34074; phi = 0.81425; rho = 0.52869;
    kappa=0.12416; epsilon=0.92278; gamma=0.07624;
    alpha=0.27295; eta=0.85791; tau=0.44533; delta=0.05442;

    dx = zeros(11,1);
    N_h = sum(x(1:6));

    dx(1) = Lambda_h - (b * beta * x(10) * x(1))/N_h - (sigma * x(1) * x(7)) - mu_h * x(1) + phi * x(2);
    dx(2) = (b * beta * x(10) * x(1))/N_h + sigma * x(1) * x(7) - (phi + mu_h) * x(2);
    dx(3) = ((1 - eta) * (c * beta * x(10) * x(2)))/N_h + (b * beta * x(10) * x(1))/N_h - (kappa * epsilon + mu_h) * x(3);
    dx(4) = kappa * epsilon * x(3) - (gamma + mu_h + mu1) * x(4);
    dx(5) = (1 - kappa) * epsilon * x(3) - (alpha + mu_h + mu2) * x(5);
    dx(6) = gamma * x(4) + alpha * x(5) - mu_h * x(6);
    dx(7) = tau * x(5) - delta * x(7);
    dx(8) = Lambda_v - (d * beta * x(8) * x(4))/N_h - (e * beta * x(8) * x(5))/N_h - mu_v * x(8);
    dx(9) = (d * beta * x(8) * x(4))/N_h + (e * beta * x(8) * x(5))/N_h - (rho + mu_v) * x(9);
    dx(10) = rho * x(9) - mu_v * x(10);
    dx(11) = (1 - kappa) * epsilon * x(3);
end
