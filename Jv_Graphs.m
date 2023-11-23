%__Gruppe 3-Simulering-3.Sem__
clear
close all
clc
%%

x = linspace(0, 500, 500);

load('step1jv.mat', 'Jv_values1');
load('step2jv.mat', 'Jv_values2');
load('step3jv.mat', 'Jv_values3');
load('step4jv.mat', 'Jv_values4');

%%

% Create a figure for the first plot
figure;
hold on;

plot(x, Jv_values1, 'LineWidth', 1.5);
plot(x, Jv_values2, 'LineWidth', 1.5);
plot(x, Jv_values3, 'LineWidth', 1.5);
plot(x, Jv_values4, 'LineWidth', 1.5);

xlabel('Time (seconds)');
ylabel('Jv (Velocity)');
title('Jv (Velocity) Over Time');
grid on;

legend('Step 1', 'Step 2', 'Step 3', 'Step 4');
%% Double Y-axis Jv graph
x = linspace(0, 500, 500);

load('step2jv.mat', 'Jv_values2');
load('step3jv.mat', 'Jv_values3');
load('step4jv.mat', 'Jv_values4');

% Remove the first value from each series
Jv_values2 = Jv_values2(2:end);
Jv_values3 = Jv_values3(2:end);
Jv_values4 = Jv_values4(2:end);

figure;

% Plot Step 2 on the left y-axis
yyaxis left;
plot(x(2:end), Jv_values2, 'LineWidth', 1.5, 'LineStyle', '-'); 
ylabel('J_{v} [m^{3} m^{-2} s^{-1}]');

hold on;

% Plot Step 3 on the left y-axis
plot(x(2:end), Jv_values3, 'LineWidth', 1.5, 'LineStyle', '--');
ylim([1.825*10^-5, 1.865*10^-5]);           % y-limit for left
ylabel('J_{v} [m^{3} m^{-2} s^{-1}]');

% Plot Step 4 on the right y-axis
yyaxis right;
ylim([1.2*10^-5, 1.24*10^-5]);               % y-limit for right

plot(x(2:end), Jv_values4, 'LineWidth', 1.5);
ylabel('J_{v} [m^{3} m^{-2} s^{-1}]');

xlabel('Tid [s]');
title('J_{v} Over Tid');
grid on;

legend('Step 2', 'Step 3', 'Step 4');


