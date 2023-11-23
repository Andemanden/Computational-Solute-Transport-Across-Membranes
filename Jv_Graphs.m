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

%%new graph
figure;
plot(x, Jv_values3, 'LineWidth', 1.5);
xlabel('Time (seconds)');
ylabel('Jv (Velocity)');
title('Jv (Velocity) Over Time 3');
grid on;
%ylim([, ])
ax = gca;
ax.YAxis.Exponent = -5;
legend('Step 3');

figure;
plot(x, Jv_values4, 'LineWidth', 1.5);
xlabel('Time (seconds)');
ylabel('Jv (Velocity)');
title('Jv (Velocity) Over Time 4');
grid on;
%ylim([, ])
ax = gca;
ax.YAxis.Exponent = -5;
legend('Step 4');

