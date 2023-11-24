%__Gruppe 3-Simulering-3.Sem__
clear
close all
clc
%% Inital Parameters
Lx = 0.1;   % Domain length [m]
domain_steps = 1000;    % Number of spatial domain_steps points
dx = Lx / (domain_steps - 1); % Position discretization
 
Lt = 500;      % Time length (seconds)
time_steps = 500;       % Number of time steps
dt = Lt / time_steps; % Temporal discretization

% DEFINED VARIABLES
D = 1.464*10^-9; % Diffusivity coefficient H2PO4- in water [m^2 s^-1]
Cf = 0.1; % Feed concentration [mol L^-1]
DeltaP = 35; %TMP: Transmembrane Pressure [bar]

Am = 0.0308; % Area of the membrane surface [m^2]
k0 = 5.7311*10^(-7); % Initial water permeability [m^3 m^-2 bar^-1 s^-1] 
mu = 0.8903*10^-9; % Water viscosity [Bar∙s]
Rm = 1/(mu*k0); % Rejection of water at the membrane (σ) [m^-1]
alpha = 1*10^14; % Specific resistance of fouling [m L mol^-1]
JuScalar = 0.5; %Percipitate Advection Coefficient

sig_i = 0.1; % Rejection of ions 

% PHYSICAL CONSTANTS
R = 8.31415*10^-2; % Gasconstant [L Bar mol^-1 K^-1]
T = 273.15+25; %Temperature [K]
 
% Anonymous functions

Udf = @(conc) 0.0011*exp(36.328*conc) + 0.263253; % concentration of udfæld in respect to phosphate increase. [mol L^-1]
Rf = @(conc) alpha*(Udf(conc)/Am); % Rate of fouling [m^-1]
k = @(conc) 1/(mu*(Rm+Rf(conc))); % Water permeability dependent on fouling [m^3 m^-2 bar^-1 s^-1]


% Create a grid for space and time
x = linspace(0, Lx, domain_steps);
t = linspace(0, Lt, time_steps);

% Initialize the concentration array AND Initial condition (1D)
C = zeros(domain_steps, time_steps)+0.1; % 0.1 molar [H2PO4] at pH 2.9


%Diffusive Stability

DS = D*dt/dx^2;
fprintf('\n Diffusivity Stability = %f', DS);
if DS>0.5 
   fprintf(2,'\n ERROR: Stabilitetsfejl');
else
    fprintf('\n Stable Diffusion Model !!');
end


%% Time-stepping loop
for j = 2:time_steps
    for i = 1:domain_steps
        Ciw = C(domain_steps, j-1);

        C(1, :) = Cf; % Set the leftmost boundary to 0.1M

        if j == 2
            Cbw = Udf(Cf); % The rate of percipitation
        else
            Cbw = Udf(Ciw) + Jkonv*Udf(Cf)*(dt/dx)*JuScalar; % Advection of Percipitate
        end

        Jkonv = (k(Ciw)*(DeltaP-(1*R*T*(Ciw))));  % Volume flux = Jkonv ,  in terms of osmotic pressure (DeltaP, R, T, C_(i_w)) and k. [m/s]

        % Calculate the second derivative in x direction
        if i == 1 %__First Cell__
            Jdiff = 0;
            Jadv = 0;

        elseif i == domain_steps %__Membrane Cell___
            % No right neighbor at the last cell and MEMBRANE
            Jdiff = (D * (C(domain_steps - 1, j-1) - C(domain_steps, j-1))) / dx^2;
            Jadv = -Jkonv * (C(domain_steps, j-1)*(1-sig_i) - C(domain_steps - 1, j-1)) / dx;

        else    %__Normal Cells__
            % Calculate the second derivative normally
            Jdiff = D * (C(i + 1, j-1) - 2 * C(i, j-1) + C(i - 1, j-1)) / dx^2;
            % Apply advection term
            Jadv = -Jkonv * (C(i, j-1) - C(i-1, j-1)) / dx;
        end
        % Apply the diffusion-advection-(electromigration) equation
        C(i, j) = C(i, j-1) + dt * (Jdiff + Jadv);
        Jkonv_values(j) = Jkonv; % Plot values
        Cbw_values(j) = Cbw; 
        AS(j) = Jkonv * dt/dx; % Stability plot values
    end
end
%% conservation and error

volprc = dx*Am*1000;  % Volume per cell in Liters

Ciw = [0,C(domain_steps, :)]; % Making the matrix line up properly and have the right size
Ciw(end) = [];

Systemdiff = [0, diff(sum(C*volprc,1))]; % Difference in system concentration per dt

inflow = Cf*(Jkonv_values)*(dt/dx)*volprc; % In - mol
outflow = Ciw.* Jkonv_values*(dt/dx)*volprc*(1-sig_i); % Out - mol
infoutdiff = (inflow - outflow);
ERROR = ((Systemdiff - infoutdiff).*Systemdiff.^(-1))*100; % The mass conservation error


%% 2D Plots

% Plot Mass Conservation Error values over time

figure;
plot(t, ERROR);
xlabel('Tid [s]');
ylabel('Afvigelse [%]');
title('Afvigelse Over Tid');
grid on;


% Plot Percipitate (DS) values over time

figure;
plot(t, Cbw_values);
xlabel('Tid [s]');
ylabel('Bundfald [mol L^{-1}]');
title('Bundfald Over Tid');
grid on;
ylim([0.3, 0.35]);

% Plot Advection Stability (DS) values over time

figure;
plot(t(2:end), AS(2:end));
xlabel('Tid [s]');
ylabel('Advection Stabilitet');
title('Advection Stabilitet Over Tid');
grid on;

% Plot Jkonv values over time
figure;
plot(t(2:end), Jkonv_values(2:end));
xlabel('Tid [s]');
ylabel('J_{konv} [m s^{-1}]');
title('J_{konv} [m s^{-1}] Over Tid');
grid on;
xlim([0, 500]);
ax = gca;
ax.YAxis.Exponent = -6;

% Define fractions of time steps you want to visualize
time_fraction = [0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];  % For example, 0.1 corresponds to 10% of time steps

% Calculate the corresponding time indices
time_instances = round(time_fraction * time_steps);

% Create a figure for the first plot
figure;
hold on;

% Plot concentration at specific time instances
for i = 1:length(time_instances)
    time_index = time_instances(i);
    plot(x, C(:, time_index));
end

xlabel('Position [m]');
ylabel('Koncentration [mol L^{-1}]');
title('Koncentration over Position ved Forskellige Tidsfraktioner');
xlim([0.099, Lx]);
ylim([0.1, 0.12]);

% Add a legend for clarity
legend(arrayfun(@(f) ['t=', num2str(f)], time_fraction, 'UniformOutput', false));

grid on;
hold off;


%% 3D Plot

% Create a 3D surface plot to visualize concentration over time and position
[T, X] = meshgrid(t, x);
figure;
h = surf(X, T, C); % Transpose removed here
xlabel('Position [m]');
ylabel('Tid [s]');
zlabel('Koncentration [mol L^{-1}]');
title('Koncentration Over Tid og Position');

% Set axis limits to start at the origin
xlim([0, Lx]);
ylim([0, Lt]);
zlim([0.1, 0.12]); % Assuming max(C(:)) is the maximum concentration in your data

set(h,'LineStyle','none')
colormap(jet)
clim([0.1, 0.11])

%%
load('step1jv.mat', 'Jv_values1');
load('step2jv.mat', 'Jv_values2');
load('step3jv.mat', 'Jv_values3');
load('step4jv.mat', 'Jv_values4');

% Remove the first value from each series
Jv_values1 = Jv_values1(2:end);
Jv_values2 = Jv_values2(2:end);
Jv_values3 = Jv_values3(2:end);
Jv_values4 = Jv_values4(2:end);

%%

% Create a figure for the first plot
figure;
hold on;
x = linspace(0, 500, 500);

plot(x(2:end), Jv_values1, 'LineWidth', 1.5);
plot(x(2:end), Jv_values2, 'LineWidth', 1.5);
plot(x(2:end), Jv_values3, 'LineWidth', 1.5);
plot(x(2:end), Jv_values4, 'LineWidth', 1.5);

xlabel('Time (seconds)');
ylabel('Jv (Velocity)');
title('Jv (Velocity) Over Time');
grid on;

legend('Step 1', 'Step 2', 'Step 3', 'Step 4');
%% Double Y-axis Jv graph
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
