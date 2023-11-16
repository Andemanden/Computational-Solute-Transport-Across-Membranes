%__Gruppe 3-Simulering-3.Sem__
clear
close all
clc
%% Inital Parameters
domain_length = 0.1;   % Domain length (meters)
domain_steps = 1000;    % Number of spatial domain_steps points
dx = domain_length / (domain_steps - 1); % Position discretization
 
time_length = 500;      % Time length (seconds)
time_steps = 500;       % Number of time steps
dt = time_length / time_steps; % Temporal discretization

% DEFINED VARIABLES
D = 0; % Diffusivity coefficient H2PO4- in water [m^2 s^-1]
feed_conc = 0.1; % Constant solute concentration at the first cell 0.1 molar [H2PO4-]
TMP = 35; %TMP: Transmembrane Pressure [bar]

area = 0.001; % Area of the membrane surface [m^2]
kw = 5.7311*10^(-7); % Initial water permeability m^3 m^-2 bar^-1 s^-1 
my = 0.8903*10^-9; % Water viscosity [Bar∙s]
Rm = 1/(my*kw); % Rejection of water at the membrane (σ) [m^-1]
InitP = 0.263253+0.0011; % Initial percipitation
alpha = 95000000000; % Specific resistance of fouling [m mol m^3]
PC = 0.5; %Percipitate Advection Coefficient

sig_m = 0.1; % Rejection of ions 

% PHYSICAL CONSTANTS
R= 8.31415*10^-2; % Gasconstant [L^3 Bar mol^-1 K^-1]
T= 273.15+25; %Temperature [K]
 
% Anonymous functions
%Lv = @(time) kw*exp(-kb*time); %Water permeability dependent on Fouling

Mp = @(conc) 0.0011*exp(36.328*conc) + 0.263253; % concentration of udfæld in respect to phosphate increase.
Rf = @(conc) alpha*(Mp(conc)/area); % Rate of fouling
Lv = @(conc) 1/(my*(Rm+Rf(conc))); % Water permeability dependent on fouling


% Create a grid for space and time
x = linspace(0, domain_length, domain_steps);
t = linspace(0, time_length, time_steps);

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
        LastC = C(domain_steps, j-1);

        C(1, :) = feed_conc; % Set the leftmost boundary to 0.1

        if j == 2
            Ptot = InitP; % The rate of percipitation WI
        else
            Ptot = Mp(LastC) + (Mp(LastC) - InitP)*Jv*(dt/dx)*PC;
        end

        Jv = (kw*(TMP-(1*R*T*(LastC))));  % Volume flux = Jv ,  in terms of osmotic pressure (TMP, R, T, delta_C) and Lv. [m/s]

        % Calculate the second derivative in x direction
        if i == 1 %_Left Cell__
            d2Cdx2 = 0;
            dCdt_advection = 0;

        elseif i == domain_steps %__Membrane Cell___
            % No right neighbor at the last cell and MEMBRANE
            d2Cdx2 = (D * (C(domain_steps - 1, j-1) - C(domain_steps, j-1)) / dx^2) + C(i,j-1)*(sig_m)*Jv/(dx);
        else    %__Normal Cells__
            % Calculate the second derivative normally
            d2Cdx2 = D * (C(i + 1, j-1) - 2 * C(i, j-1) + C(i - 1, j-1)) / dx^2;
            % Apply advection term
            dCdt_advection = -Jv * (C(i, j-1) - C(i-1, j-1)) / dx;
        end
        % Apply the diffusion-advection-(electromigration) equation
        C(i, j) = C(i, j-1) + dt * (d2Cdx2 + dCdt_advection);
        Jv_values(j) = Jv; % Plot values
        P_values(j) = Ptot;
        AS(j) = Jv * dt/dx; % Stability plot values
    end
end
%% CONSERVATION ERROR

Systemdiff = [0, diff(sum(C,1))]; % Diffrence in systemsums

inflow = Jv_values*feed_conc*(dt/dx); % Inflow array
outflow = C(domain_steps, :).* Jv_values*(1-sig_m)*(dt/dx); % Outflow array 

ERROR = Systemdiff - inflow + outflow; % The mass conservation error



%% 2D Plots

% Plot Mass Conservation Error values over time

figure;
plot(t, ERROR);
xlabel('Time (seconds)');
ylabel('Error ');
title('Mass Conservation Error Over Time');
grid on;


% Plot Percipitate (DS) values over time

figure;
plot(t, P_values);
xlabel('Time (seconds)');
ylabel('percipitate ');
title('percipitate  Over Time');
grid on;

% Plot Advection Stability (DS) values over time

figure;
plot(t, AS);
xlabel('Time (seconds)');
ylabel('Advection Stability');
title('Advection Stability Over Time');
grid on;

% Plot Jv values over time
figure;
plot(t, Jv_values);
xlabel('Time (seconds)');
ylabel('Jv (Velocity)');
title('Jv (Velocity) Over Time');
grid on;


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

xlabel('Position (meters)');
ylabel('Concentration');
title('Concentration Over Time at Different Instances');

% Add a legend for clarity
legend(arrayfun(@(f) ['t=', num2str(f)], time_fraction, 'UniformOutput', false));

grid on;
hold off;


%% 3D Plot

% Create a 3D surface plot to visualize concentration over time and position
[T, X] = meshgrid(t, x);
figure;
h = surf(X, T, C); % Transpose removed here
xlabel('Position (meter)');
ylabel('Tid (sekunder)');
zlabel('Koncentration (M)');
title('Koncentration Over Tid og Position');

% Set axis limits to start at the origin
xlim([0, domain_length]);
ylim([0, time_length]);
zlim([0, 1.5]); % Assuming max(C(:)) is the maximum concentration in your data

set(h,'LineStyle','none')
colormap(jet)
clim([0 1.5])



