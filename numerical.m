%__Gruppe 3-Simulering-3.Sem__
clear
close all
clc
%% Inital Parameters
domain_length = 0.1;   % Domain length (meters)
domain_steps = 1000;    % Number of spatial domain_steps points
dx = domain_length / (domain_steps - 1); % Position discretization

time_length = 50;      % Time length (seconds)
time_steps = 500;       % Number of time steps
dt = time_length / time_steps; % Temporal discretization

% DEFINED VARIABLES
D = 1.464*10^-9; % Diffusivity coefficient H2PO4- in water [m^2 s^-1]
feed_conc = 0.1; % Constant solute concentration at the first cell 0.1 molar [H2PO4-]
TMP = 5; %TMP: Transmembrane Pressure [bar]

area = 0.001; % Area of the membrane surface [m^2]
kw = 5.7311*10^(-7); % Initial water permeability m^3 m^-2 bar^-1 s^-1 
my = 0.891*10^-9; % Water viscosity [Bar∙s]
Rm = 1/(my*kw); % Rejection of water at the membrane (σ) [procentage]

alpha = 9.7915*10^12; % Specific resistance of fouling [m mol m^3]

sig_m = 0.1; % Rejection of ions 

% PHYSICAL CONSTANTS
R= 0.0831415; % Gasconstant [L Bar mol^-1 K^-1]
T= 273.15+25; %Temperature [K]
F= 96.485; %Faraday[C/mol]

% Anonymous functions
%Lv = @(time) kw*exp(-kb*time); %Water permeability dependent on Fouling

percip_rate = @(conc) 0.3822436157/(1+2.6885457349*exp(-13.1198866*conc)); % The rate of percipitation
Rf = @(conc) alpha*(percip_rate(conc)/area); % Rate of fouling
Lv = @(conc) 1/(my*(Rm+Rf(conc))); % Water permeability dependent on fouling

% Percipitation array - WORK IN PROGRESS
%P = zeros(domain_length,);

% Create a grid for space and time
x = linspace(0, domain_length, domain_steps);
t = linspace(0, time_length, time_steps);

% Initialize the concentration array AND Initial condition (1D)
C = zeros(domain_steps, time_steps)+0.1; % 0.1 molar [H2PO4] at pH 2.9
C(domain_steps/2) = 0.2;

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
        LastC=C(i,j-1);
        Jv = (Lv(LastC)*(TMP-(1*R*T*(C(domain_steps, j-1)))));  % Volume flux = Jv ,  in terms of osmotic pressure (TMP, R, T, delta_C) and Lv. [m/s]
        %Jv = (Lv*(TMP-(1*R*T*(C(domain_steps, j-1)))));
        C(1, :) = feed_conc; % Set the leftmost boundary to 0.1
        % Calculate the second derivative in x direction
        if i == 1 %__First Cell__
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
        AS(j) = Jv * dt/dx; % Stability plot values
    end
end

%% 2D Plots

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
xlabel('Position (meters)');
ylabel('Time (seconds)');
zlabel('Concentration');
title('Concentration Over Time and Position');

% Set axis limits to start at the origin
xlim([0, domain_length]);
ylim([0, time_length]);
zlim([0, max(C(:))]); % Assuming max(C(:)) is the maximum concentration in your data

set(h,'LineStyle','none')
%colormap(jet)




