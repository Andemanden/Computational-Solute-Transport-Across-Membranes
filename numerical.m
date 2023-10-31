%__Gruppe 3-Simulering-3.Sem__
clear
close all
clc
%% Inital Parameters 
domain_length = 0.5;   % Domain length (meters)
domain_steps = 100;    % Number of spatial domain_steps points
dx = domain_length / (domain_steps - 1);

time_length = 0.5;      % Time length (seconds)
time_steps = 1000;       % Number of time steps
dt = time_length / time_steps;

% DEFINED VARIABLES
D = 0.001;               % Diffusion coefficient
feed_conc = 1; % Constant solute concentration at the first cell
TMP=10; %TMP: Transmembrane Pressure [bar]
kb=0.0002; % Fouling Constant
kw=0.00021; % Initial water permiability

% PHYSICAL CONSTANTS
R=8.31415; % Gasconstant []
T=273.15+20; %Temperature [K]
F= 96.485;%Faraday[C/mol
% Anonymous functions
rejection_rate = @(time) 0.1; % σ_0+σ_f
Lv= @(time) kw*exp(-kb*time); %Water Permiability for all cells

% Create a grid for space and time
x = linspace(0, domain_length, domain_steps);
t = linspace(0, time_length, time_steps);

% Initialize the concentration array AND Initial condition (1D)
C = ones(domain_steps, time_steps);
C(domain_steps/2,1)=2;


%% Time-stepping loop
for j = 2:time_steps
    for i = 1:domain_steps
        Jv = (Lv(j)*(TMP-(1*R*T*(C(domain_steps-1, j-1)-C(domain_steps, j-1)))));  % Volume flux = Velocity ,  in terms of osmotic pressure (TMP, R, T, delta_C) and Lv. [m/s]
        % Calculate the second derivative in x direction
        if i == 1 %__First Cell__
             C(1,:) = feed_conc; % Fill the first cell with inflow concentration
            % No left neighbor at the first cell
            d2Cdx2 = D * (C(2, j-1) -  C(1, j-1)) / dx^2;
             % Apply advection and inflow at the first cell
            dCdt_advection = -Jv * C(i, j-1) / dx;
        elseif i == domain_steps %__Membrane Cell___
            % No right neighbor at the last cell and MEMBRANE
            d2Cdx2 = (D * (C(domain_steps - 1, j-1) - C(domain_steps, j-1)) / dx^2) + C(i,j-1)*(rejection_rate(j))*Jv/(dx);
        else    %__Normal Cells__
            % Calculate the second derivative normally
            d2Cdx2 = D * (C(i + 1, j-1) - 2 * C(i, j-1) + C(i - 1, j-1)) / dx^2;
            % Apply advection term
            dCdt_advection = -Jv * (C(i, j-1) - C(i-1, j-1)) / dx;
        end
        % Apply the diffusion-advection-(electromigration) equation
        C(i, j) = C(i, j-1) + dt * (d2Cdx2 + dCdt_advection);
        Jv_values(j) = Jv;
        AS(j) = Jv * dt/dx;
    end
end

%Diffusive Stability
DS = D*dx/dt;
fprintf(' Diffusivity Stability = %f', DS);
if DS>0.5 
   fprintf('ERROR: Stabilitetsfejl');
else
    fprintf('\n Stable Diffusion Model !!');
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




