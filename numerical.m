%__Gruppe 3-Simulering-3.Sem__
clear
close all
clc
%% Inital Parameters 
domain_length = 0.1;   % Domain length (meters)
domain_steps = 100;    % Number of spatial domain_steps points
dx = domain_length / (domain_steps - 1);

time_length = 500;      % Time length (seconds)
time_steps = 5000;       % Number of time steps
dt = time_length / time_steps;

% DEFINED VARIABLES
D = 0.846*10^-9; % Diffusivity coefficient H2PO4 in water
feed_conc = 0.1; % Constant solute concentration at the first cell 0.1 molar [H2PO4]
TMP= 5.5; %TMP: Transmembrane Pressure [bar]
kb= 0.001; % Fouling Constant
kw= 3.044001*10^-4; % Initial water permiability L m^-2 bar^-1 s^-1 
rejection_rate = 0.1; % Rejection of the membrane (σ) [procentage]


% PHYSICAL CONSTANTS
R= 0.0831415; % Gasconstant [L Bar mol^-1 K^-1]
T= 273.15+20; %Temperature [K]
F= 96.485; %Faraday[C/mol]

% Anonymous functions
Lv = @(time) kw*exp(-kb*time); %Water Permiability for all cells

% Percipitation array
%P = zeros(domain_length,);

% Create a grid for space and time
x = linspace(0, domain_length, domain_steps);
t = linspace(0, time_length, time_steps);

% Initialize the concentration array AND Initial condition (1D)
C = zeros(domain_steps, time_steps)+0.1; % 0.1 molar [H2PO4]

C(domain_steps/2,1)=0.5;

%Diffusive Stability
DS = D*dx/dt;
fprintf('\n Diffusivity Stability = %f', DS);
if DS>0.5 
   fprintf(2,'\n ERROR: Stabilitetsfejl');
else
    fprintf('\n Stable Diffusion Model !!');
end


%% Time-stepping loop
for j = 2:time_steps
    for i = 1:domain_steps
        Jv = (Lv(j)*(TMP-(1*R*T*(C(domain_steps, j-1)))));  % Volume flux = Jv ,  in terms of osmotic pressure (TMP, R, T, delta_C) and Lv. [m/s]
        C(1, :) = 0.1; % Set the leftmost boundary to 0.1
        % Calculate the second derivative in x direction
        if i == 1 %__First Cell__
            % No left neighbor at the first cell
            d2Cdx2 = D * (C(2, j-1) -  C(1, j-1)) / dx^2;
             % Apply advection and inflow at the first cell
            dCdt_advection = -Jv * C(i, j-1) / dx;
        elseif i == domain_steps %__Membrane Cell___
            % No right neighbor at the last cell and MEMBRANE
            d2Cdx2 = (D * (C(domain_steps - 1, j-1) - C(domain_steps, j-1)) / dx^2) + C(i,j-1)*(rejection_rate)*Jv/(dx);
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




