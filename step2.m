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
sig_m = 0.1; % Rejection of ions 

% PHYSICAL CONSTANTS
R= 8.31415*10^-2; % Gasconstant [L^3 Bar mol^-1 K^-1]
T= 273.15+25; %Temperature [K]

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
        Jv = (kw*(TMP-(1*R*T*(LastC))));  % Volume flux = Jv ,  in terms of osmotic pressure (TMP, R, T, delta_C) and Lv. [m/s]

        % Calculate the second derivative in x direction
        if i == 1 %_Left Cell__
            d2Cdx2 = 0;
            dCdt_advection = 0;

        elseif i == domain_steps %__Membrane Cell___
            % No right neighbor at the last cell and MEMBRANE
            d2Cdx2 = (D * (C(domain_steps - 1, j-1) - C(domain_steps, j-1))) / dx^2;
            dCdt_advection = -Jv * (C(domain_steps, j-1)*(1-sig_m) - C(domain_steps - 1, j-1)) / dx;

        else    %__Normal Cells__
            % Calculate the second derivative normally
            d2Cdx2 = D * (C(i + 1, j-1) - 2 * C(i, j-1) + C(i - 1, j-1)) / dx^2;
            % Apply advection term
            dCdt_advection = -Jv * (C(i, j-1) - C(i-1, j-1)) / dx;
        end
        % Apply the diffusion-advection-(electromigration) equation
        C(i, j) = C(i, j-1) + dt * (d2Cdx2 + dCdt_advection);
        Jv_values2(j) = Jv; % Plot values
        AS(j) = Jv * dt/dx; % Stability plot values
    end
end
%% conservation and error

volprc = dx*area*1000;  % Volume per cell in Liters

Cw = [0,C(domain_steps, :)]; % Making the matrix line up properly and have the right size
Cw(end) = [];

Systemdiff = [0, diff(sum(C*volprc,1))]; % Difference in system concentration per dt

inflow = feed_conc*(Jv_values2)*(dt/dx)*volprc; % In - mol
outflow = Cw.* Jv_values2*(dt/dx)*volprc*(1-sig_m); % Out - mol
infoutdiff = (inflow - outflow);
ERROR = ((Systemdiff - infoutdiff).*Systemdiff.^(-1))*100; % The mass conservation error


%% 2D Plots

% Plot Mass Conservation Error values over time

figure;
plot(t, ERROR);
xlabel('Time (seconds)');
ylabel('Error ');
title('Mass Conservation Error Over Time');
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
plot(t, Jv_values2);
xlabel('Time (seconds)');
ylabel('Jv (Velocity)');
title('Jv (Velocity) Over Time');
grid on;


%% 3D Plot

% 3D surface plot - concentration over time and position
[T, X] = meshgrid(t, x);
figure;
h = surf(X, T, (C/feed_conc)-1);
xlabel('Position [m]');
ylabel('Tid [s]');
zlabel('CF - 1');
title('Centreret Koncentrations Faktor Over Tid og Position');
xlim([0, domain_length]);
ylim([0, time_length]);
zlim([0, 0.12]);
set(h,'LineStyle','none')
colormap(jet)
clim([0, 0.12])

%%
save('step2jv', 'Jv_values2')

