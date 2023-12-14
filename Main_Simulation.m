%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  GITHUB - ALL RIGHTS RESERVED   %%%
%%%                                 %%%
%%% PROPERTY OF AALBORG UNIVERSITY  %%%
%%%         CREATED BY:             %%%
%%%  GROUP 3 - 3RD SEMESTER- 2023   %%%
%%%    @ANALREXIA    @ANDEMANDEN    %%%
%%%    @BONINATOR1   @PETERPASTA11  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear
close all
clc
%% Inital Parameters
Lx = 0.1;                     % Domain length [m]
domain_steps = 1000;          % Number of spatial domain steps
dx = Lx / (domain_steps - 1); % Position discretization [m]
 
Lt = 500;              % Time length [s]
time_steps = 500;      % Number of time steps
dt = Lt / time_steps;  % Temporal discretization [s]

% DEFINED VARIABLES
D = 1.464*10^-9;      % Diffusivity coefficient H2PO4- in water [m^2 s^-1]
Cf = 0.1;             % Ionconcentration of the feed [mol L^-1]
DeltaP = 35;          % Change in outer pressure [bar]
Am = 0.0308;          % Surfacearea of the membrane [m^2]
k0 = 5.7311*10^-7;    % Initial water permeability [m^3 m^-2 bar^-1 s^-1] 
mu = 0.8903*10^-9;    % Viscosity of water [bar∙s]
Rm = 1/(mu*k0);       % Membraneresitance of water [m^-1]
alpha = 1*10^14;      % Specific fouling resistance [m · kg-1] or [m L mol^-1]
JuScalar = 0.5;       % Percipitate Advection Coefficient
sig_i = 0.1;          % Rejection of ions 

% PHYSICAL CONSTANTS

R = 8.31415*10^-2;  % Gasconstant [L Bar mol^-1 K^-1]
T = 273.15+25;      % Temperature [K]
 
% Anonymous functions

Udf = @(conc) 0.0011*exp(36.328*conc) + 0.263253; % Concentration of percipitate in respect to Ciw. [mol L^-1]
Rf = @(conc) alpha*(conc/Am);                % Foulingresistance [m^-1]
k = @(conc) 1/(mu*(Rm+Rf(conc)));                 % Water permeability dependent on mu,Rm and Rf [m^3 m^-2 bar^-1 s^-1]


% Space and Time Grid
x = linspace(0, Lx, domain_steps);
t = linspace(0, Lt, time_steps);

% Initial concentration array AND Initial condition
C = zeros(domain_steps, time_steps)+0.1; % 0.1 molar [H2PO4] at pH 2.9

% Diffusive Stability

DS = D*dt/dx^2;
fprintf('\n Diffusivity Stability = %f', DS);
if DS>0.5 
   fprintf(2,'\n ERROR: Stabilitetsfejl');
else
    fprintf('\n Stable Diffusion Model !!');
end


%% Time-stepping loop
for j = 2:time_steps
    if j == 2
            Cudf = Udf(Cf); % Initial Percipitation
            Cop = 0;
            Cbw = Udf(Cf);
            Jb = 0;
        else
            Cudf = Udf(Ciw);
            Jb = Jkonv*Udf(Cf)*(dt/dx)*JuScalar;
            Cop = Cop + Jb; % Total percipitate after advection [mol L^-1]
            Cbw = Cudf + Cop;
    end
    
    
    for i = 1:domain_steps
        Ciw = C(domain_steps, j-1); % Concentration of ions at the wall [mol L^-1]

        C(1, :) = Cf; % Set the leftmost boundary to Cf [mol L^-1]

        Jkonv = (k(Cbw)*(DeltaP-(1*R*T*(Ciw))));  % Volume flux = Jkonv , in terms of osmotic pressure (DeltaP, R, T, C_(i_w)) and k. [m/s]

        if i == 1 % First Cell (no left neighbor)
            Jdiff = 0;       % Diffusive ionflux
            Jadv = 0;        % Advective flux

        elseif i == domain_steps % Membrane wall cell (no right neighbor)
            Jdiff = (D * (C(domain_steps - 1, j-1) - C(domain_steps, j-1))) / dx^2;           % Diffusive ionflux [mol · m-2 · s-1]
            Jadv = -Jkonv * (C(domain_steps, j-1)*(1-sig_i) - C(domain_steps - 1, j-1)) / dx; % Advective flux [mol · m-2 · s-1]

        else    % Bulk Cells
            % Calculate the second derivative normally
            Jdiff = D * (C(i + 1, j-1) - 2 * C(i, j-1) + C(i - 1, j-1)) / dx^2; % Diffusive ionflux at the membrane wall [mol · m-2 · s-1]
            Jadv = -Jkonv * (C(i, j-1) - C(i-1, j-1)) / dx;                     % Advective flux at the membrane wall [mol · m-2 · s-1]
        end
        % Apply the diffusion-advection equation
        C(i, j) = C(i, j-1) + dt * (Jdiff + Jadv); % Apply the diffusion-advection equation [mol L^-1]
        Jkonv_values(j) = Jkonv;                   % Storing convective flux values in respect to time
        Cbw_values(j) = Cbw;                       % Storing precipitate values in respect to time
        AS(j) = Jkonv * dt/dx;                     % Storing Stability plot values
        % Main Stability
        MainS=(1-2*DS-AS(j));                      % Calculation of specific main stability to time
        MainS_values(j)=MainS;                     % Storing main stability values

    end
end
%% conservation and error

volprc = dx*Am*1000;  % Volume per cell in Liters

Ciw = [0,C(domain_steps, :)]; % Making the matrix line up properly and have the right size
Ciw(end) = [];

Systemdiff = [0, diff(sum(C*volprc,1))]; % Difference in system concentration per dt

inflow = Cf*(Jkonv_values)*(dt/dx)*volprc;                 % Input [mol]
outflow = Ciw.* Jkonv_values*(dt/dx)*volprc*(1-sig_i);     % Output [mol]
infoutdiff = (inflow - outflow);                           % Difference in input and output
ERROR = ((Systemdiff - infoutdiff).*Systemdiff.^(-1))*100; % Mass conservation error


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
%ylim([0.3, 0.35]);

% Plot Advection Stability (DS) values over time

figure;
plot(t(2:end), AS(2:end));
xlabel('Tid [s]');
ylabel('Advection Stabilitet');
title('Advection Stabilitet Over Tid');
grid on;

% Plot Main Stability (MainS) values over time

figure;
plot(t(2:end), MainS_values(2:end));
xlabel('Tid [s]');
ylabel('Stability value');
title('Main Stabilitet Over Tid');
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

% Figure - Concentration over time sideplot
figure;
hold on;

for i = 1:length(time_instances)   % Concentration at specific time instances
    time_index = time_instances(i);
    plot(x, C(:, time_index)/Cf - 1);
end

xlabel('Position [m]');
ylabel('CF - 1');
title('Centreret CF over Position ved Forskellige Tidspunkter');
xlim([0.099, Lx]);
ylim([0, 0.12]);

% Add a legend for clarity
legend(arrayfun(@(f) ['t=', num2str(f), 's'], time_instances, 'UniformOutput', false));

grid on;
hold off;


%% 3D Plot

% 3D surface plot - concentration over time and position
[T, X] = meshgrid(t, x);
figure;
h = surf(X, T, (C/Cf)-1);
xlabel('Position [m]');
ylabel('Tid [s]');
zlabel('CF - 1');
title('Centreret Koncentrations Faktor Over Tid og Position');
xlim([0, Lx]);
ylim([0, Lt]);
zlim([0, 0.12]);
set(h,'LineStyle','none')
colormap(jet)
clim([0, 0.12])


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

%% Jv graph

% Create a figure for the first plot
figure;
hold on;
x = linspace(0, 500, 500);

plot(x(2:end), Jv_values1, 'LineWidth', 1.5);
plot(x(2:end), Jv_values2, 'LineWidth', 1.5);
plot(x(2:end), Jv_values3, 'LineWidth', 1.5);
plot(x(2:end), Jv_values4, 'LineWidth', 1.5);

xlabel('Tid [s]');
ylabel('J_{konv} [m s^{-1}]');
title('J_{konv} [m s^{-1}] Over Tid');
grid on;

legend('Advektion', 'Membran', 'Diffusion', 'Fouling');
%% Jv graph

% Create a figure for the first plot
figure;
hold on;
x = linspace(0, 500, 500);

% Subtract the last value of each series from the entire series
Jv_values2_centered = Jv_values2 - Jv_values2(end);
Jv_values3_centered = Jv_values3 - Jv_values3(end);
Jv_values4_centered = Jv_values4 - Jv_values4(end);

plot(x(2:end), Jv_values2_centered, 'LineWidth', 1.5);
plot(x(2:end), Jv_values3_centered, 'LineWidth', 1.5);
plot(x(2:end), Jv_values4_centered, 'LineWidth', 1.5);

xlabel('Tid [s]');
ylabel('J_{konv} - J_{SS} [m^{3} m^{-2} s^{-1}]');
title('Centreret J_{konv} Over Tid');
grid on;

legend('Membran', 'Diffusion', 'Fouling');
