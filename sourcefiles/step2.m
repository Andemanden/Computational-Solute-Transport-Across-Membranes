%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  GITHUB - ALL RIGHTS RESERVED   %%%
%%%                                 %%%
%%% PROPERTY OF AALBORG UNIVERSITY  %%%
%%%         CREATED BY:             %%%
%%%  GROUP 3 - 3RD SEMESTER- 2023   %%%
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
D = 0;      % Diffusivity coefficient H2PO4- in water [m^2 s^-1]
Cf = 0.1;             % Ionconcentration of the feed [mol L^-1]
DeltaP = 35;          % Change in outer pressure [bar]
Am = 0.0308;          % Surfacearea of the membrane [m^2]
k0 = 5.7311*10^-7;    % Initial water permeability [m^3 m^-2 bar^-1 s^-1] 
mu = 0.8903*10^-9;    % Viscosity of water [bar∙s]
Rm = 1/(mu*k0);       % Membraneresitance of water [m^-1]
alpha = 1*10^16;      % Specific fouling resistance [m · kg-1] or [m mol^-1]
JuScalar = 0.5;       % Percipitate Advection Coefficient
sig_i = 0.1;          % Rejection of ions 
volprc = dx*Am*1000;  % Volume per cell in Liters

% PHYSICAL CONSTANTS

R = 8.31415*10^-2;  % Gasconstant [L Bar mol^-1 K^-1]
T = 273.15+25;      % Temperature [K]
 
% Anonymous functions

Udf = @(conc) 0.0011*exp(36.328*conc) + 0.263253; % Concentration of percipitate in respect to Ciw. [mol L^-1]
Rf = @(conc) alpha*(conc*volprc/Am);                % Foulingresistance [m^-1]
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
            Jb = Jtotv*Udf(Cf)*(dt/dx)*JuScalar;
            Cop = Cop + Jb; % Total percipitate after advection [mol L^-1]
            Cbw = Cudf + Cop;
    end    
    for i = 1:domain_steps
        Ciw = C(domain_steps, j-1); % Concentration of ions at the wall [mol L^-1]

        C(1, :) = Cf; % Set the leftmost boundary to Cf [mol L^-1]

        Jtotv = (k0*(DeltaP-(1*R*T*(Ciw))));  % Volume flux = Jtotv , in terms of osmotic pressure (DeltaP, R, T, C_(i_w)) and k. [m/s]

        if i == 1 % First Cell (no left neighbor)
            Jdiff = 0;       % Diffusive ionflux
            Jadv = 0;        % Advective flux

        elseif i == domain_steps % Membrane wall cell (no right neighbor)
            Jdiff = (D * (C(domain_steps - 1, j-1) - C(domain_steps, j-1))) / dx^2;           % Diffusive ionflux [mol · m-2 · s-1]
            Jadv = -Jtotv * (C(domain_steps, j-1)*(1-sig_i) - C(domain_steps - 1, j-1)) / dx; % Advective flux [mol · m-2 · s-1]

        else    % Bulk Cells
            % Calculate the second derivative normally
            Jdiff = D * (C(i + 1, j-1) - 2 * C(i, j-1) + C(i - 1, j-1)) / dx^2; % Diffusive ionflux at the membrane wall [mol · m-2 · s-1]
            Jadv = -Jtotv * (C(i, j-1) - C(i-1, j-1)) / dx;                     % Advective flux at the membrane wall [mol · m-2 · s-1]
        end
        % Apply the diffusion-advection equation
        C(i, j) = C(i, j-1) + dt * (Jdiff + Jadv); % Apply the diffusion-advection equation [mol L^-1]
        Jv_values2(j) = Jtotv;                   % Storing convective flux values in respect to time
        Cbw_values(j) = Cbw;                       % Storing precipitate values in respect to time
        AS(j) = Jtotv * dt/dx;                     % Storing Stability plot values
        % Main Stability
        MainS=(1-2*DS-AS(j));                      % Calculation of specific main stability to time
        MainS_values(j)=MainS;                     % Storing main stability values
    end
end
%% conservation and error

Ciw = [0,C(domain_steps, :)]; % Making the matrix line up properly and have the right size
Ciw(end) = [];

Systemdiff = [0, diff(sum(C*volprc,1))]; % Difference in system concentration per dt

inflow = Cf*(Jv_values2)*(dt/dx)*volprc;                 % Input [mol]
outflow = Ciw.* Jv_values2*(dt/dx)*volprc*(1-sig_i);     % Output [mol]
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
h = surf(X, T, (C/Cf)-1);
xlabel('Position [m]');
ylabel('Tid [s]');
zlabel('CF');
title('CF over tid og position');
xlim([0, Lx]);
ylim([0, Lt]);
zlim([0, 0.12]);
set(h,'LineStyle','none')
colormap(jet)
clim([0, 0.12])

%%
save('step2jv', 'Jv_values2')
