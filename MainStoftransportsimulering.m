%% Clearing
clear %Clear all global variables
close all 
clc %Clear command window
tic %Clear elapsed time counter

%% Decleration of startvariables
D = 0.1; % Diffusionskoefficient (tilpas efter behov)
u = 1.0; % Hastighed (tilpas efter behov)
rho = 1.0; % Densitet (tilpas efter behov)
% Create domain
L = 1.0; %[m] % Domain length
N = 50; % Number of Control Volumes
dx = L / N; % Cellestørrelse
x = linspace(0, L, N + 1); % Cellegrænser
h = L/(N); % Grid Spacing (Wall thickness)

q = 1000e3;     % Viscosity
% The center of the first cell is at dx/2 & the last cell is at L-dx/2.
x = h/2 : h : L-(h/2); % Cross sectional area of the 1D domain
A = 1;          %[m^2] 
alpha = k/Row_c; % Diffusivity (Heat)
tMax = 1.0; % Simuleringsvarighed
dt = 2; % Discrete Time Steps (Descrete time increment)
t = 0:dt:tMax;
R=8.31415; % Gasconstant []
C0 = zeros(1, N+1);
C0(1) = 1.0; % Koncentrationen ved den ene ende

%% Main

value=0;

%MainStoftransportsimulering

Arr=[ 1 2 3]; %Main array

%Main Forloop of array
for i=1:length(Arr) %For looping through x-dimension
    for j=1:length(Arr(j)) %For looping through t-dimension
        
    end
end



Stoftransportsimuleringsapp %%The main script for GUI.mlapp

%% Object decleration
%obj1 = Kontrolvol;
%obj1 = Kontrolvol.MyClass(1);

%% Plotfunctions

function data=loadData()
    data=Arr;
end


