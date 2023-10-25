%% Clearing
clear %Clear all global variables
close all 
clc %Clear command window
tic %Clear elapsed time counter

%% Decleration of startvariables
N = 5; % Number of Control Volumes
L = 0.02; %[m] % Domain length
h = L/(N); % Grid Spacing (Wall thickness)
k = 10;  %[W/m-K] % Diffusivity (Thermal conductivity)
Row_c = 10e6;  %[J/m^(3)-K] % Specific Heat Capacity
q = 1000e3;     %  [W/m^(3)] % Uniform Heat Generation, A Source Term
% The center of the first cell is at dx/2 & the last cell is at L-dx/2.
x = h/2 : h : L-(h/2); % Cross sectional area of the 1D domain
A = 1;          %[m^2] 
alpha = k/Row_c; % Diffusivity (Heat)
t_Final = 8; % Simulation Time
dt = 2; % Discrete Time Steps
T_a = 0;    %[\circ c] % Left Surface Temperature
T_b = 0;    %[\circ c] % Right Surface Temperature
lambda = (alpha.*dt)/(h^2);  % Parameteric Setup
% Initializing Variable
T_Old = zeros(N,1); % Unknowns at time level n
T_New = zeros(N,1); % Unknowns at time level n+1
T_Old(:) = 200; % Initial Value (Initial Condition)

%% Main

value=0;

%MainStoftransportsimulering

Arr=[]; %Main array

%Main Forloop of array
for i=1:2 %For looping through x-dimension
    for j=1:2 %For looping through t-dimension
    
    end
end



Stoftransportsimuleringsapp %%The main script for GUI.mlapp

%% Object decleration
obj = Kontrolvol;
obj = obj.MyClass(1);



