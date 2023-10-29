clear
close all
clc

% Define parameters
domain_length = 0.5;   % Domain length (meters)
domain_steps = 100;    % Number of spatial domain_steps points
dx = domain_length / (domain_steps - 1);

time_length = 0.5;      % Time length (seconds)
time_steps = 500;       % Number of time steps
D = 0.001;

velocity = 0.2;
feed_conc = 1;
rejection_rate = 0.5;

% Create an instance of the OONumerical class
sim = OONumerical(domain_length, domain_steps, time_length, time_steps, D, velocity, feed_conc, rejection_rate);

% Run the simulation
sim.simulate();

% Define fractions of time steps you want to visualize
time_fraction = [0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9];

% Plot the results
sim.plotResults(time_fraction);