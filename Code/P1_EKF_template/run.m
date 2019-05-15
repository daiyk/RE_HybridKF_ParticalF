function trackErrorNorm=run(simConst,estConst,doplot,seed)
%
%
% Main function for Extended Kalman Filter programming exercise.
% It is used to simulate the true model, call the estimator and show the 
% results.
%
%
% Class:
% Recursive Estimation
% Spring 2019
% Programming Exercise 1
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Raffaello D'Andrea, Matthias Hofer, Carlo Sferrazza
% hofermat@ethz.ch
% csferrazza@ethz.ch
%

% clear command window, close figures
% clc;
close all;

if nargin==0
   % Define the simulation constants that are used in simulation, but not 
   % accessible to the estimator.  
   simConst = SimulationConst();
   
   % Define the constants accessible to the estimator.
   estConst = EstimatorConst();
   
   % Generate plots by default.
   doplot = true;
   
   % use random seed
   seed = 0;
end

%% Setup

% Set the random number generator state.
% Uncomment to make results reproducable. This setting was used to generate
% the plot in the problem description.
if seed ~= 0
    rand('seed',seed);
    randn('seed',seed);
end

%% Simulation
% The function 'Simulator' simulates the boat dynamics and generates
% measurements.
[tm, state, drift, input, sense] = Simulator( simConst );

% state contains rolling out of [p_x, p_y, s_x, s_y, phi]: time index from 1 to N

% input contains rolling out of [u_t, u_r]: time index from 1 to N-1

% sense contains rolling out of [z_a, z_b, z_c, z_g, z_n] (contains inf):
% time index from 1 to N, with first measurement inf (to indicate no
% measurement at initial step)

%% Run the Estimator
N = simConst.N;

% Initialize the estimator.  
estState = [];
posEst = zeros(N,2);
linVelEst = zeros(N,2);
oriEst = zeros(N,1);
driftEst = zeros(N,1);
posVar = zeros(N,2);
linVelVar = zeros(N,2);
oriVar = zeros(N,1);
driftVar = zeros(N,1);

% Call the estimator for each time step.

% Initialization
[posEst(1,:),linVelEst(1,:),oriEst(1),driftEst(1),...
 posVar(1,:),linVelVar(1,:),oriVar(1),driftVar(1),estState] = ...
     Estimator(estState,[],[],tm(1),estConst);
% Remaining steps
tstart = tic;
for n = 2:N
    [posEst(n,:),linVelEst(n,:),oriEst(n),driftEst(n),...
     posVar(n,:),linVelVar(n,:),oriVar(n),driftVar(n),estState] = ...
         Estimator(estState,input(n-1,:),sense(n,:),tm(n),estConst);
end
tEstAvg = toc(tstart)/(N-1)

%% The results

% Calculate the total tracking error as an indicator of your estimator's performance.
trackError = [state(:,1) - posEst(:,1);state(:,2) - posEst(:,2)];
trackErrorNorm = sqrt(trackError'*trackError/N);

% Plots of the results.
if doplot
    % planar position plot
    figure(1)
    hold all
    axis equal
    plot(state(:,1),state(:,2),'r.')
    plot(posEst(:,1),posEst(:,2),'b.')
    xlabel('x position [m]')
    ylabel('y position [m]')
    legend('true state', 'estimate')
    grid on
    
    % x position over time
    figure(2)
    hold all
    plot(tm, state(:,1),'r.')
    plot(tm, posEst(:,1),'b.')
    xlabel('time [s]')
    ylabel('x position [m]')
    legend('true state', 'estimate')
    grid on
    
    % y position over time
    figure(3)
    hold all
    plot(tm, state(:,2),'r.')
    plot(tm, posEst(:,2),'b.')
    xlabel('time [s]')
    ylabel('y position [m]')
    legend('true state', 'estimate')
    grid on
        
    % orientatin over time
    figure(4)
    hold all
    plot(tm, state(:,5),'r.')
    plot(tm, oriEst(:),'b.')
    xlabel('time [s]')
    ylabel('orientation [rad]')
    legend('true state', 'estimate')
    grid on
    
    % x velocity over time
    figure(5)
    hold all
    plot(tm, state(:,3),'r.')
    plot(tm, linVelEst(:,1),'b.')
    xlabel('time [s]')
    ylabel('x velocity  [m/s]')
    legend('true state', 'estimate')
    grid on
    
    % y velocity over time
    figure(6)
    hold all
    plot(tm, state(:,4),'r.')
    plot(tm, linVelEst(:,2),'b.')
    xlabel('time [s]')
    ylabel('y velocity [m/s]')
    legend('true state', 'estimate')
    grid on
    
    % drift over time
    figure(7)
    hold all
    plot(tm, drift,'r.')
    plot(tm, driftEst(:),'b.')
    xlabel('time [s]')
    ylabel('gyro drift [rad]')
    legend('true state', 'estimate')
    grid on
end
return;