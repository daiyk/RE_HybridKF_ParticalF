function const = EstimatorConst()
% 
% Define the physical constants that are available to the estimator.
%
%
% Class:
% Recursive Estimation
% Spring 2018
% Programming Exercise 1
%
% --
% ETH Zurich
% Institute for Dynamic Systems and Control
% Raffaello D'Andrea, Matthias Hofer, Carlo Sferrazza
% hofermat@ethz.ch
% csferrazza@ethz.ch
%


%% The wall contour - (x,y) coordinates of corners p1 to p10 as in Table 1
const.contour = [0.5 0; 
                 2 0;
                 2 0.5;
                 1.5 1;
                 1.5 2;
                 0.25 2;
                 0.25 1.75;
                 0 1.75;
                 0 1;
                 0.5 1];      
             
%% Initialization
const.pA = [0.9,0.4]; % Center point pA of the initial position distribution
const.pB = [1.1,1.6]; % Center point pB of the initial position distribution
const.d = 0.3;  % Radius of shaded regions for initialization

const.phi_0 = pi/4; % Initial heading is uniformly distr. in [-phi_0,phi_0]

%% Noise properties
% process noise
const.sigma_phi = 0.05; % Parameter for process noise v_phi
const.sigma_f = 0.01; % Parameter for process noise v_f

% measurement noise
const.epsilon = 0.02; % Parameter for measurement noise w
