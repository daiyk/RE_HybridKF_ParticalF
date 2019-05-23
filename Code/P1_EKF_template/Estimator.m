function [posEst,linVelEst,oriEst,driftEst,...
          posVar,linVelVar,oriVar,driftVar,estState] = ...
    Estimator(estState,actuate,sense,tm,estConst)
% [posEst,linVelEst,oriEst,driftEst,...
%    posVar,linVelVar,oriVar,driftVar,estState] = 
% Estimator(estState,actuate,sense,tm,estConst)
%
% The estimator.
%
% The function is initialized for tm == 0, otherwise the estimator does an
% iteration step (compute estimates for the time step k).
%
% Inputs:
%   estState        previous estimator state (time step k-1)
%                   May be defined by the user (for example as a struct).
%   actuate         control input u(k-1), [1x2]-vector
%                   actuate(1): u_t, thrust command
%                   actuate(2): u_r, rudder command
%   sense           sensor measurements z(k), [1x5]-vector, INF entry if no
%                   measurement
%                   sense(1): z_a, distance measurement a
%                   sense(2): z_b, distance measurement b
%                   sense(3): z_c, distance measurement c
%                   sense(4): z_g, gyro measurement
%                   sense(5): z_n, compass measurement
%   tm              time, scalar
%                   If tm==0 initialization, otherwise estimator
%                   iteration step.
%   estConst        estimator constants (as in EstimatorConst.m)
%
% Outputs:
%   posEst          position estimate (time step k), [1x2]-vector
%                   posEst(1): p_x position estimate
%                   posEst(2): p_y position estimate
%   linVelEst       velocity estimate (time step k), [1x2]-vector
%                   linVelEst(1): s_x velocity estimate
%                   linVelEst(2): s_y velocity estimate
%   oriEst          orientation estimate (time step k), scalar
%   driftEst        estimate of the gyro drift b (time step k), scalar
%   posVar          variance of position estimate (time step k), [1x2]-vector
%                   posVar(1): x position variance
%                   posVar(2): y position variance
%   linVelVar       variance of velocity estimate (time step k), [1x2]-vector
%                   linVelVar(1): x velocity variance
%                   linVelVar(2): y velocity variance
%   oriVar          variance of orientation estimate (time step k), scalar
%   driftVar        variance of gyro drift estimate (time step k), scalar
%   estState        current estimator state (time step k)
%                   Will be input to this function at the next call.
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

%% Initialization
if (tm == 0)
    % Do the initialization of your estimator here!
    % initial state mean
    const=EstimatorConst();
    
    posPol = randPos(const.StartRadiusBound); % 1x2 matrix pos est in Polar coord
    posEst =[posPol(1)*sin(posPol(2)),posPol(1)*cos(posPol(2))];
    linVelEst = [0.0, 0.0]; % 1x2 matrix
    %affine transform uniform rand in [0,1] to [-thetaB,thetaB]
    oriEst = (2.0*rand() - 1) / const.RotationStartBound; % 1x1 matrix
    driftEst = 0.0; % 1x1 matrix
    
    % initial state variance
    % uniform dist with variance
    posVarFull = uniDisk(const.StartRadiusBound,posPol);
    posVar = [posVarFull(1,1),posVarFull(2,2)];%[0.0, 0.0]; % 1x2 matrix Not 2x2 ??
    linVelVar = [0.0, 0.0]; % 1x2 matrix
    oriVar =  1.0/12*(2*const.RotationStartBound)^2; % 1x1 matrix uniform distribution variance
    driftVar = 0.0; % 1x1 matrix
    
    % estimator variance init (initial posterior variance)
    estState.Pm = diag([posVar,oriVar,linVelVar,driftVar]);
    % estimator state
    estState.xm = [posEst,oriEst,linVelEst,driftEst];
    % time of last update
    estState.tm = tm;
    return;
end

%% Estimator iteration.

% get time since last estimator update
dt = tm - estState.tm;
estState.tm = tm; % update measurement update time

% prior update
% get constant for calculation
Cd = estConst.dragCoefficient;
Cr = estConst.rudderCoefficient;
ut = actuate(1);
ur = actuate(2);
pos_radioA = estConst.pos_radioA;
pos_radioB = estConst.pos_radioB;
pos_radioC = estConst.pos_radioC;
sigmaA2 = estConst.DistNoiseA;
sigmaB2 = estConst.DistNoiseB;
sigmaC2 = estConst.DistNoiseC;
sigmaG2 = estConst.GyroNoise;
sigmaN2 = estConst.CompassNoise;
% compute process equation
tspan = [tm,tm+dt];
[txp,xp] = ode45(@(t,y) Xp(t,y,Cr,Cd,ur,ut),tspan,estState.xm); %y is the process update

% define function handle for process variance update
[tvar,var] = ode45(@(tvar,var) Pp(tvar,var,ut,Cd,txp,xp),tspan,estState.Pm); %var is the process variance update

% measurement update
% compute H(t) and M(t) for measurement update
%   Build H which is T x Xp.size
H = zeros(6,5);
H(1,:) = [(xp(end,1)-pos_radioA(1))/vecnorm(xp(end,1:2)-pos_radioA),...
          (xp(end,1)-pos_radioB(1))/vecnorm(xp(end,1:2)-pos_radioB),...
          (xp(end,1)-pos_radioC(1))/vecnorm(xp(end,1:2)-pos_radioC),0,0];%H is T x Xp

H(2,:) = [(xp(end,2)-pos_radioA(2))/vecnorm(xp(end,1:2)-pos_radioA),...
          (xp(end,2)-pos_radioB(2))/vecnorm(xp(end,1:2)-pos_radioB),...
          (xp(end,2)-pos_radioC(2))/vecnorm(xp(end,1:2)-pos_radioC),0,0];
H(3,:) = [0,0,0,1,1];
H(4,:) = [0,0,0,0,0];
H(5,:) = [0,0,0,0,0];
H(6,:) = [0,0,0,1,0];

%   initialize M
M = eye(5,5); %identity matrix

%   initialize measurement noise matrix R
R = [sigmaA2,sigmaB2,sigmaC2,sigmaG2,sigmaN2];

% compute measurement z = h(xp)
h = zeros(1,5);
h(1) = sqrt(vecnorm(xp(end,1:2)-pos_radioA));
h(2) = sqrt(vecnorm(xp(end,1:2)-pos_radioB));
h(3) = sqrt(vecnorm(xp(end,1:2)-pos_radioC));
h(4) = xp(end,3) + xp(end,end);
h(5) = xp(end,3);

z = sense;
if (isinf(sense(3))) %if measurement C does not exist modify coeffs
    H(:,3)=[];
    R(3)=[];
    M=eye(4,4);
    h(3) = [];
    z(3) = [];
end
R = diag(R);
H = transpose(H);
% compute K(h)
Ppt = reshape(var(end,:),[6,6]);
K = Ppt*transpose(H)/(H*Ppt*transpose(H)+M*R*transpose(M));

% compute estimates and variance
xpt = xp(end,:);
estState.xm = xpt+(z-h)*transpose(K);
estState.Pm = (eye(6,6)-K*H)*Ppt;
% Get resulting estimates and variances
% Output quantities
posEst = estState.xm(1:2);
linVelEst = estState.xm(4:5);
oriEst = estState.xm(3);
driftEst = estState.xm(6);

posVar = [estState.Pm(1,1),estState.Pm(1,1)];
linVelVar = [estState.Pm(4,4), estState.Pm(5,5)];
oriVar = estState.Pm(3,3);
driftVar = estState.Pm(6,6);

end

function dvardt = Pp(t,var,ut,Cd,tx,xp)
%odefcn to return ode solution of process variances

%find the interpolated time
xpt = interp1(tx,xp,t);
dvardt = zeros(6,6);
dvardt(1,:) = [0,0,0,1,0,0];
dvardt(2,:) = [0,0,0,0,1,0];
dvardt(3,:) = [0,0,0,0,0,0];
dvardt(4,:) = [0,0,-sin(xpt(3))*tanh(ut),2*Cd*xpt(4),-2*Cd*xpt(5),0];
dvardt(5,:) = [0,0,cos(xpt(3))*tanh(ut),2*Cd*xpt(4),-2*Cd*xpt(5),0];
dvardt(6,:) = [0,0,0,0,0,0];
dvardt = dvardt(:);
end

function dydt = Xp(t,y,Cr,Cd,ur,ut)
%odefcn function to return ode solution of process states

dydt=zeros(6,1);
dydt(1)=y(4);
dydt(2)=y(5);
dydt(3)=Cr*ur;
dydt(4)=cos(y(3))*(tanh(ut)-Cd*(y(4)^2+y(5)^2));
dydt(5)=sin(y(3))*(tanh(ut)-Cd*(y(4)^2+y(5)^2));
dydt(6)=0.0;
end

function pos = randPos(rMax)
%Uniform sampling of a disk with radius rMax


%generate two random numbers in (0,1) for sampling purpose
randn1=rand();
randn2=rand();

%uniform sampling start point(in polar coordinate) that with radius r
r=rMax*sqrt(randn1);
theta=2.0*pi*randn2;

%initialize pos
pos=[r,theta];
end

function var = uniDisk(rMax,posPol)
%Return variance of uniform sampling of disk with radius rMax

r = posPol(1);
theta=posPol(2);
varPol = [(rMax)^2/12.0,0.0;(2*pi)^2/12,0.0];

%variance pos in polar coordinate
%var[r]=E[r^2]-E[r]^2, E[r^2]=1/2*R^2, E[r]=2/3*R
%var[theta] = (2pi-0)^2/12
varPol = [1.0/18*r^2,0;0,(2*pi)^2/12];
%Jacobian matrix between Polar and Cartesian
Jab = [cos(theta),-r*sin(theta);sin(theta),r*cos(theta)];

%variance transformation: var(x,y)=Jab*var(r,theta)*Jab.T
var = Jab*varPol*transpose(Jab);


end

