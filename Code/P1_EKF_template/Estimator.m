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
    oriEst = (2.0*rand() - 1) * const.RotationStartBound; % 1x1 matrix
    driftEst = 0.0; % 1x1 matrix
    
    % initial state variance
    % uniform dist with variance
    posVar = [0,0];   %[posVarFull(1,1),posVarFull(2,2)]; % 1x2 matrix Not 2x2 ??
    linVelVar = [0.0, 0.0]; % 1x2 matrix
    % uniform distribution (b-a)^2*1/12
    oriVar = 0;%(1.0/12)*(2*const.RotationStartBound)^2; % 1x1 matrix
    driftVar = 0.0; % 1x1 matrix
    
    % estimator variance init (initial posterior variance)
     estState.Pm = diag([posVar,oriVar,linVelVar,driftVar]);
    % estState.Pm = [posVar,oriVar,linVelVar,driftVar];
    % estimator state
    estState.xm = [posEst,oriEst,linVelEst,driftEst];
    % time of last update
    estState.tm = tm;
    return;
end

%% Estimator iteration.

% get time since last estimator update
dt = estState.tm;
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
Qd = estConst.DragNoise; 
Qr = estConst.RudderNoise; 
Qb = estConst.GyroDriftNoise; 
% compute process equation
tspan = [dt,tm];
% compute Q
Q = diag([Qr,Qd,Qb]);

%concatenate xp and var
init = [transpose(estState.xm);estState.Pm(:)];

%build matrix for coefficient
[txp,xp] = ode45(@(t,y) ODEeq(t,y,Cr,Cd,ur,ut,Q),tspan,init); %y is the process update


% define function handle for process variance update

%recover variance
var = reshape(xp(end,7:end),[6,6]);
%[tvar,var] = ode45(@(tvar,var) Pp(tvar,var,ut,Cd,ur,Cr,Q,txp,xp),tspan,estState.Pm(:));%Debug
%Pm_old = diag(estState.Pm);
%[tvar,var] = ode45(@(tvar,var) Pp(tvar,var,ut,Cd,ur,Cr,Q,txp,xp),tspan,estState.Pm); %var is the process variance update
%[tvar,var] = ode45(@(tvar,var) Pp(tvar,var,ut,Cd,ur,Cr,Q,txp,xp),tspan,Pm_old(:));

% measurement update
% compute H(t) and M(t) for measurement update
%   Build H which is Ttimes x 6
xp_t = xp(end,1:6);
H = zeros(5,6);
H(:,1) = [(xp_t(1)-pos_radioA(1))/(vecnorm(xp_t(1:2)-pos_radioA)),...
          (xp_t(1)-pos_radioB(1))/(vecnorm(xp_t(1:2)-pos_radioB)),...
          (xp_t(1)-pos_radioC(1))/(vecnorm(xp_t(1:2)-pos_radioC)),0,0];%H is T x Xp

H(:,2) = [(xp_t(2)-pos_radioA(2))/(vecnorm(xp_t(1:2)-pos_radioA)),...
          (xp_t(2)-pos_radioB(2))/(vecnorm(xp_t(1:2)-pos_radioB)),...
          (xp_t(2)-pos_radioC(2))/(vecnorm(xp_t(1:2)-pos_radioC)),0,0];
H(:,3) = [0,0,0,1,1];
H(:,4) = [0,0,0,0,0];
H(:,5) = [0,0,0,0,0];
H(:,6) = [0,0,0,1,0];

%   initialize M
M = eye(5,5); %identity matrix

%   initialize measurement noise matrix R
R = [sigmaA2,sigmaB2,sigmaC2,sigmaG2,sigmaN2];

% compute measurement z = h(xp)
h = zeros(1,5);
h(1) = sqrt(vecnorm(xp_t(1:2)-pos_radioA));
h(2) = sqrt(vecnorm(xp_t(1:2)-pos_radioB));
h(3) = sqrt(vecnorm(xp_t(1:2)-pos_radioC));
h(4) = xp_t(3) + xp_t(end);
h(5) = xp_t(3);

z = sense;
if (isinf(sense(3))) %if measurement C does not exist modify coeffs
    H(3,:)=[];
    R(3)=[];
    M=eye(4,4);
    h(3) = [];
    z(3) = [];
end
R = diag(R);

% compute K(h)
Ppt = var; %DEBUG
K = Ppt*transpose(H)/(H*Ppt*transpose(H)+M*R*transpose(M));%+ eye(size(M)) * rand());


% measurement update
estState.xm = xp_t+(z-h)*transpose(K);
%estState.Pm = (eye(6,6)-K*H)*Ppt;%DEBUG
estState.Pm = (eye(6,6)-K*H)*Ppt;

%estState.xm(3)= mod(estState.xm(3),2*pi);

% Get resulting estimates and variances
% Output quantities
posEst = estState.xm(1:2);
linVelEst = estState.xm(4:5);
oriEst = estState.xm(3);
driftEst = estState.xm(6);

posVar = [estState.Pm(1,1),estState.Pm(2,2)];%DEBUG
linVelVar = [estState.Pm(4,4), estState.Pm(5,5)];
oriVar = estState.Pm(3,3);
driftVar = estState.Pm(6,6);
%estState.Pm = diag([posVar,oriVar,linVelVar,driftVar]);

%posVar = estState.Pm(1:2);%DEBUG
%linVelVar = estState.Pm(4:5);
%oriVar = estState.Pm(3);
%driftVar = estState.Pm(6);

end

function dvardt = Pp(t,var,ut,Cd,ur,Cr,Q,tx,xp)
%odefcn to return ode solution of process variances

%find the interpolated xp
xpt = interp1(tx,xp,t);
%reshape variance to recover variance matrix shape
Pk = reshape(var,[6,6]);
%Pk = diag(var);DEBUG
%compute A
A  =zeros(6,6);
A(1,:) = [0,0,0,1,0,0];
A(2,:) = [0,0,0,0,1,0];
A(3,:) = [0,0,0,0,0,0];
A(4,:) = [0,0,-sin(xpt(3))*(tanh(ut)-Cd*(xpt(4)^2+xpt(5)^2)),-cos(xpt(3))*2*Cd*xpt(4),-cos(xpt(3))*2*Cd*xpt(5),0];
A(5,:) = [0,0,cos(xpt(3))*(tanh(ut)-Cd*(xpt(4)^2+xpt(5)^2)),-sin(xpt(3))*2*Cd*xpt(4),-sin(xpt(3))*2*Cd*xpt(5),0];
A(6,:) = [0,0,0,0,0,0];

%compute L
L = zeros(6,3);
L(1,:) = [0,0,0];
L(2,:) = [0,0,0];
L(3,:) = [Cr*ur,0,0];
L(4,:) = [0,-cos(xpt(3))*Cd*(xpt(4)^2+xpt(5)^2),0];
L(5,:) = [0,-sin(xpt(3))*Cd*(xpt(4)^2+xpt(5)^2),0];
L(6,:) = [0,0,1];
      
%compute derivative of variance matrix at timestep t
dpdt = A * Pk + Pk * transpose(A) + L*Q*transpose(L);
%dvardt = diag(dpdt);
dvardt = dpdt(:);

end

function dydt = ODEeq(t,y,Cr,Cd,ur,ut,Q)
%odefcn function to return ode solution of process states
%state vector computation
dydt = zeros(42,1);
dydt(1)=y(4);
dydt(2)=y(5);
dydt(3)=Cr*ur;
dydt(4)=cos(y(3))*(tanh(ut)-Cd*(y(4)^2+y(5)^2));
dydt(5)=sin(y(3))*(tanh(ut)-Cd*(y(4)^2+y(5)^2));
dydt(6)=0.0;

pp = reshape(y(7:end),[6,6]);
%create matrix A
A=zeros(6,6);
A(1,:) = [0,0,0,1,0,0];
A(2,:) = [0,0,0,0,1,0];
A(3,:) = [0,0,0,0,0,0];
A(4,:) = [0,0,-sin(y(3))*(tanh(ut)-Cd*(y(4)^2+y(5)^2)),-cos(y(3))*2*Cd*y(4),-cos(y(3))*2*Cd*y(5),0];
A(5,:) = [0,0,cos(y(3))*(tanh(ut)-Cd*(y(4)^2+y(5)^2)),-sin(y(3))*2*Cd*y(4),-sin(y(3))*2*Cd*y(5),0];
A(6,:) = [0,0,0,0,0,0];

%A_sparse = blkdiag(A(1,:),A(2,:),A(3,:),A(4,:),A(5,:),A(6,:));
A = transpose(A);

%create matrix L
%compute L
L = zeros(6,3);
L(1,:) = [0,0,0];
L(2,:) = [0,0,0];
L(3,:) = [Cr*ur,0,0];
L(4,:) = [0,-cos(y(3))*Cd*(y(4)^2+y(5)^2),0];
L(5,:) = [0,-sin(y(3))*Cd*(y(4)^2+y(5)^2),0];
L(6,:) = [0,0,1];

temp = A*pp+pp*transpose(pp)+L*Q*transpose(L);
dydt(7:end) = temp(:);
%compute coefficient
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

