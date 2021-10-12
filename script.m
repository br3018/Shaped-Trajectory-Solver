%%  Script for drawing shaped low-thrust trajectories

%   Housekeeping
clear
clc
close all
fig=1;

%   Setting constants
G = 6.673e-11; %m^3/kg/s^2
g_acc = 9.81; % (ms^-2)

%   Writing spacecraft parameters
%   SPACECRAFT MANEUVERING MASSESS SET THESE DEPNENDING ON DESIRED
%   TRAJECTORY
%   Earth leg:
 SC.m_0 = 1174215; % (kg) 
%   Mars leg:
% SC.m_0 = 1013623; % (kg) (initial spacecraft mass before maneuver)

SC.ISP = 2800; % (s)
SC.spec_thrust = 18.5; %(kW/N)

%   Writing planet parameters
E.R = 6378.1363d3; % Earth radius (m)
E.M = 5.9742e24; % Earth mass (kg)
E.mu=E.M*G; % Standard gravitational parameter (m^3s^-2)
E.name='Earth';

M.R = 3389.5d3; % Mars radius (m)
M.M = 6.39*10^23; % Mars mass (kg)
M.mu=M.M*G; % Standard gravitational parameter (m^3s^-2)
M.name='Mars';

S.M= 1.989*10^30; % Sun mass (kg)
S.mu = S.M*G;
S.name='Sun';

%   Setting central body
P0=S;
%P0=E;

%   Setting departure body (P=Planet)
P1=E;
%P1=M;
%   Setting target body
%P2=E;
P2=M;

% %   Setting departure/arrival time (EARTH->MARS SAMPLE)
% t.depart = juliandate(2065,01,11);
% t.arrive = juliandate(2066,03,07);
% %   Setting departure/arrival time (MARS->EARTH SAMPLE)
% t.depart = 2476040;
% t.arrive = t.depart+440;
%   Setting departure/arrival time (EARTH->MARS OPTIMAL)
% t.depart = juliandate(2050,02,25);
% t.arrive = juliandate(2051,06,30);
%  2047 launch window EM
t.depart = juliandate(datetime(747973,'ConvertFrom','datenum'));
t.arrive = juliandate(datetime(747973+460,'ConvertFrom','datenum'));
%  2047 launch window ME
% t.depart = juliandate(datetime(748743,'ConvertFrom','datenum'));
% t.arrive = juliandate(datetime(748743+440,'ConvertFrom','datenum'));
N_rev = 0;

t.delta = (t.arrive-t.depart)*24*60*60; % (s)

%   Calculating planet state vectors at specified times
%   J2000 -> Ecliptic rotation
eta = 23.5; % rotation angle
R = [
    1,  0,              0;
    0,  cosd(eta),     sind(eta);
    0,  -sind(eta),    cosd(eta)];

[P1.rd,P1.vd]=planetEphemeris(t.depart,P0.name,P1.name,'423');
[P2.rd,P2.vd]=planetEphemeris(t.depart,P0.name,P2.name,'423');
P1.rd=R*P1.rd'.*1000; % convert to m
P1.vd=R*P1.vd'.*1000;
P2.rd=R*P2.rd'.*1000; % convert to m
P2.vd=R*P2.vd'.*1000;


[P1.ra,P1.va]=planetEphemeris(t.arrive,P0.name,P1.name,'423');
[P2.ra,P2.va]=planetEphemeris(t.arrive,P0.name,P2.name,'423');
P1.ra=R*P1.ra'.*1000; % convert to m
P1.va=R*P1.va'.*1000;
P2.ra=R*P2.ra'.*1000; % convert to m
P2.va=R*P2.va'.*1000;

P1.ra = P1.ra';
P1.va = P1.va';
P1.rd = P1.rd';
P1.vd = P1.vd';
P2.ra = P2.ra';
P2.va = P2.va';
P2.rd = P2.rd';
P2.vd = P2.vd';


%   Propagating planet orbits for specified time
[P1.r, P1.v, P1.t] = ODE45_prop(P1.rd, P1.vd, [0 t.delta], P0.mu);
[P2.r, P2.v, P2.t] = ODE45_prop(P2.rd, P2.vd, [0 t.delta], P0.mu);

theta_f = acos(dot(P1.rd,P2.ra)/(norm(P1.rd)*norm(P2.ra))); % angle between desired position vectors for shaped trajectory
n = cross(P1.rd,P2.ra)/((norm(P1.rd)*norm(P2.ra))*sin(theta_f));
%   Checking direction of theta.diff (clockwise or anticlockwise) and
%   enforcing movement in direction of planets
if n(3) < 0 
    theta_f = 2*pi - theta_f;
end

%   Calculating shaped trajectory coefficients
theta_f = theta_f + N_rev*2*pi;
q = 7; % z shaping parameter
[a, b, c, d, e, f, g, a_z, b_z, c_z, d_z, DU, TU, zDU, zTU, error] = shape_trajectory(P1.rd, P2.ra, P1.vd, P2.va, t.delta, P0.mu, N_rev, q);
%   Checking for method solution:
TOF_check = calc_TOF(a, b, c, d, e, f, g, 1, theta_f)*TU;
if ~isreal(TOF_check)
    disp('WARNING Trajectory not found ')
else
    %   NOTE UNITS ARE CONVERTED TO CANONICAL UNITS
    %   Plotting shaped trajectory
    delta_theta = 0.05; 
    theta = [0:delta_theta:theta_f];
    r = 1./(a+b*theta+c*theta.^2+d*theta.^3+e*theta.^4+f*theta.^5+g*theta.^6);
    z = a_z + b_z*theta + c_z*theta.^(q-1) + d_z*theta.^q;
    theta_dot = sqrt((1./r.^4).*(1./((1./r)+2*c+6*d*theta+12*e*theta.^2+20*f*theta.^3+30*g*theta.^4)));
    %   Calculating time scale (APPROX)
    t.array = zeros(1,length(theta));
    for i = 2:length(theta)
        t.array(i) = t.array(i-1)+1/theta_dot(i)*delta_theta;
    end
    %   Scaling from canonical units 
    t.array = t.array*TU;

    %   Calculating flight path angle
    gamma = atan(-r.*(b+2*c*theta+3*d*theta.^2+4*e*theta.^3+5*f*theta.^4+6*g*theta.^5));
    
    %   Calculating thrust acceleration
    T_a = -1./(2*r.^3.*cos(gamma)).*(6*d+24*e*theta+60*f*theta.^2+120*g*theta.^3-tan(gamma)./r)./(1./r+2*c+6*d*theta+12*e*theta.^2+20*f*theta.^3+30*g*theta.^4).^2;
    %   Converting thrust acceleration to SI units
    T_a_SI = T_a*DU/TU^2;    
    
    %   Calculating max thrust acceleration required
    T_a_SI_max = max(abs(T_a_SI));
    
    %   Calculating delta V of trip
    delta_V = trapz(theta,T_a./theta_dot); % (DU/TU)
    delta_V_SI = delta_V*DU/TU;
    
    %   Calculating propellant loss
    error = false;
    SC.m = zeros(1,length(theta));
    SC.m(1) = SC.m_0;
    T = zeros(1,length(theta));
    T(1) = T_a_SI(1)*SC.m(1); % (N) (actual thrust)
    for i = 2:length(theta)
        SC.m(i) = SC.m(i-1) - abs(T(i-1))/(SC.ISP*g_acc)*(t.array(i)-t.array(i-1));
        T(i) = T_a_SI(i)*SC.m(i);
        if SC.m(i)<=0
            error = true;
        end
    end
    %   Mass check
    if error == true
        disp('Mass becomes negative with current parameters')
    end
    
    %% Plotting trajectory
    
    figure(fig)
    fig = fig+1;
    polarplot(theta,r)
    [x, y, z] = pol2cart(theta, r, z);
    %   Calculating rotation required to put into ecliptic frame       
    
    %   Plotting thrust acceleration variation
    figure(fig)
    fig = fig+1;
    hold on
    grid on 
    grid minor
    plot(theta,T_a)
    xlabel('\theta (rad)')
    ylabel('Thrust acceleration (DU/TU^2)','interpreter','TeX')
    
    figure(fig)
    fig = fig+1;
    hold on
    grid on 
    grid minor
    plot(theta,T_a_SI)
    xlabel('\theta (rad)')
    ylabel('Thrust acceleration (m/s^2)','interpreter','TeX')
    
    %   Plotting flight path angle
    figure(fig)
    fig = fig+1;
    hold on
    grid on
    grid minor
    plot(theta,gamma)
    xlabel('\theta (rad)')
    ylabel('\gamma (rad)','interpreter','TeX')
    
    
    %   Plotting flight path angle with time
    figure(fig)
    fig = fig+1;
    hold on
    grid on
    grid minor
    plot(t.array/(60*60*24),gamma)
    xlabel('Time (days)')
    ylabel('\gamma (rad)','interpreter','TeX')
    
    %   Plotting spacecraft mass over time
    figure(fig)
    fig = fig+1;
    hold on
    grid on
    grid minor
    plot(t.array/(60*60*24),SC.m)
    xlabel('Time (days)')
    ylabel('Spacecraft mass (kg)','interpreter','TeX')
    
    %   Plotting spacecraft actual thrust over time (N)
    figure(fig)
    fig = fig+1;
    hold on
    grid on
    grid minor
    plot(t.array/(60*60*24),T)
    xlabel('Time (days)')
    ylabel('Thrust (N)','interpreter','TeX')
    
    %   Plotting propulsion system power consumption over time 
    figure(fig)
    fig = fig+1;
    hold on
    grid on 
    grid minor
    plot(t.array/(60*60*24),(abs(T)*SC.spec_thrust)/1000)
    xlabel('Time (days)')
    ylabel('Propulsion system power consumption (MW)','interpreter','TeX')
    
    %   Plotting spacecraft orbital radius over time
    figure(fig)
    fig = fig+1;
    hold on
    grid on 
    grid minor
    plot(t.array/(60*60*24),r*DU/1000)
    xlabel('Time (days)')
    ylabel('Orbital radius (km)','interpreter','TeX')
    
    %   Calculating spacecraft distance wrt Earth
    
    disp(['Trip time: ',num2str(t.delta/(60*60*24)), ' days'])
    disp(['Delta V required: ',num2str(abs(delta_V_SI)),' m/s'])
    disp(['Max acceleration required: ',num2str(T_a_SI_max),' m/s^2'])
    %   Mass sanity check
    disp(['Mass loss: ',num2str(SC.m_0*(1-exp(-abs(delta_V_SI)/(SC.ISP*g_acc)))),'kg (from rocket equation)'])

    
end

%   Rotating about z axis 
%   Angle between Pos1 and x axis
K = [1; 0; 0];
%   Calculating required theta values
theta1 = acos(dot(P1.rd',K)/norm(P1.rd)); % angle between desired position vectors
n1 = cross(P1.rd',K)/(norm(P1.rd)*sin(theta1));
%   Checking direction of theta1 (clockwise or anticlockwise) and
%   enforcing clockwise
if n1(3) < 0 % negative z implies 2*pi>theta.diff>pi
    theta1 = 2*pi-theta1;
end
%   Forming rotation matrix
R1 = [
    cos(theta1),  sin(theta1),    0;
    -sin(theta1),  cos(theta1),     0;
    0,            0,               1];
sc_pos = R1*[x; y; z];
x = sc_pos(1,:);
y = sc_pos(2,:);
z = sc_pos(3,:);


%   Setting up plot 
figure(fig)
fig=fig+1;
hold on
grid minor
plot3(P1.r(:,1),P1.r(:,2),P1.r(:,3))
plot3(P2.r(:,1),P2.r(:,2),P2.r(:,3))
plot3(x*DU,y*DU,z*zDU)
view(90,0)
% % %   Plotting planets for graphics
% % planet3D(P0.name,[0,0,0],[],[],'m',[])
% % planet3D(P1.name,P1.rd,[],[],'m',[])
% % planet3D(P2.name,P2.ra,[],[],'m',[])