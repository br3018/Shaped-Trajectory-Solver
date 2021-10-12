%%  Script for finding porkchop plots of continuous thrust trajectories

%   Housekeeping
clear
clc
close all

fig=1;
tic

%   Setting constants
G = 6.673e-11; %m^3/kg/s^2
g_acc = 9.81; % (ms^-2) gravitational acceleration 

%   Writing spacecraft parameters
SC.m_0 = 700*10^3; % (kg) (initial spacecraft mass before maneuver)
SC.ISP = 2800; % (s)

%   Writing planet parameters
E.R = 6378.1363d3; % Earth radius (m)
E.M = 5.9742e24; % Earth mass (kg)
E.mu=E.M*G; % Standard gravitational parameter (m^3s^-2)
E.T = 365*24*60*60; % Time period (s)
E.name='Earth';

M.R = 3389.5d3; % Mars radius (m)
M.M = 6.39*10^23; % Mars mass (kg)
M.mu=M.M*G; % Standard gravitational parameter (m^3s^-2)
M.T = 687*24*60*60; % Time period (s)
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

disp(['Departure from ',P1.name])
disp(['Arriving at ',P2.name])

%   Pre-allocating days over which to calculate
%   Setting departure/arrival time window
t.depart_i = juliandate(2040,01,01);
t.depart_f = juliandate(2040,02,01);
t.arrive_i = juliandate(2040,01,01)+250;
t.arrive_f = juliandate(2040,02,01)+500;
N_rev = 0;

t.depart = [t.depart_i:10:t.depart_f];
t.arrive = [t.arrive_i:5:t.arrive_f];

t.delta = zeros(length(t.depart),length(t.arrive));
delta_V = zeros(length(t.depart),length(t.arrive));
T_a_max = zeros(length(t.depart),length(t.arrive));
error = strings(length(t.depart),length(t.arrive));
theta_f = zeros(length(t.depart),length(t.arrive));

%   Beginning loop of trajectory finding and characterising
iter = 0;
tot_iter = length(t.depart)*length(t.arrive); % total number of calculations
for i = 1:length(t.depart)
    for j = 1:length(t.arrive)
        iter = iter+1; 
        disp(['Iteration ',num2str(iter),'/',num2str(tot_iter)]) % progress check
        toc
        t.delta(i,j) = (t.arrive(j)-t.depart(i))*24*60*60; % (s)
        %   Checking that TOF isnt negative
        if t.delta(i,j)<=0
            % Returning arbitrarily high value of delta_V and
            % max thrust acceleration
            delta_V(i,j) = NaN;
            T_a_max(i,j) = NaN;
            t.delta(i,j) = NaN;
            error(i,j) = "TOF negative";
            continue
        else
        %   Calculating planet state vectors at specified times
        [P1.rd,P1.vd]=planetEphemeris(t.depart(i),P0.name,P1.name);
        [P2.rd,P2.vd]=planetEphemeris(t.depart(i),P0.name,P2.name);
        P1.rd=P1.rd.*1000; % convert to m
        P1.vd=P1.vd.*1000;
        P2.rd=P2.rd.*1000; % convert to m
        P2.vd=P2.vd.*1000;


        [P1.ra,P1.va]=planetEphemeris(t.arrive(j),P0.name,P1.name);
        [P2.ra,P2.va]=planetEphemeris(t.arrive(j),P0.name,P2.name);
        P1.ra=P1.ra.*1000; % convert to m
        P1.va=P1.va.*1000;
        P2.ra=P2.ra.*1000; % convert to m
        P2.va=P2.va.*1000;

        theta_f(i,j) = acos(dot(P1.rd,P2.ra)/(norm(P1.rd)*norm(P2.ra))); % angle between desired position vectors for shaped trajectory
        n = cross(P1.rd,P2.ra)/((norm(P1.rd)*norm(P2.ra))*sin(theta_f(i,j)));
        %   Checking direction of rotation (METHOD WILL BREAKDOWN OVER HIGH
        %   INCLINATION CHANGES)
        if n(3) < 0 
            theta_f(i,j) = theta_f(i,j) + pi;
        end
        theta_f(i,j) = theta_f(i,j) + N_rev*2*pi;
        
        %   Calculating shaped trajectory coefficients
        
        [a, b, c, d, e, f, g, DU, TU, coeff_error] = shape_trajectory(P1.rd, P2.ra, P1.vd, P2.va, t.delta(i,j), P0.mu, N_rev);
        %   Checking trajectory coefficients are valid
        if coeff_error==true
            % Returning arbitrarily high value of delta_V and
            % max thrust acceleration
            delta_V(i,j) = NaN;
            T_a_max(i,j) = NaN;
            error(i,j) = "fsolve failed";
            continue
        end
        %   Characterising trajectory
        [delta_V(i,j), T_a_max(i,j), error(i,j)] = characterise_shaped_traj(a, b, c, d, e, f, g, DU, TU, theta_f(i,j));
        end
    end
end

%   Saving output data as backup 
save('output_backup.mat')

%   Setting contour plot for max acceleration values
% %   Filtering out anomalyous values (arbitraritly set at 0.1m/s^2)
% i_anom = find(T_a_max > 0.1);
% T_a_max(i_anom) = NaN;
%   Setting interesting acceleration levels
T_a_levels = [0.0002:0.0001:0.001];

%   Setting interesting delta V values 
delta_V_levels = [1000:500:7000];


%   Setting contour plot for thrust accelerations
figure(fig)
fig = fig +1;
hold on 
grid on
grid minor
contour(t.arrive,t.depart,T_a_max,T_a_levels, 'ShowText', 'on');
contour(t.arrive,t.depart,t.delta/(60*60*24),'color', 'k', 'ShowText', 'on');
legend({'Thrust acceleration (m/s^2)','TOF (days)'}, 'Interpreter', 'TeX')
title('Thrust acceleration (m/s^2)','Interpreter','TeX')
xlabel('Date of arrival')
ylabel('Date of departure')
dateaxis('x',17, datetime(t.arrive_i,'ConvertFrom','juliandate'))
dateaxis('y',17, datetime(t.depart_i,'ConvertFrom','juliandate'))

figure(fig)
fig = fig +1;
hold on 
grid on
grid minor
%   Setting contour plot for delta V values
contour(t.arrive,t.depart,abs(delta_V),delta_V_levels, 'ShowText', 'on');
contour(t.arrive,t.depart,t.delta/(60*60*24),'color', 'k', 'ShowText', 'on');
legend({'\DeltaV (m/s)','TOF (days)'}, 'Interpreter', 'TeX')
title('\DeltaV (m/s)','Interpreter','TeX')
xlabel('Date of arrival')
ylabel('Date of departure')
dateaxis('x',17, datetime(t.arrive_i,'ConvertFrom','juliandate'))
dateaxis('y',17, datetime(t.depart_i,'ConvertFrom','juliandate'))

%   Finding minimum thrust required
min_T_a.val = min(T_a_max(:));
[min_T_a.i,min_T_a.j] = find(T_a_max==min_T_a.val);
min_T_a.depart = t.depart(min_T_a.i);
min_T_a.arrive = t.arrive(min_T_a.j);
%   Converting to datetime
min_T_a.depart = datetime(min_T_a.depart, 'ConvertFrom', 'juliandate');
min_T_a.arrive = datetime(min_T_a.arrive, 'ConvertFrom', 'juliandate');
disp(['Minimum thrust acceleration: ',num2str(min_T_a.val),'m/s^2'])
disp(['Delta V required: ',num2str(abs(delta_V(min_T_a.i,min_T_a.j))),'m/s'])
disp(['Departure date: ', datestr(min_T_a.depart)])
disp(['Arrival date: ', datestr(min_T_a.arrive)])
toc