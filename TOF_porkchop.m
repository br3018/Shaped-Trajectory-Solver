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
t.depart_i = juliandate(2060,01,01);
t.depart_f = juliandate(2065,01,01);
t.depart_step = 5;
TOF.max = 600; % (days) (setting range of TOF over which to iterate)
TOF.min = 250; % (days)
TOF.step = 10; % (days)
N_rev = 0;

t.depart = [t.depart_i:t.depart_step:t.depart_f];
TOF.array = [TOF.min:TOF.step:TOF.max];

t.delta = zeros(length(t.depart),length(TOF.array));
delta_V = zeros(length(t.depart),length(TOF.array));
T_a_max = zeros(length(t.depart),length(TOF.array));
error = strings(length(t.depart),length(TOF.array));
theta_f = zeros(length(t.depart),length(TOF.array));

%   Beginning loop of trajectory finding and characterising
iter = 0;
tot_iter = length(t.depart)*length(TOF.array); % total number of calculations
for i = 1:length(t.depart)
    for j = 1:length(TOF.array)
        iter = iter+1; 
        error(i,j) = "ok";
        disp(['Iteration ',num2str(iter),'/',num2str(tot_iter)]) % progress check
        disp(['Estimated time to complete: ', num2str(((tot_iter-iter)/2)/(60*60)),' hours'])
        toc
        t.delta(i,j) = TOF.array(j)*24*60*60; % (s)
        
        %   Calculating planet state vectors at specified times
        [P1.rd,P1.vd]=planetEphemeris(t.depart(i),P0.name,P1.name);
        [P2.rd,P2.vd]=planetEphemeris(t.depart(i),P0.name,P2.name);
        P1.rd=P1.rd.*1000; % convert to m
        P1.vd=P1.vd.*1000;
        P2.rd=P2.rd.*1000; % convert to m
        P2.vd=P2.vd.*1000;


        [P1.ra,P1.va]=planetEphemeris(t.depart(i)+TOF.array(j),P0.name,P1.name);
        [P2.ra,P2.va]=planetEphemeris(t.depart(i)+TOF.array(j),P0.name,P2.name);
        P1.ra=P1.ra.*1000; % convert to m
        P1.va=P1.va.*1000;
        P2.ra=P2.ra.*1000; % convert to m
        P2.va=P2.va.*1000;

        theta_f(i,j) = acos(dot(P1.rd,P2.ra)/(norm(P1.rd)*norm(P2.ra))); % angle between desired position vectors for shaped trajectory
        n = cross(P1.rd,P2.ra)/((norm(P1.rd)*norm(P2.ra))*sin(theta_f(i,j)));
        %   Checking direction of rotation (METHOD WILL BREAKDOWN OVER HIGH
        %   INCLINATION CHANGES)
        if n(3) < 0 
            theta_f(i,j) = 2*pi - theta_f(i,j);
        end
        theta_f(i,j) = theta_f(i,j) + N_rev*2*pi;
        
        %   Calculating shaped trajectory coefficients
        
        [a, b, c, d, e, f, g, a_z, b_z, c_z, d_z, DU, TU, zDU, zTU, coeff_error] = shape_trajectory(P1.rd, P2.ra, P1.vd, P2.va, t.delta(i,j), P0.mu, N_rev);
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


%   Saving output data as backup 
save('ME_2060to2065.mat')
%   Converting to datetime array
t.depart = datetime(t.depart,'ConvertFrom','juliandate');
t.depart = datenum(t.depart);

%   Setting contour plot for max acceleration values
%   Setting interesting acceleration levels
T_a_ax.max = 0.001;
T_a_ax.min = 0.0001;
T_a_ax.step = 0.0001;
T_a_levels = [T_a_ax.min:T_a_ax.step:T_a_ax.max];

%   Setting interesting delta V values 
delta_V_ax.max = 6000;
delta_V_ax.min = 5000;
delta_V_ax.step = 100;
delta_V_levels = [delta_V_ax.min:delta_V_ax.step:delta_V_ax.max];

%   Filtering out high thrust values (arbitraritly set at 0.1m/s^2)
i_anom = find(T_a_max > T_a_ax.max);
T_a_max(i_anom) = NaN;
delta_V(i_anom) = NaN;

%   Setting contour plot for thrust accelerations
figure('Renderer', 'painters', 'Position', [10 10 900 600]*1.1)
fig = fig +1;
hold on 
grid on
grid minor
contourf(t.depart,TOF.array,T_a_max'*1000,T_a_levels'*1000);
contour(t.depart,TOF.array,t.delta'/(60*60*24),'color', 'k', 'ShowText', 'on');
legend({'Thrust acceleration (mm/s^2)','TOF (days)'}, 'Interpreter', 'TeX','location','southeast')
colormap(turbo(length(T_a_levels)-1))
caxis([T_a_ax.min T_a_ax.max]*1000)
colorbar
ylabel('TOF (days)')
xlabel('Date of departure')
datetick('x','yyyy')

figure('Renderer', 'painters', 'Position', [10 10 900 600]*1.1)
fig = fig +1;
hold on 
grid on
grid minor
%   Setting contour plot for delta V values
contourf(t.depart,TOF.array,abs(delta_V'),delta_V_levels);
contour(t.depart,TOF.array,t.delta'/(60*60*24),'color', 'k', 'ShowText', 'on');
legend({'\DeltaV (m/s)','TOF (days)'}, 'Interpreter', 'TeX','location','southeast')
colormap(turbo(length(delta_V_levels)-1))
caxis([delta_V_ax.min delta_V_ax.max])
ylabel('TOF (days)')
xlabel('Date of departure')
colorbar
datetick('x','yyyy')

%   Finding minimum thrust required
min_T_a.val = min(T_a_max(:));
[min_T_a.i,min_T_a.j] = find(T_a_max == min_T_a.val);
min_T_a.depart = t.depart(min_T_a.i);
min_T_a.arrive = t.depart(min_T_a.i) + TOF.array(min_T_a.j);
%   Converting to datetime
min_T_a.depart = datetime(min_T_a.depart,'ConvertFrom','datenum');
min_T_a.arrive = datetime(min_T_a.arrive,'ConvertFrom','datenum');
disp(['Minimum thrust acceleration: ',num2str(min_T_a.val),'m/s^2'])
disp(['Delta V required: ',num2str(abs(delta_V(min_T_a.i,min_T_a.j))),'m/s'])
disp(['Departure date: ', datestr(min_T_a.depart)])
disp(['Arrival date: ', datestr(min_T_a.arrive)])
toc