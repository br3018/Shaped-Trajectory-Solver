%%  Function for characterising given shaped orbit trajectory
function [delta_V, T_a_max, error] = characterise_shaped_traj(a, b, c, d, e, f, g, DU, TU, theta_f)
%   Initialising
error = "ok";

%   Checking for method solution:
TOF_check = calc_TOF(a, b, c, d, e, f, g, 1, theta_f)*TU;
if ~isreal(TOF_check)
    % trajectory not found, returning arbitrarily high value of delta_V and
    % max thrust acceleration
    delta_V = NaN; %(m/s)
    T_a_max = NaN; %(m/s^2)
    error = "imaginary TOF error";
    return
else
    %   NOTE UNITS ARE CONVERTED TO CANONICAL UNITS
    %   Calculating shaped trajectory
    theta = [0:0.05:theta_f];
    r = 1./(a+b*theta+c*theta.^2+d*theta.^3+e*theta.^4+f*theta.^5+g*theta.^6);
    theta_dot = sqrt((1./r.^4).*(1./((1./r)+2*c+6*d*theta+12*e*theta.^2+20*f*theta.^3+30*g*theta.^4)));
    
    %   Calculating flight path angle
    gamma = atan(-r.*(b+2*c*theta+3*d*theta.^2+4*e*theta.^3+5*f*theta.^4+6*g*theta.^5));

    %   Calculating thrust acceleration
    T_a = -1./(2*r.^3.*cos(gamma)).*(6*d+24*e*theta+60*f*theta.^2+120*g*theta.^3-tan(gamma)./r)./(1./r+2*c+6*d*theta+12*e*theta.^2+20*f*theta.^3+30*g*theta.^4).^2; % (DU/TU^2)
  
    %   Calculating max thrust acceleration required
    T_a_max = max(abs(T_a))*(DU/TU^2); % (m/s^2)
    
    %   Calculating delta V of trip
    delta_V = trapz(theta,T_a./theta_dot)*DU/TU; % (m/s)
    
end
