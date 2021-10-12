%%  Function for calculating TOF of shaped trajectory from coefficients
%   Inputs: 
%   a,b,c,d,e,f,g: shape coefficients
%   mu: standard gravitational parameter
%   theta_f: change in angle
%   Outputs:
%   TOF: time of flight

function [TOF] = calc_TOF(a, b, c, d, e, f, g, mu, theta_f)
%   Setting up anonymous functions
r = @(theta) 1./(a+b*theta+c*theta.^2+d*theta.^3+e*theta.^4+f*theta.^5+g*theta.^6);
T = @(theta) sqrt(r(theta).^4./mu.*(1./r(theta)+2*c+6*d*theta+12*e*theta.^2+20*f*theta.^3+30*g*theta.^4));
TOF = integral(T,0,theta_f);
end