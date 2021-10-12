%%  Function for calculating shape-based approach to low thrust trajectory problem
%   BIR 01493976
%   Inputs:
%   r1: position vector at point 1 (ijk) (m)
%   r2: position vector at point 2 (ijk) (m)
%   v1: velocity vector at point 1 (ijk) (m/s)
%   v2: velcoity vector at point 2 (ijk) (m/s)
%   TOF: time of flight (s)
%   mu: standard gravitational parameter of central body (m^3/s^2)
%   N_rev: number of revolutions around central body
%   q: 3D shaping parameter (default = 7)
%   
%   Outputs:
%   a,b,c,d,e,f,g: coefficients of shaped trajectory
%   DU: 1 distance unit (m)
%   TU: 1 time unit (s)
%   Assumptions:
%   Thrust is in direction of flight path angle (reasonable for optimal
%   solution)
%   Inclination change is small (<15 degrees)
%   Notes:
%   
function [a, b, c, d, e, f, g, a_z, b_z, c_z, d_z, DU, TU, zDU, zTU, error] = shape_trajectory(pos1, pos2, vel1, vel2, TOF, mu, N_rev, q)
%   Setting default q if not given
if nargin ~= 8
    q = 7;
end

error = false;

%   Setting up structures for input
r1.vect = pos1;
r1.norm = norm(r1.vect);
r2.vect = pos2;
r2.norm = norm(r2.vect);
v1.vect = vel1;
v1.norm = norm(v1.vect);
v2.vect = vel2;
v2.norm = norm(v2.vect);

%   Converting to canonical units
DU = r1.norm; % one distance unit set as initial position (m)
TU = sqrt(r1.norm^3/mu); % one time unit (s)

%   Calculating zDU and zTU
zDU = norm(r1.vect(3));
zTU = sqrt(norm(r1.vect(3))^3/mu);

%   Calculating z coords from reference plane 
z1 = r1.vect(3)/zDU;
z2 = r2.vect(3)/zDU;
v_z1 = v1.vect(3)*zTU/zDU;
v_z2 = v2.vect(3)*zTU/zDU;

mu = mu*(TU^2/DU^3);
TOF = TOF/TU;

r1.vect = r1.vect/DU;
r1.norm = r1.norm/DU;
r2.vect = r2.vect/DU;
r2.norm = r2.norm/DU;
v1.vect = v1.vect*TU/DU;
v1.norm = v1.norm*TU/DU;
v2.vect = v2.vect*TU/DU;
v2.norm = v2.norm*TU/DU;

%   Calculating flight path angles
gamma1 = asin(dot(r1.vect,v1.vect)/(r1.norm*v1.norm)); % calculate from complement of flight direction angle pg 128 Battin
gamma2 = asin(dot(r2.vect,v2.vect)/(r2.norm*v2.norm));

%   Calculating initial and final rate of change of true anomaly derivatives
theta1.dot = norm(cross(r1.vect,v1.vect))/r1.norm^2; % Using angular momentum argument in Battin
theta2.dot = norm(cross(r2.vect,v2.vect))/r2.norm^2;

%   Calculating required theta values
theta.diff = acos(dot(r1.vect,r2.vect)/(r1.norm*r2.norm)); % angle between desired position vectors
n = cross(r1.vect,r2.vect)/((r1.norm*r2.norm)*sin(theta.diff));
%   Checking direction of rotation
if n(3) < 0 % negative z implies 2pi > theta.diff > pi
    theta.diff = 2*pi - theta.diff;
end
theta.f = theta.diff + 2*pi*N_rev; % angle at r2

%   Determining coefficients for method
a = 1/r1.norm;
b = -tan(gamma1)/r1.norm;
c = 1/(2*r1.norm)*(mu/(r1.norm^3*theta1.dot^2)-1);

%   Setting up fsolve to find coefficient d
fun = @diff_TOF;
[d, ~, exitflag] = fsolve(fun, 1);
%   fsolve error catching
if exitflag <= 0
    error = true;
end

[e, f, g] = calc_coeffs(d);

%   Calculating z(theta) shaping (Wall)
a_z = z1; % z
b_z = v_z1/theta1.dot;
c_z = (1/theta.f^q)*(q*theta.f*(z2-a_z-b_z*theta.f)-theta.f^2*(v_z2/theta2.dot-b_z));
d_z = (1/theta.f^q)*((1-q)*(z2-a_z-b_z*theta.f)-theta.f*(v_z2/theta2.dot-b_z));

%   Sub-functions
%   Function for calculating value of other coefficients given a value of d
function [e, f, g] = calc_coeffs(d_given)
A = [
    30*theta.f^2,   -10*theta.f^3,  theta.f^4;
    -48*theta.f,    18*theta.f^2,   -2*theta.f^3;
    20,             -8*theta.f,     theta.f^2];
B = [
    1/r2.norm-(a+b*theta.f+c*theta.f^2+d_given*theta.f^3);
    -tan(gamma2)/r2.norm-(b+2*c*theta.f+3*d_given*theta.f^2);
    mu/(r2.norm^4*theta2.dot^2)-(1/r2.norm+2*c+6*d_given*theta.f)];
    
e = 1/(2*theta.f^6)*A(1,:)*B;
f = 1/(2*theta.f^6)*A(2,:)*B;
g = 1/(2*theta.f^6)*A(3,:)*B;
end

%   Parameterised function to pass to fzero to find d which gives a TOF
%   equivilent to the one specified
function [delta_TOF] = diff_TOF(d_guess)
[e_d,f_d,g_d] = calc_coeffs(d_guess);
TOF_d = calc_TOF(a, b, c, d_guess, e_d, f_d, g_d, mu, double(theta.f));
delta_TOF = TOF_d - TOF;
end

end