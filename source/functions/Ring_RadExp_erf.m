function [radius,velocity,t_peak] = Ring_RadExp_erf(t,ri,rf,t0,sigma)
% ===================================================================================
% [radius,velocity,t_peak] = Radius_erf_3(t,ri,rf,t0,sigma)
% Function to calculate radius as a function of time expanding as erf function
% Inputs:
% - t: time
% - ri, rf: Initial, final radii
% - t0: time at which r begins to change as erf(2*t/sigma)
% Outpus:
% - radius and velocity as function of time
% - t_peak: time of maximum |V/R|
% 2020/10/16: Written by Monica Gutierrez
% ===================================================================================

% Change of variable
tc = t0 + sigma;
beta = sigma/2.5;
x = (t-tc)/beta;

% Calculate r(t)
radius = ri + 0.5*(rf-ri)*(1+erf(x));
radius(t<t0) = ri;
radius(t>(t0+2*sigma)) = rf;

% Calculate v(t) = dR/dt
velocity = (rf-ri)/(sqrt(pi)*beta)*exp(-x.^2);
velocity(t<t0) = 0;
velocity(t>(t0+2*sigma)) = 0;

% Find time that maximizes S = abs(V/R)
S = abs(velocity./radius);
[~,idx] = max(S);
t_peak = t(idx);

end

