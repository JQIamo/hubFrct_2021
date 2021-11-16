function phmodeamp = phononAmpStaticDiffEq(tfit,amp,ph,gamma,f,c)
% =========================================================================
% This function returns the value of dchi/dt which staisfies the DE  
%   d^2chi/dt^2 + 2 gamma dchi/dt + omega^2 chi = 0
% and adds 'c' to it.
% Inputs:
%   tfit: time[s] values for which values are needed
%   amp,ph: 
%   gamma: co-efficient of the DE
%   f: frq[Hz] used to determine omega 
%   c: the constant

% Outputs:
%   phmodeamp: dchi/dt + c
% =========================================================================

    omega = 2*pi*f;
    y0 = amp/omega*sin(ph);
    y0p = amp*cos(ph);
    sol = ode45(@(t,y)diffEq(t,y,2*gamma,omega^2),[0 max(tfit)],[y0 y0p]);
    phmodeamp = deval(sol,tfit); % The solution chi and its first derivative dchi/dt
    phmodeamp = phmodeamp(2,:)+c; % dchi/dt + c
    if size(tfit,2)==1
        phmodeamp = phmodeamp';
    end

end

%% Local Functions

function [dydt] = diffEq(~,y,dydtcoeff,ycoeff)

dydt = zeros(2,1);

dydt(1) = y(2);
dydt(2) = - dydtcoeff*y(2) - ycoeff*y(1);

end