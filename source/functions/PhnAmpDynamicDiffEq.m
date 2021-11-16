function [phmodeamp,r_temp,v_temp,Omegas,Gammas,cThetas,Qs] = ...
    PhnAmpDynamicDiffEq(tfit,rfunc,ANfunc,cTheta0,Qi,Qf,amp,ph,alpha,hub)
% =========================================================================
% This function returns the value of 
% delta n_0 = dchi/dt * R^alpha 
% which staisfies the DE  
% d^2chi/dt^2 + (2 g + hub*v/r) dchi/dt + (m*cTheta/R)^2 chi = 0
% Inputs:
%   tfit: time[s]
%   rfunc: function handle to determine the radius as a function of time.
%       Should be func of time only.
%   ANfunc: function handle to determine the atom number as a function of time.
%       Should be func of time only.
%   ctheta0: Speed of sound [um/s] for 100k atoms in a 40 um radius ring
%   Qi, Qf: Quality factors in the initial and final rings
%   amp, ph: Amplitude and temporal phase of delta n
%   alpha, hub: alpha and hub friction co-efficients

% Outputs:
%   phmodeamp: delta n_0 = dchi/dt * R^alpha
% =========================================================================

    %% Find Q as a function of radius ===============================================
    % This assumes a linear dependence of Q on R. If you wish to have
    % constant Q, set variables 'Qi' and 'Qf' same
    fm_Q = fit([rfunc(0) rfunc(tfit(end))]',[Qi Qf]','poly1');

    %% Initialize ===================================================================
    % Variables 'amp' and 'phase' are Amplitude and phase of the initial phonon 
    % 'delta n_0'. In this section we find the initial values of 
    % dchi/dt = delta n_0 / R^alpha. 
    % Calculate the initial amplitudes:
    ci = SoundSpeed(cTheta0,rfunc(0),ANfunc(0),alpha);
    gi = ci./rfunc(0)/2/Qi;
    fi = ci./rfunc(0)*(1-1/(4*Qi^2))^(0.5)/2/pi;
    deltan_i = @(t)amp*sin(2*pi*fi*t+ph).*exp(-t.*gi);
    A1 = (amp/(2*1i))*exp(1i*ph)/(1i*2*pi*fi-gi).*rfunc(0).^(-alpha);
    A2 = -(amp/(2*1i))*exp(-1i*ph)/(-1i*2*pi*fi-gi).*rfunc(0).^(-alpha);
    chi_i = @(t)A1*exp(1i*2*pi*fi*t-gi*t) + A2*exp(-1i*2*pi*fi*t-gi*t);
    y0 = chi_i(0);
    y0p = deltan_i(0)./rfunc(0)^alpha;
    clear A1 A2 deltan_i chi_i

    %% Solve the differential equation ====================================
    % First define the terms as a function of time ========================
    cThetaFunc = @(t)SoundSpeed(cTheta0,rfunc(t),ANfunc(t),alpha);         % ci as a function of t
    QFunc = @(t)fm_Q.p1*rfunc(t)+fm_Q.p2;                                  % Q as a function of t
    OmegaFunc = @(t)cThetaFunc(t)./rfunc(t).*(1-1./(4*QFunc(t).^2))^(0.5); % The resultant Omega as a function of t
    gammaFunc = @(t)cThetaFunc(t)./rfunc(t)./QFunc(t)/2;                   % The resultant gamma as a function of r
    % Now solve the DE ====================================================
    sol = ode45(@(t,y)diffEq(t,y,rfunc,hub,...
            gammaFunc,cThetaFunc),[0 1.01*max(tfit)],[y0 y0p]);
    phmodeamp = deval(sol,tfit);
    phmodeamp = phmodeamp(2,:);
    if size(tfit,2)==1
        phmodeamp = phmodeamp';
    end
    phmodeamp = phmodeamp.*rfunc(tfit).^(alpha);

    if nargout>1
        [r_temp,v_temp] = rfunc(tfit);
        Omegas = OmegaFunc(r_temp);
        Gammas = gammaFunc(r_temp);
        cThetas = cThetaFunc(r_temp);
        Qs =  QFunc(r_temp);
    end

end

%% Local Functions ========================================================
% The DE function for =====================================================
function [dydt] = diffEq(t,y,rfunc,hubFrVal,gamma_func,cTheta_func)
    % Inputs: 
    %   t: time [time_units]
    %   y: a 2 element array. y(1) = chi, y(2) = dchi/dt
    %   rfunc: Function handle to radius as a function of time 't'. This
    %       function should have two outputs. 1st output should be 
    %       radius [rad_units] and second output should be rate of change of 
    %       radius [rad_units/time_units].
    %   hubFrVal: Co-efficient of hub friction. Should be positive
    %   gamma_func: Function handle to decay gamma as a function of time
    %       't'. 
    %   cTheta_func: Function handle to cTheta [rad_units/time_units] as 
    %       a function of radius. radius is defined by function rfunc.
    % Outputs:
    %   dydt: A element array. dydt(1)= dchi/dt, dydt(2) = d^2chi/dt^2 
    [r,v] = rfunc(t);
    dydt = zeros(2,1);
    dydt(1) = y(2);                                             % dchi/dt
    dydt(2) = -(2*gamma_func(t)+hubFrVal*v/r)*y(2) - ...
        (cTheta_func(t)./r).^2*y(1);                            % d^2chi/dt^2

end % =====================================================================
