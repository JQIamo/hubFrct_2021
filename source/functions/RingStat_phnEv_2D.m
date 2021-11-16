function [PHMODE, AMP] = RingStat_phnEv_2D...
    (t,phi,r,nlobes,amp,c,Q,t_phi0,ang_phi0,varargin)
% =========================================================================
% Function to generate a 2D sinusoidal pattern of ring phonon evolution 
% n1D_th  with time and angular position. This takes care of changing ring radii
% Inputs:
%   t and phi: time t[s] and angular position phi[rad]. 
%   r: radius [um]
%   nlobes: The mode# of phonon
%   amp: amplitude of phonon ossc.
%   c: speed of sound [um/s]
%   Q: quality factor 
%   t_phi0: Temporal Phase [rad] at t = 0 for the time sinusoid
%   ang_phi0: Phase[rad] for the angular sinusoid
%   vector: set to true to get a column vector [logical]

% Outputs:
% 	PHMODE: delta n_{1D} as a function of t and phi. 
%   AMP: delta n as a function of t.
% =========================================================================
   
    %% Assemble Optional Inputs =====================================================
    p = inputParser;
    p.addParameter('vector',false,@(x)islogical(x));
    p.addParameter('PersCurr',false,@(x)islogical(x));
    p.addParameter('dPhidt',0,@(x)isnumeric(x));
    p.parse(varargin{:});    
    vector = p.Results.vector;   
    if size(t)~=size(phi)
        error('RingStat_phnEv_2D: The inputs t and phi need to be of the same size.');
    end
    
    %% Evaluate Phonon ==============================================================
    omega = c/r;
    AMP = amp*exp(-omega*t/Q/2).*sin(omega*t+t_phi0);
    if p.Results.PersCurr
        ang_phi = ang_phi0 + t*p.Results.dPhidt;
        PHMODE = AMP.*sin(nlobes.*phi+ang_phi);
    else
        PHMODE = AMP.*sin(nlobes.*phi+ang_phi0); 
    end
    if vector
        PHMODE = PHMODE(:);
    end

end

