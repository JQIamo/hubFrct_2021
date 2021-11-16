function [PHMODE, AMP] = RingDy_phnEv_2D...
    (t,phi,rfunc,ANfunc,nlobes,ampi,cTheta0,Qi,Qf,t_phi0,ang_phi0,varargin)
% =========================================================================
% Function to generate a 2D sinusoidal pattern of ring phonon evolution 
% n1D_th  with time and angular position. This takes care of changing ring radii
% Inputs:
%   t and phi: time t[s] and angular position phi[rad].  
%   rfunc: function handle to determine the radius as a function of time.
%       Should be func of time only.
%   ANfunc: function handle to determine the atom number as a function of time.
%       Should be func of time only.
%   nlobes: The mode# of phonon
%   ampi: amplitude of phonon ossc.
%   ctheta0: Speed of sound [um/s] for 100k atoms in a 40 um radius ring
%   Qi: quality factor in initial ring
%   Qf: quality factor in final ring
%   t_phi0: Temporal Phase [rad] at t = 0 for the time sinusoid
%   ang_phi0: Phase[rad] for the angular sinusoid
%   vector: set to true to get a column vector [logical]
%   gamma, hub: gamma and hub friction coefficients

% Outputs:
% 	PHMODE: delta n_{1D} as a function of t and phi. 
%   AMP: delta n as a function of t.
% =========================================================================
   
    %% Assemble Optional Inputs =====================================================
    p = inputParser;
    p.addParameter('Gamma',1/2,@(x)isnumeric(x));
    p.addParameter('Hub',1/2,@(x)isnumeric(x));
    p.addParameter('PersCurr',false,@(x)islogical(x));
    p.addParameter('dPhidt',0,@(x)isnumeric(x));
    p.addParameter('vector',false,@(x)islogical(x));
    p.parse(varargin{:});    
    gamma = p.Results.Gamma;
    hub = p.Results.Hub;
    vector = p.Results.vector;
    
    if size(t)~=size(phi)
        error('RingDy_phnEv_2D: The inputs t and phi need to be of the same size.');
    end
    if (size(t,1)~=1 && size(t,2)~=1)
        error('RingDy_phnEv_2D: The inputs t and phi need to be 1D vector');
    end
        
    %% Evaluate Phonon ============================================================== 
    AMP = PhnAmpDynamicDiffEq(t,rfunc,ANfunc,cTheta0,Qi,Qf,ampi,t_phi0,gamma,hub);
    
    if p.Results.PersCurr
        ang_phi = ang_phi0 + t*p.Results.dPhidt;
        PHMODE = AMP.*sin(nlobes.*phi+ang_phi);
    else
        PHMODE = AMP.*sin(nlobes.*phi+ang_phi0);
    end
    if vector
        PHMODE = PHMODE(:);
        AMP = AMP(:);
    end

end

