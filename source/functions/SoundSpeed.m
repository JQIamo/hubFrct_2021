function cTheta = SoundSpeed(cTHETA0,rad,number,alpha)
% Function to determine the speed of sound in ring 
% Inputs:
%   cTHETA0: speed @ N_atoms = 100 k and R = 40 um
%   rad: Radius of ring [um]
%   number: number of atoms [k]
%   alpha: geometric factor of the trap

% Outputs:
% 	cTheta: speed [um/s] @ N_atoms = 'number' and R = 'rad' um
% =========================================================================
    cTheta = cTHETA0*(rad/40).^(-alpha/2).*(number/100).^(alpha/2);
end