function [N] = atomNumber(t,t_start,tau,Ni,Nf)
% ===================================================================================
% Function to determine atom number as a function of time
% Inputs:
% - t: time
% - t_start: start of decay [s]
% - tau: decay constant [s]
% - Ni and Nf: Initial and final atom numbers [k]
% Outputs:
% - N: Atom number[k]
% ===================================================================================
    N = (Ni-Nf)*exp(-(t-t_start)/tau)+Nf;
    N (t<=t_start) = Ni;
end

