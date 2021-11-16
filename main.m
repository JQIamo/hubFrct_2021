% Code to generate the SM tables ==========================================
% ©2021 Swarnav Banik <sbanik1@umd.edu>
clc, clear, tic
% ########################## INITIALIZATION ###############################
tmp = matlab.desktop.editor.getActiveFilename;
if  ispc==1, newPath = strrep(tmp,'\main.m',''); end
if ismac==1, newPath = strrep(tmp,'/main.m',''); end
restoredefaultpath;
addpath(genpath(newPath))
StartupHubble

%% ########################################################################
% ########################## INPUTS #######################################
% #########################################################################
sets = 'expansions'; % choose from 'expansions' or 'contractions'
% #########################################################################
%% Fitting ================================================================
methods = {'one','two','three','four'};
fit_val = zeros(8,6);
fit_err = zeros(8,6);
nu_fit = zeros(8,1);
% Methods I - IV
for ii = 1:length(methods)
    [nu_fit(ii,1), B, cThetai, deltani, cThetai_err, deltani_err] = ...
        globalFits('fitSets',sets,'fitMethod',methods{ii},'fitAtomNum',false,...
        'saveData',false);
    fit_val(ii,:) = [B.Qi(1) B.Qf(1) B.alpha B.gamma_H cThetai, deltani];
    fit_err(ii,:) = [B.Qi_err(1) B.Qf_err(1) B.alpha_err B.gamma_H_err ...
        cThetai_err, deltani_err];  
end
% Methods V - VIII
for ii = 1:length(methods)
    [nu_fit(ii+4,1), B, cThetai, deltani, cThetai_err, deltani_err] = ...
        globalFits('fitSets',sets,'fitMethod',methods{ii},'fitAtomNum',true,...
        'saveData',false);
    fit_val(ii+4,:) = [B.Qi(1) B.Qf(1) B.alpha B.gamma_H cThetai, deltani];
    fit_err(ii+4,:) = [B.Qi_err(1) B.Qf_err(1) B.alpha_err B.gamma_H_err ...
        cThetai_err, deltani_err];
end
%% Display Table ==========================================================
T = cell(length(nu_fit),7);
for ii = 1:length(nu_fit)
    T(ii,:) = {nu_fit(ii), sprintf('%.1f(%d)',fit_val(ii,1),ceil(fit_err(ii,1)*10)), ...
        sprintf('%.1f(%d)',fit_val(ii,2),ceil(fit_err(ii,2)*10)),...
        sprintf('%.2f(%d)',fit_val(ii,3),ceil(fit_err(ii,3)*100)),...
        sprintf('%.2f(%d)',fit_val(ii,4),ceil(fit_err(ii,4)*100)),...
        sprintf('%.2f(%d)',fit_val(ii,5),ceil(fit_err(ii,5)*100)),...
        sprintf('%.2f(%d)',fit_val(ii,6),ceil(fit_err(ii,6)*100))};
end
Tab = cell2table(T,...
    'VariableNames',{'nu','Qi','Qf','alpha','gamma_H','cThetai','deltani'},...
    'RowName',{'I','II','III','IV','V','VI','VII','VIII'}); 
clc
disp(Tab);
toc
