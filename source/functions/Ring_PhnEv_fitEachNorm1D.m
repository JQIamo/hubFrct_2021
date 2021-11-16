function [a, da, phi, dphi] = Ring_PhnEv_fitEachNorm1D(th_ax,norm_n1ds_th,nlobes,ph,...
    varargin)
% =========================================================================
% UTILITY: This code individually fit n1D_th to a 1D sinusoidal function to map the  
% phonon mode to angular position. 
% Inputs:
%   th_ax: angle[rad]
%   norm_n1ds_th: Normalized n1D
%   nlobes: mode generated
%   ph: phase[rad] of angular sinusoid 
%   YLims: Limits of plot
%   FigNum: figure to plot result
%   FiitPhase: fit the phase? [logocal]

% Outputs:
%   a, da: Fit results amplitude and error
%   phi, dphi: Fit results phase and error
% =========================================================================

    p = inputParser;
    p.addParameter('FigNum',[],@(x)isnumeric(x));
    p.addParameter('FitPhase',false,@(x)(islogical(x) || isnumeric(x)));
    p.addParameter('FitResidual',false,@(x)(islogical(x) || isnumeric(x)));
    p.addParameter('YLims',[0.5 1.5],@(x)isnumeric(x));
    p.parse(varargin{:});

    fignum = p.Results.FigNum;
    fit_phase = p.Results.FitPhase;
    fit_res = p.Results.FitResidual;
    ylims = p.Results.YLims;
    
    % Initialize
    a = zeros(size(norm_n1ds_th,2),1);
    da = zeros(size(norm_n1ds_th,2),1);
    phi = a; dphi = da;

    if ~isempty(fignum)
        if fignum<1, figure;
        else, figure(fignum);
        end
        cla;
    end

    for jj=1:size(norm_n1ds_th,2)
        if ~fit_phase
            if fit_res
                [fm_amp_vs_t,gof] = fit(th_ax,norm_n1ds_th(:,jj),...
                    ['1+a*sin(' num2str(nlobes) '*x+' num2str(ph) ')+c'],...
                    'StartPoint',[max(norm_n1ds_th(:,jj)) 1]);
            else
                [fm_amp_vs_t,gof] = fit(th_ax,norm_n1ds_th(:,jj),...
                    ['1+a*sin(' num2str(nlobes) '*x+' num2str(ph) ')'],...
                    'StartPoint',max(norm_n1ds_th(:,jj)));
            end
        else
            if fit_res
                [fm_amp_vs_t,gof] = fit(th_ax,norm_n1ds_th(:,jj),...
                    ['1+a*sin(' num2str(nlobes) '*x+phi)+c'],...
                    'StartPoint',[max(norm_n1ds_th(:,jj)) 1 ph]);
            else
                [fm_amp_vs_t,gof] = fit(th_ax,norm_n1ds_th(:,jj),...
                    ['1+a*sin(' num2str(nlobes) '*x+phi)+0'],...
                    'StartPoint',[max(norm_n1ds_th(:,jj)) ph]);
            end
        end

        if ~isempty(fignum)
            plot(th_ax,norm_n1ds_th(:,jj),'.',...
                linspace(-pi,pi,50),feval(fm_amp_vs_t,linspace(-pi,pi,50)),'-');
            grid on; xlim([-pi pi]); ylim(ylims);
            drawnow; %pause(2);
        end

        a(jj) = fm_amp_vs_t.a;
        err =  diff(confint(fm_amp_vs_t))/4;
        %da(jj) = 0;
        da(jj) = err(1);

        if fit_phase
            phi(jj) = fm_amp_vs_t.phi;
            if fit_res
                dphi(jj) = err(3);
            else
                dphi(jj) = err(2);
            end
        else
            phi(jj) = ph;
            dphi(jj) = 0;
        end
    end

end

