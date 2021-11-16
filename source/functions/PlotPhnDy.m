% Function to plot Dynamic Phonon Data
% Inputs:
%   fignum: figure number
%   N_timeSlice: Array of time slices for each Dy data-set
%   (X_data,data1D_avg,data1D_err): 1D array of data to plot
%   (X_fit,Fit): 1D array of fit to plot
%   Ri, Rf, t_phn1, cTheta0, AtomNum, alpha: details to determine phase,
%       speed etc.
% Output:
%   The plotted figure.

function [] = PlotPhnDy(fignum,N_timeSlice,X_data,data1D_avg,...
        data1D_err,X_fit,Fit,R_i,R_f,t_phn1,t_exp,cTHETA0,...
        ANfit_tau,ANfit_Nis,ANfit_Nfs,...
        alpha,varargin)
    %% Assemble variable inputs and check ==================================
    p = inputParser;
    p.addParameter('plotXtra',[],@(x)(isnumeric(x) && length(x)==length(X_fit)));
    p.addParameter('plotErrors1D',true,@(x)(islogical(x)));
    p.addParameter('xlimits',[0 250],@(x)(isnumeric(x) && length(x)==2));
    p.addParameter('ylimitsLeft',[-12 12],@(x)(isnumeric(x) && length(x)==2));
    p.addParameter('ylimitsRight',[0 50],@(x)(isnumeric(x) && length(x)==2));
    p.addParameter('RdotByR_min',0,@(x)isnumeric(x));
    p.addParameter('RdotByR_max',400,@(x)isnumeric(x));
    p.parse(varargin{:});    
    RdotByR_min = p.Results.RdotByR_min;
    RdotByR_max = p.Results.RdotByR_max;
    if (size(N_timeSlice) ~= size(R_i) | size(R_i) ~= size(R_f) | size(R_f) ~= size(t_phn1)...
            | size(t_phn1) ~= size(t_exp))
        error('PlotPhnDy: R_i,R_f,t_phn1,t_exp and N_timeSlice should be the same size.');
    end
    if (size(X_data) ~= size(data1D_avg) | size(data1D_avg) ~= size(data1D_err)...
            | length(data1D_err) ~= sum(N_timeSlice))
        error('PlotPhnDy: Mismatch in deltan0 and X_data sizes.');
    end
    if (size(X_fit) ~= size(Fit))
        error('PlotPhnDy: Mismatch in Fit and X_fit sizes.');
    end
    %%  Make the Speed of exp/ contraction bar at top =====================
    fig = figure(fignum); clf
    subplot(2*length(N_timeSlice)+1,1,1), cla
    axis off
    colormap(flipud(gray));
    caxis([RdotByR_min RdotByR_max]);
    h = colorbar('location','north');
    title(h,'$|\dot{R}/R|$ (Hz)','Interpreter','latex');
    %% loop through and plot all data sets.
    for ii = 1:1:length(N_timeSlice)
        if ii == 1 
            idx_start_data = 1;
            idx_start_fit = 1;
        else 
            idx_start_data = sum(N_timeSlice(1:ii-1))+1; 
            idx_start_fit = (ii-1)*length(Fit)/length(N_timeSlice)+1;
        end
        idx_end_data = sum(N_timeSlice(1:ii)); 
        idx_end_fit = ii*length(Fit)/length(N_timeSlice);
        [Rad_plot, vel, t_peak] = Ring_RadExp_erf(X_fit(idx_start_fit:idx_end_fit),R_i(ii),...
            R_f(ii),t_phn1(ii),t_exp(ii)/2);
        RdotByR = abs(vel./Rad_plot);
        AN = atomNumber(X_fit(idx_start_fit:idx_end_fit),t_phn1(ii),ANfit_tau(ii),ANfit_Nis(ii),ANfit_Nfs(ii));
        omega = SoundSpeed(cTHETA0,Rad_plot,AN,alpha)./Rad_plot;
        phase = sum(omega(X_fit(idx_start_fit:idx_end_fit)<=t_peak))...
            *mean(diff(X_fit(idx_start_fit:idx_end_fit)))/pi; 
        if RdotByR_max == 0
            temp = RdotByR/abs(RdotByR_min);
        elseif RdotByR_min == 0
            temp = RdotByR/RdotByR_max;
        end
        bkg = repmat(temp',1000,1); clear temp
               
        subplot(2*length(N_timeSlice)+1,1,[2*(ii-1)+2 2*(ii-1)+3])  
        yyaxis left
        set(gca, 'YDir','reverse');
        imagesc(1e3*X_fit(idx_start_fit:idx_end_fit),...
            linspace(p.Results.ylimitsLeft(1),p.Results.ylimitsLeft(2),1000),...
            bkg,[0 1]),hold on;  
        if p.Results.plotErrors1D
            errorbar(1e3*X_data(idx_start_data:idx_end_data),1e-3*data1D_avg(idx_start_data:idx_end_data),...
            1e-3*data1D_err(idx_start_data:idx_end_data),'.','MarkerSize',10); hold on;
        else
            plot(1e3*X_data(idx_start_data:idx_end_data),1e-3*data1D_avg(idx_start_data:idx_end_data),...
                '.','MarkerSize',30); hold on;
        end
        plot(1e3*X_fit(idx_start_fit:idx_end_fit),1e-3*Fit(idx_start_fit:idx_end_fit),...
            '-r','LineWidth',2); 
        if ~isempty(p.Results.plotXtra)
            hold on;
            plot(1e3*X_fit(idx_start_fit:idx_end_fit),1e-3*p.Results.plotXtra(idx_start_fit:idx_end_fit),...
            '--k','LineWidth',2); 
        end
        hold off;

        xlim(p.Results.xlimits); ylim(p.Results.ylimitsLeft); 
        yyaxis right
        ax = gca;
        ax.YColor = 'k'; 
        plot(1e3*X_fit(idx_start_fit:idx_end_fit),Rad_plot,'-k','LineWidth',2);
        ylim(p.Results.ylimitsRight);
        if ii~=length(N_timeSlice)
            set(ax, 'XTickLabels', []);
        end
        box on; grid on;
        txt = [sprintf('$t_{peak}$ = %0.2f ms',t_peak*1e3) newline ...
            sprintf('$\\phi$($t_{peak})-\\phi (0)$ = %0.1f $\\pi$',phase)];
        text(abs(diff(p.Results.xlimits))*0.7,...
            abs(diff(p.Results.ylimitsRight))*0.25,txt,...
            'fontweight','bold','FontSize',12) ;
    end
    ax = axes(fig);
    yyaxis(ax, 'left');
    han = gca;
    han.Visible = 'off';
    han.XLabel.Visible = 'on';
    han.YLabel.Visible = 'on';
    ylabel(sprintf('$\\delta n_0$ ($10^3$ atoms/rad)'));
    xlabel(sprintf('$t$ (ms)'));
    yyaxis(ax, 'right');
    ylabel(sprintf('$R$ ($\\mu m$)'));
    han.YLabel.Visible = 'on';
    ax.YAxis(2).Color = 'k'; 
    
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [N] = atomNumber(t,t_start,tau,Ni,Nf)
    N = (Ni-Nf)*exp(-(t-t_start)/tau)+Nf;
    N (t<=t_start) = Ni;
end