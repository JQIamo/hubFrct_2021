% Function to plot Static Phonon Data
% Inputs:
%   fignum: figure number
%   N_timeSlice: Array of time slices for each Dy data-set
%   (X_data,data1D_avg,data1D_err): 1D array of data to plot
%   (X_fit,Fit): 1D array of fit to plot
%   Ri, Rf, t_phn1, cTheta0, AtomNum, alpha: details to determine phase,
%       speed etc.
% Output:
%   The plotted figure.

function [] = PlotPhnStat(fignum,N_timeSlice,X_data,data1D_avg,...
        data1D_err,X_fit,Fit,Rad,t_phn1,t_exp,cTHETA0,AtomNum,alpha,Q_val,Q_err,varargin)
    %% Assemble variable inputs and check ==================================
    p = inputParser;
    p.addParameter('plotErrors1D',true,@(x)(islogical(x)));
    p.addParameter('xlimits',[0 250],@(x)(isnumeric(x) && length(x)==2));
    p.addParameter('ylimits',[-12 12],@(x)(isnumeric(x) && length(x)==2));
    p.parse(varargin{:});    
    
    if (size(N_timeSlice) ~= size(Rad) | size(Rad) ~= size(t_phn1)...
            | size(t_phn1) ~= size(t_exp))
        error('PlotPhnStat: R_i,t_phn1,t_exp and N_timeSlice should be the same size.');
    end
    if (size(X_data) ~= size(data1D_avg) | size(data1D_avg) ~= size(data1D_err)...
            | length(data1D_err) ~= sum(N_timeSlice))
        error('PlotPhnStat: Mismatch in deltan0 and X_data sizes.');
    end
    if (size(X_fit) ~= size(Fit))
        error('PlotPhnStat: Mismatch in Fit and X_fit sizes.');
    end
        
    fig = figure(fignum); clf
    %% Loop through and plot each static ring =============================
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
        frq = SoundSpeed(cTHETA0,Rad(ii),AtomNum(ii),alpha)./Rad(ii)/2/pi;
               
        subplot(length(N_timeSlice),1,ii)  
        if p.Results.plotErrors1D
            errorbar(1e3*X_data(idx_start_data:idx_end_data),1e-3*data1D_avg(idx_start_data:idx_end_data),...
            1e-3*data1D_err(idx_start_data:idx_end_data),'.','MarkerSize',10); hold on;
        else
            plot(1e3*X_data(idx_start_data:idx_end_data),1e-3*data1D_avg(idx_start_data:idx_end_data),...
                '.','MarkerSize',30); hold on;
        end
        plot(1e3*X_fit(idx_start_fit:idx_end_fit),1e-3*Fit(idx_start_fit:idx_end_fit),...
            '-r','LineWidth',2); 
        xlim(p.Results.xlimits); ylim(p.Results.ylimits); 

        if ii~=length(N_timeSlice)
            ax = gca;
            set(ax, 'XTickLabels', []);
        end
        box on; grid on;
%         txt = [sprintf('$R$ = %0.2f $\\mu m$, $N_{atoms}$ = %0.1f $k$',Rad(ii),AtomNum(ii)) newline ...
%             sprintf('$\\Omega/2\\pi$ = %0.1f Hz, $Q$ = %0.2f $\\pm$ %0.2f',frq,Q_val(ii),Q_err(ii))];
%         text(abs(diff(p.Results.xlimits))*0.5,...
%             abs(diff(p.Results.ylimits))*(-0.3),txt,...
%             'fontweight','bold','FontSize',12) ;
    end
    ax = axes(fig);
    yyaxis(ax, 'left');
    han = gca;
    han.Visible = 'off';
    han.XLabel.Visible = 'on';
    han.YLabel.Visible = 'on';
    ylabel(sprintf('$\\delta n_0$ ($10^3$ atoms/rad)'));
    xlabel(sprintf('$t$ (ms)'));  
end