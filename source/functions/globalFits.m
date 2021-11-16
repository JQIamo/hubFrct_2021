function [nu_fit, B, cThetai, deltani, cThetai_err, deltani_err] = globalFits(varargin)
    % =========================================================================
    % #########################################################################
    % #################### DESCRIPTION: FIT METHODS ###########################
    % #########################################################################
    % zero  == plot with guess
    % one   == same ph_t, ph_az, Q, A                       through Dy datasets
    % two   == same ph_az, Q, A      - diff ph_t            through Dy datasets
    % three == same ph_t, Q and A    - diff ph_az           through Dy datasets
    % four  == same Q, A             - diff ph_t and ph_az  through Dy datasets
    % fitAtomNum: Do you want to fit with variable atom number? [logical]
    % #########################################################################
    % Inputs:
    %   fitSets: The data set to be fitted ['expansions' or 'contractions']
    %   fitMethod: fitting method ['one' or 'two' or 'three' or 'four']
    %   fitAtomNum: Do you want to fit with variable atom number? [logical]
    %   saveData: Do you want to save the results? [logical]
    % Outputs:
    %   nu_fit, Qi,Qf, alpha, hub, cThetai, deltani: SM table values
    % =========================================================================
    p = inputParser;
    p.addParameter('fitSets','expansions',@(x)(strcmp(x,'expansions') || ...
        strcmp(x,'contractions')));
    p.addParameter('fitMethod','one',@(x)(strcmp(x,'one') || ...
        strcmp(x,'two') || strcmp(x,'three') || strcmp(x,'four')));
    p.addParameter('fitAtomNum',false,@(x)islogical(x));
    p.addParameter('saveData',true,@(x)islogical(x));
    p.parse(varargin{:});
    
    fitSets = p.Results.fitSets;
    saveData = p.Results.saveData;
    fitAtomNum = p.Results.fitAtomNum;
    fitMethod = p.Results.fitMethod;
    
    %% Save and Plot ==========================================================
    ANbase = 100;               % Base atom number [k] for quoting amplitudes.
    yLimsDeltan = [-12 12];     % Y Limits for delta n plots
    xLimitsDy = [0 190];        % X Limits for delta n plots for Dy datasets
    yLimsAN = [60 140];         % Display limits for AN plot
    N_fitPlot = 500;            % Number of points in each fitted curve
    % Initial values Phonon parameters ========================================
    nlobes = 1;                 % Number of lobes
    tau_guess = 12e-3;          % Atom number decay time [s] guess
    alpha_guess = 0.5;          % Geometric Constant Guess
    hub_guess = 0.5;            % Hub friction co-eff Guess
    cTheta0_guess = 4.33e3;     % Speed of sound for R = 40um, N = 100k [um/s]
    %% Generate file names and sort them acc to t_phn1 ========================
    % Generate list of all phonon mat files ===================================
    currentDir = pwd;
    dataDir = fullfile(pwd,'data',fitSets);
    cd(dataDir );
    file_list_start = {dir('PhnRuns_*.mat').name}';
    % Find all Dy files and sort them acc to t_phn1 ===========================
    kk = 1;
    for jj = 1:1:numel(file_list_start)
        load(file_list_start{jj});
        if strcmp(A.dataRingType,'All')
            file_name_Dy{kk} = file_list_start{jj};
            temp(kk) = A.t_phn1;
            kk = kk + 1;
        end
    end
    [~,I] = sort(temp); 
    file_name_Dy = file_name_Dy(I);
    clear I temp
    % Find all Stat files and sort them acc to mean Rad =======================
    kk = 1;
    for jj = 1:1:numel(file_list_start)
        load(file_list_start{jj});
        if strcmp(getfield(A,'dataRingType'),'Initial')
            file_name_Stat{kk} = file_list_start{jj};
            temp(kk) = mean(A.MeanRad_avg);
            kk = kk + 1;
        end
    end
    [~,I] = sort(temp); 
    file_name_Stat = file_name_Stat(I);
    clear I temp 
    % Now sort the file_name_Stat to have the most relevant at the top ========
    kk = 0;
    temp_rad = 0;
    for jj = 1:1:numel(file_name_Dy)
        load(file_name_Dy{jj});
        temp_rad = temp_rad + mean(A.MeanRad_avg(A.EvTime<A.t_phn1));
        kk = kk+1;
    end
    temp_rad = temp_rad/kk;
    temp_comp = 5;
    stop_idx = 0;
    for jj = 1:1:numel(file_name_Stat)
        load(file_name_Stat{jj});
        if abs(mean(A.MeanRad_avg)-temp_rad)<temp_comp
            stop_idx = jj;
            temp_comp = abs(mean(A.MeanRad_avg)-temp_rad);
        end    
    end
    file_name_Stat = file_name_Stat([stop_idx,1:stop_idx-1, stop_idx+1:numel(file_name_Stat)]);
    clear temp_rad temp_comp
    %% Initialize =============================================================
    file_names = [file_name_Dy, file_name_Stat];
    Rad_i = zeros(length(file_names),1);
    Rad_f = zeros(length(file_names),1);
    t_exp = zeros(length(file_names),1);
    t_phn1 = zeros(length(file_names),1);
    Amp = zeros(length(file_names),1);
    Ph_t = zeros(length(file_names),1);
    Ph_az = zeros(length(file_names),1);
    dPhiaz_dt = zeros(length(file_names),1);
    Qi = zeros(length(file_names),1);
    Qi_err = zeros(length(file_names),1);
    MeanAtomNum = zeros(length(file_names),1);
    N_timeSlices = zeros(length(file_names),1);
    %% Import the data ========================================================
    for ii = 1:1:length(file_names)
        load(strcat(file_names{ii}), 'A');
        Rad_i(ii) = mean(A.MeanRad_avg(A.EvTime<A.t_phn1));
        Rad_f(ii) = mean(A.MeanRad_avg(A.EvTime>A.t_phn1+A.t_exp));
        t_exp(ii) = A.t_exp;
        t_phn1(ii) = A.t_phn1;
        Amp(ii) = A.Ampi;
        Ph_t(ii) = A.TemporalPhasei;
        Ph_az(ii) = A.AzPhase;
        dPhiaz_dt(ii) = A.dAzPhase_dt;
        Qi(ii) = A.Qi;
        Qi_err(ii) = A.Qi_err;
        MeanAtomNum(ii) = mean(A.AN_phnEv_avg);
        N_timeSlices(ii) = length(A.EvTime);    
        if ii == 1
            Phi = A.AzAngle;
        else
            if (A.AzAngle ~= Phi)
                error('main: Anglular positions are not same across datasets.');
            end
        end
        clear A
    end
    n1D_asmbld_avg = zeros(length(Phi),sum(N_timeSlices));
    n1D_asmbld_dev = zeros(length(Phi),sum(N_timeSlices));
    deltan0_asmbld_avg = zeros(sum(N_timeSlices),1);
    deltan0_asmbld_dev = zeros(sum(N_timeSlices),1);
    AN_asmbld_avg = zeros(sum(N_timeSlices),1);
    AN_asmbld_err = zeros(sum(N_timeSlices),1);
    EvTime_asmbld = zeros(sum(N_timeSlices),1);
    DSidx = zeros(sum(N_timeSlices),1);
    for ii = 1:1:length(file_names)
        load(strcat(file_names{ii}), 'A');
        if ii == 1 
            idx_start_data = 1;
        else
            idx_start_data = sum(N_timeSlices(1:ii-1))+1;
        end
        idx_end_data = sum(N_timeSlices(1:ii));
        n1D_asmbld_avg(:,idx_start_data:idx_end_data) = A.n1D_phnEv_avg;
        n1D_asmbld_dev(:,idx_start_data:idx_end_data) = A.n1D_phnEv_dev;
        deltan0_asmbld_avg(idx_start_data:idx_end_data) = A.deltan_avg;
        deltan0_asmbld_dev(idx_start_data:idx_end_data) = A.deltan_dev;
        AN_asmbld_avg(idx_start_data:idx_end_data) = A.AN_phnEv_avg;
        AN_asmbld_err(idx_start_data:idx_end_data) = A.AN_phnEv_dev;
        EvTime_asmbld(idx_start_data:idx_end_data) = A.EvTime;
        DSidx(idx_start_data:idx_end_data) = ii*ones(idx_end_data-idx_start_data+1,1);
        clear A
    end
    Rad_f(isnan(Rad_f)) = [];
    Qf = 7*ones(length(file_name_Dy),1);
    %% Generate the [xData, yData] for curve fitting and fit plotting =========
    xVal = [EvTime_asmbld DSidx];
    xVal_fitPlot = zeros(N_fitPlot*length(file_names),2);
    for ii = 1:1:length(file_names)
        idx_start = N_fitPlot*(ii-1)+1 ;
        idx_end = N_fitPlot*ii;
        xVal_fitPlot(idx_start:idx_end,1) = ...
            linspace(0,max(EvTime_asmbld(sum(N_timeSlices(1:ii-1))+1:sum(N_timeSlices(1:ii)))),...
            N_fitPlot);
        xVal_fitPlot(idx_start:idx_end,2) = ii;
        clear idx_start idx_end
    end

    %% Determine the start of atom number decay ===============================
    t_starts = zeros(length(file_name_Dy),1);
    for ii = 1:1:length(file_name_Dy)
        if ii == 1 
            idx_start_data = 1;
        else
            idx_start_data = sum(N_timeSlices(1:ii-1))+1;
        end
        idx_end_data = sum(N_timeSlices(1:ii));
        [~,~,t_starts(ii)] = Ring_RadExp_erf(xVal(idx_start_data:idx_end_data,1),...
            Rad_i(ii),Rad_f(ii),t_phn1(ii),t_exp(ii,1)/2);
    end
    %% Fit the atom number decay ==============================================
    if fitAtomNum
        fitFunc = @(params,x)atomNumber_all(x,params(1),t_starts(1:length(file_name_Dy)),...
            params(2:1+length(file_name_Dy))',params(2+length(file_name_Dy):1+2*length(file_name_Dy))');
        Initial_guess = [tau_guess,mean(AN_asmbld_avg)*ones(1,2*length(file_name_Dy))];
        lower = [0,0*ones(1,2*length(file_name_Dy))];
        upper = [1,300*ones(1,2*length(file_name_Dy))];
        options = optimoptions('lsqcurvefit','Display','off');
        [beta,~,residual,~,~,~,jacobian]  = ...
            lsqcurvefit(fitFunc,Initial_guess,xVal(1:sum(N_timeSlices(1:length(file_name_Dy))),:),...
            AN_asmbld_avg(1:sum(N_timeSlices(1:length(file_name_Dy)))),lower,upper,options);

        ANfit_tau = beta(1);
        ANfit_Nis = [beta(2:1+length(file_name_Dy))'; MeanAtomNum(length(file_name_Dy)+1:end)] ;
        ANfit_Nfs = beta(2+length(file_name_Dy):1+2*length(file_name_Dy))';
        fprintf('====================================================== \n');
        fprintf('AN Fit results: ====================================== \n');
        fprintf('tau = %0.2f ms\n',ANfit_tau*1e3);
        fprintf('avg. Ni = %0.2f k\n',mean(ANfit_Nis));
        fprintf('avg. Nf = %0.2f k\n',mean(ANfit_Nfs));
        fprintf('====================================================== \n');
        nu_fit_AN = length(beta);
    else
        ANfit_tau = 1;
        ANfit_Nis = MeanAtomNum;
        ANfit_Nfs = MeanAtomNum(1:length(file_name_Dy));
        nu_fit_AN = 0;
    end
    %% Start the fitting ======================================================
    if strcmp(fitMethod,'one')
        %% Fitting method 'one' ===============================================
        B = phnFits_sim2Dfits_struct(file_name_Dy,file_name_Stat,'varType_Dy_ph_t','S_Dy+',...
            'varType_Dy_ph_az','S_Dy+','varType_Dy_amp','S_Dy_scaled+',...
            'varType_Dy_Qi','S_Dy+','varType_Dy_Qf','S_Dy','PersCurr',false,...
            'varType_Stat_dPhidt','I','varType_Dy_dPhidt','None');
        fitFunc = @(params,x)Phn_AllRadDy_2D(x,Phi,nlobes,params(2),params(1),params(3),...
            Rad_i,t_exp,t_phn1,...
            [params(4+length(file_name_Stat)*3)*ones(length(file_name_Dy),1);...                % Qi
            params(4+length(file_name_Stat)*3);...                                              % Qi
            params(5+length(file_name_Stat)*3:3+length(file_name_Stat)*4)'],...                 % Qi
            [params(4+length(file_name_Stat)*2)*ANfit_Nis(1:length(file_name_Dy))/ANbase;...      % amplitude
            params(4+length(file_name_Stat)*2)*ANfit_Nis(length(file_name_Dy)+1)/ANbase;...       % amplitude
            params(5+length(file_name_Stat)*2:3+length(file_name_Stat)*3)'],...                 % amplitude
            [params(4)*ones(length(file_name_Dy),1);...                                         % ph_t
            params(4);...                                                                       % ph_t
            params(5:3+length(file_name_Stat)*1)'],...                                          % ph_t
            [params(4+length(file_name_Stat)*1)*ones(length(file_name_Dy),1);...                % ph_az
            params(4+length(file_name_Stat)*1);...                                              % ph_az
            params(5+length(file_name_Stat)*1:3+length(file_name_Stat)*2)'],...                 % ph_az
            ANfit_Nis,Rad_f,...
            params(4+length(file_name_Stat)*4)*ones(length(file_name_Dy),1),...                 % Qf
            ANfit_Nfs,ANfit_tau*ones(length(file_name_Dy),1),...
            'vector',true,'N_Dy',length(file_name_Dy),'N_Stat',length(file_name_Stat),...
            'PersCurr',false);
        % Order: hub, alpha, cTheta0, Ph_t, Ph_az, Amp, Qf, Qi
        Initial_guess = [hub_guess, alpha_guess,cTheta0_guess,...
            mean(Ph_t(1:length(file_name_Dy))), Ph_t(2+length(file_name_Dy):length(file_names))',...
            mean(Ph_az(1:length(file_name_Dy))), Ph_az(2+length(file_name_Dy):length(file_names))',...
            mean(Amp(1:length(file_name_Dy))), Amp(2+length(file_name_Dy):length(file_names))',...
            mean(Qi(1:length(file_name_Dy))), Qi(2+length(file_name_Dy):length(file_names))',7];
        lower = [-Inf,0,0,...
            -2*pi*ones(1,length(file_name_Stat)),...
            -2*pi*ones(1,length(file_name_Stat)),...
            zeros(1,length(file_name_Stat)),...
            zeros(1,length(file_name_Stat)),0];
        upper = [Inf,1,Inf,...
            +2*pi*ones(1,length(file_name_Stat)),...
            +2*pi*ones(1,length(file_name_Stat)),...
            20e3*ones(1,length(file_name_Stat)),...
            50*ones(1,length(file_name_Stat)),50];
        options = optimoptions('lsqcurvefit','Display','off');
        [beta,~,residual,~,~,~,jacobian]  = ...
            lsqcurvefit(fitFunc,Initial_guess,xVal,n1D_asmbld_avg(:),lower,upper,options);
        ci = nlparci(beta,residual,'jacobian',jacobian,'alpha',1-0.6827); ci = diff(ci,1,2)/2;
        % Get useful data =====================================================
        hub_fit_val = beta(1,1); hub_fit_err = ci(1);
        alpha_fit_val = beta(1,2); alpha_fit_err = ci(2);
        cTheta0_fit_val = beta(1,3); cTheta0_fit_err = ci(3);
        [~, Amp2D_fit] = fitFunc(beta,xVal_fitPlot);
        [~, Amp2D_noHub] = fitFunc([0 beta(2:end)],xVal_fitPlot);
        MSE = sum(residual.^2)/(length(residual)-length(beta));
        % Generate the 1D fits ===============================================
        % First normalize =====================================================
        norm_constant = (MeanAtomNum*1e3/2/pi);
        deltan0_asmbld_avg = zeros(size(EvTime_asmbld));
        deltan0_asmbld_dev = zeros(size(EvTime_asmbld));
        for ii=1:length(file_names)
            if ii == 1
                idx_start_data = 1;
            else
                idx_start_data = sum(N_timeSlices(1:ii-1))+1;
            end
            idx_end_data = sum(N_timeSlices(1:ii));
            n1D_norm = n1D_asmbld_avg(:,idx_start_data:idx_end_data)/norm_constant(ii);
            time_norm = xVal(idx_start_data:idx_end_data,1);
            xtradPhidt = 0;

            % Now fit each normalized n1D to a 1D sine ========================           
            for jj=1:size(n1D_norm,2)
                if ii>length(file_name_Dy)
                    [a(jj), da(jj),~,~] = Ring_PhnEv_fitEachNorm1D(Phi,n1D_norm(:,jj)+1,nlobes,...
                        beta(3+length(file_name_Stat)*1+ii-length(file_name_Dy))+...
                        xtradPhidt*(time_norm(jj)-0),...
                        'FitPhase',false,'YLims',[1-0.5 1+0.5]);
                else
                    [a(jj), da(jj),~,~] = Ring_PhnEv_fitEachNorm1D(Phi,n1D_norm(:,jj)+1,nlobes,...
                        beta(4+length(file_name_Stat)*1)+0,...
                        'FitPhase',false,'YLims',[1-0.5 1+0.5]);
                end
            end
            % Now de-normalize and save n1D amplitudes in real units ==========
            a = a*norm_constant(ii); da = da*norm_constant(ii);
            deltan0_asmbld_avg(idx_start_data:idx_end_data) = a;
            deltan0_asmbld_dev(idx_start_data:idx_end_data) = da;
            clear n1D_norm time_norm a da
        end
        % Plot the dynamic rings =============================================
        idx_start_data = 1; idx_end_data = sum(N_timeSlices(1:length(file_name_Dy)));
        idx_start_fit = 1; idx_end_fit = length(file_name_Dy)*N_fitPlot;
        PlotPhnDy(1,N_timeSlices(1:length(file_name_Dy)),EvTime_asmbld(idx_start_data:idx_end_data),...
            deltan0_asmbld_avg(idx_start_data:idx_end_data),deltan0_asmbld_dev(idx_start_data:idx_end_data),...
            xVal_fitPlot(idx_start_fit:idx_end_fit,1),Amp2D_fit(idx_start_fit:idx_end_fit),...
            Rad_i(1:length(file_name_Dy)),Rad_f(1:length(file_name_Dy)),...
            t_phn1(1:length(file_name_Dy)),t_exp(1:length(file_name_Dy)),...
            cTheta0_fit_val,...
            ANfit_tau*ones(length(file_name_Dy),1),ANfit_Nis,ANfit_Nfs,...
            alpha_fit_val,...
            'xlimits',xLimitsDy,'ylimitsLeft',yLimsDeltan);
        % Plot the static rings ===============================================
        idx_start_data = sum(N_timeSlices(1:length(file_name_Dy)))+1; 
        idx_end_data = sum(N_timeSlices);
        idx_start_fit = length(file_name_Dy)*N_fitPlot + 1; 
        idx_end_fit = length(file_names)*N_fitPlot;
        PlotPhnStat(2,N_timeSlices(length(file_name_Dy)+1:end),EvTime_asmbld(idx_start_data:idx_end_data),...
            deltan0_asmbld_avg(idx_start_data:idx_end_data),deltan0_asmbld_dev(idx_start_data:idx_end_data),...
            xVal_fitPlot(idx_start_fit:idx_end_fit,1),Amp2D_fit(idx_start_fit:idx_end_fit),...
            Rad_i(length(file_name_Dy)+1:end),t_phn1(length(file_name_Dy)+1:end),...
            t_exp(length(file_name_Dy)+1:end),...
            cTheta0_fit_val,MeanAtomNum,alpha_fit_val,...
            B.Qi(length(file_name_Dy)+1:end),B.Qi_err(length(file_name_Dy)+1:end),...
            'xlimits',[0 250],'ylimits',yLimsDeltan);
        % Print the fit values ================================================
        fprintf('====================================================== \n');
        fprintf('Fit results: one ===================================== \n');
        fprintf('alpha = %0.3f \x00B1 %0.3f \n',alpha_fit_val,alpha_fit_err);
        fprintf('hub = %0.3f \x00B1 %0.3f \n',hub_fit_val,hub_fit_err);
        fprintf('cTheta (N = 100k, R = 40 um) = %0.2f \x00B1 %0.2f mm/s\n',cTheta0_fit_val*1e-3,cTheta0_fit_err*1e-3);
        fprintf('cTheta (N = %0.0f k, Ri) = %0.2f \x00B1 %0.2f mm/s\n',...
            mean(MeanAtomNum(1:length(file_name_Dy))),...
            SoundSpeed(cTheta0_fit_val,mean(Rad_i(1:length(file_name_Dy))),...
            mean(MeanAtomNum(1:length(file_name_Dy))),alpha_fit_val)*1e-3,cTheta0_fit_err*1e-3);
        fprintf('cTheta (N = %0.0f k, Rf) = %0.2f \x00B1 %0.2f mm/s\n',...
            mean(MeanAtomNum(1:length(file_name_Dy))),...
            SoundSpeed(cTheta0_fit_val,mean(Rad_f(1:length(file_name_Dy))),...
            mean(MeanAtomNum(1:length(file_name_Dy))),alpha_fit_val)*1e-3,cTheta0_fit_err*1e-3);
        fprintf('deltani = %0.3f \x00B1 %0.3f k\n',beta(4+length(file_name_Stat)*2)...
            *mean(ANfit_Nis(1:length(file_name_Dy)))/ANbase*1e-3,...
            ci(4+length(file_name_Stat)*2)...
            *mean(ANfit_Nis(1:length(file_name_Dy)))/ANbase*1e-3);    
        fprintf('sqrt(MSE) = %0.2f k\n',sqrt(MSE)*1E-3);
        fprintf('DOF = %d - %d\n',length(residual),length(beta));
        fprintf('====================================================== \n');
        cThetai = SoundSpeed(cTheta0_fit_val,mean(Rad_i(1:length(file_name_Dy))),...
            mean(MeanAtomNum(1:length(file_name_Dy))),alpha_fit_val)*1e-3;
        cThetai_err = cTheta0_fit_err*1e-3;
        deltani = beta(4+length(file_name_Stat)*2)...
            *mean(ANfit_Nis(1:length(file_name_Dy)))/ANbase*1e-3;
        deltani_err = ci(4+length(file_name_Stat)*2)...
            *mean(ANfit_Nis(1:length(file_name_Dy)))/ANbase*1e-3;
    elseif strcmp(fitMethod,'two')
        %% Fitting method 'two' ===============================================
        B = phnFits_sim2Dfits_struct(file_name_Dy,file_name_Stat,'varType_Dy_ph_t','I',...
            'varType_Dy_ph_az','S_Dy+','varType_Dy_amp','S_Dy_scaled+',...
            'varType_Dy_Qi','S_Dy+','varType_Dy_Qf','S_Dy','PersCurr',false,...
            'varType_Stat_dPhidt','I','varType_Dy_dPhidt','None');
        fitFunc = @(params,x)Phn_AllRadDy_2D(x,Phi,nlobes,params(2),params(1),params(3),...
            Rad_i,t_exp,t_phn1,...
            [params(4+length(file_names)*1+length(file_name_Stat)*2)*ones(length(file_name_Dy),1);...           % Qi
            params(4+length(file_names)*1+length(file_name_Stat)*2);...                                         % Qi
            params(5+length(file_names)*1+length(file_name_Stat)*2:...                                          % Qi
            3+length(file_names)*1+length(file_name_Stat)*3)'],...                                              % Qi
            [params(4+length(file_name_Stat)*1+length(file_names)*1)*ANfit_Nis(1:length(file_name_Dy))/ANbase;... % amplitude
            params(4+length(file_name_Stat)*1+length(file_names)*1)*ANfit_Nis(length(file_name_Dy)+1)/ANbase;...  % amplitude
            params(5+length(file_name_Stat)*1+length(file_names)*1:...                                          % amplitude
            3+length(file_names)*1+length(file_name_Stat)*2)'],...                                              % amplitude
            params(4:3+length(file_names))',...                                                                 % ph_t
            [params(4+length(file_names)*1+length(file_name_Stat)*0)*ones(length(file_name_Dy),1);...           % ph_az
            params(4+length(file_names)*1+length(file_name_Stat)*0);...                                         % ph_az
            params(5+length(file_names)*1+length(file_name_Stat)*0:...                                          % ph_az
            3+length(file_names)*1+length(file_name_Stat)*1)'],...                                              % ph_az
            ANfit_Nis,Rad_f,...
            params(4+length(file_names)*1+length(file_name_Stat)*3)*ones(length(file_name_Dy),1),...            % Qf
            ANfit_Nfs,ANfit_tau*ones(length(file_name_Dy),1),...
            'vector',true,'N_Dy',length(file_name_Dy),'N_Stat',length(file_name_Stat),...
            'PersCurr',false);
        % Order: hub, alpha, cTheta0, Ph_t, Ph_az, Amp, Qf, Qi
        Initial_guess = [hub_guess, alpha_guess,cTheta0_guess,Ph_t',...
            mean(Ph_az(1:length(file_name_Dy))), Ph_az(2+length(file_name_Dy):length(file_names))',...
            mean(Amp(1:length(file_name_Dy))), Amp(2+length(file_name_Dy):length(file_names))',...
            mean(Qi(1:length(file_name_Dy))), Qi(2+length(file_name_Dy):length(file_names))',7];
        lower = [-Inf,0,0,...
            -2*pi*ones(1,length(file_names)),...
            -2*pi*ones(1,length(file_name_Stat)),...
            zeros(1,length(file_name_Stat)),...
            zeros(1,length(file_name_Stat)),0];
        upper = [Inf,1,Inf,...
            +2*pi*ones(1,length(file_names)),...
            +2*pi*ones(1,length(file_name_Stat)),...
            20e3*ones(1,length(file_name_Stat)),...
            50*ones(1,length(file_name_Stat)),50];

        options = optimoptions('lsqcurvefit','Display','off');
        [beta,~,residual,~,~,~,jacobian]  = ...
            lsqcurvefit(fitFunc,Initial_guess,xVal,n1D_asmbld_avg(:),lower,upper,options);
        ci = nlparci(beta,residual,'jacobian',jacobian,'alpha',1-0.6827); ci = diff(ci,1,2)/2;
        % Get useful data =====================================================
        hub_fit_val = beta(1,1); hub_fit_err = ci(1);
        alpha_fit_val = beta(1,2); alpha_fit_err = ci(2);
        cTheta0_fit_val = beta(1,3); cTheta0_fit_err = ci(3);
        [~, Amp2D_fit] = fitFunc(beta,xVal_fitPlot);
        [~, Amp2D_noHub] = fitFunc([0 beta(2:end)],xVal_fitPlot);
        MSE = sum(residual.^2)/(length(residual)-length(beta));
        % Generate the 1D fits ===============================================
        % First normalize =====================================================
        norm_constant = (MeanAtomNum*1e3/2/pi);
        deltan0_asmbld_avg = zeros(size(EvTime_asmbld));
        deltan0_asmbld_dev = zeros(size(EvTime_asmbld));
        for ii=1:length(file_names)
            if ii == 1 
                idx_start_data = 1;
            else
                idx_start_data = sum(N_timeSlices(1:ii-1))+1;
            end
            idx_end_data = sum(N_timeSlices(1:ii));
            n1D_norm = n1D_asmbld_avg(:,idx_start_data:idx_end_data)/norm_constant(ii);
            time_norm = xVal(idx_start_data:idx_end_data,1);
            xtradPhidt = 0;

            % Now fit each normalized n1D to a 1D sine ========================           
            for jj=1:size(n1D_norm,2)
                if ii>length(file_name_Dy)
                    [a(jj), da(jj),~,~] = Ring_PhnEv_fitEachNorm1D(Phi,n1D_norm(:,jj)+1,nlobes,...
                        beta(3+length(file_names)*1+ii-length(file_name_Dy))+...
                        xtradPhidt*(time_norm(jj)-0),...
                        'FitPhase',false,'YLims',[1-0.5 1+0.5]);
                else
                    [a(jj), da(jj),~,~] = Ring_PhnEv_fitEachNorm1D(Phi,n1D_norm(:,jj)+1,nlobes,...
                        beta(4+length(file_names)*1)+0,...
                        'FitPhase',false,'YLims',[1-0.5 1+0.5]);
                end
            end
            % Now de-normalize and save n1D amplitudes in real units ==========
            a = a*norm_constant(ii); da = da*norm_constant(ii);
            deltan0_asmbld_avg(idx_start_data:idx_end_data) = a;
            deltan0_asmbld_dev(idx_start_data:idx_end_data) = da;
            clear n1D_norm time_norm a da
        end
        % Plot the dynamic rings ==================================================
        idx_start_data = 1; idx_end_data = sum(N_timeSlices(1:length(file_name_Dy)));
        idx_start_fit = 1; idx_end_fit = length(file_name_Dy)*N_fitPlot;
        PlotPhnDy(1,N_timeSlices(1:length(file_name_Dy)),EvTime_asmbld(idx_start_data:idx_end_data),...
            deltan0_asmbld_avg(idx_start_data:idx_end_data),deltan0_asmbld_dev(idx_start_data:idx_end_data),...
            xVal_fitPlot(idx_start_fit:idx_end_fit,1),Amp2D_fit(idx_start_fit:idx_end_fit),...
            Rad_i(1:length(file_name_Dy)),Rad_f(1:length(file_name_Dy)),...
            t_phn1(1:length(file_name_Dy)),t_exp(1:length(file_name_Dy)),...
            cTheta0_fit_val,...
            ANfit_tau*ones(length(file_name_Dy),1),ANfit_Nis,ANfit_Nfs,...
            alpha_fit_val,...
            'xlimits',xLimitsDy,'ylimitsLeft',yLimsDeltan);
        % plot the static rings ===============================================
        idx_start_data = sum(N_timeSlices(1:length(file_name_Dy)))+1; 
        idx_end_data = sum(N_timeSlices);
        idx_start_fit = length(file_name_Dy)*N_fitPlot + 1; 
        idx_end_fit = length(file_names)*N_fitPlot;
        PlotPhnStat(2,N_timeSlices(length(file_name_Dy)+1:end),EvTime_asmbld(idx_start_data:idx_end_data),...
            deltan0_asmbld_avg(idx_start_data:idx_end_data),deltan0_asmbld_dev(idx_start_data:idx_end_data),...
            xVal_fitPlot(idx_start_fit:idx_end_fit,1),Amp2D_fit(idx_start_fit:idx_end_fit),...
            Rad_i(length(file_name_Dy)+1:end),t_phn1(length(file_name_Dy)+1:end),...
            t_exp(length(file_name_Dy)+1:end),...
            cTheta0_guess,MeanAtomNum,alpha_guess,...
            B.Qi(length(file_name_Dy)+1:end),B.Qi_err(length(file_name_Dy)+1:end),...
            'xlimits',[0 250],'ylimits',yLimsDeltan);
        % Print the fit values ================================================
        fprintf('====================================================== \n');
        fprintf('Fit results: two ===================================== \n');
        fprintf('alpha = %0.3f \x00B1 %0.3f \n',alpha_fit_val,alpha_fit_err);
        fprintf('hub = %0.3f \x00B1 %0.3f \n',hub_fit_val,hub_fit_err);
        fprintf('cTheta (N = 100k, R = 40 um) = %0.2f \x00B1 %0.2f mm/s\n',cTheta0_fit_val*1e-3,cTheta0_fit_err*1e-3);
        fprintf('cTheta (N = %0.0f k, Ri) = %0.2f \x00B1 %0.2f mm/s\n',...
            mean(MeanAtomNum(1:length(file_name_Dy))),...
            SoundSpeed(cTheta0_fit_val,mean(Rad_i(1:length(file_name_Dy))),...
            mean(MeanAtomNum(1:length(file_name_Dy))),alpha_fit_val)*1e-3,cTheta0_fit_err*1e-3);
        fprintf('cTheta (N = %0.0f k, Rf) = %0.2f \x00B1 %0.2f mm/s\n',...
            mean(MeanAtomNum(1:length(file_name_Dy))),...
            SoundSpeed(cTheta0_fit_val,mean(Rad_f(1:length(file_name_Dy))),...
            mean(MeanAtomNum(1:length(file_name_Dy))),alpha_fit_val)*1e-3,cTheta0_fit_err*1e-3);
        fprintf('deltani = %0.3f \x00B1 %0.3f k\n',...
            beta(4+length(file_name_Stat)*1+length(file_names)*1)...
            *mean(ANfit_Nis(1:length(file_name_Dy)))/ANbase*1e-3,...
            ci(4+length(file_name_Stat)*1+length(file_names)*1)...
            *mean(ANfit_Nis(1:length(file_name_Dy)))/ANbase*1e-3);   
        fprintf('sqrt(MSE) = %0.2f k\n',sqrt(MSE)*1E-3);
        fprintf('DOF = %d - %d\n',length(residual),length(beta));
        fprintf('====================================================== \n');
        cThetai = SoundSpeed(cTheta0_fit_val,mean(Rad_i(1:length(file_name_Dy))),...
            mean(MeanAtomNum(1:length(file_name_Dy))),alpha_fit_val)*1e-3;
        cThetai_err = cTheta0_fit_err*1e-3;
        deltani = beta(4+length(file_name_Stat)*1+length(file_names)*1)...
            *mean(ANfit_Nis(1:length(file_name_Dy)))/ANbase*1e-3;
        deltani_err = ci(4+length(file_name_Stat)*1+length(file_names)*1)...
            *mean(ANfit_Nis(1:length(file_name_Dy)))/ANbase*1e-3;

    elseif strcmp(fitMethod,'three')
        %% Fitting method 'three' ===============================================
        B = phnFits_sim2Dfits_struct(file_name_Dy,file_name_Stat,'varType_Dy_ph_t','S_Dy+',...
            'varType_Dy_ph_az','I','varType_Dy_amp','S_Dy_scaled+',...
            'varType_Dy_Qi','S_Dy+','varType_Dy_Qf','S_Dy','PersCurr',false,...
            'varType_Stat_dPhidt','I','varType_Dy_dPhidt','None');
        fitFunc = @(params,x)Phn_AllRadDy_2D(x,Phi,nlobes,params(2),params(1),params(3),...
            Rad_i,t_exp,t_phn1,...
            [params(4+length(file_names)*1+length(file_name_Stat)*2)*ones(length(file_name_Dy),1);...           % Qi
            params(4+length(file_names)*1+length(file_name_Stat)*2);...                                         % Qi
            params(5+length(file_names)*1+length(file_name_Stat)*2:...                                          % Qi
            3+length(file_names)*1+length(file_name_Stat)*3)'],...                                              % Qi
            [params(4+length(file_name_Stat)*1+length(file_names)*1)*ANfit_Nis(1:length(file_name_Dy))/ANbase;... % amplitude
            params(4+length(file_name_Stat)*1+length(file_names)*1)*ANfit_Nis(length(file_name_Dy)+1)/ANbase;...  % amplitude
            params(5+length(file_name_Stat)*1+length(file_names)*1:...                                          % amplitude
            3+length(file_names)*1+length(file_name_Stat)*2)'],...                                              % amplitude
            [params(4)*ones(length(file_name_Dy),1);...                                                         % ph_t
            params(4);...                                                                                       % ph_t
            params(5:3+length(file_name_Stat)*1)'],...                                                          % ph_t
            params(4+length(file_name_Stat)*1:3+length(file_name_Stat)*1+length(file_names))',...               % ph_az
            ANfit_Nis,Rad_f,...
            params(4+length(file_names)*1+length(file_name_Stat)*3)*ones(length(file_name_Dy),1),...            % Qf
            ANfit_Nfs,ANfit_tau*ones(length(file_name_Dy),1),...
            'vector',true,'N_Dy',length(file_name_Dy),'N_Stat',length(file_name_Stat),...
            'PersCurr',false);
        % Order: hub, alpha, cTheta0, Ph_t, Ph_az, Amp, Qf, Qi
        Initial_guess = [hub_guess, alpha_guess,cTheta0_guess,...
            mean(Ph_t(1:length(file_name_Dy))), Ph_t(2+length(file_name_Dy):length(file_names))',...
            Ph_az',...
            mean(Amp(1:length(file_name_Dy))), Amp(2+length(file_name_Dy):length(file_names))',...
            mean(Qi(1:length(file_name_Dy))), Qi(2+length(file_name_Dy):length(file_names))',7];
        lower = [-Inf,0,0,...
            -2*pi*ones(1,length(file_name_Stat)),...
            -2*pi*ones(1,length(file_names)),...
            zeros(1,length(file_name_Stat)),...
            zeros(1,length(file_name_Stat)),0];
        upper = [Inf,1,Inf,...
            +2*pi*ones(1,length(file_name_Stat)),...
            +2*pi*ones(1,length(file_names)),...
            20e3*ones(1,length(file_name_Stat)),...
            50*ones(1,length(file_name_Stat)),50];

        options = optimoptions('lsqcurvefit','Display','off');
        [beta,~,residual,~,~,~,jacobian]  = ...
            lsqcurvefit(fitFunc,Initial_guess,xVal,n1D_asmbld_avg(:),lower,upper,options);
        ci = nlparci(beta,residual,'jacobian',jacobian,'alpha',1-0.6827); ci = diff(ci,1,2)/2;
        % Get useful data =====================================================
        hub_fit_val = beta(1,1); hub_fit_err = ci(1);
        alpha_fit_val = beta(1,2); alpha_fit_err = ci(2);
        cTheta0_fit_val = beta(1,3); cTheta0_fit_err = ci(3);
        [~, Amp2D_fit] = fitFunc(beta,xVal_fitPlot);
        [~, Amp2D_noHub] = fitFunc([0 beta(2:end)],xVal_fitPlot);
        MSE = sum(residual.^2)/(length(residual)-length(beta));
        % Generate the 1D fits ===============================================
        % First normalize =====================================================
        norm_constant = (MeanAtomNum*1e3/2/pi);
        deltan0_asmbld_avg = zeros(size(EvTime_asmbld));
        deltan0_asmbld_dev = zeros(size(EvTime_asmbld));
        for ii=1:length(file_names)
            if ii == 1 
                idx_start_data = 1;
            else
                idx_start_data = sum(N_timeSlices(1:ii-1))+1;
            end
            idx_end_data = sum(N_timeSlices(1:ii));
            n1D_norm = n1D_asmbld_avg(:,idx_start_data:idx_end_data)/norm_constant(ii);
            time_norm = xVal(idx_start_data:idx_end_data,1);
            xtradPhidt = 0;

            % Now fit each normalized n1D to a 1D sine ========================           
            for jj=1:size(n1D_norm,2)
                [a(jj), da(jj),~,~] = Ring_PhnEv_fitEachNorm1D(Phi,n1D_norm(:,jj)+1,nlobes,...
                        beta(3+length(file_name_Stat)*1+ii)+...
                        xtradPhidt*(time_norm(jj)-0),...
                        'FitPhase',false,'YLims',[1-0.5 1+0.5]);
            end
            % Now de-normalize and save n1D amplitudes in real units ==========
            a = a*norm_constant(ii); da = da*norm_constant(ii);
            deltan0_asmbld_avg(idx_start_data:idx_end_data) = a;
            deltan0_asmbld_dev(idx_start_data:idx_end_data) = da;
            clear n1D_norm time_norm a da
        end
        % plot the dynamic rings ==================================================
        idx_start_data = 1; idx_end_data = sum(N_timeSlices(1:length(file_name_Dy)));
        idx_start_fit = 1; idx_end_fit = length(file_name_Dy)*N_fitPlot;
        PlotPhnDy(1,N_timeSlices(1:length(file_name_Dy)),EvTime_asmbld(idx_start_data:idx_end_data),...
            deltan0_asmbld_avg(idx_start_data:idx_end_data),deltan0_asmbld_dev(idx_start_data:idx_end_data),...
            xVal_fitPlot(idx_start_fit:idx_end_fit,1),Amp2D_fit(idx_start_fit:idx_end_fit),...
            Rad_i(1:length(file_name_Dy)),Rad_f(1:length(file_name_Dy)),...
            t_phn1(1:length(file_name_Dy)),t_exp(1:length(file_name_Dy)),...
            cTheta0_guess,...
            ANfit_tau*ones(length(file_name_Dy),1),ANfit_Nis,ANfit_Nfs,...
            alpha_guess,...
            'xlimits',xLimitsDy,'ylimitsLeft',yLimsDeltan);
        % plot the static rings ===============================================
        idx_start_data = sum(N_timeSlices(1:length(file_name_Dy)))+1; 
        idx_end_data = sum(N_timeSlices);
        idx_start_fit = length(file_name_Dy)*N_fitPlot + 1; 
        idx_end_fit = length(file_names)*N_fitPlot;
        PlotPhnStat(2,N_timeSlices(length(file_name_Dy)+1:end),EvTime_asmbld(idx_start_data:idx_end_data),...
            deltan0_asmbld_avg(idx_start_data:idx_end_data),deltan0_asmbld_dev(idx_start_data:idx_end_data),...
            xVal_fitPlot(idx_start_fit:idx_end_fit,1),Amp2D_fit(idx_start_fit:idx_end_fit),...
            Rad_i(length(file_name_Dy)+1:end),t_phn1(length(file_name_Dy)+1:end),...
            t_exp(length(file_name_Dy)+1:end),...
            cTheta0_guess,MeanAtomNum,alpha_guess,...
            B.Qi(length(file_name_Dy)+1:end),B.Qi_err(length(file_name_Dy)+1:end),...
            'xlimits',[0 250],'ylimits',yLimsDeltan);
        % Print the fit values ================================================
        fprintf('====================================================== \n');
        fprintf('Fit results: three ===================================== \n');
        fprintf('alpha = %0.3f \x00B1 %0.3f \n',alpha_fit_val,alpha_fit_err);
        fprintf('hub = %0.3f \x00B1 %0.3f \n',hub_fit_val,hub_fit_err);
        fprintf('cTheta (N = 100k, R = 40 um) = %0.2f \x00B1 %0.2f mm/s\n',cTheta0_fit_val*1e-3,cTheta0_fit_err*1e-3);
        fprintf('cTheta (N = %0.0f k, Ri) = %0.2f \x00B1 %0.2f mm/s\n',...
            mean(MeanAtomNum(1:length(file_name_Dy))),...
            SoundSpeed(cTheta0_fit_val,mean(Rad_i(1:length(file_name_Dy))),...
            mean(MeanAtomNum(1:length(file_name_Dy))),alpha_fit_val)*1e-3,cTheta0_fit_err*1e-3);
        fprintf('cTheta (N = %0.0f k, Rf) = %0.2f \x00B1 %0.2f mm/s\n',...
            mean(MeanAtomNum(1:length(file_name_Dy))),...
            SoundSpeed(cTheta0_fit_val,mean(Rad_f(1:length(file_name_Dy))),...
            mean(MeanAtomNum(1:length(file_name_Dy))),alpha_fit_val)*1e-3,cTheta0_fit_err*1e-3);
        fprintf('deltani = %0.3f \x00B1 %0.3f k\n',...
            beta(4+length(file_name_Stat)*1+length(file_names)*1)*...
            mean(ANfit_Nis(1:length(file_name_Dy)))/ANbase*1e-3,...
            ci(4+length(file_name_Stat)*1+length(file_names)*1)*...
            mean(ANfit_Nis(1:length(file_name_Dy)))/ANbase*1e-3);

        fprintf('sqrt(MSE) = %0.2f k\n',sqrt(MSE)*1E-3);
        fprintf('DOF = %d - %d\n',length(residual),length(beta));
        fprintf('====================================================== \n');
        cThetai = SoundSpeed(cTheta0_fit_val,mean(Rad_i(1:length(file_name_Dy))),...
            mean(MeanAtomNum(1:length(file_name_Dy))),alpha_fit_val)*1e-3;
        cThetai_err = cTheta0_fit_err*1e-3;
        deltani = beta(4+length(file_name_Stat)*1+length(file_names)*1)*...
            mean(ANfit_Nis(1:length(file_name_Dy)))/ANbase*1e-3;
        deltani_err = ci(4+length(file_name_Stat)*1+length(file_names)*1)*...
            mean(ANfit_Nis(1:length(file_name_Dy)))/ANbase*1e-3;
    elseif strcmp(fitMethod,'four')
        %% Fitting method 'four' ===============================================
        B = phnFits_sim2Dfits_struct(file_name_Dy,file_name_Stat,'varType_Dy_ph_t','I',...
            'varType_Dy_ph_az','I','varType_Dy_amp','S_Dy_scaled+',...
            'varType_Dy_Qi','S_Dy+','varType_Dy_Qf','S_Dy','PersCurr',false,...
            'varType_Stat_dPhidt','I','varType_Dy_dPhidt','None');
        fitFunc = @(params,x)Phn_AllRadDy_2D(x,Phi,nlobes,params(2),params(1),params(3),...
            Rad_i,t_exp,t_phn1,...
            [params(4+length(file_names)*2+length(file_name_Stat)*1)*ones(length(file_name_Dy),1);...           % Qi
            params(4+length(file_names)*2+length(file_name_Stat)*1);...                                         % Qi
            params(5+length(file_names)*2+length(file_name_Stat)*1:...                                          % Qi
            3+length(file_names)*2+length(file_name_Stat)*2)'],...                                              % Qi
            [params(4+length(file_name_Stat)*0+length(file_names)*2)*ANfit_Nis(1:length(file_name_Dy))/ANbase;... % amplitude
            params(4+length(file_name_Stat)*0+length(file_names)*2)*ANfit_Nis(length(file_name_Dy)+1)/ANbase;...  % amplitude
            params(5+length(file_name_Stat)*0+length(file_names)*2:...                                          % amplitude
            3+length(file_names)*2+length(file_name_Stat)*1)'],...                                              % amplitude
            params(4:3+length(file_names))',...                                                                 % ph_t
            params(4+length(file_names):3+length(file_names)*2)',...                                            % ph_az
            ANfit_Nis,Rad_f,...
            params(4+length(file_names)*2+length(file_name_Stat)*2)*ones(length(file_name_Dy),1),...            % Qf
            ANfit_Nfs,ANfit_tau*ones(length(file_name_Dy),1),...
            'vector',true,'N_Dy',length(file_name_Dy),'N_Stat',length(file_name_Stat),...
            'PersCurr',false);
        % Order: hub, alpha, cTheta0, Ph_t, Ph_az, Amp, Qf, Qi
        Initial_guess = [hub_guess, alpha_guess,cTheta0_guess,Ph_t', Ph_az',...
            mean(Amp(1:length(file_name_Dy))), Amp(2+length(file_name_Dy):length(file_names))',...
            mean(Qi(1:length(file_name_Dy))), Qi(2+length(file_name_Dy):length(file_names))',7];
        lower = [-Inf,0,0,...
            -2*pi*ones(1,length(file_names)),...
            -2*pi*ones(1,length(file_names)),...
            zeros(1,length(file_name_Stat)),...
            zeros(1,length(file_name_Stat)),0];
        upper = [Inf,1,Inf,...
            +2*pi*ones(1,length(file_names)),...
            +2*pi*ones(1,length(file_names)),...
            20e3*ones(1,length(file_name_Stat)),...
            50*ones(1,length(file_name_Stat)),50];

        options = optimoptions('lsqcurvefit','Display','off');
        [beta,~,residual,~,~,~,jacobian]  = ...
            lsqcurvefit(fitFunc,Initial_guess,xVal,n1D_asmbld_avg(:),lower,upper,options);
        ci = nlparci(beta,residual,'jacobian',jacobian,'alpha',1-0.6827); ci = diff(ci,1,2)/2;
        % Get useful data =====================================================
        hub_fit_val = beta(1,1); hub_fit_err = ci(1);
        alpha_fit_val = beta(1,2); alpha_fit_err = ci(2);
        cTheta0_fit_val = beta(1,3); cTheta0_fit_err = ci(3);
        [~, Amp2D_fit] = fitFunc(beta,xVal_fitPlot);
        [~, Amp2D_noHub] = fitFunc([0 beta(2:end)],xVal_fitPlot);
        MSE = sum(residual.^2)/(length(residual)-length(beta));
        % Generate the 1D fits ===============================================
        % First normalize =====================================================
        norm_constant = (MeanAtomNum*1e3/2/pi);
        deltan0_asmbld_avg = zeros(size(EvTime_asmbld));
        deltan0_asmbld_dev = zeros(size(EvTime_asmbld));
        for ii=1:length(file_names)
            if ii == 1
                idx_start_data = 1;
            else
                idx_start_data = sum(N_timeSlices(1:ii-1))+1;
            end
            idx_end_data = sum(N_timeSlices(1:ii));
            n1D_norm = n1D_asmbld_avg(:,idx_start_data:idx_end_data)/norm_constant(ii);
            time_norm = xVal(idx_start_data:idx_end_data,1);
            xtradPhidt = 0;

            % Now fit each normalized n1D to a 1D sine ========================
            for jj=1:size(n1D_norm,2)
                [a(jj), da(jj),~,~] = Ring_PhnEv_fitEachNorm1D(Phi,n1D_norm(:,jj)+1,nlobes,...
                    beta(3+length(file_names)*1+ii)+...
                    xtradPhidt*(time_norm(jj)-0),...
                    'FitPhase',false,'YLims',[1-0.5 1+0.5]);
            end
            % Now de-normalize and save n1D amplitudes in real units ==========
            a = a*norm_constant(ii); da = da*norm_constant(ii);
            deltan0_asmbld_avg(idx_start_data:idx_end_data) = a;
            deltan0_asmbld_dev(idx_start_data:idx_end_data) = da;
            clear n1D_norm time_norm a da
        end
        % Plot the dynamic rings ==================================================
        idx_start_data = 1; idx_end_data = sum(N_timeSlices(1:length(file_name_Dy)));
        idx_start_fit = 1; idx_end_fit = length(file_name_Dy)*N_fitPlot;
        PlotPhnDy(1,N_timeSlices(1:length(file_name_Dy)),EvTime_asmbld(idx_start_data:idx_end_data),...
            deltan0_asmbld_avg(idx_start_data:idx_end_data),deltan0_asmbld_dev(idx_start_data:idx_end_data),...
            xVal_fitPlot(idx_start_fit:idx_end_fit,1),Amp2D_fit(idx_start_fit:idx_end_fit),...
            Rad_i(1:length(file_name_Dy)),Rad_f(1:length(file_name_Dy)),...
            t_phn1(1:length(file_name_Dy)),t_exp(1:length(file_name_Dy)),...
            cTheta0_guess,...
            ANfit_tau*ones(length(file_name_Dy),1),ANfit_Nis,ANfit_Nfs,...
            alpha_guess,...
            'xlimits',xLimitsDy,'ylimitsLeft',yLimsDeltan);
        % plot the static rings ===============================================
        idx_start_data = sum(N_timeSlices(1:length(file_name_Dy)))+1;
        idx_end_data = sum(N_timeSlices);
        idx_start_fit = length(file_name_Dy)*N_fitPlot + 1;
        idx_end_fit = length(file_names)*N_fitPlot;
        PlotPhnStat(2,N_timeSlices(length(file_name_Dy)+1:end),EvTime_asmbld(idx_start_data:idx_end_data),...
            deltan0_asmbld_avg(idx_start_data:idx_end_data),deltan0_asmbld_dev(idx_start_data:idx_end_data),...
            xVal_fitPlot(idx_start_fit:idx_end_fit,1),Amp2D_fit(idx_start_fit:idx_end_fit),...
            Rad_i(length(file_name_Dy)+1:end),t_phn1(length(file_name_Dy)+1:end),...
            t_exp(length(file_name_Dy)+1:end),...
            cTheta0_guess,MeanAtomNum,alpha_guess,...
            B.Qi(length(file_name_Dy)+1:end),B.Qi_err(length(file_name_Dy)+1:end),...
            'xlimits',[0 250],'ylimits',yLimsDeltan);
        % Print the fit values ================================================
        fprintf('====================================================== \n');
        fprintf('Fit results: four ===================================== \n');
        fprintf('alpha = %0.3f \x00B1 %0.3f \n',alpha_fit_val,alpha_fit_err);
        fprintf('hub = %0.3f \x00B1 %0.3f \n',hub_fit_val,hub_fit_err);
        fprintf('cTheta (N = 100k, R = 40 um) = %0.2f \x00B1 %0.2f mm/s\n',cTheta0_fit_val*1e-3,cTheta0_fit_err*1e-3);
        fprintf('cTheta (N = %0.0f k, Ri) = %0.2f \x00B1 %0.2f mm/s\n',...
            mean(MeanAtomNum(1:length(file_name_Dy))),...
            SoundSpeed(cTheta0_fit_val,mean(Rad_i(1:length(file_name_Dy))),...
            mean(MeanAtomNum(1:length(file_name_Dy))),alpha_fit_val)*1e-3,cTheta0_fit_err*1e-3);
        fprintf('cTheta (N = %0.0f k, Rf) = %0.2f \x00B1 %0.2f mm/s\n',...
            mean(MeanAtomNum(1:length(file_name_Dy))),...
            SoundSpeed(cTheta0_fit_val,mean(Rad_f(1:length(file_name_Dy))),...
            mean(MeanAtomNum(1:length(file_name_Dy))),alpha_fit_val)*1e-3,cTheta0_fit_err*1e-3);
        fprintf('deltani = %0.3f \x00B1 %0.3f k\n',...
            beta(4+length(file_names)*2)*...
            mean(ANfit_Nis(1:length(file_name_Dy)))/ANbase*1e-3,...
            ci(4+length(file_names)*2)*...
            mean(ANfit_Nis(1:length(file_name_Dy)))/ANbase*1e-3);
        fprintf('sqrt(MSE) = %0.2f k\n',sqrt(MSE)*1E-3);
        fprintf('DOF = %d - %d\n',length(residual),length(beta));
        fprintf('====================================================== \n');
        cThetai = SoundSpeed(cTheta0_fit_val,mean(Rad_i(1:length(file_name_Dy))),...
            mean(MeanAtomNum(1:length(file_name_Dy))),alpha_fit_val)*1e-3;
        cThetai_err = cTheta0_fit_err*1e-3;
        deltani = beta(4+length(file_names)*2)*...
            mean(ANfit_Nis(1:length(file_name_Dy)))/ANbase*1e-3;
        deltani_err = ci(4+length(file_names)*2)*...
            mean(ANfit_Nis(1:length(file_name_Dy)))/ANbase*1e-3;
    end

    %% Organize the fit results ============================================
    if exist('B','var')
        B = assblFitVals(B,beta,MeanAtomNum,ANbase);
        B = assblFitErrs(B,ci,MeanAtomNum,ANbase);
        nu_fit = length(beta) + nu_fit_AN;
        fprintf('Qi = %0.3f \x00B1 %0.3f \n',B.Qi(1),B.Qi_err(1));
        fprintf('Qf = %0.3f \x00B1 %0.3f \n',B.Qf(1),B.Qf_err(1));
    end
    %% Save the plots =========================================================
    cd(currentDir);
    if saveData
        figure(1)
        set(gcf, 'PaperUnits', 'inches','PaperPosition', [0 0 10 2*length(file_name_Dy)]);
        print(['ovrlFits_' fitSets '_Dy'],'-dpng','-r200');
        figure(2)
        set(gcf, 'PaperUnits', 'inches','PaperPosition', [0 0 10 2*length(file_name_Stat)]);
        print(['ovrlFits_' fitSets '_Stat'],'-dpng','-r200');
        save(['ovrlFits_' fitSets '.mat']);
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Local Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [results2D, amp] = Phn_AllRadDy_2D(X,phi,nlobes,...
    ALPHA,HUB,cTHETA0,...
    Ri,t_exp,t_phn1,Qi,A,ph_t,ph_az,Ni,...
    Rf,Qf,Nf,tau, varargin)
    % Inputs:
    %   X: 2D array, 1st col is time [s], 2nd col is index referring to
    %       dataset.
    %   phi: 1D array of azimuthal angles [-pi to pi]
    %   nlobes: mode of phonon
    %   ALPHA: single value of alpha for all datasets
    %   HUB: single value of alpha for all datasets
    %   cTHETA0: single value of cTheta [um/s] for N=100k and R = 40um 
    %   Ri,t_exp,t_phn1,Qi,A,ph_t,ph_az,Ni: 1D arrays with a value 
    %       for each dataset
    %   Rf, Qf, Nf, tau: 1D arrays with a value for each Dy dataset
    % Outputs:
    %   results2D: the Phn oscc [k atoms/ rad] as a 2D array.
    %   amp: amp [k atoms/ rad] as a function of time
    % Get optional inputs =================================================
    p = inputParser;
    p.addParameter('vector',false,@(x)islogical(x));
    p.addParameter('N_Stat',7,@(x)(rem(x,1) == 0 & x>0));
    p.addParameter('N_Dy',2,@(x)(rem(x,1) == 0 & x>0));
    p.addParameter('PersCurr',false,@(x)islogical(x));
    p.addParameter('dPhidt',[],@(x)isnumeric(x));
    p.parse(varargin{:});    
    vector = p.Results.vector;
    N_tot = p.Results.N_Stat + p.Results.N_Dy;
    dPhidt = p.Results.dPhidt;
    if isempty(dPhidt)
        dPhidt = zeros(N_tot,1);
    elseif size(dPhidt,1) == p.Results.N_Dy
        dPhidt = [dPhidt; zeros(p.Results.N_Stat,1)];
    elseif size(dPhidt,1) == p.Results.N_Stat
        dPhidt = [zeros(p.Results.N_Dy,1); dPhidt];
    elseif ~size(dPhidt,1) == N_tot
        error('Phn_AllRadDy_2D: dPhidt dimensions is wrong.');
    end
    
    % Initialize ===========================================================
    Time = X(:,1);
    Idx = X(:,2); 
    results2D = zeros(size(phi,1),size(Time,1));
    amp = zeros(size(Time,1),1);
    % Check Inputs ========================================================
    if ( size(t_phn1,1)~=N_tot | size(t_exp,1)~=N_tot | size(Ri,1)~=N_tot | ...
            size(Qi,1)~=N_tot | size(A,1)~=N_tot | ...
            size(ph_t,1)~=N_tot |  size(ph_az,1)~=N_tot |  size(Ni,1)~=N_tot )
        error('Phn_AllRadDy_2D: Size of t_phn1, t_exp, Ri, Qi, A, ph_t, ph_az and AN should be same as total number of files.');
    end
    if (length(unique(Idx))~=N_tot)
        error('Phn_AllRadDy_2D: The number of unique indices should be same as the size of AN.');
    end
    if size(Rf,1)~=p.Results.N_Dy
        error('Phn_AllRadDy_2D: Size of Rf should be same as N_Dy.');
    end
    if size(Qf,1)~=p.Results.N_Dy
        error('Phn_AllRadDy_2D: Size of Qf should be same as N_Dy.');
    end
    if size(Nf,1)~=p.Results.N_Dy
        error('Phn_AllRadDy_2D: Size of Nf should be same as N_Dy.');
    end
    if size(tau,1)~=p.Results.N_Dy
        error('Phn_AllRadDy_2D: Size of tau should be same as N_Dy.');
    end
    % Arrange the data for analysis =======================================
    [Idx, sort_idx] = sort(Idx);
    Time = Time(sort_idx);
    StartPoints = zeros(N_tot,1);
    StopPoints = zeros(N_tot,1);       
    StartPoints(1) = 1; StopPoints(end) = size(Time,1);
    test = Idx(1,1); kk = 2;        
    for ii = 1:1:size(Idx,1)
        if  Idx(ii,1)~= test
            test = Idx(ii,1);
            StartPoints(kk) = ii;
            StopPoints(kk-1) = ii-1;
            kk = kk + 1;
        end
    end
    % Call the 2D Phn functions ==========================================
    for ii = 1:1:N_tot
        % Make the theta-time 2D grid =====================================
        [TS_FIT,THS_FIT] = meshgrid(Time(StartPoints(ii):StopPoints(ii)),phi); 
        % Call the appropriate 2D function  ===============================
        if ii>p.Results.N_Dy
            cTheta_i = SoundSpeed(cTHETA0,Ri(ii),Ni(ii),ALPHA);
            [set_temp, amp_temp] =  RingStat_phnEv_2D(TS_FIT(:),THS_FIT(:),Ri(ii),nlobes,...
                A(ii),cTheta_i,Qi(ii),ph_t(ii),ph_az(ii),...
                'vector',true,'PersCurr',p.Results.PersCurr,'dPhidt',dPhidt(ii));
        else
            rfuncpot = @(t)Ring_RadExp_erf(t,Ri(ii),Rf(ii),t_phn1(ii),t_exp(ii)/2);
            ANfunc = @(t)atomNumber(t,t_phn1(ii),tau(ii),Ni(ii),Nf(ii));
            [set_temp, amp_temp] =  RingDy_phnEv_2D(TS_FIT(:),THS_FIT(:),rfuncpot,ANfunc,nlobes,...
                A(ii),cTHETA0,Qi(ii),Qf(ii),ph_t(ii),ph_az(ii),...
                'vector',true,'Gamma',ALPHA,'Hub',HUB,...
                'PersCurr',false,'dPhidt',0);
        end
        results2D(:,StartPoints(ii):StopPoints(ii)) = ...
            reshape(set_temp,[size(phi,1),StopPoints(ii)-StartPoints(ii)+1]); 
        amp_temp = reshape(amp_temp,[size(phi,1),StopPoints(ii)-StartPoints(ii)+1]);
        amp_temp = amp_temp(1,:)';
        amp(StartPoints(ii):StopPoints(ii)) = amp_temp;
        clear amp_temp set_temp
    end
    if vector, results2D = results2D(:); end
end

function [results] = atomNumber_all(X,tau,t_start,Ni,Nf)
    % Inputs:
    %   X: 2D array, 1st col is time [s], 2nd col is index referring to
    %       dataset.
    %   tau: decay time [s]
    %   t_start, Ni and Nf: 1D arrays with a value for each
    %       dataset
    % Outputs:
    %   results: the AN [k] as a fucntion of time.
    
    % Initialize ===========================================================
    Time = X(:,1);
    Idx = X(:,2); 
    results = zeros(size(Time,1),1);
    % Check Inputs ========================================================
    if ( size(t_start,1)~=size(Ni,1) || size(Ni,1)~=size(Nf,1))
        error('atomNumber_all: Size of t_start, Ni and Nf should be same as total number of files.');
    else
        N_tot = size(t_start,1);
    end
    if (length(unique(Idx))~=N_tot)
        error('atomNumber_all: The number of unique indices should be same as the size of Ni.');
    end
    % Arrange the data for analysis =======================================
    [Idx, sort_idx] = sort(Idx);
    Time = Time(sort_idx);
    StartPoints = zeros(N_tot,1);
    StopPoints = zeros(N_tot,1);       
    StartPoints(1) = 1; StopPoints(end) = size(Time,1);
    test = Idx(1,1); kk = 2;        
    for ii = 1:1:size(Idx,1)
        if  Idx(ii,1)~= test
            test = Idx(ii,1);
            StartPoints(kk) = ii;
            StopPoints(kk-1) = ii-1;
            kk = kk + 1;
        end
    end
    % Call the AN fit function ==========================================
    for ii = 1:1:N_tot   
        results(StartPoints(ii):StopPoints(ii)) = ...
            atomNumber(Time(StartPoints(ii):StopPoints(ii)),t_start(ii),tau,Ni(ii),Nf(ii));
    end
end

    