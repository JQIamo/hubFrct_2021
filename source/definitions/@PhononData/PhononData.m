classdef PhononData
    properties
      % Raw Data Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      dataRingType          % Type of rings ['Initial' or 'Final' or 'Intermediate' or 'All']
      RunNos_ring_I         % 1D array of normalization initial ring run numbers, e.g., '1:100'
      RunNos_ring_Exp       % 1D array of normalization intermediate ring run numbers, e.g., '1:100'
      RunNos_ring_F         % 1D array of normalization final ring run numbers, e.g., '1:100'
      RunNos_phnEv          % 1D array of phonon run numbers, e.g., '1:100'
      dateFile              % date of images
      basepath              % basepath
      
      % Data Analysis Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      CropX                 % ROI along the column
      CropY                 % ROI along the row
      PrRe                  % Probe Reconst. Type ['None' or 'PCA' or 'GS']
      PCAbasisSet           % a 'PCAset' object containing PCA basis
      GSbasisSet            % a 'PCAset' object containing PCA basis
      N_ev                  % Number of PCs used to evaluate OD
      Mask4PrRe             % Mask for Probe Reconst. (0 in ROI)
      Mask4Ring             % Mask for Ring (1 in ROI)
      I_sat                 % I_sat [counts/pix]
      fitMethod             % Data Fitting method for 1D fits [Choose from 'FitOfAvg' and 'AvgOfFit']
      ChangingVar           % The Variable that changes in each run
      t_extra               % Extra time to be added manually to Evolution time [s]
      n_r                   % The exponent for radial potential
      
      % Save/ Display Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SaveFitRing           % Save the fitted ring image? [logical]
      Saven1D               % Save the n1D image? [logical]

      % Extracted Info %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      AzAngle               % 1D row vec of Azimuthal Angle [radians]
      EvTime                % 1D row vec evolution time [s]
      t_phn1                % The phonon evolution time in the initial ring [s]
      t_exp                 % The duration of expansion of the ring [s]
      AN_phnEv_avg          % 1D array of Atom Num of rings with phonon - Average for each time slice
      AN_phnEv_dev          % 1D array of Atom Num of rings with phonon - Standard dev for each time slice
      MaxOD_phnEv_avg       % 1D array of max OD of rings with phonon - Average for each time slice
      MaxOD_phnEv_dev       % 1D array of max OD of rings with phonon - Err for each time slice
      AN_ringI_avg          % Single value of average Atom Num of normalization rings - INITIAL ring  
      AN_ringI_dev          % Single value of the std dev of Atom Num of normalization rings - INITIAL ring
      AN_ringF_avg          % Single value of average Atom Num of normalization rings - INITIAL ring  
      AN_ringF_dev          % Single value of the std dev of Atom Num of normalization rings - INITIAL ring
      MeanRad_avg           % 1D array of mean radius of rings with phonon - Average for each time slice
      MeanRad_dev           % 1D array of mean radius of rings with phonon - Standard dev for each time slice
      MeanThk_avg           % 1D array of mean radius of rings with phonon - Average for each time slice
      MeanThk_dev           % 1D array of mean radius of rings with phonon - Standard dev for each time slice
      n1D_phnEv_avg         % 2D array of n1D(\theta,t) of phonons - Average for each time slice 
      n1D_phnEv_dev         % 2D array of n1D(\theta,t) of phonons - Standard dev for each time slice
      n1D_ringI_avg         % 1D array of n1D(\theta) of initial bare ring - Average for each time slice 
      n1D_ringI_dev         % 1D array of n1D(\theta) of initial bare ring - Standard dev for each time slice
      n1D_ringF_avg         % 1D array of n1D(\theta) of final bare ring - Average for each time slice 
      n1D_ringF_dev         % 1D array of n1D(\theta) of final bare ring - Standard dev for each time slice
      n1D_ringExp_avg       % 1D array of n1D(\theta,t) of the expanding ring - Average for each time slice 
      n1D_ringExp_dev       % 1D array of n1D(\theta,t) of the expanding ring - Standard dev for each time slice
      n1D_ringExpTimes      % 1D array of times corresponding to the expanding ring

      % Fit Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      AzPhase               % Azimuthal Phase [rad]
      dAzPhase_dt           % Azimuthal Phase Shift Gradient [rad] - val
      dAzPhase_dt_err       % Azimuthal Phase Shift Gradient [rad] - err
      TemporalPhasei        % Azimuthal Phase [rad]
      TemporalPhasef        % Azimuthal Phase [rad]
      deltan_avg            % 1D array of amplitude for each time slice [um^-3] 
      deltan_dev            % 1D array of amplitude conf int for each time slice [um^-3]
      Ampi                  % Amplitude of decaying sinusoid [um^-3]
      Ampi_err              % Amplitude of decaying sinusoid [um^-3]
      ci                    % Initial Speed of sound [um/s] - val
      ci_err                % Initial Speed of sound [um/s] - err
      fi                    % Frequency of decaying sinusoid [Hz]
      fi_err                % Frequency of decaying sinusoid [Hz]
      Qi                    % Q of decaying sinusoid
      Qi_err                % Q of decaying sinusoid
      Ampf                  % Amplitude of decaying sinusoid [um^-3]
      Ampf_err              % Amplitude of decaying sinusoid [um^-3]
      ff                    % Frequency of decaying sinusoid [Hz]
      ff_err                % Frequency of decaying sinusoid [Hz]
      Qf                    % Q of decaying sinusoid
      Qf_err                % Q of decaying sinusoid
      alpha                 % geometric factor
      hub                   % hubble friction co-effcient
      alpha_err                 % geometric factor
      hub_err                   % hubble friction co-effcient

    end

    %% Methods
    methods
        %% ================================================================
        % Constructor =====================================================
        function obj = PhononData(dateFile,RunNos_phnEv,RunNos_ring_I,RunNos_ring_Exp,...
                RunNos_ring_F,varargin)
            % dateFile == YYYY/MM/DD
            % RunNo == 1D array of run numbers, e.g., '1:23'
            % Required Inputs =============================================
            obj.dateFile = dateFile;
            obj.RunNos_phnEv = RunNos_phnEv;
            obj.RunNos_ring_I = RunNos_ring_I;
            obj.RunNos_ring_Exp = RunNos_ring_Exp;
            obj.RunNos_ring_F = RunNos_ring_F;

            % Optional Inputs =============================================            
            p = inputParser;
            p.addParameter('basepath','\\systemadministr\data\raw',@(x)ischar(x));
            p.addParameter('CropX',[1,512],@(x)isnumeric(x));
            p.addParameter('CropY',[1,512],@(x)isnumeric(x));
            addParameter(p,'PrRe','None',@(x)(ischar(x) && ...
                (strcmp(x,'PCA') || strcmp(x,'GS') || strcmp(x,'None'))));
            p.addParameter('PCAbasisSet',[],@(x)isa(x,'PCAset'));
            p.addParameter('GSbasisSet',[],@(x)isa(x,'GSset'));
            p.addParameter('N_ev',[],@isnumeric);
            p.addParameter('Mask4PrRe',[],@isnumeric);
            p.addParameter('Mask4Ring',[],@isnumeric);
            p.addParameter('I_sat',2410,@(x)isnumeric(x));
            p.addParameter('fitMethod','FitOfAvg',@(x)ischar(x));
            p.addParameter('n_r',2,@(x)isnumeric(x));
            p.addParameter('t_extra',0,@(x)isnumeric(x));
            p.addParameter('SaveFitRing',true,@(x)islogical(x));
            p.addParameter('Saven1D',true,@(x)islogical(x));
            p.addParameter('t_phn1_override',[],@isnumeric);
            p.addParameter('t_exp_override',[],@isnumeric);
            parse(p,varargin{:}); 
           
            obj.basepath = p.Results.basepath;            
            obj.CropX = p.Results.CropX;
            obj.CropY = p.Results.CropY;
            obj.PrRe = p.Results.PrRe;
            obj.PCAbasisSet = p.Results.PCAbasisSet;
            obj.GSbasisSet = p.Results.GSbasisSet;
            obj.N_ev = p.Results.N_ev;
            obj.I_sat = p.Results.I_sat;
            obj.fitMethod = p.Results.fitMethod;
            obj.n_r = p.Results.n_r;
            obj.t_extra = p.Results.t_extra;
            obj.Mask4PrRe = p.Results.Mask4PrRe;
            obj.Mask4Ring = p.Results.Mask4Ring;
            obj.SaveFitRing = p.Results.SaveFitRing;
            obj.Saven1D = p.Results.Saven1D;
            obj.t_phn1 = p.Results.t_phn1_override;
            obj.t_exp = p.Results.t_exp_override;
                     
        end % =============================================================
        %% ================================================================
        % Function to determine Az angle ==================================
        % Inputs:
        %   obj: The 'PhnData' class object
        function obj = DetAzAngle(obj) % ==================================             
            A = DataExp(obj.dateFile,obj.RunNos_phnEv(1),'Andor',obj.basepath);
            if strcmp(obj.PrRe,'None')
                [~,~,~,~,~,~,obj.AzAngle,~] = Ring_ExtractData_FromRuns...
                    (A,obj.CropX,obj.CropY,'IntFlcCorr',2,'I_sat',obj.I_sat,...
                    'PrRe',obj.PrRe,'Mask4PrRe',obj.Mask4PrRe,...
                    'Mask4Ring',obj.Mask4Ring,...
                    'n_r',obj.n_r,...
                    'PlotFitRing',false,'Plotn1D',false,...
                    'SaveFitRing',false,'Saven1D',false);
            elseif strcmp(obj.PrRe,'PCA')
                [~,~,~,~,~,~,obj.AzAngle,~] = Ring_ExtractData_FromRuns...
                    (A,obj.CropX,obj.CropY,'I_sat',obj.I_sat,...
                    'PrRe',obj.PrRe,'PCAsetValue',obj.PCAbasisSet,'N_ev',obj.N_ev,'Mask4PrRe',obj.Mask4PrRe,...
                    'Mask4Ring',obj.Mask4Ring,...
                    'n_r',obj.n_r,...
                    'PlotFitRing',false,'Plotn1D',false,...
                    'SaveFitRing',false,'Saven1D',false);
            elseif strcmp(obj.PrRe,'GS')
                [~,~,~,~,~,~,obj.AzAngle,~] = Ring_ExtractData_FromRuns...
                    (A,obj.CropX,obj.CropY,'I_sat',obj.I_sat,...
                    'PrRe',obj.PrRe,'GSsetValue',obj.GSbasisSet,'N_ev',obj.N_ev,'Mask4PrRe',obj.Mask4PrRe,...
                    'Mask4Ring',obj.Mask4Ring,...
                    'n_r',obj.n_r,...
                    'PlotFitRing',false,'Plotn1D',false,...
                    'SaveFitRing',false,'Saven1D',false);
            end  
          if size(obj.AzAngle,2) > 1, obj.AzAngle = obj.AzAngle'; end
        end % =============================================================
        %% ================================================================
        % Function to determine the evolution time ========================
        % Inputs:
        %   obj: The 'PhnData' class object
        function obj = DetEvTime(obj) % ===================================
            params = extract_parameters_erna(obj.dateFile,obj.RunNos_phnEv,...
                {'t_ev','t_phn1','t_exp'},'BasePath',obj.basepath);
            R1 = range(params(:,1)); R2 = range(params(:,2)); R3 = range(params(:,3));
            % Find t_phn1 and t_exp if overrides are empty ========================
            if (isempty(obj.t_phn1) && R2 ~= 0)
                error('PhnData:: DetEvTime: All Runs do not have the same t_phn1.');
            elseif isempty(obj.t_phn1)
                obj.t_phn1 = params(1,2);
            end
            if (isempty(obj.t_exp) && R3 ~= 0)
                error('PhnData:: DetEvTime: All Runs do not have the same t_exp.');
            elseif isempty(obj.t_exp)
                obj.t_exp = params(1,3);
            end
            if (R1 == 0)
                error('PhnData:: DetEvTime: All Runs have the same evolution time t_ev.');
            end
            evTime = params(:,1)+obj.t_extra;
            obj.ChangingVar = 't_ev';
            obj.EvTime = unique(evTime);
            if size(obj.EvTime,2) > 1, obj.EvTime = obj.EvTime'; end
            
            if (max(obj.EvTime)<obj.t_phn1)
                obj.dataRingType = 'Initial';
            elseif (min(obj.EvTime)>obj.t_phn1+obj.t_exp)
                obj.dataRingType = 'Final';
            elseif (max(obj.EvTime)>obj.t_phn1+obj.t_exp && min(obj.EvTime)<obj.t_phn1)
                obj.dataRingType = 'All';
            end   
            
        end % =============================================================
        %% ================================================================
        % Function to determine the ODs and n1Ds of all runs ==============
        % Inputs:
        %   obj: The 'PhnData' class object
        function [Runs,time,AN,MeanRad,Thk,n1Ds,ODs_masked,ODs_ringFit,MaxOD] = Det_n1Ds(obj)            
            Runs = [obj.RunNos_phnEv, obj.RunNos_ring_I, obj.RunNos_ring_Exp, obj.RunNos_ring_F];
            if isempty(obj.AzAngle)
                obj = DetAzAngle(obj);          
            end
            % Find MeanRad, Thk, ODs, n1Ds and AN =========================
            A = DataExp(obj.dateFile,Runs,'Andor',obj.basepath);
            if strcmp(obj.PrRe,'None')
                [ODs_masked,~,~,MeanRad,Thk,ODs_ringFit,~,n1Ds,AN] = Ring_ExtractData_FromRuns...
                    (A,obj.CropX,obj.CropY,'IntFlcCorr',2,'I_sat',obj.I_sat,...
                    'PrRe',obj.PrRe,'Mask4PrRe',obj.Mask4PrRe,...
                    'Mask4Ring',obj.Mask4Ring,...
                    'n_r',obj.n_r,...
                    'PlotFitRing',false,'Plotn1D',false,...
                    'SaveFitRing',obj.SaveFitRing,'Saven1D',obj.Saven1D);
            elseif strcmp(obj.PrRe,'PCA')
                [ODs_masked,~,~,MeanRad,Thk,ODs_ringFit,~,n1Ds,AN] = Ring_ExtractData_FromRuns...
                    (A,obj.CropX,obj.CropY,'I_sat',obj.I_sat,...
                    'PrRe',obj.PrRe,'PCAsetValue',obj.PCAbasisSet,'N_ev',obj.N_ev,'Mask4PrRe',obj.Mask4PrRe,...
                    'Mask4Ring',obj.Mask4Ring,...
                    'n_r',obj.n_r,...
                    'PlotFitRing',false,'Plotn1D',false,...
                    'SaveFitRing',obj.SaveFitRing,'Saven1D',obj.Saven1D);
            elseif strcmp(obj.PrRe,'GS')
                [ODs_masked,~,~,MeanRad,Thk,ODs_ringFit,~,n1Ds,AN] = Ring_ExtractData_FromRuns...
                    (A,obj.CropX,obj.CropY,'I_sat',obj.I_sat,...
                    'PrRe',obj.PrRe,'GSsetValue',obj.GSbasisSet,'N_ev',obj.N_ev,'Mask4PrRe',obj.Mask4PrRe,...
                    'Mask4Ring',obj.Mask4Ring,...
                    'n_r',obj.n_r,...
                    'PlotFitRing',false,'Plotn1D',false,...
                    'SaveFitRing',obj.SaveFitRing,'Saven1D',obj.Saven1D);
            end
            % Find Max OD =================================================
            MaxOD = permute(max(max(ODs_ringFit)),[3 2 1]);
            % Find t_ev ===================================================
            if isempty(obj.ChangingVar)
                obj = DetEvTime(obj);
                time = extract_parameters_erna(obj.dateFile,Runs,...
                  {obj.ChangingVar},'BasePath',obj.basepath)+obj.t_extra;
            elseif strcmp(obj.ChangingVar,'DMD_F2')
                t_expTime = extract_parameters_erna(obj.dateFile,obj.RunNos_phnEv,...
              {'t_exp'},'BasePath',obj.basepath);
                time = extract_parameters_erna(obj.dateFile,Runs,...
                  {obj.ChangingVar},'BasePath',obj.basepath).*t_expTime(1)/90+obj.t_extra;
            else
                time = extract_parameters_erna(obj.dateFile,Runs,...
                  {obj.ChangingVar},'BasePath',obj.basepath)+obj.t_extra;
            end           
            Runs = Runs';
        end % =============================================================
        %% ================================================================
        % Function to sort n1Ds to extract n1D_ring for normalization =====
        % Inputs:
        %   obj: The 'PhnData' class object
        %   Runs,times,AN,n1D: The outputs from Det_n1Ds
        % Outputs:
        %   obj: return the object with n1D_ringI, n1D_ringF, n1D_ringExp and AN_ringI
        %       values saved into the object.
        function [obj] = Sort_n1D_ring_fromRuns(obj,Runs,time,AN,n1D) 
                % Sort the initial ring =======================================
                if ~isempty(obj.RunNos_ring_I)
                    Runs_ring = obj.RunNos_ring_I';
                    kk = 1;
                    for ii = 1:1:size(Runs,1)
                        for jj = 1:1:size(Runs_ring,1)
                            if Runs(ii) == Runs_ring(jj)
                                AN_temp(kk) = AN(ii);
                                n1D_temp(:,kk) = n1D(:,ii);
                                kk = kk+1;
                            end
                        end
                    end
                    obj.n1D_ringI_avg = nanmean(n1D_temp,2);
                    obj.AN_ringI_avg = nanmean(AN_temp,2);
                    if kk == 2
                        obj.n1D_ringI_dev = nanstd(n1D_temp,0,2);
                        obj.AN_ringI_dev = nanstd(AN_temp,0,2);
                    else
                        obj.n1D_ringI_dev = nanstd(n1D_temp,0,2)/sqrt(kk-2);
                        obj.AN_ringI_dev = nanstd(AN_temp,0,2)/sqrt(kk-2);
                    end
                else
                    obj.n1D_ringI_avg = mean(obj.n1D_phnEv_avg(:,obj.EvTime<=obj.t_phn1),2);
                    obj.n1D_ringI_dev = std(obj.n1D_phnEv_avg(:,obj.EvTime<=obj.t_phn1),0,2);
                    obj.AN_ringI_avg = mean(obj.AN_phnEv_avg(obj.EvTime<=obj.t_phn1));
                    obj.AN_ringI_dev = std(obj.AN_phnEv_avg(obj.EvTime<=obj.t_phn1))...
                        /sqrt(size(obj.AN_phnEv_avg(obj.EvTime<=obj.t_phn1),1)-1);
                end                    
                clear AN_temp n1D_temp % ======================================
                % Sort the final Ring =========================================
                if ~isempty(obj.RunNos_ring_F)
                    Runs_ring = obj.RunNos_ring_F';
                    kk = 1;
                    for ii = 1:1:size(Runs,1)
                        for jj = 1:1:size(Runs_ring,1)
                            if Runs(ii) == Runs_ring(jj)
                                AN_temp(kk) = AN(ii);
                                n1D_temp(:,kk) = n1D(:,ii);
                                kk = kk+1;
                            end
                        end
                    end
                    obj.n1D_ringF_avg = nanmean(n1D_temp,2);
                    obj.AN_ringF_avg = nanmean(AN_temp,2);
                    if kk == 2
                        obj.n1D_ringF_dev = nanstd(n1D_temp,0,2);
                        obj.AN_ringF_dev = nanstd(AN_temp,0,2);
                    else
                        obj.n1D_ringF_dev = nanstd(n1D_temp,0,2)/sqrt(kk-2);
                        obj.AN_ringF_dev = nanstd(AN_temp,0,2)/sqrt(kk-2);
                    end
                else
                    obj.n1D_ringF_avg = mean(obj.n1D_phnEv_avg(:,obj.EvTime>=obj.t_phn1+obj.t_exp),2);
                    obj.n1D_ringF_dev = std(obj.n1D_phnEv_avg(:,obj.EvTime>=obj.t_phn1+obj.t_exp),0,2);
                    obj.AN_ringF_avg = mean(obj.AN_phnEv_avg(obj.EvTime>=obj.t_phn1+obj.t_exp));
                    obj.AN_ringF_dev = std(obj.AN_phnEv_avg(obj.EvTime>=obj.t_phn1+obj.t_exp))...
                        /sqrt(size(obj.AN_phnEv_avg(obj.EvTime>=obj.t_phn1+obj.t_exp),1)-1);
                end
                clear AN_temp n1D_temp % ======================================
                % Sort the expanding Ring =====================================
                if ~isempty(obj.RunNos_ring_Exp)
                    Runs_ring = obj.RunNos_ring_Exp';  
                    ExpTimes = obj.EvTime(obj.EvTime>obj.t_phn1*1.01 & obj.EvTime<(obj.t_phn1+obj.t_exp)*0.99);
                    obj.n1D_ringExp_avg = zeros(size(obj.AzAngle,1),size(ExpTimes,1));
                    obj.n1D_ringExp_dev = zeros(size(obj.AzAngle,1),size(ExpTimes,1));
                    for ll = 1:1:size(ExpTimes,1)
                        kk = 1;
                        for ii = 1:1:size(Runs,1)
                            for jj = 1:1:size(Runs_ring,1)
                                if (Runs(ii) == Runs_ring(jj) && time(ii)== ExpTimes(ll))
                                    n1D_temp(:,kk) = n1D(:,ii);
                                    kk = kk+1;
                                end
                            end
                        end
                        if ~exist('n1D_temp','var')
                            error('PhononData::Sort_n1D_ring: Norm ring for t_ev = %0.2f ms missing.',ExpTimes(ll)*1e3);
                        end
                        obj.n1D_ringExp_avg(:,ll) = nanmean(n1D_temp,2);
                        obj.n1D_ringExpTimes(ll) = ExpTimes(ll);
                        if kk == 2
                            obj.n1D_ringExp_dev(:,ll) = nanstd(n1D_temp,0,2);
                        else
                            obj.n1D_ringExp_dev(:,ll) = nanstd(n1D_temp,0,2)/sqrt(kk-2);
                        end
                        clear n1D_temp
                    end
                end               
        end % =============================================================
        %% Function to sort n1Ds to extract n1D_phnEv =====================
        % Inputs:
        %   obj: The 'PhnData' class object
        %   Runs,AN,MeanRad,MeanThk,n1D: The outputs from Det_n1Ds
        % Outputs:
        %   obj: The object with n1D_phnEv_avg and n1D_phnEv_dev saved.
        function [obj] = Sort_n1D_phnEv(obj,Runs,time,AN,MeanRad,MeanThk,n1D,MaxODs)   
            if isempty(obj.EvTime)  
                obj = DetEvTime(obj);           
            end
            % Initialize ==================================================
            obj.n1D_phnEv_avg = zeros(size(n1D,1),size(obj.EvTime,1));
            obj.AN_phnEv_avg = zeros(size(obj.EvTime,1),1);
            obj.AN_phnEv_dev = zeros(size(obj.EvTime,1),1);
            obj.MeanRad_avg = zeros(size(obj.EvTime,1),1);
            obj.MeanRad_dev = zeros(size(obj.EvTime,1),1);  
            obj.MeanThk_avg = zeros(size(obj.EvTime,1),1);
            obj.MeanThk_dev = zeros(size(obj.EvTime,1),1);
            obj.MaxOD_phnEv_avg = zeros(size(obj.EvTime,1),1);
            obj.MaxOD_phnEv_dev = zeros(size(obj.EvTime,1),1);

            Runs_phnEv = obj.RunNos_phnEv';        
            for ll = 1:1:size(obj.EvTime,1)
                kk = 1;
                for ii = 1:1:size(Runs,1)
                    for jj = 1:1:size(Runs_phnEv,1)                   
                        if (time(ii) == obj.EvTime(ll) && Runs(ii) == Runs_phnEv(jj)) 
                            AN_temp(kk) = AN(ii);
                            MeanRad_temp(kk) = MeanRad(ii);
                            MeanThk_temp(kk) = MeanThk(ii);
                            n1D_temp(:,kk) = n1D(:,ii);
                            MaxOD_temp(kk) = MaxODs(ii);
                            kk = kk+1;
                        end
                    end
                end
                obj.n1D_phnEv_avg(:,ll) = nanmean(n1D_temp,2);
                obj.AN_phnEv_avg(ll) = nanmean(AN_temp,2);
                obj.MeanRad_avg(ll) = nanmean(MeanRad_temp,2);
                obj.MeanThk_avg(ll) = nanmean(MeanThk_temp,2);
                obj.MaxOD_phnEv_avg(ll) = nanmean(MaxOD_temp,2);
                if kk == 2
                    obj.n1D_phnEv_dev(:,ll) = nanstd(n1D_temp,0,2);
                    obj.AN_phnEv_dev(ll) = nanstd(AN_temp,0,2);
                    obj.MeanRad_dev(ll) = nanstd(MeanRad_temp,0,2);
                    obj.MeanThk_dev(ll) = nanstd(MeanThk_temp,0,2);
                    obj.MaxOD_phnEv_dev(ll) = nanstd(MaxOD_temp,0,2);
                else
                    obj.n1D_phnEv_dev(:,ll) = nanstd(n1D_temp,0,2)/sqrt(kk-2);
                    obj.AN_phnEv_dev(ll) = nanstd(AN_temp,0,2)/sqrt(kk-2);
                    obj.MeanRad_dev(ll) = nanstd(MeanRad_temp,0,2)/sqrt(kk-2);
                    obj.MeanThk_dev(ll) = nanstd(MeanThk_temp,0,2)/sqrt(kk-2);
                    obj.MaxOD_phnEv_dev(ll) = nanstd(MaxOD_temp,0,2)/sqrt(kk-2);
                end               
                clear AN_temp MeanRad_temp n1D_temp
            end    
        end % =============================================================
        %% Function to normalize the n1D_phnEv ============================
        % Inputs:
        %   obj: The 'PhononData' class object
        %   Runs,AN,MeanRad,n1D: The outputs from Det_n1Ds
        function [obj] = Norm_n1D_phnEv(obj) % ============================
                for ii = 1:1:size(obj.n1D_phnEv_avg,2)
                    if obj.EvTime(ii) <= obj.t_phn1
                        n1D_ring_avg = obj.n1D_ringI_avg;
                        n1D_ring_dev = obj.n1D_ringI_dev;
                    elseif obj.EvTime(ii) >= (obj.t_phn1 + obj.t_exp)
                        n1D_ring_avg = obj.n1D_ringF_avg;
                        n1D_ring_dev = obj.n1D_ringF_dev;
                    else
                        for jj = 1:1:size(obj.n1D_ringExpTimes,2)
                            if obj.EvTime(ii) == obj.n1D_ringExpTimes(jj)
                                n1D_ring_avg = obj.n1D_ringExp_avg(:,jj);
                                n1D_ring_dev = obj.n1D_ringExp_dev(:,jj);
                            end
                        end
                    end
                    if (isempty(n1D_ring_avg) || isempty(n1D_ring_dev))
                        error('PhononData:Norm_n1D_phnEv:: Normalization Bare ring is missing.');
                    end  
                    % Normalization: Subtract the bare ring
                    alpha =  sum(obj.n1D_phnEv_avg(:,ii))/sum(n1D_ring_avg);
                    obj.n1D_phnEv_avg(:,ii) = obj.n1D_phnEv_avg(:,ii) - alpha.*n1D_ring_avg;
                    obj.n1D_phnEv_dev = sqrt(obj.n1D_phnEv_dev.^2 + (alpha.*n1D_ring_dev).^2);
                end             
        end % =============================================================
        %% Function to determine n1D_phnEv ================================
        % Inputs:
        %   obj: The 'PhononData' class object
        %   Runs,AN,MeanRad,n1D: The outputs from Det_n1Ds
        function [obj] = Det_n1D_phnEv(obj,Runs,time,AN,MeanRad,Thk,n1D,MaxODs)
            obj = Sort_n1D_phnEv(obj,Runs,time,AN,MeanRad,Thk,n1D,MaxODs);
            obj = Sort_n1D_ring_fromRuns(obj,Runs,time,AN,n1D);
            obj = Norm_n1D_phnEv(obj);
        end % =======================================================================        
        %% Function to normalize all n1Ds =================================
        % Inputs:
        %   obj: The 'PhnData' class object
        %   n1Ds: 2D array of n1Ds to be normalized with ring n1D
        %   time: the corresponding values of time
        %       Both n1Ds and time are obtained from the function Det_n1Ds
        % Outputs
        %   n1Ds: normalized n1Ds
        function [n1Ds] = Norm_n1Ds(obj,n1Ds,time) % ======================           
            % Check Inputs
            if size(n1Ds,2)~=size(time,1)
                error('PhononData::Norm_n1Ds: size(n1Ds,2) needs to be same as size(time,1).');
            end
            if size(n1Ds,1)~= size(obj.AzAngle,1)
                error('PhononData::Norm_n1Ds:size(n1Ds,1) needs to be same as size(obj.AzAngle,1)');
            end                                
            for ii = 1:1:size(n1Ds,2)
                    if time(ii) <= obj.t_phn1
                        n1D_ring_avg = obj.n1D_ringI_avg;
                        n1D_ring_dev = obj.n1D_ringI_dev;
                    elseif time(ii) >= (obj.t_phn1 + obj.t_exp)
                        n1D_ring_avg = obj.n1D_ringF_avg;
                        n1D_ring_dev = obj.n1D_ringF_dev;
                    else
                        for jj = 1:1:size(obj.n1D_ringExpTimes,2)
                            if time(ii) == obj.n1D_ringExpTimes(jj)
                                n1D_ring_avg = obj.n1D_ringExp_avg(:,jj);
                                n1D_ring_dev = obj.n1D_ringExp_dev(:,jj);
                            end
                        end
                    end                   
                if (isempty(n1D_ring_avg) || isempty(n1D_ring_dev))                    
                    error('PhononData:Norm_n1Ds:: Ev. time = %0.2f ms encountered. Normalization data is missing.',time(ii));
                end
                % Normalization: Subtraction
                alpha =  sum(n1Ds(:,ii))/sum(n1D_ring_avg);
                n1Ds(:,ii) = n1Ds(:,ii) - alpha.*n1D_ring_avg;
            end                      
        end % =============================================================
        %% Function to sort through all n1Ds to extract deltan values ======
        % Inputs:
        %   obj: The 'PhnData' class object
        %   Runs,time: The outputs from Det_n1Ds
        %   a: The corresponding amplitudes of 1D sine fit to n1Ds
        % Outputs:
        %   obj: The object is returned with deltan_avg and deltan_dev
        %       values saved.
        function [obj] = Sort_n1DfitAmp(obj,Runs,time,a) % =============       
            if isempty(obj.EvTime)  
                obj = DetEvTime(obj);           
            end  
            obj.deltan_avg = zeros(size(obj.EvTime,1),1);
            obj.deltan_dev = zeros(size(obj.EvTime,1),1); 
            Runs_phnEv = obj.RunNos_phnEv';        
            for ll = 1:1:size(obj.EvTime,1)
                kk = 1;
                for ii = 1:1:size(Runs,1)
                    for jj = 1:1:size(Runs_phnEv,1)                   
                        if (time(ii) == obj.EvTime(ll) && Runs(ii) == Runs_phnEv(jj)) 
                            a_temp(kk) = a(ii);
                            kk = kk+1;
                        end
                    end
                end
                obj.deltan_avg(ll) = nanmean(a_temp,2);
                if kk == 2
                    obj.deltan_dev(ll) = nanstd(a_temp,0,2);
                else
                    obj.deltan_dev(ll) = nanstd(a_temp,0,2)/sqrt(kk-2);
                end
                clear a_temp 
            end    
        end % =======================================================================
    end
end
