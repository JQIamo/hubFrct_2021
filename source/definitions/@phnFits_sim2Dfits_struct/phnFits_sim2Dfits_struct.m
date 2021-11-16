classdef phnFits_sim2Dfits_struct
    properties
      % Raw Data Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      file_names            % file names of all time-traces
      file_names_Dy         % file names of all dynamic time-traces
      file_names_Stat       % file names of all static time-traces
      % Fitting Scheme Employed %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % The following properties can either of the following values.
      % 'I'                 Independent
      % 'S_all'             Shared across all time-traces
      % 'S_all_scaled'      Shared across all time-traces, but scaled with atom number.
      % 'S_Dy'      	    Shared across all dynamic ring time-traces.
      % 'S_Dy+'      	    Shared across all dynamic ring time-traces + the appropriate static ring time-trace
      % 'S_Dy_scaled'	    Shared across all dynamic ring time-traces, but scaled with atom number.
      % 'S_Dy_scaled+'	    Shared across all dynamic ring time-traces + the appropriate static ring time-trace, but scaled with atom number.
      % 'S_Stat'          	Shared across all static ring time-traces.
      % 'S_Stat_scaled'     Shared across all static ring time-traces, but scaled with atom number.
      varType_Dy_ph_t       % 
      varType_Dy_ph_az      %
      varType_Dy_amp        % 
      varType_Dy_Qi         %
      varType_Dy_Qf         %
      varType_Dy_dPhidt     %
      varType_Stat_ph_t     %
      varType_Stat_ph_az    %
      varType_Stat_amp      %
      varType_Stat_Qi       %
      varType_Stat_dPhidt   %
      PersCurr              % true of false
      % Data Analysis Properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      alpha                 % geometric factor
      gamma_H               % hub friction coeff
      cTheta0               % cTheta (N = 100k, R = 40 um)
      ph_t                  % temporal phase
      ph_az                 % azimuthal phase
      amp                   % amplitude of oscillation
      Qi                    % initial quality factor
      Qf                    % final quality factor
      dPhidt                % rate of change of az. phase
      alpha_err                 % ERR geometric factor 
      gamma_H_err               % ERR hub friction coeff 
      cTheta0_err               % ERR cTheta (N = 100k, R = 40 um) 
      ph_t_err                  % ERR temporal phase 
      ph_az_err                 % ERR azimuthal phase 
      amp_err                   % ERR amplitude of oscillation 
      Qi_err                    % ERR initial quality factor 
      Qf_err                    % ERR final quality factor 
      dPhidt_err                % ERR rate of change of az. phase 
      % Number of fit parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      Nfit_Dy_ph_t       % 
      Nfit_Dy_ph_az      %
      Nfit_Dy_amp        % 
      Nfit_Dy_Qi         %
      Nfit_Dy_Qf         %
      Nfit_Dy_dPhidt     %
      Nfit_Stat_ph_t     %
      Nfit_Stat_ph_az    %
      Nfit_Stat_amp      %
      Nfit_Stat_Qi       %
      Nfit_Stat_dPhidt   %
      Ntraces_Dy
      Ntraces_Stat
    end

    %% Methods
    methods
        %% ================================================================
        % Constructor =====================================================
        function obj = phnFits_sim2Dfits_struct(file_names_Dy,file_names_Stat,varargin)
            % file_names_Dy == str of filenames for dynamic rings
            % file_names_Stat == str of filenames for static rings
            % Required Inputs =============================================
            obj.file_names_Dy = file_names_Dy;
            obj.file_names_Stat = file_names_Stat;
            obj.file_names = [obj.file_names_Dy, obj.file_names_Stat];
            % Optional Inputs =============================================
            p = inputParser;
            p.addParameter('varType_Dy_ph_t','I',@(x)(ischar(x) && (strcmp(x,'I')...
                || strcmp(x,'S_Dy') || strcmp(x,'S_Dy+'))));
            p.addParameter('varType_Dy_ph_az','I',@(x)(ischar(x) && (strcmp(x,'I')...
                || strcmp(x,'S_Dy') || strcmp(x,'S_Dy+'))));
            p.addParameter('varType_Dy_amp','I',@(x)(ischar(x) && (strcmp(x,'I')...
                || strcmp(x,'S_Dy') || strcmp(x,'S_Dy+') || ...
                strcmp(x,'S_Dy_scaled') || strcmp(x,'S_Dy_scaled+'))));
            p.addParameter('varType_Dy_Qi','I',@(x)(ischar(x) && (strcmp(x,'I')...
                || strcmp(x,'S_Dy') || strcmp(x,'S_Dy+'))));
            p.addParameter('varType_Dy_Qf','I',@(x)(ischar(x) && (strcmp(x,'I')...
                || strcmp(x,'S_Dy'))));  
            p.addParameter('varType_Dy_dPhidt','I',@(x)(ischar(x) && (strcmp(x,'I')...
                || strcmp(x,'S_Dy') || strcmp(x,'None'))));
            p.addParameter('varType_Stat_ph_t','I',@(x)(ischar(x) && (strcmp(x,'I'))));
            p.addParameter('varType_Stat_ph_az','I',@(x)(ischar(x) && (strcmp(x,'I'))));
            p.addParameter('varType_Stat_amp','I',@(x)(ischar(x) && (strcmp(x,'I'))));
            p.addParameter('varType_Stat_Qi','I',@(x)(ischar(x) && (strcmp(x,'I'))));
            p.addParameter('PersCurr',false,@(x)islogical(x));
            p.addParameter('varType_Stat_dPhidt','I',@(x)(ischar(x) && (strcmp(x,'I'))));
            parse(p,varargin{:});            
            % Set the variable types ======================================
            obj.varType_Dy_ph_t = p.Results.varType_Dy_ph_t;
            obj.varType_Dy_ph_az = p.Results.varType_Dy_ph_az;
            obj.varType_Dy_amp = p.Results.varType_Dy_amp;
            obj.varType_Dy_Qi = p.Results.varType_Dy_Qi;
            obj.varType_Dy_Qf = p.Results.varType_Dy_Qf;
            obj.varType_Stat_ph_t = p.Results.varType_Stat_ph_t;
            obj.varType_Stat_ph_az = p.Results.varType_Stat_ph_az;
            obj.varType_Stat_amp = p.Results.varType_Stat_amp;
            obj.varType_Stat_Qi = p.Results.varType_Stat_Qi;
            obj.PersCurr = p.Results.PersCurr;
            obj.varType_Dy_dPhidt = p.Results.varType_Dy_dPhidt;
            obj.varType_Stat_dPhidt = p.Results.varType_Stat_dPhidt;
            % Evaluate the number of fit parameters =======================
            obj.Ntraces_Dy = length(obj.file_names_Dy);
            obj.Ntraces_Stat = length(obj.file_names_Stat);
            if strcmp(obj.varType_Dy_ph_t,'I')
                obj.Nfit_Dy_ph_t = obj.Ntraces_Dy;
            elseif (strcmp(obj.varType_Dy_ph_t,'S_Dy') || strcmp(obj.varType_Dy_ph_t,'S_Dy+'))
                obj.Nfit_Dy_ph_t = 1;
            end
            if strcmp(obj.varType_Dy_ph_az,'I')
                obj.Nfit_Dy_ph_az = obj.Ntraces_Dy;
            elseif (strcmp(obj.varType_Dy_ph_az,'S_Dy') || strcmp(obj.varType_Dy_ph_az,'S_Dy+'))
                obj.Nfit_Dy_ph_az = 1;
            end
            if strcmp(obj.varType_Dy_amp,'I')
                obj.Nfit_Dy_amp = obj.Ntraces_Dy;
            elseif (strcmp(obj.varType_Dy_amp,'S_Dy') || strcmp(obj.varType_Dy_amp,'S_Dy+')...
                    || strcmp(obj.varType_Dy_amp,'S_Dy_scaled') || strcmp(obj.varType_Dy_amp,'S_Dy_scaled+'))
                obj.Nfit_Dy_amp = 1;
            end
            if strcmp(obj.varType_Dy_Qi,'I')
                obj.Nfit_Dy_Qi = obj.Ntraces_Dy;
            elseif (strcmp(obj.varType_Dy_Qi,'S_Dy') || strcmp(obj.varType_Dy_Qi,'S_Dy+'))
                obj.Nfit_Dy_Qi = 1;
            end
            if strcmp(obj.varType_Dy_Qf,'I')
                obj.Nfit_Dy_Qf = obj.Ntraces_Dy;
            elseif (strcmp(obj.varType_Dy_Qf,'S_Dy') || strcmp(obj.varType_Dy_Qf,'S_Dy+'))
                obj.Nfit_Dy_Qf = 1;
            end            
            if strcmp(obj.varType_Dy_ph_t,'S_Dy+')
                obj.Nfit_Stat_ph_t = obj.Ntraces_Stat-1;
            else 
                obj.Nfit_Stat_ph_t = obj.Ntraces_Stat;
            end
            if strcmp(obj.varType_Dy_ph_az,'S_Dy+')
                obj.Nfit_Stat_ph_az = obj.Ntraces_Stat-1;
            else 
                obj.Nfit_Stat_ph_az = obj.Ntraces_Stat;
            end
            if (strcmp(obj.varType_Dy_amp,'S_Dy+') || strcmp(obj.varType_Dy_amp,'S_Dy_scaled+'))
                obj.Nfit_Stat_amp = obj.Ntraces_Stat-1;
            else 
                obj.Nfit_Stat_amp = obj.Ntraces_Stat;
            end
            if strcmp(obj.varType_Dy_Qi,'S_Dy+')
                obj.Nfit_Stat_Qi = obj.Ntraces_Stat-1;
            else 
                obj.Nfit_Stat_Qi = obj.Ntraces_Stat;
            end 
            if obj.PersCurr
                if (strcmp(obj.varType_Dy_dPhidt,'S_Dy') || strcmp(obj.varType_Dy_dPhidt,'S_Dy+'))
                    obj.Nfit_Dy_dPhidt = 1;
                elseif strcmp(obj.varType_Dy_dPhidt,'None')
                    obj.Nfit_Dy_dPhidt = 0;
                else
                    obj.Nfit_Dy_dPhidt = obj.Ntraces_Dy;
                end
                if strcmp(obj.varType_Dy_dPhidt,'S_Dy+')
                    obj.Nfit_Stat_dPhidt = obj.Ntraces_Stat-1;
                else
                    obj.Nfit_Stat_dPhidt = obj.Ntraces_Stat;
                end
            else
                obj.Nfit_Dy_dPhidt = 0;
                obj.Nfit_Stat_dPhidt = 0;
            end
                
        end % =============================================================
        % Assemble fit values =============================================
        function obj = assblFitVals(obj,beta,AtomNum,ANbase,varargin)
            % Optional Inputs =============================================
            p = inputParser;
            p.addParameter('alphaEqGammaH',false,@(x)islogical(x));
            parse(p,varargin{:}); 
            if p.Results.alphaEqGammaH
                baseParams = 2;
            else
                baseParams = 3;
            end
            % check data ==================================================
            NfitVals = baseParams+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+...
                obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az+ obj.Nfit_Dy_amp+obj.Nfit_Stat_amp+...
                obj.Nfit_Dy_Qi+obj.Nfit_Stat_Qi+obj.Nfit_Dy_Qf+...
                obj.Nfit_Dy_dPhidt+obj.Nfit_Stat_dPhidt;
            if length(beta) ~= NfitVals
                error('phnFits_sim2Dfits_struct::assblFitVals: beta not the right size');
            end
            % initialize ==================================================
            obj.ph_t = zeros(obj.Ntraces_Dy+obj.Ntraces_Stat,1);
            obj.ph_az = zeros(obj.Ntraces_Dy+obj.Ntraces_Stat,1);
            obj.amp = zeros(obj.Ntraces_Dy+obj.Ntraces_Stat,1);
            obj.Qi = zeros(obj.Ntraces_Dy+obj.Ntraces_Stat,1);
            obj.Qf = zeros(obj.Ntraces_Dy,1);
            % assign ======================================================
            obj.gamma_H = beta(1);
            if p.Results.alphaEqGammaH
                obj.alpha = beta(1);
            else
                obj.alpha = beta(2);
            end
            obj.cTheta0 = beta(baseParams);
            if obj.Nfit_Dy_ph_t == 1
                obj.ph_t(1:obj.Ntraces_Dy) = beta(baseParams+1)*ones(obj.Ntraces_Dy,1);
                obj.ph_t(obj.Ntraces_Dy+1) = beta(baseParams+1);
                obj.ph_t(obj.Ntraces_Dy+2:obj.Ntraces_Dy+obj.Ntraces_Stat)...
                    = beta(baseParams+2:baseParams+obj.Ntraces_Stat)';
            else
                obj.ph_t = beta(baseParams+1:baseParams+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t)';
            end
            if obj.Nfit_Dy_ph_az == 1
                obj.ph_az(1:obj.Ntraces_Dy) = beta(baseParams+1+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t)*ones(obj.Ntraces_Dy,1);
                obj.ph_az(obj.Ntraces_Dy+1) = beta(baseParams+1+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t);
                obj.ph_az(obj.Ntraces_Dy+2:obj.Ntraces_Dy+obj.Ntraces_Stat)...
                    = beta(baseParams+2+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t:baseParams+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+obj.Ntraces_Stat)';
            else
                obj.ph_az = beta(baseParams+1+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t:...
                    baseParams+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az)';
            end
            if obj.Nfit_Dy_amp == 1
                obj.amp(1:obj.Ntraces_Dy) = beta(baseParams+1+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+...
                    obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az)*AtomNum(1:obj.Ntraces_Dy)/ANbase;
                obj.amp(obj.Ntraces_Dy+1) = beta(baseParams+1+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+...
                    obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az)*AtomNum(obj.Ntraces_Dy+1)/ANbase;
                obj.amp(obj.Ntraces_Dy+2:obj.Ntraces_Dy+obj.Ntraces_Stat)...
                    = beta(baseParams+2+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az:...
                    baseParams+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az+obj.Ntraces_Stat)';
            else
                obj.amp = beta(baseParams+1+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az:...
                    baseParams+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az+...
                    obj.Nfit_Dy_amp+obj.Nfit_Stat_amp)';
            end
            if obj.Nfit_Dy_Qi == 1
                obj.Qi(1:obj.Ntraces_Dy) = beta(baseParams+1+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+...
                    obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az+obj.Nfit_Dy_amp+obj.Nfit_Stat_amp)...
                    *ones(obj.Ntraces_Dy,1);
                obj.Qi(obj.Ntraces_Dy+1) = beta(baseParams+1+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+...
                    obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az+obj.Nfit_Dy_amp+obj.Nfit_Stat_amp);
                obj.Qi(obj.Ntraces_Dy+2:obj.Ntraces_Dy+obj.Ntraces_Stat)...
                    = beta(baseParams+2+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az+...
                    obj.Nfit_Dy_amp+obj.Nfit_Stat_amp:...
                    baseParams+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az+...
                    obj.Nfit_Dy_amp+obj.Nfit_Stat_amp+obj.Ntraces_Stat)';
            else
                obj.Qi = beta(baseParams+1+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az+...
                    obj.Nfit_Dy_amp+obj.Nfit_Stat_amp:...
                    baseParams+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az+...
                    obj.Nfit_Dy_amp+obj.Nfit_Stat_amp+obj.Nfit_Dy_Qi+obj.Nfit_Stat_Qi)';
            end
            if obj.Nfit_Dy_Qf == 1
                obj.Qf(1:obj.Ntraces_Dy) = beta(baseParams+1+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+...
                    obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az+obj.Nfit_Dy_amp+obj.Nfit_Stat_amp+...
                    obj.Nfit_Dy_Qi+obj.Nfit_Stat_Qi)...
                    *ones(obj.Ntraces_Dy,1);
            else
                obj.Qf = beta(baseParams+1+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az+...
                    obj.Nfit_Dy_amp+obj.Nfit_Stat_amp+obj.Nfit_Dy_Qi+obj.Nfit_Stat_Qi:...
                    baseParams+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az+...
                    obj.Nfit_Dy_amp+obj.Nfit_Stat_amp+obj.Nfit_Dy_Qi+obj.Nfit_Stat_Qi+obj.Nfit_Dy_Qf)';
            end
         
            if (obj.Nfit_Dy_dPhidt+obj.Nfit_Stat_dPhidt) ~= 0
                obj.dPhidt = beta(baseParams+1+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az+...
                    obj.Nfit_Dy_amp+obj.Nfit_Stat_amp+obj.Nfit_Dy_Qi+obj.Nfit_Stat_Qi+obj.Nfit_Dy_Qf:...
                    baseParams+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az+...
                    obj.Nfit_Dy_amp+obj.Nfit_Stat_amp+obj.Nfit_Dy_Qi+obj.Nfit_Stat_Qi+obj.Nfit_Dy_Qf+...
                    obj.Nfit_Stat_dPhidt+obj.Nfit_Dy_dPhidt)';
            end
        end % =============================================================
        % Assemble fit error values =======================================
        function obj = assblFitErrs(obj,beta,AtomNum,ANbase,varargin)
            % Optional Inputs =============================================
            p = inputParser;
            p.addParameter('alphaEqGammaH',false,@(x)islogical(x));
            parse(p,varargin{:}); 
            if p.Results.alphaEqGammaH
                baseParams = 2;
            else
                baseParams = 3;
            end
            % check data ==================================================
            NfitVals = baseParams+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+...
                obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az+ obj.Nfit_Dy_amp+obj.Nfit_Stat_amp+...
                obj.Nfit_Dy_Qi+obj.Nfit_Stat_Qi+obj.Nfit_Dy_Qf+...
                obj.Nfit_Dy_dPhidt+obj.Nfit_Stat_dPhidt;
            if length(beta) ~= NfitVals
                error('phnFits_sim2Dfits_struct::assblFitErrs: ci not the right size');
            end
            % initialize ==================================================
            obj.ph_t_err = zeros(obj.Ntraces_Dy+obj.Ntraces_Stat,1);
            obj.ph_az_err = zeros(obj.Ntraces_Dy+obj.Ntraces_Stat,1);
            obj.amp_err = zeros(obj.Ntraces_Dy+obj.Ntraces_Stat,1);
            obj.Qi_err = zeros(obj.Ntraces_Dy+obj.Ntraces_Stat,1);
            obj.Qf_err = zeros(obj.Ntraces_Dy,1);
            % assign ======================================================
            obj.gamma_H_err = beta(1);
            if p.Results.alphaEqGammaH
                obj.alpha_err = beta(1);
            else
                obj.alpha_err = beta(2);
            end
            obj.cTheta0_err = beta(baseParams);
            if obj.Nfit_Dy_ph_t == 1
                obj.ph_t_err(1:obj.Ntraces_Dy) = beta(baseParams+1)*ones(obj.Ntraces_Dy,1);
                obj.ph_t_err(obj.Ntraces_Dy+1) = beta(baseParams+1);
                obj.ph_t_err(obj.Ntraces_Dy+2:obj.Ntraces_Dy+obj.Ntraces_Stat)...
                    = beta(baseParams+2:baseParams+obj.Ntraces_Stat)';
            else
                obj.ph_t_err = beta(baseParams+1:baseParams+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t)';
            end
            if obj.Nfit_Dy_ph_az == 1
                obj.ph_az_err(1:obj.Ntraces_Dy) = beta(baseParams+1+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t)*ones(obj.Ntraces_Dy,1);
                obj.ph_az_err(obj.Ntraces_Dy+1) = beta(baseParams+1+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t);
                obj.ph_az_err(obj.Ntraces_Dy+2:obj.Ntraces_Dy+obj.Ntraces_Stat)...
                    = beta(baseParams+2+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t:baseParams+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+obj.Ntraces_Stat)';
            else
                obj.ph_az_err = beta(baseParams+1+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t:...
                    baseParams+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az)';
            end
            if obj.Nfit_Dy_amp == 1
                obj.amp_err(1:obj.Ntraces_Dy) = beta(baseParams+1+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+...
                    obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az)*AtomNum(1:obj.Ntraces_Dy)/ANbase;
                obj.amp_err(obj.Ntraces_Dy+1) = beta(baseParams+1+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+...
                    obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az)*AtomNum(obj.Ntraces_Dy+1)/ANbase;
                obj.amp_err(obj.Ntraces_Dy+2:obj.Ntraces_Dy+obj.Ntraces_Stat)...
                    = beta(baseParams+2+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az:...
                    baseParams+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az+obj.Ntraces_Stat)';
            else
                obj.amp_err = beta(baseParams+1+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az:...
                    baseParams+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az+...
                    obj.Nfit_Dy_amp+obj.Nfit_Stat_amp)';
            end
            if obj.Nfit_Dy_Qi == 1
                obj.Qi_err(1:obj.Ntraces_Dy) = beta(baseParams+1+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+...
                    obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az+obj.Nfit_Dy_amp+obj.Nfit_Stat_amp)...
                    *ones(obj.Ntraces_Dy,1);
                obj.Qi_err(obj.Ntraces_Dy+1) = beta(baseParams+1+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+...
                    obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az+obj.Nfit_Dy_amp+obj.Nfit_Stat_amp);
                obj.Qi_err(obj.Ntraces_Dy+2:obj.Ntraces_Dy+obj.Ntraces_Stat)...
                    = beta(baseParams+2+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az+...
                    obj.Nfit_Dy_amp+obj.Nfit_Stat_amp:...
                    baseParams+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az+...
                    obj.Nfit_Dy_amp+obj.Nfit_Stat_amp+obj.Ntraces_Stat)';
            else
                obj.Qi_err = beta(baseParams+1+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az+...
                    obj.Nfit_Dy_amp+obj.Nfit_Stat_amp:...
                    baseParams+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az+...
                    obj.Nfit_Dy_amp+obj.Nfit_Stat_amp+obj.Nfit_Dy_Qi+obj.Nfit_Stat_Qi)';
            end
            if obj.Nfit_Dy_Qf == 1
                obj.Qf_err(1:obj.Ntraces_Dy) = beta(baseParams+1+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+...
                    obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az+obj.Nfit_Dy_amp+obj.Nfit_Stat_amp+...
                    obj.Nfit_Dy_Qi+obj.Nfit_Stat_Qi)...
                    *ones(obj.Ntraces_Dy,1);
            else
                obj.Qf_err = beta(baseParams+1+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az+...
                    obj.Nfit_Dy_amp+obj.Nfit_Stat_amp+obj.Nfit_Dy_Qi+obj.Nfit_Stat_Qi:...
                    baseParams+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az+...
                    obj.Nfit_Dy_amp+obj.Nfit_Stat_amp+obj.Nfit_Dy_Qi+obj.Nfit_Stat_Qi+obj.Nfit_Dy_Qf)';
            end
         
            if (obj.Nfit_Dy_dPhidt+obj.Nfit_Stat_dPhidt) ~= 0
                obj.dPhidt_err = beta(baseParams+1+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az+...
                    obj.Nfit_Dy_amp+obj.Nfit_Stat_amp+obj.Nfit_Dy_Qi+obj.Nfit_Stat_Qi+obj.Nfit_Dy_Qf:...
                    baseParams+obj.Nfit_Dy_ph_t+obj.Nfit_Stat_ph_t+obj.Nfit_Dy_ph_az+obj.Nfit_Stat_ph_az+...
                    obj.Nfit_Dy_amp+obj.Nfit_Stat_amp+obj.Nfit_Dy_Qi+obj.Nfit_Stat_Qi+obj.Nfit_Dy_Qf+...
                    obj.Nfit_Stat_dPhidt+obj.Nfit_Dy_dPhidt)';
            end
        end % =============================================================   
                      
    end
end
