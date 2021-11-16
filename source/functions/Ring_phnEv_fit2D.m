function [bfp,ConfIntv,FittedData,AMP] = ...
    Ring_phnEv_fit2D(t,th,norm_n1ds,rfunc,ANfunc,nlobes,varargin)
% =========================================================================
% Function to fit n1D_th to a 2D sinusoidal function to map the phonon
% evolution with time and angular position. This takes care of changing ring radii
% Inputs:
%   t: time [s]
%   th: angle [rad]
%   norm_n1Ds: Normalized n1D
%   rfunc: Function handle for radius. Should be a fucntion of time only
%   nlobes: mode generated
%   ExpectedGamma: Coefficient such that c_theta \propto r^{-\gamma/2}
%   CLims: Limits of color plot
%   FigNum: figure to plot result
%   InitialGuess: Initial guess to fit (amp,c,Q,tphi0,ang_phi)
%   amp: amplitude[DL such that n_0 = 1 is mean value] of phonon ossc.

% Outputs:
%   bfp: Fit results (amp,ci,Qi,Qf, t_phi0,ang_phi0,)
%   ConfIntv: The error bars (95% confidence interval) on above values
% =========================================================================
%% Assemble optional inputs =====================================================
p = inputParser;
p.addParameter('startOffset',0,@(x)isnumeric(x));
p.addParameter('ExpectedGamma',1/2,@(x)isnumeric(x));
p.addParameter('ExpectedHub',1/2,@(x)isnumeric(x));
p.addParameter('FitGammaHub',false,@(x)islogical(x));
p.addParameter('FitVaryingQ',false,@(x)islogical(x));
p.addParameter('PersCurr',false,@(x)islogical(x));
p.addParameter('dPhidt',0,@(x)isnumeric(x));
p.addParameter('CLims',[0.5 1.5],@(x)(isnumeric(x) && length(x)==2));
p.addParameter('FigNum',2,@(x)isnumeric(x));
p.addParameter('InitialGuess',[0.5 4000 10 pi pi/2],...
    @(x)(isnumeric(x) && length(x)==5));
p.addParameter('lower',[0 0 0 -0*pi -0*pi],...
    @(x)(isnumeric(x) && length(x)==6));
p.addParameter('upper',[Inf Inf Inf 2*pi 2*pi],...
    @(x)(isnumeric(x) && length(x)==6));
p.addParameter('PlotResiduals',false,@(x)islogical(x));
p.parse(varargin{:});

t_startOffset = p.Results.startOffset;
gamma_exp = p.Results.ExpectedGamma;
hub_exp = p.Results.ExpectedHub;
fignum = p.Results.FigNum;
clims = p.Results.CLims;
ig = p.Results.InitialGuess;
lb = p.Results.lower;
ub = p.Results.upper;
% Extract from the rfunc the values of r ========================================
r = rfunc(t);
if ( p.Results.PersCurr)
    ig = [ig p.Results.dPhidt] ;
    lb = [lb -Inf];
    ub = [ub +Inf];
end
if (p.Results.FitVaryingQ && (r(1) ~= r(end)))
    ig = [ig ig(3)];
    lb = [lb lb(3)];
    ub = [ub ub(3)];
end
if (p.Results.FitGammaHub && (r(1) ~= r(end)))
    ig = [ig gamma_exp hub_exp];
    lb = [lb 0 -Inf];
    ub = [ub 1 Inf];
end

%% Arrange Inputs ===============================================================
% Let's focus on ts greater than perttime =======================================
inds = t>t_startOffset;
[TS,THS] = meshgrid(t,th);
TS_FIT = TS(:,inds);
THS_FIT = THS(:,inds);
norm_n1ds_FIT = norm_n1ds(:,inds);

% Plot the raw data =============================================================
if p.Results.PlotResiduals
    figure(fignum),clf, subplot(2,2,1),
    imagesc(1e3*t(inds),th,norm_n1ds(:,inds),clims);colorbar;
    title('Data');
    xlabel('${ t} {\rm(ms)}$','Interpreter','Latex');
    ylabel('$\theta {\rm (radians)}$','Interpreter','Latex');
    drawnow
end

%% Now do the fitting ===========================================================
if r(1) == r(end)
    if p.Results.PersCurr
        fitFunc =  @(x)RingStat_phnEv_2D(TS_FIT(:),THS_FIT(:),r(1),nlobes,...
            x(1),x(2),x(3),x(4),x(5),'vector',true,...
            'PersCurr',p.Results.PersCurr,'dPhidt',x(6));
    else
        fitFunc =  @(x)RingStat_phnEv_2D(TS_FIT(:),THS_FIT(:),r(1),nlobes,...
            x(1),x(2),x(3),x(4),x(5),'vector',true,...
            'PersCurr',p.Results.PersCurr,'dPhidt',0);
    end
    figure(fignum), subplot(2,2,2);
    imagesc(1e3*t(inds),th,reshape(fitFunc(ig),size(norm_n1ds_FIT)),clims);
    title('Initial Guess');colorbar;
    xlabel('${ t} {\rm(ms)}$','Interpreter','Latex');
    [bfp,~,Red4CI,~,~,~,jacobian4CI] = lsqnonlin(@(x)(norm_n1ds_FIT(:)-fitFunc(x)),...
        ig,lb,ub);
    ConfIntv = nlparci(bfp,Red4CI,'jacobian',jacobian4CI);
    [FittedData, AMP] = fitFunc(bfp);
    FittedData = reshape(FittedData,size(norm_n1ds_FIT));
else
    if p.Results.FitVaryingQ
        if p.Results.FitGammaHub
            if p.Results.PersCurr
                fitFunc =  @(x)RingDy_phnEv_2D(TS_FIT(:),THS_FIT(:),rfunc,ANfunc,nlobes,...
                    x(1),x(2),x(3),x(7),x(4),x(5),'vector',true,...
                    'PersCurr',p.Results.PersCurr,'dPhidt',x(6),...
                    'Gamma',x(8),'Hub',x(9));
            else
                fitFunc =  @(x)RingDy_phnEv_2D(TS_FIT(:),THS_FIT(:),rfunc,ANfunc,nlobes,...
                    x(1),x(2),x(3),x(6),x(4),x(5),'vector',true,...
                    'PersCurr',p.Results.PersCurr,'dPhidt',0,...
                    'Gamma',x(7),'Hub',x(8));
            end
            figure(fignum), subplot(2,2,2);
            imagesc(1e3*t(inds),th,reshape(fitFunc(ig),size(norm_n1ds_FIT)),clims);
            title('Initial Guess');colorbar;
            xlabel('${ t} {\rm(ms)}$','Interpreter','Latex');
            [bfp,~,Red4CI,~,~,~,jacobian4CI] = lsqnonlin(@(x)(norm_n1ds_FIT(:)-fitFunc(x)),...
                ig,lb,ub);
            ConfIntv = nlparci(bfp,Red4CI,'jacobian',jacobian4CI);
            [FittedData, AMP] = fitFunc(bfp);
            FittedData = reshape(FittedData,size(norm_n1ds_FIT));
        else
            if p.Results.PersCurr
                fitFunc =  @(x)RingDy_phnEv_2D(TS_FIT(:),THS_FIT(:),rfunc,ANfunc,nlobes,...
                    x(1),x(2),x(3),x(7),x(4),x(5),'vector',true,...
                    'PersCurr',p.Results.PersCurr,'dPhidt',x(6),...
                    'Gamma',gamma_exp,'Hub',hub_exp);
            else
                fitFunc =  @(x)RingDy_phnEv_2D(TS_FIT(:),THS_FIT(:),rfunc,ANfunc,nlobes,...
                    x(1),x(2),x(3),x(6),x(4),x(5),'vector',true,...
                    'PersCurr',p.Results.PersCurr,'dPhidt',0,...
                    'Gamma',gamma_exp,'Hub',hub_exp);
            end
            figure(fignum), subplot(2,2,2);
            imagesc(1e3*t(inds),th,reshape(fitFunc(ig),size(norm_n1ds_FIT)),clims);
            title('Initial Guess');colorbar;
            xlabel('${ t} {\rm(ms)}$','Interpreter','Latex');
            [bfp,~,Red4CI,~,~,~,jacobian4CI] = lsqnonlin(@(x)(norm_n1ds_FIT(:)-fitFunc(x)),...
                ig,lb,ub);
            ConfIntv = nlparci(bfp,Red4CI,'jacobian',jacobian4CI);
            [FittedData, AMP] = fitFunc(bfp);
            FittedData = reshape(FittedData,size(norm_n1ds_FIT));
        end
    else
        if p.Results.FitGammaHub
            if p.Results.PersCurr
                fitFunc =  @(x)RingDy_phnEv_2D(TS_FIT(:),THS_FIT(:),rfunc,ANfunc,nlobes,...
                    x(1),x(2),x(3),x(3),x(4),x(5),'vector',true,...
                    'PersCurr',p.Results.PersCurr,'dPhidt',x(6),'Gamma',x(7),'Hub',x(8));
            else
                fitFunc =  @(x)RingDy_phnEv_2D(TS_FIT(:),THS_FIT(:),rfunc,ANfunc,nlobes,...
                    x(1),x(2),x(3),x(3),x(4),x(5),'vector',true,...
                    'PersCurr',p.Results.PersCurr,'dPhidt',0,'Gamma',x(6),'Hub',x(7));
            end
            figure(fignum), subplot(2,2,2);
            imagesc(1e3*t(inds),th,reshape(fitFunc(ig),size(norm_n1ds_FIT)),clims);
            title('Initial Guess');colorbar;
            xlabel('${ t} {\rm(ms)}$','Interpreter','Latex');
            [bfp,~,Red4CI,~,~,~,jacobian4CI] = lsqnonlin(@(x)(norm_n1ds_FIT(:)-fitFunc(x)),...
                ig,lb,ub);
            ConfIntv = nlparci(bfp,Red4CI,'jacobian',jacobian4CI);
            [FittedData, AMP] = fitFunc(bfp);
            FittedData = reshape(FittedData,size(norm_n1ds_FIT));
        else
            if p.Results.PersCurr
                fitFunc =  @(x)RingDy_phnEv_2D(TS_FIT(:),THS_FIT(:),rfunc,ANfunc,nlobes,...
                    x(1),x(2),x(3),x(3),x(4),x(5),'vector',true,...
                    'PersCurr',p.Results.PersCurr,'dPhidt',x(6),...
                    'Gamma',gamma_exp,'Hub',hub_exp);
            else
                fitFunc =  @(x)RingDy_phnEv_2D(TS_FIT(:),THS_FIT(:),rfunc,ANfunc,nlobes,...
                    x(1),x(2),x(3),x(3),x(4),x(5),'vector',true,...
                    'PersCurr',p.Results.PersCurr,'dPhidt',0,...
                    'Gamma',gamma_exp,'Hub',hub_exp);
            end
            figure(fignum), subplot(2,2,2);
            imagesc(1e3*t(inds),th,reshape(fitFunc(ig),size(norm_n1ds_FIT)),clims);
            title('Initial Guess');colorbar;
            xlabel('${ t} {\rm(ms)}$','Interpreter','Latex');
            [bfp,~,Red4CI,~,~,~,jacobian4CI] = lsqnonlin(@(x)(norm_n1ds_FIT(:)-fitFunc(x)),...
                ig,lb,ub);
            ConfIntv = nlparci(bfp,Red4CI,'jacobian',jacobian4CI);
            [FittedData, AMP] = fitFunc(bfp);
            FittedData = reshape(FittedData,size(norm_n1ds_FIT));
        end
    end
end
AMP = reshape(AMP,size(norm_n1ds_FIT));
AMP = AMP(1,:);
if p.Results.PlotResiduals
    figure(fignum),subplot(2,2,3)
    imagesc(1e3*t(inds),th,FittedData,clims);
    title('Fit');colorbar;
    ylabel('$\theta {\rm (radians)}$','Interpreter','Latex');
    xlabel('${ t} {\rm(ms)}$','Interpreter','Latex');
    figure(fignum), subplot(2,2,4),
    imagesc(1e3*t(inds),th,norm_n1ds(:,inds) - ...
        FittedData,clims);
    title('Residuals');colorbar;
    xlabel('${ t} {\rm(ms)}$','Interpreter','Latex');
end

end
