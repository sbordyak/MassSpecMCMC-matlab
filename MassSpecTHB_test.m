clear
addpath("matcodes/");

addpath(genpath('./matcodes'))


%% Define input dataset and collector parameters

%runname = 'test_prop';
runname = 'RealData22';

if strcmp(runname,'test_prop')
    iset = 1; % Dataset index number
    Isotopes = [206 208]';
    Iso_Name = {'Pb206','Pb208'};
    
    % Faraday Index (row # corresponds to isotope measured on Daly
    % This is based on the setup of the collectors, the piece of
    % information lacking in the typical data files.
    F_ind = [0 0 0 0 0    2 0 0 0 ;...
        0 0 0 0 1    0 0 0 0 ];
    
elseif strcmp(runname,'RealData22')
    iset = 1; % Dataset index number
    Isotopes = [204 205 206 207 208]';
    Iso_Name = {'Pb204','Pb205','Pb206','Pb207','Pb208'};
    
    % Faraday Index (row # corresponds to isotope measured on Daly
    % This is based on the setup of the collectors, the piece of
    % information lacking in the typical data files.
    F_ind = [0 0 0 0 0    2 3 4 5;...
         0 0 0 0 1    3 4 5 0;...
         0 0 0 1 2    4 5 0 0;...
         0 0 1 2 3    5 0 0 0;...
         0 1 2 3 4    0 0 0 0];
end



foldername = ['./results/' runname];

outfolder = ['./results/' runname '/output/'];
if ~exist(outfolder,'dir')
    mkdir(outfolder)
end


datafile = sprintf( '%s/data/SyntheticDataset_%02d.txt',foldername,iset);
datamat  = sprintf( '%s/data/SyntheticDataset_%02d.mat',foldername,iset);



%% Load Data and create data information struct d0

d0 = LoadMSdata_synth(datafile,Isotopes,F_ind);

% Matrix to project spline knots for intensity function to time samples
InterpMat = d0.InterpMat; 



%% CREATE INITIAL MODEL AND MODEL DATA
[x0,d,Intensity] = InitializeModel_synth(d0);

% Hardcode values for some test I did
%x0.BL = [3e4; -1e4];
%x0.DFgain = 0.9;


%% INVERSION WITH rj-MCMC
% Hardcoded here, presumably set by user in file or GUI?


% MCMC Parameters
maxcnt = 2000;  % Maximum number of models to save
hier = 1;  % Hierachical?
datsav=100;  % Save model every this many steps

burn = 10;  % Burn-in, start doing stats after this many saved models

temp=1; % Unused parameter for parallel tempering algorithm

% Baseline multiplier - weight Daly more strongly (I think)
blmult = ones(size(d0.data));
blmult(d0.axflag)=0.1;


Ndata=d0.Ndata; % Number of picks
Nsig = d0.Nsig; % Number of noise variables


% Range for ratios and intensity parameters
prior.BL = [-1 1]*1e6;  % Faraday baseline
prior.BLdaly = [0 0];   % Daly baseline (no baseline uncertainty)
prior.lograt = [-20 20]; % Log ratio
prior.I = [0 1.5*max([x0.I{:}])];  % Intensity
prior.DFgain = [0.8 1.0];  % Daly-Faraday gain

prior.sig = [0 1e6];  % Noise hyperparameter for Faraday
prior.sigdaly = [0 0]; % Gaussian noise on Daly
prior.sigpois = [0 10]; % Poisson noise on Daly


% "Proposal Sigmas"
% Standard deviations for proposing changes to model
psig.BL = max(x0.BLstd)/10*1;  % Faraday Baseline
psig.BLdaly = 1e-1*1;  % Daly Baseline 
psig.lograt = 0.0005*.2;  % Log Ratio
psig.I = max(max([x0.I{:}])-min([x0.I{:}]))/100*1 ; % Intensity
psig.DFgain = 0.001; % Daly-Faraday gain

psig.sig = max(x0.BLstd); % Noise hyperparameter for Faraday
psig.sigpois = 0.5; % Poisson noise on Daly
psig.sigdaly = 0;  % Gaussian noise on Daly


% Index for data covariance term (change if multiple types of data)
% Can use more complex model for day - correlated or time variable.
% Not implemented here
%sig_ind = d0.sig_ind;

% data covariance vector
%Dsig = x.sig(sig_ind);

% Assign initial values for model x
x=x0;


%% Forward model data from initial model

% Forward model baseline measurements
for mm=1:d0.Nfar%+1  % Iterate over Faradays
    d(d0.blflag & d0.det_ind(:,mm),1) = x0.BL(mm); % Faraday Baseline
    dnobl(d0.blflag & d0.det_ind(:,mm),1) = 0; % Data with No Baseline
end

% Forward model isotope measurements
for n = 1:d0.Nblock  % Iterate over blocks
    
    % Calculate block intensity from intensity variables
    Intensity{n} = InterpMat{n}*x0.I{n};  
    Intensity2{n} = Intensity{n};
    
    %Iterate over Isotopes
    for mm=1:d0.Niso;
        % Calculate Daly data
        itmp = d0.iso_ind(:,mm) & d0.axflag & d0.block(:,n); % If isotope and axial and block number
        d(itmp) = exp(x0.lograt(mm))*Intensity{n}(d0.time_ind(itmp));
        dnobl(itmp) = d(itmp);
        
        % Calculate Faraday datas
        itmp = d0.iso_ind(:,mm) & ~d0.axflag & d0.block(:,n);
        dnobl(itmp) = exp(x0.lograt(mm))*x0.DFgain^-1 *Intensity{n}(d0.time_ind(itmp)); % Data w/o baseline
        d(itmp) = dnobl(itmp) + x0.BL(d0.det_vec(itmp)); % Add baseline
        
    end
end

% New data covariance vector
Dsig = sqrt(x0.sig(d0.det_vec).^2 + x0.sig(d0.iso_vec+d0.Ndet).*dnobl); 

% Initialize data residual vectors
restmp=zeros(size(Dsig));
restmp2=zeros(size(Dsig));

% Calculate data residuals from starting model
restmp = (d0.data-d).^2;


% Calculate error function
E=sum(restmp.*blmult./Dsig);  % Weighted by noise variance (for acceptance)
E0=sum(restmp);  % Unweighted (for tracking convergence)


%     ensname = sprintf('%s/Ensemble_run%02d.mat',outfolder,iset);
%     load(ensname,'ensemble')
%     [xmean,xcov] = UpdateMeanCovMS(x,0,0,ensemble,10000,0);
%     xcov0=xcov;

%% Initialize MCMC loop variables
cnt=0; % Counter
kept=zeros(5,4); % For displaying how many updates are accepted

clear ens*
ensemble=[]; % Make sure to start with new ensemble

Ndf = 1; % Number of DF gains = 1

% Size of model: # isotopes + # intensity knots + # baselines + # df gain
Nmod = d0.Niso + sum(d0.Ncycle) + d0.Nfar + Ndf ;

% Data and data covariance vectors
xmean = zeros(Nmod,1);
xcov = zeros(Nmod,Nmod);

% Adaptive MCMC proposal term
delx_adapt=zeros(Nmod,datsav);

%%
d0.iso_vec(d0.iso_vec==0)=d0.Niso; %Set BL to denominator iso

% Find beginning and end of each block for Faraday and Daly (ax)
for ii=1:d0.Nblock
    block0(ii,1) = find(d0.block(:,ii)&~d0.axflag,1,'first');
    blockf(ii,1) = find(d0.block(:,ii)&~d0.axflag,1,'last');
    blockax0(ii,1) = find(d0.block(:,ii)&d0.axflag,1,'first');
    blockaxf(ii,1) = find(d0.block(:,ii)&d0.axflag,1,'last');
end




%% MCMC Iterations
tic


for m = 1:maxcnt*datsav
    %%
    % Choose an operation for updating model
    oper = RandomOperMS(hier);

    
    clear delx
    
    
    if cnt<100000 
        adaptflag = 0;  % For first NN iterations, do regular MCMC
    else
        adaptflag = 1;  % After, switch to Adaptive
        temp = 1;
    end
    
    allflag = adaptflag;
    
    % Update model and save proposed update values (delx)
    [x2,delx] = UpdateMSv2(oper,x,psig,prior,ensemble,xcov,delx_adapt(:,mod(m,datsav)+1),adaptflag,allflag);
    
      
    %% Create updated data based on new model
    % I was working on making this more compact and some of the details
    % elude me.
    tmpBLind = [x2.BL; 0]; tmpBL = tmpBLind(d0.det_vec);
    tmpDF = ones(d0.Ndata,1); tmpDF(~d0.axflag) = x2.DFgain^-1;
    tmpLR = exp(x2.lograt(d0.iso_vec)); 
    tmpI = zeros(d0.Ndata,1);
    for n=1:d0.Nblock
        Intensity2{n} = InterpMat{n}*x2.I{n};
        tmpI(block0(n):blockf(n)) = Intensity2{n}(d0.time_ind(block0(n):blockf(n)));
        tmpI(blockax0(n):blockaxf(n)) = Intensity2{n}(d0.time_ind(blockax0(n):blockaxf(n)));
    end
    
    dnobl2 = tmpDF.*tmpLR.*tmpI;
    
    % New data vector
    d2 = dnobl2 + tmpBL;
    
    
    % New data covariance vector
    Dsig2 = x2.sig(d0.det_vec).^2 + x2.sig(d0.iso_vec+d0.Ndet).*dnobl2;
    
    % Calculate residuals for current and new model
    restmp = (d0.data-d).^2;
    restmp2 = (d0.data-d2).^2;

    
    E02=sum(restmp2);  % Unweighted error func (for visualization)
    
    
    if strcmp(oper,'noise')  
        % If noise operation
        E=sum(restmp./Dsig);
        E2=sum(restmp2./Dsig2);
        dE=E2-E; % Change in misfit
    else                     
        % If any other model update
        E=sum(restmp.*blmult./Dsig);
        E2=sum(restmp2.*blmult./Dsig2);
        dE=temp^-1*(E2-E); % Change in misfit
    end
    
    %%
    % Decide whether to accept or reject model
    keep = AcceptItMS(oper,dE,psig,delx,prior,Dsig,Dsig2,d0);
    
    % Update kept variables for display
    kept(OpNumMS(oper),2) = kept(OpNumMS(oper),2)+1;
    kept(OpNumMS(oper),4) = kept(OpNumMS(oper),4)+1;
    
    % If we accept the new model update values
    if keep>=rand(1)
        
        E=E2; % Misfit
        E0=E02; % Unweighted misfit
        d=d2; % Data
        x=x2; % Model
        Dsig=Dsig2;  % Model variance
        dnobl=dnobl2;  % Data without baseline
        Intensity=Intensity2;  % Intensity
        
        % Display info
        kept(OpNumMS(oper),1) = kept(OpNumMS(oper),1)+1;
        kept(OpNumMS(oper),3) = kept(OpNumMS(oper),3)+1;
        
    end
        
    
    % Save model values to ensemble at regular intervals
    if  mod(m,datsav)==0
        
        cnt=cnt+1; % Increment counter
        
        ensemble(cnt).lograt=x.lograt; % Log ratios
        for mm=1:d0.Nblock
            ensemble(cnt).I{mm}=x.I{mm}; % Intensity by block
        end
        ensemble(cnt).BL=x.BL;  % Baselines
        ensemble(cnt).DFgain=x.DFgain;  %Daly-Faraday gain
        ensemble(cnt).sig=x.sig;  % Noise hyperparameter
        ensemble(cnt).E=E;  % Misfit
        ensemble(cnt).E0=E0; % Unweighted misfit
        
        covstart = 50;  % After this many iterations, begin calculating covariance iteratively
        if cnt>=covstart+1
            % Iterative covariance
            [xmean,xcov] = UpdateMeanCovMS(x,xcov,xmean,ensemble,cnt-covstart,0);
            
            % Draw random numbers based on covariance for next update
            delx_adapt = mvnrnd(zeros(Nmod,1),2.38^2*xcov/Nmod,datsav)';
        end
             
        % Display update info to screen 
        if  mod(m,10*datsav)==0
            DisplayInfo
            tic
            kept(:,1:2) = 0;
        end
        
        
        % If you want, calculate stats and plot regularly during iterations
        if  mod(m,100*datsav)==0
            %%
            burn = min(1000,cnt-50);
            ens_rat =[ensemble.lograt];
            ens_sig =[ensemble.sig];
            ens_DF =[ensemble.DFgain];
            ens_BL =[ensemble.BL];
            ens_E =[ensemble.E];
            ens_E0 =[ensemble.E0];
            ratmean = mean(ens_rat(:,burn:cnt),2);
            ratstd = std(ens_rat(:,burn:cnt),[],2);
            
            BLmean = mean(ens_BL(:,burn:cnt),2);
            BLstd = std(ens_BL(:,burn:cnt),[],2);
            
            sigmean = mean(ens_sig(:,burn:cnt),2);
            sigstd = std(ens_sig(:,burn:cnt),[],2);
            
            DFmean = mean(ens_DF(:,burn:cnt),2);
            DFstd = std(ens_DF(:,burn:cnt),[],2);
            
            for m=1:d0.Nblock
                for n = 1:cnt;
                    ens_I{m}(:,n) =[ensemble(n).I{m}];
                end
                Imean{m} = mean(ens_I{m}(:,burn:cnt),2);
                Istd{m} = std(ens_I{m}(:,burn:cnt),[],2);
            end
            
            %PlotEnsemble_synth;  % Uncomment for plotting
            %PlotEnsemble_MS;
            %PlotData_MS
        end
        
        
    end % End saving/displaying/plotting
    
    
end  % End of MCMC iterations




%% Analysis and Plotting

burn = 1000; % Number of models to discard
ens_rat =[ensemble.lograt];
ens_sig =[ensemble.sig];
ens_DF =[ensemble.DFgain];
ens_BL =[ensemble.BL];
ens_E =[ensemble.E];
ens_E0 =[ensemble.E0];


for m=1:d0.Nblock
    for n = 1:cnt
        ens_I{m}(:,n) =[ensemble(n).I{m}];
    end
    Imean{m} = mean(ens_I{m}(:,burn:cnt),2);  % Mean intensity knots
    Istd{m} = std(ens_I{m}(:,burn:cnt),[],2); % Std dev intensity knots
end



% Calculate mean and st dev of ratios after burn in time
ratmean = mean(ens_rat(:,burn:cnt),2);  % Log ratios
ratstd = std(ens_rat(:,burn:cnt),[],2);

BLmean = mean(ens_BL(:,burn:cnt),2);  % Baselines
BLstd = std(ens_BL(:,burn:cnt),[],2);

sigmean = mean(ens_sig(:,burn:cnt),2);   % Noise hyperparams
sigstd = std(ens_sig(:,burn:cnt),[],2);

DFmean = mean(ens_DF(:,burn:cnt),2);   % Daly-Far gain
DFstd = std(ens_DF(:,burn:cnt),[],2);


%%



ensname = sprintf('%s/Ensemble_run%02d.mat',outfolder,iset);
save(ensname,'ensemble','d0','datafile','Isotopes','Iso_Name','prior',...
    'psig','maxcnt','datsav','burn','cnt','F_ind','*mean','*std','InterpMat','x')

if strcmp(runname,'test_prop')
    PlotEnsemble_synth   
elseif strcmp(runname,'RealData22')
    PlotEnsemble_MS
end

PlotProgress_MS

PlotData_MS
