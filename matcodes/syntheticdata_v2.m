%% Synthetic data for beam interpolation
%runname = 'test_Ifunc';
runname = 'test_prop';

foldername = ['./results/' runname '/data/'];
if ~exist(foldername,'dir')
    mkdir(foldername)
end


Nset = 8;
%rat =  [1/50000 1/500 1/5 1 20 200 2000 20000];
rat =  [1/5 1/3 1/2 1 2 5 20 50];

intens = [1e5 2e5 5e5 1e6 2e6 5e6 1e7 2e7 5e7 1e8];
ncyc = 5:5:40;
Ifunc = 0:0.2:1.4;
Nblock = 1:8;


for ii = 1:Nset
    
    
    
    
    % set up analytical parameters
    % one block, 15 cycles, 30-second baseline, only two isotopes (206Pb and 208Pb)
    % true ratio of 208Pb/206Pb = 0.02
    
    integrationTime.Pb206 = 5; % measure 206Pb for 5 seconds every cycle
    integrationTime.Pb208 = 5; % measure 208Pb for 5 seconds every cycle
    integrationTime.settle = 2; % leave 2 seconds between isotope measurements for the magnet to 'settle'
    integrationTime.baseline = 30; % seconds of baseline at the beginning of the block
    integrationTime.reportinterval = 0.1; % an intensity reported every 0.1 seconds
    massspec.numberOfCycles = 15; %ncyc(ii); % measure: {settle, 206Pb, settle, 208Pb} fifteen times before the next baseline
    truedata.ratio208206 = rat(ii); %200; %rat(ii);%20; % true ratio
    
    maxintensity = 1e5; %intens(ii); %1e7;
    
    if truedata.ratio208206>=1
        intensity.Pb208start = maxintensity; % counts per second of 208Pb at start of run
        intensity.Pb206start = intensity.Pb208start/truedata.ratio208206; % counts per second of 206Pb at start of run
    else
        intensity.Pb206start = maxintensity; % counts per second of 208Pb at start of run
        intensity.Pb208start = intensity.Pb206start*truedata.ratio208206; % counts per second of 206Pb at start of run
    end
    
    intensity.fun = Ifunc(4); %Ifunc(ii);
    
    % set up function of intensity variation with time
    % normalized intensity, domain is totalTimeSpan (in seconds)
    % intensityFunction = @(t) 1 + 0*t;
    % intensityFunction = @(t) 1 - 0.001*t;
    intensityFunction = @(t) exp(-t.^2/(integrationTime.baseline + massspec.numberOfCycles*...
        (integrationTime.Pb206+integrationTime.Pb206+2*integrationTime.settle))^2*intensity.fun);
    
    % set up mass spectrometer intrinsic parameters
    massspec.dalyFaradayGain = 0.9; % Daly gain is 90% of Faraday signal
    massspec.resistor = 1e11; % ohms resistance
    massspec.L1baseline = 3e4; % cps baseline
    massspec.H1baseline = -1e4; % cps baseline
    massspec.noiseVarianceV = 1.3806e-23*(18+273)*massspec.resistor; % Johnson-Nyquist noise on Faradays, bandwidth = 1 second, in volts^2
    massspec.voltsToCPS = 6.24150934e18/massspec.resistor; % convert voltage to counts per second
    massspec.FnoiseVarianceCPS = massspec.noiseVarianceV * massspec.voltsToCPS^2; % Johnson-Nyquist noise on Faradays, bandwidth = 1 s, in cps^2
    massspec.dalyBaseline = 6/60; % counts per second baseline (dark noise).  6 counts per minute/60 = cps
    
    
    %% make a home for the data
    
    integrationTime.totaltime = integrationTime.baseline + ...
        massspec.numberOfCycles*(integrationTime.Pb206 + integrationTime.Pb208 + 2*integrationTime.settle);
    
    integrationTime.totalIntervals = integrationTime.totaltime/integrationTime.reportinterval;
    integrationTime.totalMeasurements = (integrationTime.baseline+massspec.numberOfCycles*(integrationTime.Pb206+integrationTime.Pb208))...
        /integrationTime.reportinterval;
    data = zeros(integrationTime.totalMeasurements,4); % [timestamp, L1, Daly, H1]
        
    
    for jj=1:Nblock(1)

    
    %% Sort out when things are happening
    
    % When is data (on peak or baseline) actually being recorded? indices when on-peak data is collected.
    time.MeasurementTimeIndices = [ones(integrationTime.baseline/integrationTime.reportinterval,1);...    %baseline
        repmat([zeros(integrationTime.settle/integrationTime.reportinterval,1); ... %settle
        ones(integrationTime.Pb206/integrationTime.reportinterval,1);...    %206Pb
        zeros(integrationTime.settle/integrationTime.reportinterval,1); ... %settle
        ones(integrationTime.Pb208/integrationTime.reportinterval,1)],...   %208Pb
        massspec.numberOfCycles, 1)];
    time.MeasurementTimeIndices = logical(time.MeasurementTimeIndices);
    time.totalTimeSpan = (0:integrationTime.reportinterval:(integrationTime.totaltime-integrationTime.reportinterval))';
    
    % data is only recorded when a measurement is happening
    data(:,1) = time.totalTimeSpan(time.MeasurementTimeIndices);% + (jj-1)*1.1*max(time.totalTimeSpan(time.MeasurementTimeIndices));
    
    % When is 206 in the axial (Daly) collector?
    time.Pb206AxialTimeIndices = [zeros(integrationTime.baseline/integrationTime.reportinterval,1);...
        repmat([ones(integrationTime.Pb206/integrationTime.reportinterval,1);... %206Pb
        zeros(integrationTime.Pb208/integrationTime.reportinterval,1)],...
        massspec.numberOfCycles, 1)]; % indices when on-peak data is collected
    time.Pb206AxialTimeIndices = logical(time.Pb206AxialTimeIndices);
    
    % When is 208 in the axial (Daly) collector?
    time.Pb208AxialTimeIndices = [zeros(integrationTime.baseline/integrationTime.reportinterval,1);...
        repmat([zeros(integrationTime.Pb206/integrationTime.reportinterval,1);...
        ones(integrationTime.Pb208/integrationTime.reportinterval,1)],... %208Pb
        massspec.numberOfCycles, 1)]; % indices when on-peak data is collected
    time.Pb208AxialTimeIndices = logical(time.Pb208AxialTimeIndices);
    
    
        
        
        %% Generate true isotope beam intensities
        
        truedata.iPb208 = intensity.Pb208start * intensityFunction(data(:,1));
        truedata.iPb206 = truedata.iPb208 / truedata.ratio208206;
        
        
        %% Make background collector noise and for Faradays, add 'true' baselines
        
        % background on Faradays is normally distributed Johnson noise, also depends on measurement interval
        % also, add the 'background', the reference voltage against which each signal is measured against
        measdata.backgroundFaradays = randn(integrationTime.totalMeasurements,2)*sqrt(massspec.FnoiseVarianceCPS/integrationTime.reportinterval)+...
            + repmat([massspec.L1baseline massspec.H1baseline], integrationTime.totalMeasurements, 1);
        % background on the Daly is just dark noise
        measdata.backgroundDaly     = poissrnd(massspec.dalyBaseline*integrationTime.reportinterval, ...
            [integrationTime.totalMeasurements,1])/integrationTime.reportinterval;
        % both Faradays and the Daly always experience a background
        data(:,2:4) = [measdata.backgroundFaradays(:,1) measdata.backgroundDaly measdata.backgroundFaradays(:,2)];
        
        % Each isotope beam has 'shot noise', Poisson-distributed.  Also, the Daly detects <100% of ion arrivals, refected in its gain
        measdata.iPb208IC = poissrnd(truedata.iPb208(time.Pb208AxialTimeIndices)*integrationTime.reportinterval)...
            /integrationTime.reportinterval * massspec.dalyFaradayGain;
        measdata.iPb206IC = poissrnd(truedata.iPb206(time.Pb206AxialTimeIndices)*integrationTime.reportinterval)...
            /integrationTime.reportinterval * massspec.dalyFaradayGain;
        measdata.iPb208F  = poissrnd(truedata.iPb208(time.Pb206AxialTimeIndices)*integrationTime.reportinterval)/integrationTime.reportinterval;
        measdata.iPb206F  = poissrnd(truedata.iPb206(time.Pb208AxialTimeIndices)*integrationTime.reportinterval)/integrationTime.reportinterval;
        
        % put the measured isotope beams in the right place
        data(time.Pb206AxialTimeIndices,3) = data(time.Pb206AxialTimeIndices,3) + measdata.iPb206IC;
        data(time.Pb206AxialTimeIndices,4) = data(time.Pb206AxialTimeIndices,4) + measdata.iPb208F;
        data(time.Pb208AxialTimeIndices,3) = data(time.Pb208AxialTimeIndices,3) + measdata.iPb208IC;
        data(time.Pb208AxialTimeIndices,2) = data(time.Pb208AxialTimeIndices,2) + measdata.iPb206F;
        
        return
        %% Reproduce the comma-separated .txt file produced by Ionvantage
        
        
        %save(sprintf( '%sSyntheticDataset_%02d.mat',foldername,ii),'data','measdata','truedata','integrationTime','intensity','intensityFunction','time','massspec')
        
        if jj==1
            filenamestr = sprintf( '%sSyntheticDataset_%02d.txt',foldername,ii);
            fid = fopen(filenamestr, 'w');
            
            % Header
            fprintf(fid, 'Version,1.00\r\n');
            fprintf(fid, 'Filename,SyntheticData_v2\r\n');
            fprintf(fid, 'Sample No,7\r\n');
            fprintf(fid, 'Sample ID,SyntheticDataSet1\r\n');
            fprintf(fid, 'Sample Type,\r\n');
            fprintf(fid, 'Analysis name,\r\n');
            fprintf(fid, 'User name,\r\n');
            fprintf(fid, 'Method name,syntheticdata_v2.m\r\n');
            fprintf(fid, 'MS Tune file,\r\n');
            fprintf(fid, 'Inlet file,\r\n');
            fprintf(fid, 'Serial Number,2019\r\n');
            fprintf(fid, 'Centre Channel,RM\r\n');
            fprintf(fid, 'DataDump Level,CSV-RAW\r\n');
            fprintf(fid, ['Acquire Date,' '\r\n']);% char(datetime('now','Format','eeee MMMM d y HH:mm:ss')) '\r\n']);
            fprintf(fid, ['Cycles per Block,' num2str(massspec.numberOfCycles) '\r\n']);
            fprintf(fid, 'Num. Blocks,100\r\n');
            fprintf(fid, 'Num. Sequences,2\r\n');
            fprintf(fid, '#START\r\n');
            fprintf(fid, 'ID,Block,Cycle,Integ,Time,Mass,Low5,Low4,Low3,Low2,Ax Fara,Axial,High1,High2,High3,High4\r\n');
            
            
           end 

        blocktime = (jj-1)*1.1*max(time.totalTimeSpan(time.MeasurementTimeIndices));
        
        % write baseline data
        integrationTime.baselineIntegrations = integrationTime.baseline/integrationTime.reportinterval;
        for ndxBL = 1:integrationTime.baselineIntegrations
            
            fprintf(fid, ['BL1,' num2str(jj) ',0,' num2str(ndxBL) ',' num2str(data(ndxBL,1)+blocktime) ',205.500,0,0,0,0,' ...
                num2str(data(ndxBL,2)) ',' num2str(data(ndxBL,3)) ',' num2str(data(ndxBL,4)) ',0,0,0\r\n']);
        end
        
        
        
        
        % write on-peak data
        ndxData = integrationTime.baselineIntegrations; % indexes into data matrix created above
        for ndxCycle = 1:massspec.numberOfCycles
            
            % for the sequence where 206 is on the Daly
            for ndxInteg = 1:integrationTime.Pb206/integrationTime.reportinterval
                ndxData = ndxData + 1;
                fprintf(fid, ['S1,' num2str(jj) ',' num2str(ndxCycle) ',' num2str(ndxInteg) ',' num2str(data(ndxData,1)+blocktime) ',206.000,0,0,0,0,' ...
                    num2str(data(ndxData,2)) ',' num2str(data(ndxData,3)) ',' num2str(data(ndxData,4)) ',0,0,0\r\n']);
            end
            
            % for the sequence where 208 is on the Daly
            for ndxInteg = 1:integrationTime.Pb208/integrationTime.reportinterval
                ndxData = ndxData + 1;
                fprintf(fid, ['S2,' num2str(jj) ',' num2str(ndxCycle) ',' num2str(ndxInteg) ',' num2str(data(ndxData,1)+blocktime) ',208.000,0,0,0,0,' ...
                    num2str(data(ndxData,2)) ',' num2str(data(ndxData,3)) ',' num2str(data(ndxData,4)) ',0,0,0\r\n']);
            end
            
        end
        
        
    end
    
    
    fclose(fid);
    
    
end