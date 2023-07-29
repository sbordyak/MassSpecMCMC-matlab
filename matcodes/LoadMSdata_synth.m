function d0 = LoadMSdata_synth(datafile,Isotopes,F_ind)

%%

if 1
    
    
    % Add subdirectory to matlab path
    addpath(genpath('./matcodes/'))
    
    
%     filen = './data/SyntheticDataset_001.txt';
    %filen = './data/pb linearity 2Mcps 1a-1475.txt';
    %dtmp=textread(filen,'%s','delimiter',','); %#ok<DTXTRD>
        
    fid=fopen(datafile,'r');
    dtmp=textscan(fid,'%s','delimiter',',','Headerlines',19); %#ok<DTXTRD>
    dall=reshape(dtmp{1},16,[])';
    fclose(fid);
    
    for ii=[2:16];
        d(:,ii) = str2num(char(dall(:,ii))); %#ok<*ST2NM>
    end
    
    Niso = length(Isotopes);
    
    Mass = d(:,6);
    Time = d(:,5);
    
    Axial = d(:,12);
    
    Block = d(:,2);
    Cycle = d(:,3);  
    
    if strcmp(datafile,'./results/RealData22/data/SyntheticDataset_12.txt')
        Cycle(1111:end) = Cycle(1111:end)+3; %sb629  JUST FOR ONE JANKY SYNTH DATA
    end
    Faraday = d(:,[7:11 13:16]);
    
    %clear d dtmp
end


% F_ind = [0 0 0 0 0    2 0 0 0 ;...
%          0 0 0 0 1    0 0 0 0 ];
  %%  
  
FarsUsed = find(sum(F_ind));    
 
F_ind = F_ind(:,FarsUsed) ;
Faraday = Faraday(:,FarsUsed);
     
     
Nfar = size(Faraday,2);
Ndet = Nfar + 1;
Nsamptot = length(Axial);
Nblock = max(1,length(unique(Block))-1);
%Nblock = 1; %length(unique(Block))-1;
for mm=1:Nblock
    Ncycle(mm) = length(unique(Cycle(Block==mm)));
end


% Find Faradays used



Cycle_Far = repmat(Cycle,1,Nfar);
Block_Far = repmat(Block,1,Nfar);
Time_Far = repmat(Time,1,Nfar);
Detector_Far = repmat(1:Nfar,Nsamptot,1);


Ax_ind = zeros(Niso,1);
Far_ind = zeros(Niso,Nfar);
for m = 1:Niso
    mmass = abs(Mass-Isotopes(m))<0.25;
    
    Ax_ind(mmass,1) = m;
    
    Far_ind(mmass,:) = repmat(F_ind(m,:),sum(mmass),1);
end

%%Far_ind(:,5) = 0; % Axial Faraday data is missing


for m=1:Nfar
    BL(m) = mean(Faraday(Ax_ind==0,m));
end
BL_Far = repmat(BL,Nsamptot,1);

%figure;plot(Time_Far(Far_ind==4),Faraday(Far_ind==4)-BL_Far(Far_ind==4),'.')

for m = 1:Nblock
    for n = 1:Ncycle(m)-1
        medCycleTime(m,n) = median(Time(Block==m & Cycle==n));
        minCycleTime(m,n) = min(Time(Block==m & Cycle==n));
        maxCycleTime(m,n) = max(Time(Block==m & Cycle==n));
        
        iminCT(m,n) = find(Time==minCycleTime(m,n));
        imaxCT(m,n) = find(Time==maxCycleTime(m,n));
        
    end
end

Tknots0 = [minCycleTime maxCycleTime(:,end)];
iTknots0 =  [iminCT imaxCT(:,end)];



%%


for m=1:Nblock
    Block_Time{m} = Time(iTknots0(m,1):iTknots0(m,end));
    InterpMat{m} = interp1(Tknots0(m,:),eye(length(Tknots0(m,:))),Block_Time{m},'spline');
    
    Nknots(m) = length(Tknots0(m,:));
    Ntb(m) = length(Block_Time{m});
    
    %sb726 Added this to include cycle number for InterpMat/time index
    CycleMat{m} = zeros(Ntb,1);
    for n = 1:Ncycle(m)-1
        CycleMat{m}((iminCT(n):imaxCT(n))-iminCT(1)+1) = n;
    end
        
    
end

ftimeind = repmat([1:Nsamptot]',1,Nfar);

for m = 1:Nblock
    for n = 1:Niso
        k = Far_ind==n & Block_Far==m;
        Data_Faraday{m,n} = Faraday(k);
        Time_Faraday{m,n} = Time_Far(k);
        Timeind_Faraday{m,n} = ftimeind(k)-iTknots0(m,1)+1;
        Detind_Faraday{m,n} = Detector_Far(k);
        Cycle_Faraday{m,n} = Cycle_Far(k);  %sb629  Adding cycle to d0

        
        k = Ax_ind==n & Block==m;
        Data_Axial{m,n} = Axial(k);
        Time_Axial{m,n} = Time(k);
        Timeind_Axial{m,n} = find(k)-iTknots0(m,1)+1;
        Cycle_Axial{m,n} = Cycle(k); %sb629  
    end
    
    
end

%Data_BL(1:Nfar) are baseline measurements for Faradays
for n = 1:Nfar
    Data_BL{n} = Faraday(Ax_ind==0,n);
    x0.BL(n,1) = mean(Data_BL{n});
end

% %Data_BL(Nfar+1) is baseline measurement for Daly
% Data_BL{Nfar+1} = Axial(Ax_ind==0,1);
% x0.BL(Nfar+1,1) = mean(Data_BL{Nfar+1});



%% SET UP DATA VECTOR
d0.data = [];
d0.iso_vec = [];
d0.blflag = [];
d0.sig_ind = [];
d0.det_vec = [];
d0.time = [];
d0.time_ind =[];
d0.axflag = [];
d0.cycle = [];
dblocktmp = [];
d0.Include = [];

d0.Nfar = Nfar;
d0.Nt = Nsamptot;
d0.Niso = Niso;
d0.Nblock = Nblock;
d0.Ndet = Ndet;
d0.Nknots = Nknots;
d0.Ncycle = Ncycle; %sb629  
d0.Ntb = Ntb;

d0.InterpMat = InterpMat;
d0.CycleMat = CycleMat;  %sb726
for m = 1:Nblock
    d0.IncludeMat{m} = ones(size(d0.CycleMat{m})); %sb726 Include that cycle?
end

d0.ReportInterval = 0.1; %sb629  This should be initialized based on header information
d0.Isotopes = Isotopes; %sb630 Add isotope field so it can be used for Mass Bias later

% Add Faraday baseline data
for n = 1:d0.Nfar
    
    Ndata = length(Data_BL{n});
    
    d0.data = [d0.data; Data_BL{n}];
    d0.iso_vec =[d0.iso_vec; zeros(Ndata,1)];
    d0.blflag = [d0.blflag; ones(Ndata,1)];
    d0.sig_ind = [d0.sig_ind; n*ones(Ndata,1)];
    d0.det_vec = [d0.det_vec; n*ones(Ndata,1)];
    d0.time = [d0.time; zeros(Ndata,1)];
    d0.time_ind =[d0.time_ind; zeros(Ndata,1)];
    d0.axflag = [d0.axflag; zeros(Ndata,1)];
    d0.cycle = [d0.cycle; zeros(Ndata,1)]; %sb629  
    dblocktmp = [dblocktmp; zeros(Ndata,1)];
            
    %sb726 Vector collects which entries to include in initialization
    %and MCMC evalutation - to be updated by user data exclusion in
    %Tripoli
    d0.Include = [d0.Include; true(Ndata,1)];
end


% % Add Axial baseline data
% for n = d0.Nfar+1
%     
%     Ndata = length(Data_BL{n});
%     
%     d0.data = [d0.data; Data_BL{n}];
%     d0.iso_vec =[d0.iso_vec; zeros(Ndata,1)];
%     d0.blflag = [d0.blflag; ones(Ndata,1)];
%     d0.sig_ind = [d0.sig_ind; n*ones(Ndata,1)];
%     d0.det_vec = [d0.det_vec; n*ones(Ndata,1)];
%     d0.time = [d0.time; zeros(Ndata,1)];
%     d0.time_ind =[d0.time_ind; zeros(Ndata,1)];
%     d0.axflag = [d0.axflag; zeros(Ndata,1)];
%     dblocktmp = [dblocktmp; zeros(Ndata,1)];
%     
% end

isig_ind = max(d0.sig_ind);

% Add Faraday data
for m = 1:Nblock
    for n = 1:Niso
        
        Ndata = length(Data_Faraday{m,n});
        
        si = n + isig_ind;
        
        
        d0.data = [d0.data; Data_Faraday{m,n}];
        d0.iso_vec =[d0.iso_vec; n*ones(Ndata,1)];
        d0.blflag = [d0.blflag; zeros(Ndata,1)];
        d0.sig_ind = [d0.sig_ind; si*ones(Ndata,1)];
        d0.det_vec = [d0.det_vec; Detind_Faraday{m,n}];
        d0.time = [d0.time; Time_Faraday{m,n}];
        d0.time_ind =[d0.time_ind; Timeind_Faraday{m,n}];
        d0.axflag = [d0.axflag; zeros(Ndata,1)];
        d0.cycle = [d0.cycle; Cycle_Faraday{m,n}]; %sb629  
        dblocktmp = [dblocktmp; m*ones(Ndata,1)];
                
        d0.Include = [d0.Include; true(Ndata,1)]; 
    end
end

isig_ind = max(d0.sig_ind);

% Add axial data
for m = 1:Nblock
    for n = 1:Niso
        
        Ndata = length(Data_Axial{m,n});
        
        si = n + isig_ind;
        
        
        d0.data = [d0.data; Data_Axial{m,n}];
        d0.iso_vec =[d0.iso_vec; n*ones(Ndata,1)];
        d0.blflag = [d0.blflag; zeros(Ndata,1)];
        d0.sig_ind = [d0.sig_ind; si*ones(Ndata,1)];
        d0.det_vec = [d0.det_vec; (Nfar+1)*ones(Ndata,1)];
        d0.time = [d0.time; Time_Axial{m,n}];
        d0.time_ind =[d0.time_ind; Timeind_Axial{m,n}];
        d0.axflag = [d0.axflag; ones(Ndata,1)];
        d0.cycle = [d0.cycle; Cycle_Axial{m,n}]; %sb629  
        dblocktmp = [dblocktmp; m*ones(Ndata,1)];
        
        d0.Include = [d0.Include; true(Ndata,1)]; 
    end
end


d0.Nsig = max(d0.sig_ind);
d0.Ndata = length(d0.data);


d0.block = false(d0.Ndata,d0.Nblock);
for m=1:d0.Nblock
    d0.block(dblocktmp==m,m) = true;
end

d0.iso_ind = false(d0.Ndata,d0.Niso);
for m=1:d0.Niso
    d0.iso_ind(d0.iso_vec==m,m) = true;
end

d0.det_ind = false(d0.Ndata,d0.Nfar);
for m=1:d0.Nfar+1
    d0.det_ind(d0.det_vec==m,m) = true;
end

d0.axflag = logical(d0.axflag);
d0.blflag = logical(d0.blflag);
d0.Include  = logical(d0.Include);
for m = 1:Nblock
d0.IncludeMat{m} = logical(d0.IncludeMat{m});
end

