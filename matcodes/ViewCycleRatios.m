
function ViewCycleRatios(x,d0,Iso_Name);
% Function to compute mean and standard error for isotope ratios for each
% cycle after correcting for Daly-Faraday gain, baselines, and varying 
% intensity. To be used in workflow for user data rejection

figure

for m=1:d0.Nblock % Loop over block
    II = d0.InterpMat{m};
    IntensityFn = II*x.I{m};
    
    for ii = 1:d0.Niso %[1 3 4 5];%  Loop over ratio
        dind = (d0.iso_ind(:,ii) &   d0.block(:,m));
        dd=d0.data(dind);
        tmpdetvec = d0.det_vec(dind);
        tmpcycvec = d0.cycle(dind); 
        
        % Correct raw data for baseline and DF gain
        dd(~d0.axflag(dind)) = (dd(~d0.axflag(dind)) - x.BL(tmpdetvec(~d0.axflag(dind))))*x.DFgain ;
        
        % Sort by time
        tmptime = d0.time_ind(dind);
        [tsort,dsort]=sort(tmptime);
        
        dd=dd(dsort);
        cyc = tmpcycvec(dsort); % Cycle number vector for isotope
       
        % Ratio between corrected data and intensity function
        RawRatio{ii} = (dd./IntensityFn(tsort)); 
        
        for jj = 1:d0.Ncycle-1
            % Mean of ratio for each cycle
            MeanCycleIso(jj,ii) = (mean((RawRatio{ii}(cyc==jj))));   %sb1121
            
            %var_data = (mean(x.BLstd).^2+IntensityFn(cyc==jj));
            var_data = mean(x.BLstd).^2;
            
            %[mean_of_logs, std_of_logs] = MeanofLogs(dd(cyc==jj),IntensityFn(cyc==jj),var_data);
            %MeanCycleIso(jj,ii) = exp(mean_of_logs);         
            %StdCycleIso(jj,ii) = std_of_logs; % This isn't quite working %sb726
            
          
            
            % Standard deviation and standard error for each cycle
            Numcyc = sum(cyc==jj);
            StdCycleIso(jj,ii) = (std((RawRatio{ii}(cyc==jj))));  %sb1121
            SterrCycleIso(jj,ii) = StdCycleIso(jj,ii)/sqrt(Numcyc);      %sb1121
            
            %MeanCycleIso(jj,ii) = mean(log(RawRatio{ii}(cyc==jj))); 
            %SterrCycleIso(jj,ii) = std(log(RawRatio{ii}(cyc==jj)))/sqrt(Numcyc); 
            
            
        end
    end
    
    for ii = 1:d0.Niso-1 %[1 3 4] %
        % Ratio between isotopes and denominator isotopes by cycle
        MeanCycleRatio(:,ii) = MeanCycleIso(:,ii)./MeanCycleIso(:,d0.Niso); 
        %MeanCycleRatio(:,ii) = exp(MeanCycleIso(:,ii)-MeanCycleIso(:,d0.Niso)); 
        
        SterrCycleRatio(:,ii) = sqrt(SterrCycleIso(:,ii).^2+SterrCycleIso(:,d0.Niso).^2);  %sb1121
        %tmps = sqrt(SterrCycleIso(:,ii).^2+SterrCycleIso(:,d0.Niso).^2);
        %tmpr = MeanCycleIso(:,ii)-MeanCycleIso(:,d0.Niso);
        %SterrCycleRatio(:,ii) = exp(tmpr+tmps) - exp(tmpr);
        
        % Plot them up with error bars
        subplot(d0.Niso-1,d0.Nblock,m+d0.Nblock*(ii-1))
        
        errorbar(1:d0.Ncycle-1,MeanCycleRatio(:,ii),SterrCycleRatio(:,ii),'Marker','o')
        set(gca,'XTick',1:d0.Ncycle)
        if ii == d0.Niso-1
            xlabel('Cycle number')
        elseif ii == 1
            title(sprintf('Block %d',m'))
            
        end
        if m == 1
            ylabel(sprintf('%s/%s',Iso_Name{ii},Iso_Name{d0.Niso}))
        end  
    end

    
end

return



