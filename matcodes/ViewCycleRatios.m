
function ViewCycleRatios(x,d0,Iso_Name);
% Function to compute mean and standard error for isotope ratios for each
% cycle after correcting for Daly-Faraday gain, baselines, and varying 
% intensity. To be used in workflow for user data rejection

figure

for m=1:d0.Nblock % Loop over block
    II = d0.InterpMat{m};
    IntensityFn = II*x.I{m};
    
    for ii = 1:d0.Niso % Loop over ratio
        dind = (d0.iso_ind(:,ii) &   d0.block(:,m));
        dd=d0.data(dind);
        tmpdetvec = d0.det_vec(dind);
        tmpcycvec = d0.cycle(dind); 
        
        % Correct raw data for baseline and DF gain
        dd(~d0.axflag(dind)) = (dd(~d0.axflag(dind)) - x.BL(tmpdetvec(~d0.axflag(dind))))*x.DFgain ;
        
        % Sort by time
        tmptime = d0.time_ind(dind);
        [~,dsort]=sort(tmptime);
        
        dd=dd(dsort);
        cyc = tmpcycvec(dsort); % Cycle number vector for isotope
       
        % Ratio between corrected data and intensity function
        RawRatio(:,ii) = (dd./IntensityFn); 
        
        for jj = 1:d0.Ncycle
            % Mean of ratio for each cycle
            MeanCycleIso(jj,ii) = (mean((RawRatio(cyc==jj,ii)))); % 
            
            % Standard deviation and standard error for each cycle
            Numcyc = sum(cyc==jj);
            StdCycleIso(jj,ii) = (std((RawRatio(cyc==jj,ii))));
            SterrCycleIso(jj,ii) = StdCycleIso(jj,ii)/sqrt(Numcyc);         
        end
    end
    
    for ii = 1:d0.Niso-1
        % Ratio between isotopes and denominator isotopes by cycle
        MeanCycleRatio(:,ii) = MeanCycleIso(:,ii)./MeanCycleIso(:,d0.Niso); 
        SterrCycleRatio(:,ii) = sqrt(SterrCycleIso(:,ii).^2+SterrCycleIso(:,d0.Niso).^2);
        
        % Plot them up with error bars
        subplot(d0.Niso-1,d0.Nblock,m+d0.Nblock*(ii-1))
        
        errorbar(1:d0.Ncycle,MeanCycleRatio(:,ii),SterrCycleRatio(:,ii),'Marker','o')
        set(gca,'XTick',1:d0.Ncycle)
        if ii == d0.Niso-1
            xlabel('Cycle number')
        end
        if m == 1
            ylabel(sprintf('%s/%s',Iso_Name{ii},Iso_Name{d0.Niso}))
        end  
    end

    
end





