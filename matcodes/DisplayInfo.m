%function DisplayInfoP(aniso,hier,iproc,countr,m,ttoc,dispnum,Nvor,EE,kept,temp)

if ~allflag
    display(sprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'));
    %display(sprintf('Processor %d, %d models accepted out of %d',iproc,countr,m));
    display(sprintf('Elapsed time = %0.2f Seconds for %d realizations (%d total)',toc,10*datsav,m));
    %display(sprintf('Cell count = %d, Temperature = %.2f',Nvor,temp));
    display(sprintf('Error function = %.0f',sqrt(E0/Ndata)));
    display(sprintf('Change Log Ratio:   %d of %d accepted (%.1f%% total)',kept(1,1:2),100*kept(1,3)/kept(1,4)));
    display(sprintf('Change Intensity:   %d of %d accepted (%.1f%% total)',kept(2,1:2),100*kept(2,3)/kept(2,4)));
    display(sprintf('Change DF Gain:     %d of %d accepted (%.1f%% total)',kept(3,1:2),100*kept(3,3)/kept(3,4)));
    display(sprintf('Change Baseline:    %d of %d accepted (%.1f%% total)',kept(4,1:2),100*kept(4,3)/kept(4,4)));
    if hier==1;
        display(sprintf('Noise:              %d of %d accepted (%.1f%% total)',kept(5,1:2),100*kept(5,3)/kept(5,4)));
    end
    display(sprintf(' '));
    
else
    
    display(sprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'));
    %display(sprintf('Processor %d, %d models accepted out of %d',iproc,countr,m));
    display(sprintf('Elapsed time = %0.2f Seconds for %d realizations (%d total)',toc,10*datsav,m));
    display(sprintf('Error function = %.0f',sqrt(E0/d0.Ndata)));
    display(sprintf('Change all variables:   %d of %d accepted (%.1f%% total)',sum(kept(1:4,1:2)),100*sum(kept(1:4,3))/sum(kept(1:4,4))));
    
    if hier==1;
        display(sprintf('Noise:              %d of %d accepted (%.1f%% total)',kept(5,1:2),100*kept(5,3)/kept(5,4)));
    end
    display(sprintf(' '));
    
end