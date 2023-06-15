
function d = ModelInitData(x0,d0);


%% MODEL DATA WITH INITIAL MODEL
II = d0.InterpMat;

for m=1:d0.Nfar%+1
    d(d0.blflag & d0.det_ind(:,m),1) = x0.BL(m);
end

for n = 1:d0.Nblock
    Intensity{n} = II{n}*x0.I{n};
    for m=1:d0.Niso;
        itmp = d0.iso_ind(:,m) & d0.axflag & d0.block(:,n);
        d(itmp) = exp(x0.lograt(m))*Intensity{n}(d0.time_ind(itmp));
        
        itmp = d0.iso_ind(:,m) & ~d0.axflag & d0.block(:,n);
        d(itmp) = exp(x0.lograt(m))*x0.DFgain^-1 *Intensity{n}(d0.time_ind(itmp)) + x0.BL(d0.det_vec(itmp));
    end
end

