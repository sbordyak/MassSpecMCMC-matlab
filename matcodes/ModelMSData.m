
function d = ModelMSData(x,d0);


%% MODEL DATA WITH INITIAL MODEL
II = d0.InterpMat;

d = zeros(size(d0.data));

for m=1:d0.Nfar%+1
    d(d0.blflag & d0.det_ind(:,m),1) = x.BL(m);
end


LogRatios = [x.lograt; 0];

for n = 1:d0.Nblock
    Intensity = II{n}*x.I{n};
    for m=1:d0.Niso;
        itmp = d0.iso_ind(:,m) & d0.axflag & d0.block(:,n);
        d(itmp) = exp(LogRatios(m))*Intensity(d0.time_ind(itmp));
        
        itmp = d0.iso_ind(:,m) & ~d0.axflag & d0.block(:,n);
        d(itmp) = exp(LogRatios(m))*x.DFgain^-1 *Intensity(d0.time_ind(itmp)) + x.BL(d0.det_vec(itmp));
    end
end

