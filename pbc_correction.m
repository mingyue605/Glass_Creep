function delta = pbc_correction(delta,Lbox)

dim = size(delta,2);

for ii=1:dim
    delta(delta(:,ii)>Lbox(ii)/2,ii) =delta(delta(:,ii)>Lbox(ii)/2,ii) - Lbox(ii);
    delta(delta(:,ii)<-Lbox(ii)/2,ii) =delta(delta(:,ii)<-Lbox(ii)/2,ii) + Lbox(ii);
end

end