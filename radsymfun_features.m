function  Fall = radsymfun_features(coords,atom_type,Lbox)

% parameters
rcut = 5;
incr = 0.1;
rbins = (0.1:incr:rcut);
nbins = numel(rbins);
delta = 0.1;

natoms = size(coords,1);

GxA = zeros(natoms,nbins);
GxB = zeros(natoms,nbins);

for ii=1:natoms

    coords_atom = coords(ii,:);
    coords_other = coords(setdiff(1:end,ii),:);
    type_other = atom_type(setdiff(1:end,ii),:);
    
    coordsA = coords_other(type_other==1,:);
    coordsB = coords_other(type_other==2,:);

    %% 1. Distance symmetry functions

    % xA RDFs
    delta_xA = pbc_correction(coordsA - coords_atom,Lbox);
    dist_xA = sqrt(sum(delta_xA.^2,2));
    for kk=1:nbins
        GxA(ii,kk) = sum( exp( -(rbins(kk)-dist_xA).^2 / (2*delta^2) ) );
    end

    % xB RDFs
    delta_xB = pbc_correction(coordsB - coords_atom,Lbox);
    dist_xB = sqrt(sum(delta_xB.^2,2));
    for kk=1:nbins
        GxB(ii,kk) = sum( exp( -(rbins(kk)-dist_xB).^2 / (2*delta^2) ) );
    end

end

Fall = [GxA GxB]; % a total of 100 features

end
