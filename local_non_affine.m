function [Dsquare_min,Dsquare_min_norm] = local_non_affine(coords_ref,coords_curr,Lbox_ref,Lbox_curr)

% Calculate non-affine displacements according to Falk & Lagner:
% Dynamics of viscoplastic deformation in amorphous solids, PRE, 1998

% cutoff radius of an atoms' environment
Rcut = 2.5;

% number of atoms in the system
natoms = size(coords_ref,1);

Dsquare_min = zeros(natoms,1);
Dsquare_min_norm = zeros(natoms,1);

deltaij = eye(3);

for ii=1:natoms

    % 1. Find the nearest neighbors in the old configuration
    delta = pbc_correction(coords_ref - coords_ref(ii,:),Lbox_ref);
    [idx,~] = find(sqrt(sum(delta.^2,2))<=Rcut);

    % remove self from the mix
    [self,~] = find(idx==ii);
    idx(self) = [];

    % 2. Calculate the strain tensor: Formula 2.12-2.14
    Xij = zeros(3,3);
    Yij = zeros(3,3);
    for jj=1:numel(idx)
        rcurr = pbc_correction(coords_curr(idx(jj),:) - coords_curr(ii,:),Lbox_curr);
        rref = pbc_correction(coords_ref(idx(jj),:) - coords_ref(ii,:),Lbox_ref);
        % Formula 2.12
        Xij = Xij + rcurr'*rref;
        % Formula 2.13
        Yij = Yij + rref'*rref;
    end
    % Formula 2.14
    Eij = Xij/Yij - deltaij;

    % 3. Calcualte the non-affine displacements
    for jj=1:numel(idx)
        rcurr = pbc_correction(coords_curr(idx(jj),:) - coords_curr(ii,:),Lbox_curr);
        rref = pbc_correction(coords_ref(idx(jj),:) - coords_ref(ii,:),Lbox_ref);
        % Formula 2.11
        Dsquare_min(ii) = Dsquare_min(ii) + sum((rcurr - ((Eij + deltaij)*rref')').^2);
    end

    Dsquare_min_norm(ii) = Dsquare_min(ii)/numel(idx);

end

end