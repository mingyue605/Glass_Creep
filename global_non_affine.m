function [disp_components, disp_total] = global_non_affine(coords_ref,coords_curr,Lbox_ref,Lbox_curr)

% NOTE: the non-affine correction only works for orthorhombic deformations,
%       not for shear or torsion

% recenter the reference and current coordinates
coords_ref = coords_ref - mean(coords_ref);
coords_curr = coords_curr - mean(coords_curr);

% remove the affine deformations from the current coordinates
coords_curr_nonaffine = coords_curr.*(Lbox_ref./Lbox_curr);

% calculate displacements between the current coordinates corrected for
% non-affine deformations and the reference coordinates
disp_components = pbc_correction(coords_curr_nonaffine-coords_ref,Lbox_ref);

% calculate the total displacements
disp_total = sqrt(sum(disp_components.^2,2));

end