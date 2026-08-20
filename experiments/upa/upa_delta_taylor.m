function delta = upa_delta_taylor(Nx, Ny, d, r, Wx, Wy, include_cross)
% UPA_DELTA_TAYLOR  Second-order Taylor of the UPA path-length difference.
%
% Expanding upa_delta_exact to 2nd order:
%
%   r_{nx,ny}-r ≈ -(nx d Wx + ny d Wy)                                 (linear: 2 angles)
%                 + [ (nx d)^2 (1-Wx^2) + (ny d)^2 (1-Wy^2) ] / (2r)   (diagonal curvature)
%                 - (nx d)(ny d) Wx Wy / r                             (CROSS-COUPLING term)
%
% include_cross = false  -> the phase a SEPARABLE TD architecture can realize
%                           (per-axis delays only: N_x + N_y TDs, cheap). It matches
%                           the two linear terms + the two diagonal curvatures, but
%                           CANNOT produce the non-separable cross term nx*ny.
% include_cross = true   -> the phase a FULLY-CONNECTED TD (N_x*N_y TDs, expensive)
%                           can realize. Used as the reference upper bound.
%
% Returns delta as (Nx*Ny) x 1 column, matching upa_delta_exact's ordering.

nx = (-(Nx-1)/2:(Nx-1)/2);
ny = (-(Ny-1)/2:(Ny-1)/2);
[NX, NY] = ndgrid(nx, ny);
px = NX*d; py = NY*d;

delta = -(px*Wx + py*Wy) ...
        + ( px.^2*(1-Wx^2) + py.^2*(1-Wy^2) ) / (2*r);
if nargin >= 7 && include_cross
    delta = delta - (px .* py) * (Wx*Wy) / r;
end
delta = delta(:);
end
