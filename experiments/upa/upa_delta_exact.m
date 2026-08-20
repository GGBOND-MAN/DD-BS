function delta = upa_delta_exact(Nx, Ny, d, r, Wx, Wy)
% UPA_DELTA_EXACT  Exact near-field path-length difference (r_{nx,ny} - r) for a UPA.
%
% Antennas sit on a plane at (nx*d, ny*d, 0); the user is at distance r in the
% direction with direction cosines (Wx, Wy) = (sinϑcosφ, sinϑsinφ), Wx^2+Wy^2<=1.
%
%   r_{nx,ny} = || r*u - p_{nx,ny} ||,   u=(Wx,Wy,Wz),  p=(nx d, ny d, 0)
%             = sqrt( r^2 - 2 r (nx d Wx + ny d Wy) + (nx d)^2 + (ny d)^2 )
%
% Returns delta = r_{nx,ny} - r as an (Nx*Ny) x 1 column (column-major flatten).
% This is the ground truth the beam must match to focus at (Wx,Wy,r).

nx = (-(Nx-1)/2:(Nx-1)/2);
ny = (-(Ny-1)/2:(Ny-1)/2);
[NX, NY] = ndgrid(nx, ny);
px = NX*d; py = NY*d;
rn = sqrt(r.^2 - 2*r*(px*Wx + py*Wy) + px.^2 + py.^2);
delta = rn(:) - r;
end
