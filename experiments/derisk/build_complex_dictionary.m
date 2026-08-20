function [D, grid_sin, grid_r, colnorm] = build_complex_dictionary(p, w_2D_near, Gtheta, Gr)
% BUILD_COMPLEX_DICTIONARY  Complex (phase-preserving) DDBS sensing dictionary.
%
% Unlike the baseline match filter (which correlates POWER |y|^2 and discards
% cross-subcarrier phase), this builds the COMPLEX atom each candidate path
% produces at the receiver across all pilots and subcarriers:
%
%     s_{theta,r}(k,m) = h^{(theta,r)}_m * w_DDBS(:,k,m)
%
% where h^{(theta,r)} is the single-path (unit-gain) near-field channel
% (near_field_channel.m, which includes the exp(-j2pi f r/c) common phase).
% Stacked over (k,m) this gives a length-(k*M) complex atom. The full multipath
% observation y is (approximately) a sparse complex combination of these atoms,
% so complex OMP on D is the phase-aware test of path separability.
%
%   Gtheta, Gr : # grid points in physical angle and distance
% Outputs:
%   D        : (k*M) x (Gtheta*Gr) complex, UNIT-NORM columns
%   grid_sin : 1 x G  sine-angle of each column
%   grid_r   : 1 x G  distance [m] of each column
%   colnorm  : 1 x G  original column norms (to rescale recovered gains)

Nt = p.Nt; M = p.M; k = p.k;
th_list = linspace(p.theta_min, p.theta_max, Gtheta);   % physical angle [rad]
r_list  = linspace(p.Rmin,      p.Rmax,      Gr);
G = Gtheta*Gr;

D = zeros(k*M, G);
grid_sin = zeros(1,G);
grid_r   = zeros(1,G);

col = 0;
for ia = 1:Gtheta
    for ir = 1:Gr
        col = col + 1;
        H = near_field_channel(Nt, p.d, p.fc, p.B, M, r_list(ir), th_list(ia)); % M x Nt
        atom = zeros(k*M, 1);
        for t = 1:k
            wt = squeeze(w_2D_near(:, t, :)).';     % M x Nt
            atom((t-1)*M + (1:M)) = sum(H .* wt, 2);% M x 1
        end
        D(:, col) = atom;
        grid_sin(col) = sin(th_list(ia));
        grid_r(col)   = r_list(ir);
    end
end

colnorm = sqrt(sum(abs(D).^2, 1));
colnorm(colnorm==0) = 1;
D = D ./ colnorm;                                    % unit-norm columns
end
