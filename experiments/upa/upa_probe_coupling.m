% ========================================================================
% UPA / 3D-DDBS feasibility probe -- CORE: the near-field cross-coupling wall
% ------------------------------------------------------------------------
% The whole 3D-DDBS idea needs a UPA to FOCUS a near-field beam at an arbitrary
% (Wx, Wy, r). A cheap SEPARABLE true-time-delay architecture (per-axis delays,
% N_x + N_y TDs) can realize the two angular linear terms and the two diagonal
% curvatures -- but NOT the non-separable near-field cross term nx*ny*Wx*Wy/r.
% Only a fully-connected TD array (N_x*N_y TDs, expensive) can.
%
% This probe quantifies the focusing loss that ignoring the cross term costs,
% over the angular region and across distances. It decides the paper's story:
%
%   * loss SMALL (< ~1-3 dB over most of the +/-60 deg region at practical r)
%     -> the cheap SEPARABLE 3D-DDBS works; paper = separable design + coverage
%        / pilot-overhead analysis for UPA beam training.  GREEN.
%   * loss LARGE over much of the region
%     -> near-field cross-coupling forces a joint / (semi-)fully-connected design;
%        paper = "the near-field 3D coupling problem + a low-cost correction".
%        Harder but MORE novel.  Either outcome is publishable -- the probe just
%        tells you WHICH narrative, and de-risks the core physics.
%
% "full-Taylor" (cross included) is plotted too: it should be ~0 dB loss, which
% (a) validates the 2nd-order model and (b) isolates the cross term as the sole
% culprit behind the separable-architecture loss.
% ========================================================================
clear; clc;

% ---- UPA system -- a REAL XL-MIMO planar aperture -----------------------
% NOTE: the near-field cross term only matters when BOTH axes have a large
% aperture, i.e. r < ~Rayleigh distance. Splitting a small Nt into a square grid
% (e.g. 16x16 -> 8 cm aperture, Rayleigh ~1.3 m) is FAR-FIELD at metre ranges and
% makes this probe trivially "green" for the wrong reason. Use a large UPA and
% parameterise distance by the Rayleigh distance.
Nx = 64; Ny = 64;                 % Nt = 4096 (per-axis aperture ~32 cm)
fc = 30e9; c = 3e8; d = (c/fc)/2;
f_eval = fc;                      % evaluate focusing at centre freq (pure geometry)
k = 2*pi*f_eval/c;
Nt = Nx*Ny;
aperture = (Nx-1)*d;
Rayleigh = 2*aperture^2/(c/fc);

smax = sin(pi/3);                 % +/-60 deg angular coverage  -> |W| <= 0.866
G = 41;                           % grid per axis
Wg = linspace(-smax, smax, G);
frac_list = [0.05 0.12 0.30];     % r as a fraction of Rayleigh (deep / mod / shallow NF)
r_list = frac_list * Rayleigh;

fprintf('UPA %dx%d (Nt=%d), aperture=%.0f cm, Rayleigh=%.1f m, region |W|<=%.2f\n', ...
        Nx, Ny, Nt, aperture*100, Rayleigh, smax);

loss_sep  = nan(G, G, numel(r_list));   % separable architecture (no cross term)
loss_full = nan(G, G, numel(r_list));   % fully-connected (cross included) -> ~0

for ir = 1:numel(r_list)
    r = r_list(ir);
    for ix = 1:G
        for iy = 1:G
            Wx = Wg(ix); Wy = Wg(iy);
            if Wx^2 + Wy^2 > smax^2, continue; end   % keep the angular disk
            de   = upa_delta_exact (Nx, Ny, d, r, Wx, Wy);
            dsep = upa_delta_taylor(Nx, Ny, d, r, Wx, Wy, false);
            dful = upa_delta_taylor(Nx, Ny, d, r, Wx, Wy, true);
            % array gain of a beam that pre-compensates 'd*' against the exact steering:
            g_sep = abs(mean(exp(-1j*k*(de - dsep))))^2;
            g_ful = abs(mean(exp(-1j*k*(de - dful))))^2;
            loss_sep (ix,iy,ir) = 10*log10(g_sep);
            loss_full(ix,iy,ir) = 10*log10(g_ful);
        end
    end
    L = loss_sep(:,:,ir); L = L(~isnan(L));
    Lf = loss_full(:,:,ir); Lf = Lf(~isnan(Lf));
    fprintf(['r=%6.1f m (%.2f Rayl) | SEP loss: median %6.2f dB, worst %6.2f dB, ' ...
             'frac<3dB %3.0f%%  || full-Taylor(cross) median %5.2f dB\n'], ...
            r, frac_list(ir), median(L), min(L), 100*mean(L>-3), median(Lf));
end

save('upa_probe_coupling_results.mat','Wg','r_list','loss_sep','loss_full','Nx','Ny','fc');

% ---- heatmaps: separable-architecture focusing loss (dB) ----------------
figure;
for ir = 1:numel(r_list)
    subplot(1,numel(r_list),ir);
    imagesc(Wg, Wg, loss_sep(:,:,ir).'); axis xy equal tight;
    caxis([-12 0]); colormap(flipud(hot)); colorbar;
    xlabel('\Omega_x = sin\vartheta cos\phi'); ylabel('\Omega_y = sin\vartheta sin\phi');
    title(sprintf('r = %.1f m (%.2f Rayleigh)', r_list(ir), frac_list(ir)));
end
sgtitle('UPA 3D-DDBS separable-TD focusing loss (dB): the deep-near-field, double-off-broadside corner is the wall');
