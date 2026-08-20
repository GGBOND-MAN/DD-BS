% ========================================================================
% GATE C -- Phase-aware separability probe (the real tie-breaker)
% ------------------------------------------------------------------------
% Gate B used the baseline POWER-domain match filter and could only recover the
% LoS path (detection saturated at 1/L, flat in SNR). But the power domain
% discards the cross-subcarrier phase that the "path = structured trajectory"
% hook relies on. Gate C keeps the COMPLEX observation y = h_m * w_DDBS and runs
% complex OMP against a complex DDBS dictionary. This is the fair test of
% whether the DDBS structure carries enough information to separate multipath
% with K pilots.
%
% GO/NO-GO (decides theory-vs-pivot):
%   * If complex OMP reliably recovers >=2 of L paths at moderate SNR (and
%     detection improves with SNR) -> the separability is REAL and was hidden by
%     the power-domain MF. Proceed to design the (complex / super-resolution)
%     estimator; this becomes the paper's core.
%   * If complex OMP ALSO saturates at ~1/L (only LoS) even at high SNR
%     -> the DDBS observation genuinely cannot separate these paths with K=3.
%     Abandon the separation hook (revisit the whole direction).
%
% OMP here is given the ORACLE number of paths (best case) -- a generous probe.
% ========================================================================
clear; clc;

BASELINE_DIR = fullfile('..','..','baseline_code','code_nf_distance_dependent_rainbow');
addpath(genpath(BASELINE_DIR));
addpath(genpath(fileparts(mfilename('fullpath'))));

p = derisk_params();

% ---- knobs -------------------------------------------------------------
Gtheta = 121;   Gr = 15;      % complex-dictionary grid (angle x distance)
L_probe = 3;    kappa_dB = -3;
minsep  = 0.08;               % resolvable-case angular separation (sine)
SNR_list = 0:5:20;
N_iter = 50;
tol_sin  = 0.03;  tol_rrel = 0.25;
rng(11);

focus_loc = cal_loc(p.fc, p.f, p.theta1, p.theta2, p.alpha1, p.alpha2, p.M, p.k); %#ok
w_2D_near = delay_polar_2d(p.Nt, p.B, p.fc, p.M, p.d, ...
                           p.theta1, p.theta2, p.alpha1, p.alpha2, p.k);

fprintf('Building complex dictionary (%d x %d = %d atoms)...\n', Gtheta, Gr, Gtheta*Gr);
[D, grid_sin, grid_r] = build_complex_dictionary(p, w_2D_near, Gtheta, Gr);

det_prob = zeros(numel(SNR_list),1);
rmse_sin = zeros(numel(SNR_list),1);
rmse_r   = zeros(numel(SNR_list),1);

for is = 1:numel(SNR_list)
    SNR = 10^(SNR_list(is)/10);
    n_det=0; n_tot=0; se_sin=0; se_r=0; nm=0;
    for it = 1:N_iter
        s = gen_multipath_scenario(p, L_probe, kappa_dB, minsep);
        [h,~] = near_field_channel_multipath(p.Nt,p.d,p.fc,p.B,p.M,s.r,s.theta,s.g);

        % complex received vector y (k*M), no magnitude taken
        y = zeros(p.k*p.M,1);
        for t = 1:p.k
            wt = squeeze(w_2D_near(:,t,:)).';        % M x Nt
            y((t-1)*p.M+(1:p.M)) = sum(h .* wt, 2);
        end
        % complex AWGN at the specified measurement SNR
        Ps = mean(abs(y).^2);
        n  = sqrt(Ps/SNR/2)*(randn(size(y))+1j*randn(size(y)));
        yn = y + n;

        % complex OMP with oracle sparsity L_probe
        est_idx = omp_complex(D, yn, L_probe);
        est_sin = grid_sin(est_idx);
        est_r   = grid_r(est_idx);

        % greedy nearest assignment of true paths to estimates
        used = false(numel(est_idx),1);
        for l = 1:L_probe
            n_tot=n_tot+1; sin_t=s.theta_sin(l); r_t=s.r(l);
            best=inf; bi=0;
            for e=1:numel(est_idx)
                if used(e), continue; end
                dd = abs(est_sin(e)-sin_t) + abs(est_r(e)-r_t)/r_t;
                if dd<best, best=dd; bi=e; end
            end
            if bi>0
                used(bi)=true;
                if abs(est_sin(bi)-sin_t)<=tol_sin && abs(est_r(bi)-r_t)/r_t<=tol_rrel
                    n_det=n_det+1; se_sin=se_sin+(est_sin(bi)-sin_t)^2;
                    se_r=se_r+(est_r(bi)-r_t)^2; nm=nm+1;
                end
            end
        end
    end
    det_prob(is)=n_det/n_tot; rmse_sin(is)=sqrt(se_sin/max(nm,1)); rmse_r(is)=sqrt(se_r/max(nm,1));
    fprintf('(C) SNR=%2d dB : det.prob=%.2f | RMSE sin=%.4f | RMSE r=%.2f m\n', ...
            SNR_list(is), det_prob(is), rmse_sin(is), rmse_r(is));
end

save('derisk_gateC_results.mat','SNR_list','det_prob','rmse_sin','rmse_r', ...
     'L_probe','kappa_dB','N_iter','Gtheta','Gr');

figure; hold on; box on; grid on;
plot(SNR_list, det_prob, '-o','LineWidth',1.6);
yline(1/L_probe,'--','1/L (LoS only)');
ylim([0 1]); xlabel('SNR (dB)'); ylabel('Detection probability');
title(sprintf('(C) Complex OMP multipath separability (L=%d, NLoS/LoS=%d dB, K=%d)', ...
      L_probe, kappa_dB, p.k));

% ------------------------------------------------------------------------
function idx = omp_complex(D, y, K)
% Plain complex Orthogonal Matching Pursuit. Returns K selected column indices.
r = y; idx = [];
for kk = 1:K
    c = abs(D' * r);
    c(idx) = -inf;                 % don't re-pick
    [~, j] = max(c);
    idx = [idx, j];                %#ok<AGROW>
    Ds = D(:, idx);
    x  = Ds \ y;                   % least squares on current support
    r  = y - Ds * x;
end
end
