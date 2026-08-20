% ========================================================================
% GATE B -- Does the DDBS structure make multipath components SEPARABLE?
% ------------------------------------------------------------------------
% This is the test of the actual novelty premise. It has two parts:
%
%   (B1) Separability picture: for one well-separated L-path realisation, plot
%        the match-filter correlation surface over the (theta, alpha) grid and
%        overlay the true path locations. If the surface shows L distinct maxima
%        near the true paths, the paths ARE separable through the DDBS
%        observation with only K pilots.
%
%   (B2) Quantitative probe: a crude greedy multi-peak extractor on that same
%        surface. Sweep SNR and report, over many realisations,
%          - detection probability (fraction of true paths localised within tol)
%          - localisation RMSE in sine-angle and in distance [m].
%
% GO/NO-GO:
%   * If >=2 paths are reliably separated at moderate SNR with K=3 pilots
%     -> the "path = structured trajectory, resolvable by beam-split" premise
%        holds. PROCEED to the theoretical estimator design.
%   * If the surface fuses the paths into one blob / greedy extraction fails
%     even at high SNR -> the structural hook needs more pilots or a smarter
%     estimator; revisit the framing before investing in theory.
%
% NOTE: the greedy extractor here is deliberately dumb (peak-pick + null). It is
% NOT the proposed algorithm -- it only checks whether the *information* to
% separate paths is present in the DDBS observation.
% ========================================================================
clear; clc;

BASELINE_DIR = fullfile('..','..','baseline_code','code_nf_distance_dependent_rainbow');
addpath(genpath(BASELINE_DIR));
addpath(genpath(fileparts(mfilename('fullpath'))));

p = derisk_params();

% ---- knobs -------------------------------------------------------------
g1d   = 256;             % angle grid (sine)   -- controls memory of dictionary
g2d   = 40;              % alpha grid (finer than Gate A to resolve distance)
minsep = 0.08;           % min sine-angle separation between paths (resolvable case)
rng(7);

focus_loc = cal_loc(p.fc, p.f, p.theta1, p.theta2, p.alpha1, p.alpha2, p.M, p.k);
w_2D_near = delay_polar_2d(p.Nt, p.B, p.fc, p.M, p.d, ...
                           p.theta1, p.theta2, p.alpha1, p.alpha2, p.k);
fprintf('Building fine dictionary (%d x %d)...\n', g1d, g2d);
[match_y, g1_list, g2_list] = build_match_dictionary(p, w_2D_near, g1d, g2d);
match_y = squeeze(match_y);            % g1 x g2 x k x M   (Q=1 dropped)

% ================= (B1) single-realisation separability picture =========
L_show = 3;  kappa_show_dB = -3;  SNR_show_dB = 20;
s = gen_multipath_scenario(p, L_show, kappa_show_dB, minsep);
[h,~] = near_field_channel_multipath(p.Nt,p.d,p.fc,p.B,p.M,s.r,s.theta,s.g);
y_all = ddbs_received_power(h, w_2D_near, p, SNR_show_dB);
surf_corr = correlation_surface(y_all, match_y, p.k, p.M);

figure; imagesc(g2_list, g1_list, surf_corr); axis xy; hold on;
xlabel('\alpha  (distance ring)'); ylabel('sin\theta');
title(sprintf('(B1) MF correlation surface, L=%d, NLoS/LoS=%d dB, SNR=%d dB', ...
      L_show, kappa_show_dB, SNR_show_dB));
colorbar;
for l = 1:L_show
    a_true = (1 - s.theta_sin(l)^2) / (2*s.r(l));      % true alpha of path l
    plot(a_true, s.theta_sin(l), 'rx', 'MarkerSize', 12, 'LineWidth', 2);
end
legend('true path locations','Location','northeast');

% ================= (B2) detection / RMSE vs SNR =========================
SNR_list = 0:5:20;
L_probe  = 3;  kappa_probe_dB = -3;  N_iter = 50;
det_prob = zeros(numel(SNR_list),1);
rmse_sin = zeros(numel(SNR_list),1);
rmse_r   = zeros(numel(SNR_list),1);

tol_sin  = 0.03;         % angle match tolerance (sine)
tol_rrel = 0.25;         % distance match tolerance (relative)

for is = 1:numel(SNR_list)
    SNR_dB = SNR_list(is);
    n_det = 0; n_tot = 0; se_sin = 0; se_r = 0; n_match = 0;
    for it = 1:N_iter
        s = gen_multipath_scenario(p, L_probe, kappa_probe_dB, minsep);
        [h,~] = near_field_channel_multipath(p.Nt,p.d,p.fc,p.B,p.M,s.r,s.theta,s.g);
        y_all = ddbs_received_power(h, w_2D_near, p, SNR_dB);
        surf_corr = correlation_surface(y_all, match_y, p.k, p.M);

        est = greedy_peaks(surf_corr, g1_list, g2_list, L_probe);  % L_probe x 2 [sin, alpha]
        est_r = (1 - est(:,1).^2) ./ (2*est(:,2));                 % -> distance

        % greedy nearest assignment of each TRUE path to an estimated peak
        used = false(size(est,1),1);
        for l = 1:L_probe
            n_tot = n_tot + 1;
            r_true = s.r(l); sin_true = s.theta_sin(l);
            best = inf; bi = 0;
            for e = 1:size(est,1)
                if used(e), continue; end
                dd = abs(est(e,1)-sin_true) + abs(est_r(e)-r_true)/r_true;
                if dd < best, best = dd; bi = e; end
            end
            if bi>0
                used(bi) = true;
                if abs(est(bi,1)-sin_true)<=tol_sin && ...
                   abs(est_r(bi)-r_true)/r_true<=tol_rrel
                    n_det = n_det + 1;
                    se_sin = se_sin + (est(bi,1)-sin_true)^2;
                    se_r   = se_r   + (est_r(bi)-r_true)^2;
                    n_match = n_match + 1;
                end
            end
        end
    end
    det_prob(is) = n_det / n_tot;
    rmse_sin(is) = sqrt(se_sin / max(n_match,1));
    rmse_r(is)   = sqrt(se_r   / max(n_match,1));
    fprintf('(B2) SNR=%2d dB : det.prob=%.2f | RMSE sin=%.4f | RMSE r=%.2f m\n', ...
            SNR_dB, det_prob(is), rmse_sin(is), rmse_r(is));
end

save('derisk_gateB_results.mat','SNR_list','det_prob','rmse_sin','rmse_r', ...
     'L_probe','kappa_probe_dB','N_iter');

figure; hold on; box on; grid on;
yyaxis left;  plot(SNR_list, det_prob, '-o','LineWidth',1.6);
ylabel('Detection probability'); ylim([0 1]);
yyaxis right; plot(SNR_list, rmse_r, '-s','LineWidth',1.6);
ylabel('Distance RMSE (m)');
xlabel('SNR (dB)');
title(sprintf('(B2) Multipath separability via DDBS  (L=%d, NLoS/LoS=%d dB, K=%d)', ...
      L_probe, kappa_probe_dB, p.k));

% ======================================================================
function y_all = ddbs_received_power(h, w_2D_near, p, SNR_dB)
% Noisy received power over pilots x subcarriers, mirroring the baseline
% training functions (awgn(..., SNR_dB*2/sqrt(3)), sum over Q, abs^2).
y_all = zeros(p.k, p.M);
for idx = 1:p.k
    yk = zeros(p.M, p.Q);
    for m = 1:p.M
        yk(m,:) = awgn(repmat(h(m,:)*w_2D_near(:,idx,m),[1,p.Q]), SNR_dB*2/sqrt(3));
    end
    y_all(idx,:) = abs(sum(yk,2)).^2;
end
end

function S = correlation_surface(y_all, match_y, k, M)
% S(theta,alpha) = sum_{k,m} y_all(k,m) * dictionary_power(theta,alpha,k,m)
[g1,g2,~,~] = size(match_y);
S = zeros(g1,g2);
for t = 1:k
    for m = 1:M
        S = S + y_all(t,m) * match_y(:,:,t,m);
    end
end
end

function est = greedy_peaks(S, g1_list, g2_list, L)
% Greedy peak pick + local null. Returns L x 2 rows [sin_theta, alpha].
% Deliberately simple -- an information-presence probe, not the final estimator.
[g1,g2] = size(S);
est = zeros(L,2);
wtheta = max(3, round(0.02/ mean(diff(g1_list))));   % null half-window (angle)
walpha = max(1, round(g2/12));                        % null half-window (alpha)
Swork = S;
for l = 1:L
    [~,ind] = max(Swork(:));
    [i,j] = ind2sub([g1,g2], ind);
    est(l,:) = [g1_list(i), g2_list(j)];
    i0=max(1,i-wtheta); i1=min(g1,i+wtheta);
    j0=max(1,j-walpha); j1=min(g2,j+walpha);
    Swork(i0:i1, j0:j1) = -inf;
end
end
