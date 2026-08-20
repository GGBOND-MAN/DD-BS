% ========================================================================
% GATE A -- Does the baseline actually break under multipath?
% ------------------------------------------------------------------------
% Reproduces the baseline DDBS beam training (on-grid + single-location match
% filter) but drives it with an L-path near-field channel. Produces two plots:
%
%   (A1) Average rate vs number of paths L      (fixed NLoS/LoS ratio, fixed SNR)
%   (A2) Average rate vs NLoS/LoS power ratio    (fixed L,               fixed SNR)
%
% GO/NO-GO for the whole research direction:
%   * If DDBS on-grid and DDBS+match-filter fall CLEARLY below Perfect CSI as L
%     grows / as the NLoS ratio rises  -> the motivation is real. PROCEED.
%   * If they barely move (baseline stays near-optimal even with strong NLoS)
%     -> the multipath weakness is not consequential in this regime. STOP and
%     re-scope (e.g. blockage / multi-beam-combining scenario, or UPA).
%
% The training functions training_near_rainbow_2d.m and
% training_near_rainbow_match.m are the UNMODIFIED baseline functions: they
% accept a generic per-subcarrier channel h (M x Nt), so feeding them a
% multipath h is a faithful, minimal change.
% ========================================================================
clear; clc;

% ---- point this at the extracted baseline code folder ------------------
BASELINE_DIR = fullfile('..','..','baseline_code','code_nf_distance_dependent_rainbow');
addpath(genpath(BASELINE_DIR));
addpath(genpath(fileparts(mfilename('fullpath'))));   % this de-risk folder

p = derisk_params();

% ---- experiment knobs (tune for speed vs fidelity) ---------------------
N_iter   = 100;          % Monte-Carlo realisations per point
SNR_dB   = 15;           % fixed operating SNR
g1d      = 256;          % match-filter dictionary angle grid  (baseline: 1024)
g2d      = 10;           % match-filter dictionary alpha grid
L_sweep      = 1:5;      % (A1)
kappa_A1_dB  = -6;       % (A1) fixed NLoS/LoS power ratio
kappa_sweep  = [-20 -15 -10 -6 -3 0];   % (A2)
L_A2         = 3;        % (A2) fixed number of paths
rng(2024);               % reproducibility

SNR_t = 10^(SNR_dB/10);

% ---- build DDBS training beams + focus grid + MF dictionary once --------
focus_loc = cal_loc(p.fc, p.f, p.theta1, p.theta2, p.alpha1, p.alpha2, p.M, p.k);
w_2D_near = delay_polar_2d(p.Nt, p.B, p.fc, p.M, p.d, ...
                           p.theta1, p.theta2, p.alpha1, p.alpha2, p.k);
fprintf('Building match-filter dictionary (%d x %d)...\n', g1d, g2d);
[match_y, g1_list, g2_list] = build_match_dictionary(p, w_2D_near, g1d, g2d);

% ======================= (A1) rate vs L =================================
rate_ongrid_L = zeros(numel(L_sweep),1);
rate_mf_L     = zeros(numel(L_sweep),1);
rate_perf_L   = zeros(numel(L_sweep),1);
for iL = 1:numel(L_sweep)
    L = L_sweep(iL);
    acc = zeros(1,3);
    for it = 1:N_iter
        s = gen_multipath_scenario(p, L, kappa_A1_dB);
        [h, hc] = near_field_channel_multipath(p.Nt, p.d, p.fc, p.B, p.M, ...
                                               s.r, s.theta, s.g);
        r_og = training_near_rainbow_2d(p.Nt,p.B,p.fc,p.f,p.M,p.d,h, ...
                    w_2D_near,focus_loc,SNR_t,SNR_dB,p.Q,p.k);
        r_mf = training_near_rainbow_match(p.Nt,p.B,p.fc,p.f,p.M,p.d,h, ...
                    w_2D_near,focus_loc,SNR_t,SNR_dB,p.Q,p.k,match_y,g1_list,g2_list);
        acc(1) = acc(1) + r_og(p.k);
        acc(2) = acc(2) + r_mf(p.k);
        acc(3) = acc(3) + perfect_csi_rate(h, SNR_t, p.Nt, p.M);
    end
    acc = acc / N_iter;
    rate_ongrid_L(iL)=acc(1); rate_mf_L(iL)=acc(2); rate_perf_L(iL)=acc(3);
    fprintf('(A1) L=%d : on-grid %.3f | MF %.3f | perfect %.3f\n', ...
            L, acc(1), acc(2), acc(3));
end

% ======================= (A2) rate vs kappa =============================
rate_ongrid_k = zeros(numel(kappa_sweep),1);
rate_mf_k     = zeros(numel(kappa_sweep),1);
rate_perf_k   = zeros(numel(kappa_sweep),1);
for ik = 1:numel(kappa_sweep)
    kappa_dB = kappa_sweep(ik);
    acc = zeros(1,3);
    for it = 1:N_iter
        s = gen_multipath_scenario(p, L_A2, kappa_dB);
        [h, hc] = near_field_channel_multipath(p.Nt, p.d, p.fc, p.B, p.M, ...
                                               s.r, s.theta, s.g);
        r_og = training_near_rainbow_2d(p.Nt,p.B,p.fc,p.f,p.M,p.d,h, ...
                    w_2D_near,focus_loc,SNR_t,SNR_dB,p.Q,p.k);
        r_mf = training_near_rainbow_match(p.Nt,p.B,p.fc,p.f,p.M,p.d,h, ...
                    w_2D_near,focus_loc,SNR_t,SNR_dB,p.Q,p.k,match_y,g1_list,g2_list);
        acc(1) = acc(1) + r_og(p.k);
        acc(2) = acc(2) + r_mf(p.k);
        acc(3) = acc(3) + perfect_csi_rate(h, SNR_t, p.Nt, p.M);
    end
    acc = acc / N_iter;
    rate_ongrid_k(ik)=acc(1); rate_mf_k(ik)=acc(2); rate_perf_k(ik)=acc(3);
    fprintf('(A2) kappa=%d dB : on-grid %.3f | MF %.3f | perfect %.3f\n', ...
            kappa_dB, acc(1), acc(2), acc(3));
end

save('derisk_gateA_results.mat', ...
     'L_sweep','rate_ongrid_L','rate_mf_L','rate_perf_L','kappa_A1_dB', ...
     'kappa_sweep','rate_ongrid_k','rate_mf_k','rate_perf_k','L_A2', ...
     'SNR_dB','N_iter');

% ============================== plots ===================================
figure; hold on; box on; grid on;
plot(L_sweep, rate_perf_L,   '--k','LineWidth',1.6);
plot(L_sweep, rate_ongrid_L, '-s','LineWidth',1.6);
plot(L_sweep, rate_mf_L,     '->','LineWidth',1.6);
xlabel('Number of paths L'); ylabel('Average rate (bit/s/Hz)');
title(sprintf('(A1) Rate vs L   (NLoS/LoS = %d dB, SNR = %d dB)', kappa_A1_dB, SNR_dB));
legend('Perfect CSI','DDBS on-grid','DDBS + match filter','Location','southwest');

figure; hold on; box on; grid on;
plot(kappa_sweep, rate_perf_k,   '--k','LineWidth',1.6);
plot(kappa_sweep, rate_ongrid_k, '-s','LineWidth',1.6);
plot(kappa_sweep, rate_mf_k,     '->','LineWidth',1.6);
xlabel('NLoS/LoS power ratio (dB)'); ylabel('Average rate (bit/s/Hz)');
title(sprintf('(A2) Rate vs NLoS ratio   (L = %d, SNR = %d dB)', L_A2, SNR_dB));
legend('Perfect CSI','DDBS on-grid','DDBS + match filter','Location','southwest');

% ------------------------------------------------------------------------
function R = perfect_csi_rate(h, SNR_t, Nt, M)
% Per-subcarrier analog (phase-only) MRT upper bound, matching the baseline's
% Perfect-CSI convention (wc_opt = exp(1j*angle(h'))/sqrt(Nt)).
R = 0;
for m = 1:M
    w  = exp(1j*angle(h(m,:)'))/sqrt(Nt);
    ag = abs(h(m,:)*w)^2;
    R  = R + log2(1 + SNR_t*ag)/M;
end
end
