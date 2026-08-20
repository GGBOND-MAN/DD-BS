% ========================================================================
% GATE D -- Does phase-aware path recovery translate into RATE? (currency test)
% ------------------------------------------------------------------------
% Gate C showed complex OMP recovers LoS + ~1 dominant NLoS (but hits a
% structural K=3 ceiling on recovering ALL paths). The question that decides
% whether the direction is worth a paper is not "how many paths" but "how much
% of the rate gap does recovering the dominant paths close".
%
% Pipeline (a real receiver could do exactly this):
%   1. K=3 DDBS pilots -> complex measurement y.
%   2. complex OMP on the DDBS dictionary -> support + gains.
%   3. reconstruct the wideband channel hhat from the recovered paths.
%   4. per-subcarrier analog MRT beam from hhat  (interpolates single-beam LoS
%      <-> Perfect CSI as more paths are recovered).
%
% Reported: closure = (rate_D - rate_MF) / (rate_perf - rate_MF), i.e. the
% fraction of the baseline-to-optimal gap that phase-aware recovery recovers.
%
% DECISION RULE (my proposal):
%   * closure >= ~60% at a realistic kappa (-10..-6 dB) -> a viable (if modest)
%     rate-oriented contribution: "phase-aware multipath beamforming from K=3
%     DDBS pilots". Proceed, but position carefully vs 2505.08267.
%   * closure small / only appears at kappa ~ 0 dB -> multipath is not worth it
%     here (structural ceiling + off-brand regime). Pivot (UPA/3D DDBS).
% ========================================================================
clear; clc;

BASELINE_DIR = fullfile('..','..','baseline_code','code_nf_distance_dependent_rainbow');
addpath(genpath(BASELINE_DIR));
addpath(genpath(fileparts(mfilename('fullpath'))));

p = derisk_params();

% ---- knobs -------------------------------------------------------------
Gtheta = 121;  Gr = 15;         % complex dictionary grid
g1d = 256;     g2d = 10;        % power-domain MF dictionary (baseline single-beam)
L_true  = 3;
Lest    = 3;                    % # atoms OMP extracts (oracle-generous)
kappa_sweep = [-10 -6 -3 0];
SNR_dB  = 15;  SNR_t = 10^(SNR_dB/10);
N_iter  = 100;
rng(21);

focus_loc = cal_loc(p.fc, p.f, p.theta1, p.theta2, p.alpha1, p.alpha2, p.M, p.k);
w_2D_near = delay_polar_2d(p.Nt, p.B, p.fc, p.M, p.d, ...
                           p.theta1, p.theta2, p.alpha1, p.alpha2, p.k);
fprintf('Building complex dictionary + MF dictionary...\n');
[D, grid_sin, grid_r, colnorm] = build_complex_dictionary(p, w_2D_near, Gtheta, Gr);
[match_y, g1_list, g2_list]    = build_match_dictionary(p, w_2D_near, g1d, g2d);

rate_MF   = zeros(numel(kappa_sweep),1);
rate_D    = zeros(numel(kappa_sweep),1);
rate_perf = zeros(numel(kappa_sweep),1);

for ik = 1:numel(kappa_sweep)
    kap = kappa_sweep(ik);
    acc = zeros(1,3);
    for it = 1:N_iter
        s = gen_multipath_scenario(p, L_true, kap);      % random geometry (honest)
        [h,~] = near_field_channel_multipath(p.Nt,p.d,p.fc,p.B,p.M,s.r,s.theta,s.g);

        % --- baseline single-beam match filter rate ---
        r_mf = training_near_rainbow_match(p.Nt,p.B,p.fc,p.f,p.M,p.d,h, ...
                    w_2D_near,focus_loc,SNR_t,SNR_dB,p.Q,p.k,match_y,g1_list,g2_list);
        acc(1) = acc(1) + r_mf(p.k);

        % --- phase-aware recovery + reconstructed-channel beamforming ---
        y = zeros(p.k*p.M,1);
        for t = 1:p.k
            wt = squeeze(w_2D_near(:,t,:)).';
            y((t-1)*p.M+(1:p.M)) = sum(h .* wt, 2);
        end
        Ps = mean(abs(y).^2);
        yn = y + sqrt(Ps/SNR_t/2)*(randn(size(y))+1j*randn(size(y)));
        [idx, x] = omp_complex(D, yn, Lest);
        g = x(:) ./ colnorm(idx).';                       % channel-space gains
        hhat = zeros(p.M, p.Nt);
        for e = 1:numel(idx)
            Hp = near_field_channel(p.Nt,p.d,p.fc,p.B,p.M, grid_r(idx(e)), asin(grid_sin(idx(e))));
            hhat = hhat + g(e)*Hp;
        end
        rD = 0;
        for m = 1:p.M
            w = exp(1j*angle(hhat(m,:)'))/sqrt(p.Nt);
            rD = rD + log2(1 + SNR_t*abs(h(m,:)*w)^2)/p.M;
        end
        acc(2) = acc(2) + rD;

        % --- perfect CSI ---
        acc(3) = acc(3) + perfect_csi_rate(h, SNR_t, p.Nt, p.M);
    end
    acc = acc/N_iter;
    rate_MF(ik)=acc(1); rate_D(ik)=acc(2); rate_perf(ik)=acc(3);
    closure = (acc(2)-acc(1))/max(acc(3)-acc(1),1e-9);
    fprintf('kappa=%3d dB : MF %.3f | Recovered %.3f | Perfect %.3f | closure %.0f%%\n', ...
            kap, acc(1), acc(2), acc(3), 100*closure);
end

save('derisk_gateD_results.mat','kappa_sweep','rate_MF','rate_D','rate_perf','SNR_dB','N_iter');

figure; hold on; box on; grid on;
plot(kappa_sweep, rate_perf,'--k','LineWidth',1.6);
plot(kappa_sweep, rate_D,   '-o','LineWidth',1.6);
plot(kappa_sweep, rate_MF,  '-s','LineWidth',1.6);
xlabel('NLoS/LoS power ratio (dB)'); ylabel('Average rate (bit/s/Hz)');
title(sprintf('(D) Rate closure via phase-aware recovery (L=%d, SNR=%d dB, K=%d)', L_true, SNR_dB, p.k));
legend('Perfect CSI','Recovered multipath beamforming','Baseline single-beam MF','Location','northwest');

% ------------------------------------------------------------------------
function [idx, x] = omp_complex(D, y, K)
r = y; idx = [];
for kk = 1:K
    c = abs(D'*r); c(idx) = -inf;
    [~,j] = max(c); idx = [idx, j]; %#ok<AGROW>
    Ds = D(:,idx); x = Ds \ y; r = y - Ds*x;
end
end

function R = perfect_csi_rate(h, SNR_t, Nt, M)
R = 0;
for m = 1:M
    w = exp(1j*angle(h(m,:)'))/sqrt(Nt);
    R = R + log2(1 + SNR_t*abs(h(m,:)*w)^2)/M;
end
end
