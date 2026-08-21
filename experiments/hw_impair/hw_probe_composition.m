% ========================================================================
% (B) COMPOSITION: contribution A applied ON TOP of the AoSA architecture
% ------------------------------------------------------------------------
% THE paper figure. It supports the only claim that is defensible without
% reproducing the OJ-COMS scheme:
%
%   "On an AoSA TTD-PS array, the proposed reparameterization lowers the TTD
%    resolution needed for a given performance by ~4 bits."
%
% This is a COMPOSITION result, not a competition result. Both arms use the same
% sub-array-shared TTD architecture at the OJ-COMS operating point (P = 2, i.e.
% TTDs halved to 128 for Nt = 256, which is the setting they report). The arms
% differ ONLY in how the per-antenna phase is split between TTD and PS:
%
%   arch='ojcoms' : TTD carries the full 2*pi*f_m*tau (their eqs (13)-(14));
%                   the PS is the plain per-antenna DDBS phase term.
%   arch='shared' : identical TTD sharing, but the PS carries the exact
%                   frequency-independent 2*pi*fc*tau_n and the TTD only the
%                   residual 2*pi*(f_m-fc)*tau_n  <- contribution A.
%
% Neither arm reproduces their full scheme (parameter redesign, sector-shift,
% distance interleaving, ~12 pilots), so NEITHER may be reported as their
% performance. The claim is strictly the DELTA BETWEEN THE TWO ARMS at a fixed
% architecture -- which is exactly the contribution.
%
% Expected from the focusing-gain surface at P=2:
%   B_td   :   Inf     14      13      12      11      10
%   ojcoms :  0.478  0.431   0.301   0.055   0.040   0.054   <- collapses at ~12
%   +A     :  0.668  0.668   0.667   0.665   0.657   0.624   <- flat to ~10
% i.e. the same performance is held with ~4 fewer bits.
% ========================================================================
clear; clc;

BASELINE_DIR = fullfile('..','..','baseline_code','code_nf_distance_dependent_rainbow');
addpath(genpath(BASELINE_DIR));
addpath(genpath(fileparts(mfilename('fullpath'))));

Nt=256; fc=30e9; B=5e9; M=1024; c=3e8; d=(c/fc)/2; Q=1; K=3;
f = fc + B/M*((1:M)-1-(M-1)/2);
Rmin=5; Rmax=30; th_min=-pi/3; th_max=pi/3;
theta1=349/11; theta2=-30; alpha1=-427/800; alpha2=1859/3200;
focus_loc = cal_loc(fc,f,theta1,theta2,alpha1,alpha2,M,K);
nn=(-(Nt-1)/2:(Nt-1)/2)';
tau_rng = max((nn*d*theta1-(nn*d).^2*alpha1)/c) - min((nn*d*theta1-(nn*d).^2*alpha1)/c);

P = 2;                       % OJ-COMS operating point: TTDs halved (128 for Nt=256)
SNR_dB=15; SNR_t=10^(SNR_dB/10); N_iter=100; rng(11);
CH = cell(N_iter,1);
for it=1:N_iter
    CH{it} = near_field_channel(Nt,d,fc,B,M, Rmin+rand*(Rmax-Rmin), th_min+rand*(th_max-th_min));
end

w_id = ddbs_beam_arch(Nt,B,fc,M,d,theta1,theta2,alpha1,alpha2,K,'full',1,Inf,Inf,[]);
r_ideal = mean(cellfun(@(h) run1(h,w_id,'full',1,[],Nt,B,fc,f,M,d,focus_loc,SNR_t,SNR_dB,Q,K), CH));
fprintf('Nt=%d, P=%d (%d shared TTDs), SNR=%d dB, K=%d\n', Nt, P, Nt/P, SNR_dB, K);
fprintf('ideal full-TTD DDBS = %.3f bit/s/Hz\n\n', r_ideal);

Btd_list = [Inf 15 14 13 12 11 10];
fprintf('%-26s |','architecture'); fprintf('%8s', string(Btd_list)); fprintf('   <- B_td\n');
res = zeros(2,numel(Btd_list));
labs = {'AoSA form (their split)','AoSA + reparam. (ours)'};
archs = {'ojcoms','shared'};
for a = 1:2
    fprintf('%-26s |', labs{a});
    for j = 1:numel(Btd_list)
        Btd = Btd_list(j);
        lsb = tau_rng/2^Btd;  if ~isfinite(Btd), lsb = []; end
        w = ddbs_beam_arch(Nt,B,fc,M,d,theta1,theta2,alpha1,alpha2,K,archs{a},P,Inf,Inf,lsb);
        r = mean(cellfun(@(h) run1(h,w,archs{a},P,lsb,Nt,B,fc,f,M,d,focus_loc,SNR_t,SNR_dB,Q,K), CH));
        res(a,j) = r; fprintf('%8.3f', r);
    end
    fprintf('\n');
end

% bits needed to stay within 5% of each arm's own infinite-resolution rate
fprintf('\nBits needed to stay within 5%% of the same arm at infinite resolution:\n');
for a = 1:2
    ref = res(a,1); need = NaN;
    for j = numel(Btd_list):-1:2
        if res(a,j) >= 0.95*ref, need = Btd_list(j); end
    end
    fprintf('  %-26s : %s bits\n', labs{a}, string(need));
end
fprintf('\n-> the gap between those two numbers is the contribution.\n');
save('hw_probe_composition_results.mat','Btd_list','res','P','r_ideal');

figure; hold on; box on; grid on;
plot(Btd_list(2:end), res(1,2:end), '-s','LineWidth',1.7);
plot(Btd_list(2:end), res(2,2:end), '-o','LineWidth',1.7);
set(gca,'XDir','reverse'); xlabel('TD resolution (bits)'); ylabel('Average rate (bit/s/Hz)');
title(sprintf('Composition on the AoSA architecture (P=%d, %d shared TTDs)', P, Nt/P));
legend(labs,'Location','southwest');

% ------------------------------------------------------------------------
function r = run1(h,w_train,arch,P,lsb,Nt,B,fc,f,M,d,focus_loc,SNR_t,SNR_dB,Q,K)
best = -inf;
for idx = 1:K
    y = zeros(M,Q);
    for m = 1:M
        y(m,:) = awgn( repmat( h(m,:)*w_train(:,idx,m), [1,Q] ), SNR_dB*2/sqrt(3) );
    end
    [~,i] = max(abs(sum(y,2)).^2);
    ws = ddbs_beam_arch(Nt,B,fc,M,d, focus_loc(i,1,idx), 0, focus_loc(i,2,idx), 0, 1, arch,P,Inf,Inf,lsb);
    t = 0;
    for m = 1:M
        t = t + log2(1 + SNR_t*abs(h(m,:)*ws(:,1,m))^2)/M;
    end
    best = max(best,t);
end
r = best;
end
