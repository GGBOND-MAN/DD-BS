% ========================================================================
% (B) Hardware-impairment robustness of DDBS -- PROBE 2: beam-training RATE
% ------------------------------------------------------------------------
% Ties impairments to the paper's own metric (average rate), reusing the baseline
% K=3 DDBS beam training UNCHANGED, but driving it with impaired training beams.
% Single-path LoS channel (isolates the hardware effect from multipath).
%
% The system still MAPS the argmax subcarrier to its INTENDED focus (focus_loc
% from the ideal design); impairments move where the beam actually points, so the
% mismatch shows up directly as rate loss -- exactly the robustness gap a robust
% design must close.
% ========================================================================
clear; clc;

BASELINE_DIR = fullfile('..','..','baseline_code','code_nf_distance_dependent_rainbow');
addpath(genpath(BASELINE_DIR));
addpath(genpath(fileparts(mfilename('fullpath'))));

Nt=256; fc=30e9; B=5e9; M=1024; c=3e8; d=(c/fc)/2; Q=1;
f = fc + B/M*((1:M)-1-(M-1)/2);
Rmin=5; Rmax=30; th_min=-pi/3; th_max=pi/3;
K=3;
% baseline K=3 DDBS params (exactly Rate_snr.m)
theta1=349/11; theta2=-30; alpha1=-427/800; alpha2=1859/3200;
focus_loc = cal_loc(fc, f, theta1, theta2, alpha1, alpha2, M, K);

SNR_dB=15; SNR_t=10^(SNR_dB/10); N_iter=100; rng(3);

runcase = @(imp) montecarlo(imp, Nt,B,fc,f,M,d,theta1,theta2,alpha1,alpha2,K, ...
                             focus_loc,SNR_t,SNR_dB,Q,Rmin,Rmax,th_min,th_max,N_iter);

fprintf('Beam-training rate vs impairment (SNR=%d dB, K=%d, single-path)\n', SNR_dB, K);
[r0,rp] = runcase(struct());
fprintf('  ideal DDBS = %.3f  | Perfect CSI = %.3f  bit/s/Hz\n\n', r0, rp);

fprintf('-- TD bits --\n');
for Btd=[Inf 12 10 8]
    r = runcase(struct('Btd',Btd));
    fprintf('  B_td=%4s : rate %.3f  (%.0f%% of ideal)\n', num2str(Btd), r, 100*r/r0);
end
fprintf('-- TD jitter sigma_tau [ps] --\n');
for st=[0 1 2 5]
    r = runcase(struct('sigma_tau',st*1e-12));
    fprintf('  sigma_tau=%4.1f ps : rate %.3f  (%.0f%% of ideal)\n', st, r, 100*r/r0);
end
fprintf('-- PS bits --\n');
for Bps=[Inf 3 2]
    r = runcase(struct('Bps',Bps));
    fprintf('  B_ps=%4s : rate %.3f  (%.0f%% of ideal)\n', num2str(Bps), r, 100*r/r0);
end

% ------------------------------------------------------------------------
function [rate, rate_perf] = montecarlo(imp, Nt,B,fc,f,M,d,t1,t2,a1,a2,K, ...
                             focus_loc,SNR_t,SNR_dB,Q,Rmin,Rmax,th_min,th_max,N_iter)
racc=0; pacc=0;
for it=1:N_iter
    r = Rmin + rand*(Rmax-Rmin);
    theta = th_min + rand*(th_max-th_min);
    [h,~] = near_field_channel(Nt, d, fc, B, M, r, theta);   % M x Nt (single path)
    w_imp = ddbs_beam_impaired(Nt,B,fc,M,d,t1,t2,a1,a2,K,imp);
    rr = training_near_rainbow_2d(Nt,B,fc,f,M,d,h,w_imp,focus_loc,SNR_t,SNR_dB,Q,K);
    racc = racc + rr(K);
    % Perfect-CSI per-subcarrier analog MRT bound
    rp=0;
    for m=1:M
        wopt=exp(1j*angle(h(m,:)'))/sqrt(Nt);
        rp = rp + log2(1+SNR_t*abs(h(m,:)*wopt)^2)/M;
    end
    pacc = pacc + rp;
end
rate = racc/N_iter; rate_perf = pacc/N_iter;
end
