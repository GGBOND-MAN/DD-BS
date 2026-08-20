% ========================================================================
% (B) END-TO-END rate: naive vs hybrid, BOTH stages impaired  [money figure]
% ------------------------------------------------------------------------
% Fixes the modeling gap in hw_probe_rate.m (which served with an IDEAL beam).
% Here training beams AND the serving beam share one hardware delay LSB, so the
% comparison reflects a real low-resolution-TD implementation.
%
% Sweeps the hardware TD resolution and reports average rate for:
%   - naive  DDBS (all delay in the TTD)   = baseline architecture
%   - hybrid DDBS (PS absorbs the fc term) = proposed
% against ideal DDBS and Perfect CSI. This is the figure the paper is built on.
% ========================================================================
clear; clc;

BASELINE_DIR = fullfile('..','..','baseline_code','code_nf_distance_dependent_rainbow');
addpath(genpath(BASELINE_DIR));
addpath(genpath(fileparts(mfilename('fullpath'))));

Nt=256; fc=30e9; B=5e9; M=1024; c=3e8; d=(c/fc)/2; Q=1; K=3;
f = fc + B/M*((1:M)-1-(M-1)/2);
Rmin=5; Rmax=30; th_min=-pi/3; th_max=pi/3;
theta1=349/11; theta2=-30; alpha1=-427/800; alpha2=1859/3200;
focus_loc = cal_loc(fc, f, theta1, theta2, alpha1, alpha2, M, K);
nn=(-(Nt-1)/2:(Nt-1)/2)';
tau_rng = max((nn*d*theta1-(nn*d).^2*alpha1)/c) - min((nn*d*theta1-(nn*d).^2*alpha1)/c);

SNR_dB=15; SNR_t=10^(SNR_dB/10); N_iter=100; rng(5);

% pre-draw channels so all schemes see identical realisations
CH = cell(N_iter,1);
for it=1:N_iter
    r = Rmin+rand*(Rmax-Rmin); th = th_min+rand*(th_max-th_min);
    CH{it} = near_field_channel(Nt,d,fc,B,M,r,th);
end

runavg = @(wfun,lsb) mean(cellfun(@(h) ...
    max(training_ddbs_e2e(Nt,B,fc,f,M,d,h,wfun(),focus_loc,SNR_t,SNR_dB,Q,K,lsb,0)), CH));

wl_ideal = @() delay_polar_2d(Nt,B,fc,M,d,theta1,theta2,alpha1,alpha2,K);
r_ideal = runavg(wl_ideal, Inf);
rp = mean(cellfun(@(h) sum(arrayfun(@(m) ...
        log2(1+SNR_t*abs(h(m,:)*(exp(1j*angle(h(m,:)'))/sqrt(Nt)))^2)/M, 1:M)), CH));
fprintf('ideal DDBS %.3f | Perfect CSI %.3f  (SNR %d dB, K=%d)\n\n', r_ideal, rp, SNR_dB, K);

fprintf('%5s %9s | %9s %9s\n','B_td','LSB[ps]','naive','HYBRID');
for Btd = [16 14 13 12 11 10 9 8]
    lsb = tau_rng/2^Btd; imp = struct('Btd',Btd);
    rn = runavg(@() ddbs_beam_impaired(Nt,B,fc,M,d,theta1,theta2,alpha1,alpha2,K,imp), lsb);
    rh = runavg(@() ddbs_beam_hybrid  (Nt,B,fc,M,d,theta1,theta2,alpha1,alpha2,K,imp), lsb);
    fprintf('%5d %9.1f | %9.3f %9.3f\n', Btd, lsb*1e12, rn, rh);
end
