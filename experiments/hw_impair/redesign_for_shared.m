% ========================================================================
% TASK 1: redesign the DD-BS parameters FOR a shared-TTD architecture
% ------------------------------------------------------------------------
% The baseline picks its parameters assuming one TTD per antenna. Under shared
% TTDs that choice is no longer optimal, because the sharing error is
%       2*pi*(f_m - fc) * (tau_g - tau_n),
% whose intra-group delay spread scales as P*d*|theta_t|.
%
% MECHANISM (worth stating carefully -- it is the opposite of the naive reading):
%   theta_m = theta_t + (fc/f_m)*(theta_p + 2 p_m)
%   * the TTD term steers FREQUENCY-INDEPENDENTLY -> theta_t is a constant offset
%   * the PS term squints with frequency          -> theta_p drives the SWEEP
% So the beam split comes from the PHASE SHIFTER, and |theta_t| is large only to
% cancel the large offset that the fast sweep forces. Slowing the sweep by a
% factor s therefore shrinks |theta_t| (better under sharing) but lays down fewer
% strips per pilot, so more pilots are needed: K ~ K0/s.
%
% Verified scaling from the baseline's own Algorithm 1:
%   gamma   theta_p+2pM   theta_t   delay range   K (eq 39)
%   0.95       36.78      -32.95      140.6 ns        3
%   0.50       19.36      -16.87       72.0 ns        4->6
%   0.25        9.68       -7.94       33.9 ns        8->12
%
% This probe sweeps the design point s against the TTD budget and reports the
% rate achieved with the RECALIBRATED lookup (the free algorithmic fix), so the
% question it answers is:
%     for a given TTD budget, which sweep rate / pilot count is best?
% ========================================================================
clear; clc;

ensure_baseline();          % locate the baseline code wherever it lives
addpath(genpath(fileparts(mfilename('fullpath'))));

Nt=256; fc=30e9; B=5e9; M=1024; c=3e8; d=(c/fc)/2; Q=1;
f = fc + B/M*((1:M)-1-(M-1)/2);
Rmin=5; Rmax=30; th_min=-pi/3; th_max=pi/3;
% baseline design point (Rate_snr.m): s = 1 reproduces it exactly
TH1_0=349/11; TH2_0=-30; AL1=-427/800; AL2=1859/3200; K0=3;

SNR_dB=15; SNR_t=10^(SNR_dB/10); N_iter=60; rng(23);
CH=cell(N_iter,1);
for it=1:N_iter
    CH{it}=near_field_channel(Nt,d,fc,B,M, Rmin+rand*(Rmax-Rmin), th_min+rand*(th_max-th_min));
end

s_list = [1 0.5 0.25];            % sweep-rate scaling (1 = baseline design)
P_list = [1 2 4 8];               % 1 = full per-antenna TTD

fprintf('Rate with the RECALIBRATED lookup. Angular sweep scaled by s; the\n');
fprintf('distance parameters are left alone (their coverage requirement is\n');
fprintf('independent of s). K is the resulting pilot overhead.\n\n');
fprintf('%6s %5s %8s %9s |', 's', 'K', 'theta_t', 'range');
fprintf('%12s', string(compose("P=%d(%dTTD)", P_list(:), Nt./P_list(:)))); fprintf('\n');

for s = s_list
    TH1 = TH1_0*s; TH2 = TH2_0*s;      % scale the angular pair together
    K   = ceil(K0/s);                  % eq (39): K ~ 1/(theta_p + 2 p_1)
    nn=(-(Nt-1)/2:(Nt-1)/2)';
    rng_tau = (max(nn*d*TH1-(nn*d).^2*AL1)-min(nn*d*TH1-(nn*d).^2*AL1))/c;
    focus_loc = cal_loc(fc,f,TH1,TH2,AL1,AL2,M,K);
    fprintf('%6.2f %5d %8.2f %8.1fns |', s, K, TH1, rng_tau*1e9);
    for P = P_list
        arch = 'full'; if P>1, arch='shared'; end
        w  = ddbs_beam_arch(Nt,B,fc,M,d,TH1,TH2,AL1,AL2,K,arch,P,Inf,Inf,[]);
        fl = actual_focus(w,focus_loc,Nt,fc,B,M,d,K);
        r  = mean(cellfun(@(h) rate_ongrid(h,w,fl,Nt,B,fc,M,d,SNR_t,SNR_dB,Q,K), CH));
        fprintf('%12.3f', r);
    end
    fprintf('\n');
end
fprintf(['\nRead each column (a fixed TTD budget) downward: does giving up pilot\n' ...
         'efficiency (larger K) buy back enough sharing robustness to be worth it?\n' ...
         'Compare rates at EQUAL K across columns to judge the hardware saving.\n']);

% ------------------------------------------------------------------------
function r = rate_ongrid(h,w,focus_loc,Nt,B,fc,M,d,SNR_t,SNR_dB,Q,K)
best=-inf;
for idx=1:K
    y=zeros(M,Q);
    for m=1:M, y(m,:)=awgn(repmat(h(m,:)*w(:,idx,m),[1,Q]),SNR_dB*2/sqrt(3)); end
    [~,i]=max(abs(sum(y,2)).^2);
    ws = TTD_beam(Nt,B,fc,M,d,focus_loc(i,1,idx),focus_loc(i,2,idx));
    t=0; for m=1:M, t=t+log2(1+SNR_t*abs(h(m,:)*ws(:,m))^2)/M; end
    best=max(best,t);
end
r=best;
end
