% ========================================================================
% TASK 2: break the K-vs-TTD tradeoff by DECOUPLING theta_t from theta_p
% ------------------------------------------------------------------------
% frontier_ttd_pilot.m measured  K_min = ceil(1.12*P), i.e. the sharing limit is
%       P * |theta_t| <= ~85      -- and it is INDEPENDENT OF B.
%
% The B-independence rules out the worst-case-subcarrier phase bound (which is
% proportional to B) and points at COVERAGE TRUNCATION instead:
%   * intra-group residual = 2*pi*(f_m-fc)*(tau_g-tau_n) is LINEAR in the
%     intra-group index -> it TILTS each sub-array pattern by (f_m-fc)*theta_t/f_m;
%   * the tilt leaves the sub-array beamwidth (~2/P) once
%         |f_m - fc| > eta * c/(P*d*|theta_t|) = eta * 2*fc/(P*|theta_t|),
%     so only that slice of the band still produces a usable beam;
%   * DDBS sweeps at d(theta_m)/df ~ -theta_p/fc, so the USABLE angular coverage
%     per pilot is
%         2 * (theta_p/fc) * eta*2*fc/(P*|theta_t|) = 4*eta*theta_p/(P*|theta_t|)
%     -- fc cancels and, crucially, SO DOES B. That is exactly what PART 2 saw.
%
% Covering Theta_req with K pilots then needs
%       K >= Theta_req * P * |theta_t| / (4*eta*|theta_p|)                (*)
% and with the 3 dB usable width (eta = 0.443) and the baseline's own numbers
% (Theta_req = 2 sin 60 = 1.732, |theta_t0| = 31.73, |theta_p,eff| ~ 30..37)
% this gives K >= (0.84 .. 1.04) * P against the measured 1.12 * P.
%
% THE CONSEQUENCE THIS SCRIPT TESTS. In (*) theta_t and theta_p enter as a RATIO.
% redesign_for_shared scaled them TOGETHER (s), so it paid for sharing robustness
% in pilots. But the two have different jobs:
%       theta_p  = sweep rate       -> wants to be LARGE (few pilots)
%       theta_t  = DC offset that re-centres the sweep, ~ -theta_p
%                + per-pilot increments 2(s-1)/K, which span only ~2
% Only the DC role of theta_t is large, and only the DC role causes the sharing
% failure. If the sweep can be re-centred WITHOUT a large TTD ramp, then (*) is
% satisfied at K = 3 and P = 8 simultaneously -- the tradeoff disappears.
%
% This scans (theta_t, theta_p) INDEPENDENTLY and reports rate AND the realized
% angular coverage of the recalibrated focus map, so a high rate that is really
% just "the beams moved off the user region" cannot be mistaken for a win.
% ========================================================================
clear; clc;

ensure_baseline();
addpath(genpath(fileparts(mfilename('fullpath'))));

Nt=256; fc=30e9; B=5e9; M=1024; c=3e8; d=(c/fc)/2; Q=1;
f = fc + B/M*((1:M)-1-(M-1)/2);
Rmin=5; Rmax=30; th_min=-pi/3; th_max=pi/3;
TH1_0=349/11; TH2_0=-30; AL1=-427/800; AL2=1859/3200;
SNR_dB=15; SNR_t=10^(SNR_dB/10); N_iter=30;

rng(23); CH=cell(N_iter,1);
for it=1:N_iter
    CH{it}=near_field_channel(Nt,d,fc,B,M, Rmin+rand*(Rmax-Rmin), th_min+rand*(th_max-th_min));
end

P = 8;  K = 3;                       % 32 TTDs at the BASELINE pilot budget
% reference points
r_ideal = runpt(TH1_0,TH2_0,1,3, Nt,B,fc,M,d,AL1,AL2,f,CH,SNR_t,SNR_dB,Q);
fprintf('reference: full TTD, baseline params, K=3  -> %.3f\n', r_ideal(1));
r_base8 = runpt(TH1_0,TH2_0,P,K, Nt,B,fc,M,d,AL1,AL2,f,CH,SNR_t,SNR_dB,Q);
fprintf('           P=8 (32 TTD), baseline params, K=3 -> %.3f (%.0f%%)\n\n', ...
        r_base8(1), 100*r_base8(1)/r_ideal(1));

fprintf('Independent scan at P=%d, K=%d.  "cov" = span of the recalibrated focus\n', P, K);
fprintf('angles as a fraction of the user region sin(theta) in [-0.87, 0.87].\n');
fprintf('A win needs HIGH rate AND cov ~ 1.\n\n');
fprintf('%7s %7s %9s | %8s %6s %8s %9s\n','s_t','s_p','P*|th_t|','rate','%ideal','cov','maxgap/BW');
for s_t = [1 0.5 0.25 0.125 0.0625]
    for s_p = [1 0.5]
        TH1 = TH1_0*s_t;  TH2 = TH2_0*s_p;
        out = runpt(TH1,TH2,P,K, Nt,B,fc,M,d,AL1,AL2,f,CH,SNR_t,SNR_dB,Q);
        fprintf('%7.4f %7.2f %9.1f | %8.3f %5.0f%% %8.2f %9.1f\n', ...
                s_t, s_p, P*abs(TH1), out(1), 100*out(1)/r_ideal(1), out(2), out(3));
    end
end
fprintf(['\nIf a small-theta_t / large-theta_p cell reaches ~100%% of ideal WITH\n' ...
         'cov ~ 1 at K=3, the K-vs-TTD tradeoff is broken and the remaining job is\n' ...
         'an ARCHITECTURE that re-centres the sweep without a large TTD ramp.\n' ...
         'If instead every high-rate cell has low cov, the offset is load-bearing\n' ...
         'and the tradeoff K_min = 1.12 P is fundamental to this beam family.\n']);

% ------------------------------------------------------------------------
function out = runpt(TH1,TH2,P,K, Nt,B,fc,M,d,AL1,AL2,f,CH,SNR_t,SNR_dB,Q)
focus_loc = cal_loc(fc,f,TH1,TH2,AL1,AL2,M,K);
arch='full'; if P>1, arch='shared'; end
w  = ddbs_beam_arch(Nt,B,fc,M,d,TH1,TH2,AL1,AL2,K,arch,P,Inf,Inf,[]);
fl = actual_focus(w,focus_loc,Nt,fc,B,M,d,K);
r  = mean(cellfun(@(h) rate_ongrid(h,w,fl,Nt,B,fc,M,d,SNR_t,SNR_dB,Q,K), CH));
[gp, cov] = cov_gap(fl, Nt);      % boundary-aware -- see cov_gap.m
out = [r, cov, gp];
end

function r = rate_ongrid(h,w,focus_loc,Nt,B,fc,M,d,SNR_t,SNR_dB,Q,K)
best=-inf;
for idx=1:K
    y=zeros(M,Q);
    for m=1:M, y(m,:)=repmat(h(m,:)*w(:,idx,m),[1,Q]); end
    y = add_awgn(y, SNR_dB*2/sqrt(3));
    [~,i]=max(abs(sum(y,2)).^2);
    ws = TTD_beam(Nt,B,fc,M,d,focus_loc(i,1,idx),focus_loc(i,2,idx));
    t=0; for m=1:M, t=t+log2(1+SNR_t*abs(h(m,:)*ws(:,m))^2)/M; end
    best=max(best,t);
end
r=best;
end
