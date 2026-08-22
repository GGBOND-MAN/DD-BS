% ========================================================================
% TASK 1b: the TTD-count x PILOT-count DESIGN FRONTIER, and its scaling law
% ------------------------------------------------------------------------
% redesign_for_shared.m showed that slowing the angular sweep (s<1, hence more
% pilots K) buys back the rate that shared TTDs destroy. Reading that table
% column-by-column gives, for each TTD budget, the SMALLEST K that still reaches
% the ideal rate:
%
%       TTD elements   256   128    64    32
%       K_min            3     3     6    12          ->  K_min ~ 1.5 * P
%
% i.e. K * N_TTD ~ 1.5 * Nt  on the frontier. This script (a) fills that frontier
% in at a finer K resolution and extends it to P=16, and (b) TESTS the analytic
% prediction for where the constant 1.5 comes from.
%
% WHERE THE LAW COMES FROM (sharing-error bound):
%   antennas in a group of P share one delay tau_g, so the residual phase is
%       2*pi*(f_m - fc)*(tau_g - tau_n),  |f_m - fc| <= B/2,
%   and the intra-group delay spread is dominated by the LINEAR DDBS term,
%       d_tau ~ P*d*|theta_t|/c.
%   Requiring the worst-case residual phase to stay below a fixed phi_max:
%       pi*(B/2)*P*d*|theta_t|/c <= phi_max,   d = c/(2 fc)
%    => P*|theta_t| <= 4*fc*phi_max/(pi*B)                       ... (INVARIANT)
%   The redesign ties theta_t to K through the sweep rate, |theta_t| = |theta_t0|*K0/K,
%   so the invariant becomes
%       K / P >= (pi*B/(4*fc*phi_max)) * |theta_t0| * K0   ~  proportional to B/fc.
%
%   At the measured frontier (fc=30 GHz, B=5 GHz) K/P = 1.5, which pins
%   phi_max ~ 8.6 rad and predicts   K/P = 9 * (B/fc).
%   PART 2 tests exactly that: at P=8, K_min should be 6 / 12 / 24 for
%   B = 2.5 / 5 / 10 GHz. If it does, the constant is a device-independent
%   property of the DDBS beam and the law is publishable as a design rule.
%
% Runtime warning: actual_focus() is an M*K grid search, so large K is slow.
% ========================================================================
clear; clc;

ensure_baseline();
addpath(genpath(fileparts(mfilename('fullpath'))));

Nt=256; fc0=30e9; B0=5e9; M=1024; c=3e8; d=(c/fc0)/2; Q=1;
Rmin=5; Rmax=30; th_min=-pi/3; th_max=pi/3;
TH1_0=349/11; TH2_0=-30; AL1=-427/800; AL2=1859/3200; K0=3;
SNR_dB=15; SNR_t=10^(SNR_dB/10); N_iter=30; PASS=0.95;

% ---------------- PART 1: fill in the frontier at (fc=30, B=5) ----------------
rng(23); CH=cell(N_iter,1);
for it=1:N_iter
    CH{it}=near_field_channel(Nt,d,fc0,B0,M, Rmin+rand*(Rmax-Rmin), th_min+rand*(th_max-th_min));
end
r_ideal = eval_point(1, 3, Nt,B0,fc0,M,d,TH1_0,TH2_0,AL1,AL2,K0,CH,SNR_t,SNR_dB,Q);

fprintf('PART 1  frontier at fc=%.0f GHz, B=%.0f GHz (ideal rate %.3f, pass = %.0f%%)\n', ...
        fc0/1e9, B0/1e9, r_ideal, 100*PASS);
fprintf('%4s %8s | %6s %6s %8s %9s %10s\n','P','N_TTD','K_min','K/P','theta_t','range','P*|th_t|');
for P = [2 4 8 16]
    ladder = unique(max(3, ceil([0.375 0.75 1.125 1.5 2 3]*P)));
    Kmin = NaN;
    for K = ladder
        r = eval_point(P, K, Nt,B0,fc0,M,d,TH1_0,TH2_0,AL1,AL2,K0,CH,SNR_t,SNR_dB,Q);
        fprintf('   (P=%2d K=%2d -> %.3f  %3.0f%%)\n', P, K, r, 100*r/r_ideal);
        if r >= PASS*r_ideal, Kmin = K; break; end
    end
    if isnan(Kmin), fprintf('%4d %8d | %6s\n', P, Nt/P, 'none'); continue; end
    TH1 = TH1_0*K0/Kmin; nn=(-(Nt-1)/2:(Nt-1)/2)';
    rg = (max(nn*d*TH1-(nn*d).^2*AL1)-min(nn*d*TH1-(nn*d).^2*AL1))/c;
    fprintf('%4d %8d | %6d %6.2f %8.2f %7.1fns %10.1f\n', P, Nt/P, Kmin, Kmin/P, TH1, rg*1e9, P*abs(TH1));
end

% ---------------- PART 2: does the constant scale as B/fc? ----------------
fprintf('\nPART 2  P=8 (32 TTDs), sweep B. Prediction: K_min/P = 9*(B/fc)\n');
fprintf('%8s | %8s %8s %8s\n','B [GHz]','K_min','K/P','pred K/P');
for B = [2.5e9 5e9 10e9]
    rng(23); CHb=cell(N_iter,1);
    for it=1:N_iter
        CHb{it}=near_field_channel(Nt,d,fc0,B,M, Rmin+rand*(Rmax-Rmin), th_min+rand*(th_max-th_min));
    end
    r_id = eval_point(1, 3, Nt,B,fc0,M,d,TH1_0,TH2_0,AL1,AL2,K0,CHb,SNR_t,SNR_dB,Q);
    Kmin = NaN;
    for K = [3 6 9 12 18 24 32]
        r = eval_point(8, K, Nt,B,fc0,M,d,TH1_0,TH2_0,AL1,AL2,K0,CHb,SNR_t,SNR_dB,Q);
        fprintf('   (B=%.1f GHz K=%2d -> %.3f  %3.0f%% of %.3f)\n', B/1e9, K, r, 100*r/r_id, r_id);
        if r >= PASS*r_id, Kmin = K; break; end
    end
    fprintf('%8.1f | %8s %8.2f %8.2f\n', B/1e9, num2str(Kmin), Kmin/8, 9*(B/fc0));
end

% ------------------------------------------------------------------------
function r = eval_point(P, K, Nt,B,fc,M,d,TH1_0,TH2_0,AL1,AL2,K0,CH,SNR_t,SNR_dB,Q)
% Rate at the design point (sweep scaled so that K pilots are used) under P-way
% TTD sharing, decoded with the RECALIBRATED lookup (actual_focus).
s = K0/K;  TH1 = TH1_0*s;  TH2 = TH2_0*s;
f = fc + B/M*((1:M)-1-(M-1)/2);
focus_loc = cal_loc(fc,f,TH1,TH2,AL1,AL2,M,K);
arch='full'; if P>1, arch='shared'; end
w  = ddbs_beam_arch(Nt,B,fc,M,d,TH1,TH2,AL1,AL2,K,arch,P,Inf,Inf,[]);
fl = actual_focus(w,focus_loc,Nt,fc,B,M,d,K);
r  = mean(cellfun(@(h) rate_ongrid(h,w,fl,Nt,B,fc,M,d,SNR_t,SNR_dB,Q,K), CH));
end

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
