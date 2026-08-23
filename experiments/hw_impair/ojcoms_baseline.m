% ========================================================================
% THEIR SCHEME AS A REAL BASELINE -- staged, with a validation gate first
% ------------------------------------------------------------------------
% THEORY.md sec 22 established that OJ-COMS v7 2026 already contains most of
% sec 12-21. Their own Fig. 7 also shows where they stop: with the pilot budget
% their formulas require, the rate collapses as the sub-array grows --
%
%      L        1     2     3     5     7     9
%      K(L)     3    12    21    33    45    57      (their (73))
%      rate   ~6.5  ~5.0  ~3.7  ~1.9  ~1.3  ~0.5     (their Fig. 7)
%
% -- and they conclude "L = 2 provides a favorable tradeoff". Sharing beyond 2x
% is, on their evidence, impractical. The compensated arm measured 97-99% of the
% ungrouped bound at L = 4, 8 and 16 (sec 21.4). If that survives a like-for-like
% comparison at THEIR operating point, the contribution is removing their L <= 2
% ceiling -- large, specific, and citing them for the framework it builds on.
%
% "If" is doing real work there: sec 21.4 was measured at MY operating point
% (r in [5,30] m, baseline DDBS parameters, uniform comb) and theirs is different
% (r in [5,200] m, alpha in [0.0025,0.1], sector-shifted, SNR 20 dB). Nothing
% transfers without being re-measured.
%
% STAGE 1 (this script by default) is a VALIDATION GATE, not a result. It runs
% only their scheme, at their published operating point, and checks two numbers
% read off their Fig. 7:
%      L=1, K=3   ->  ~6.5 bit/s/Hz
%      L=2, K=12  ->  ~5.0 bit/s/Hz
% Their reported sector-1 parameters (theta't = -31.797, theta'p = 28.4) are
% checked too. If the gate fails, the reconstruction is wrong and NOTHING may be
% compared -- the fix is to obtain their Table 3 rather than to tune until it
% matches. Set RUN_STAGE2 = true only after the gate passes.
%
% Three arms, so that they get the benefit of every doubt:
%   O-cf  their beams + their closed-form LS lookup (29)-(30)
%   O-af  their beams + the grid-search focus table (>= O-cf by construction)
%   C-af  fc-compensated beams + the same grid-search table   <- ours
% ========================================================================
clear; clc;
RUN_STAGE2 = false;          % flip only after the gate passes

ensure_baseline();
addpath(genpath(fileparts(mfilename('fullpath'))));

Nt=256; fc=30e9; B=5e9; M=1024; c=3e8; d=(c/fc)/2; Q=1;
f = fc + B/M*((1:M)-1-(M-1)/2);
thlim=[-1 1]; alim=[0.0025 0.1]; gamma=0.9;      % their Sec VII-A settings
Rmin=5; Rmax=200;                                 % their 5-200 m
SNR_dB=20; SNR_t=10^(SNR_dB/10); N_iter=20;

rng(23); CH=cell(N_iter,1);
for it=1:N_iter
    r = Rmin + rand*(Rmax-Rmin);  th = asin(-0.9 + 1.8*rand);
    CH{it} = near_field_channel(Nt,d,fc,B,M,r,th);
end

% ---------------- parameter reconstruction + unit test ----------------
[~, rho2, dbg2] = ojcoms_algorithm1(Nt,fc,B,M,d,2,thlim,alim,gamma,4,3);
fprintf('=== Algorithm 1 reconstruction, L=2 ===\n');
fprintf('rho_theta = %.4f, rho_alpha = %.4f   (both -> 1 as L -> 1)\n', rho2.theta, rho2.alpha);
fprintf('S bound (59)-(60) = %.3f | S boundary-focusing = %.3f | p_M = %d\n', ...
        dbg2.Sbound, dbg2.Sbf, dbg2.pM_bf);
fprintf('sector 1: theta_t = %8.3f   theta_p = %8.3f\n', dbg2.theta_t1, dbg2.theta_p1);
fprintf('  PAPER  : theta_t =  -31.797   theta_p =   28.400   <-- unit test\n');
fprintf('Nsec by (54) = %d (paper uses 4).  alpha_p = %.4f, alpha_t range = [%.4f %.4f]\n\n', ...
        dbg2.Nsec_formula, dbg2.alpha_p, dbg2.at_range(1), dbg2.at_range(2));

% ---------------- STAGE 1: the gate ----------------
fprintf('=== STAGE 1 gate (their scheme only, SNR=%d dB, N_iter=%d) ===\n', SNR_dB, N_iter);
fprintf('%3s %5s | %8s %8s | %s\n','L','K','O-cf','O-af','their Fig. 7');
tgt = [1 3 6.5; 2 12 5.0];
for i = 1:size(tgt,1)
    L = tgt(i,1); Ktot = tgt(i,2);
    if L==1, Nsec=1; Ka=Ktot; else, Nsec=4; Ka=Ktot/Nsec; end
    [sec,rho] = ojcoms_algorithm1(Nt,fc,B,M,d,L,thlim,alim,gamma,Nsec,Ka);
    [rcf, raf] = runscheme(sec, rho, 'ojcoms', L, Nt,B,fc,M,d,f,CH,SNR_t,SNR_dB,Q);
    fprintf('%3d %5d | %8.3f %8.3f | ~%.1f\n', L, numel(sec), rcf, raf, tgt(i,3));
end
fprintf(['\nGATE: O-cf (or O-af) must land near their Fig. 7 column. If it does\n' ...
         'not, STOP -- get their Table 3 (the sector-wise parameters); do not\n' ...
         'tune parameters until the number matches.\n']);

if ~RUN_STAGE2, return; end

% ---------------- STAGE 2: the like-for-like comparison ----------------
fprintf('\n=== STAGE 2 (only valid if the gate passed) ===\n');
fprintf('%3s %7s | %6s %8s | %6s %8s | %s\n', ...
        'L','N_TTD','K_them','O-af','K_ours','C-af','their Fig. 7 rate');
fig7 = containers.Map({2,4,8},{5.0,2.6,0.9});     % read off their Fig. 7
for L = [2 4 8]
    Kthem = round(6.4*L);                          % their (73), ~6.4L at Nt=256
    Nsec = max(1,round(Kthem/3));
    [secT,rhoT] = ojcoms_algorithm1(Nt,fc,B,M,d,L,thlim,alim,gamma,Nsec,3);
    [~, rafT]   = runscheme(secT, rhoT, 'ojcoms', L, Nt,B,fc,M,d,f,CH,SNR_t,SNR_dB,Q);

    Kours = max(1,ceil(L/1.19));                   % sec 21.2
    [secO,rhoO] = ojcoms_algorithm1(Nt,fc,B,M,d,L,thlim,alim,gamma,Kours,1);
    [~, rafO]   = runscheme(secO, rhoO, 'shared', L, Nt,B,fc,M,d,f,CH,SNR_t,SNR_dB,Q);

    fprintf('%3d %7d | %6d %8.3f | %6d %8.3f | ~%.1f\n', ...
            L, Nt/L, numel(secT), rafT, numel(secO), rafO, fig7(L));
end

% ------------------------------------------------------------------------
function [r_cf, r_af] = runscheme(sec, rho, arch, L, Nt,B,fc,M,d,f,CH,SNR_t,SNR_dB,Q)
K = numel(sec); c=3e8; kc=2*pi*fc/c; km=2*pi*f/c;
w = zeros(Nt,K,M); cf = zeros(M,2,K);
for s = 1:K
    w(:,s,:) = ddbs_beam_arch(Nt,B,fc,M,d, sec(s).theta_t, sec(s).theta_p, ...
                              sec(s).alpha_t, sec(s).alpha_p, 1, arch, L, Inf, Inf, []);
    % (29)-(30): the closed-form LS effective focus
    cf(:,1,s) = (kc./km).'*sec(s).theta_p + rho.theta*sec(s).theta_t;
    cf(:,2,s) = (kc./km).'*sec(s).alpha_p + rho.alpha*sec(s).alpha_t;
end
af = actual_focus(w, cf, Nt, fc, B, M, d, K);
r_cf = mean(cellfun(@(h) rate_ongrid(h,w,cf,Nt,B,fc,M,d,SNR_t,SNR_dB,Q,K), CH));
r_af = mean(cellfun(@(h) rate_ongrid(h,w,af,Nt,B,fc,M,d,SNR_t,SNR_dB,Q,K), CH));
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
