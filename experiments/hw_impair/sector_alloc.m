% ========================================================================
% SECOND CONTRIBUTION: under compensation, THEIR pilot allocation is wrong
% ------------------------------------------------------------------------
% ablation_ps.m left a puzzle worth reading carefully. The compensated arm hits
% 99% of the ungrouped bound at L=4 but only 69% at L=8 and 34% at L=16 -- yet
% sec 21.2 measured K_min = ceil(P/1.19), which at L=8 is 7, comfortably inside
% the K=12 budget those runs were given. More pilots than the law requires, and
% still short. Why?
%
% Because the law is about the ANGULAR comb, and their allocation spends the
% budget elsewhere. With N_sec sector shifts over theta_tot = 2 the comb spacing
% is Delta = 2/N_sec, so the master-curve variable is
%
%       Delta * P / 2  =  (2/N_sec) * L / 2  =  L / N_sec
%
% and their fixed N_sec = 4 gives
%
%       L        1     2     4     8    16
%       L/N_sec  1.0   1.0   1.0   2.0   4.0
%       predicted 100%  100%  100%  ~55%  <<55%      (sec 20.2 master curve)
%       measured  100%   85%   99%   69%   34%
%
% -- the curve predicts the shape, including exactly where it breaks. Their (54)
% computes N_sec from the sector width of the UNCOMPENSATED architecture; under
% compensation the sector is wider, so the optimal split between sector shifts
% and distance interleaving moves, and their (71) allocation over-spends on
% distance.
%
% PREDICTION, registered before the run: reallocating the SAME or FEWER pilots
% from the distance dimension to the angular dimension should beat their split.
%       L=8  : N_sec=7,  K_alpha=1  ->  7 pilots, predicted ~98%  (vs 12 -> 69%)
%       L=16 : N_sec=14, K_alpha=1  -> 14 pilots, predicted ~98%  (vs 12 -> 34%)
% If it holds, the contribution is two-part and the second half costs NEGATIVE
% overhead: fewer pilots AND higher rate, purely by allocating them correctly.
% ========================================================================
clear; clc;

ensure_baseline();
addpath(genpath(fileparts(mfilename('fullpath'))));

Nt=256; fc=30e9; B=5e9; M=1024; c=3e8; d=(c/fc)/2; Q=1;
f = fc + B/M*((1:M)-1-(M-1)/2);
thlim=[-1 1]; alim=[0.0025 0.1]; gamma=0.9;
Rmin=5; Rmax=200; SNR_dB=20; SNR_t=10^(SNR_dB/10); N_iter=200;                     % final-figure sample size; sec 33.3

rng(23); CH=cell(N_iter,1);
for it=1:N_iter
    CH{it} = near_field_channel(Nt,d,fc,B,M, Rmin+rand*(Rmax-Rmin), asin(-0.9+1.8*rand));
end

% ungrouped reference at the same pilot budget
[sec0, rho0] = ojcoms_algorithm1(Nt,fc,B,M,d,1,thlim,alim,gamma,1,12);
r_ref = armrate(sec0, rho0, 'full', 1, Nt,B,fc,M,d,f,CH,SNR_t,SNR_dB,Q,'ideal');
fprintf('ungrouped reference (L=1, 12 pilots) = %.3f\n', r_ref);
fprintf('compensated PS throughout; only the (N_sec, K_alpha) split changes.\n\n');
fprintf('serving beam: "ideal" = full per-antenna TTD (old); "shared" = same Nt/L delays\n');
fprintf('rates as mean +/- 95%% CI over %d channels\n', N_iter);
fprintf('%3s %7s | %5s %7s %6s | %6s | %13s %5s | %14s %5s | %s\n', ...
        'L','N_TTD','Nsec','Kalpha','K_tot','L/Nsec','ideal serving','%ref','shared serving','%ref','note');

cells = { 8,[4 3],'their split'; 8,[7 1],'ours: 7 pilots'; 8,[8 1],''; 8,[12 1],'equal K, all angular';
         16,[4 3],'their split'; 16,[12 1],'equal K, all angular'; 16,[14 1],'ours: 14 pilots';
          4,[4 3],'their split'; 4,[4 1],'ours: 4 pilots' };
for i = 1:size(cells,1)
    L = cells{i,1}; Nsec = cells{i,2}(1); Ka = cells{i,2}(2);
    [sec, rho] = ojcoms_algorithm1(Nt,fc,B,M,d,L,thlim,alim,gamma,Nsec,Ka);
    [ri,ci] = armrate(sec, rho, 'shared', L, Nt,B,fc,M,d,f,CH,SNR_t,SNR_dB,Q,'ideal');
    [rs,cs] = armrate(sec, rho, 'shared', L, Nt,B,fc,M,d,f,CH,SNR_t,SNR_dB,Q,'shared');
    fprintf('%3d %7d | %5d %7d %6d | %6.2f | %6.3f+-%-5.3f %4.0f%% | %6.3f+-%-5.3f %4.0f%% | %s\n', ...
            L, Nt/L, Nsec, Ka, numel(sec), L/Nsec, ri,ci,100*ri/r_ref, rs,cs,100*rs/r_ref, cells{i,3});
end
fprintf(['\nRead: compare rows at equal L. If the "ours" row beats the "their\n' ...
         'split" row with FEWER pilots, the allocation rule is a contribution in\n' ...
         'its own right and it costs negative overhead. L/Nsec is the master-curve\n' ...
         'variable -- rows with L/Nsec near 1 should land near 100%% regardless of L.\n']);

% ------------------------------------------------------------------------
function [r, gp] = armrate(sec, rho, arch, L, Nt,B,fc,M,d,f,CH,SNR_t,SNR_dB,Q,SERVE)
K = numel(sec); c=3e8; kc=2*pi*fc/c; km=2*pi*f/c;
w = zeros(Nt,K,M); cf = zeros(M,2,K);
for s = 1:K
    w(:,s,:) = ddbs_beam_arch(Nt,B,fc,M,d, sec(s).theta_t, sec(s).theta_p, ...
                              sec(s).alpha_t, sec(s).alpha_p, 1, arch, L, Inf, Inf, []);
    cf(:,1,s) = (kc./km).'*sec(s).theta_p + rho.theta*sec(s).theta_t;
    cf(:,2,s) = (kc./km).'*sec(s).alpha_p + rho.alpha*sec(s).alpha_t;
end
fl = actual_focus(w, cf, Nt, fc, B, M, d, K);
v  = cellfun(@(h) rate_ongrid(h,w,fl,Nt,B,fc,M,d,SNR_t,SNR_dB,Q,K,SERVE,L), CH);
r  = mean(v);  gp = 1.96*std(v)/sqrt(numel(v));   % 2nd output is now the 95% CI
end

function r = rate_ongrid(h,w,focus_loc,Nt,B,fc,M,d,SNR_t,SNR_dB,Q,K,SERVE,P)
best=-inf;
for idx=1:K
    y=zeros(M,Q);
    for m=1:M, y(m,:)=repmat(h(m,:)*w(:,idx,m),[1,Q]); end
    y = add_awgn(y, SNR_dB*2/sqrt(3));
    [~,i]=max(abs(sum(y,2)).^2);
    ws = serve_beam(Nt,B,fc,M,d,focus_loc(i,1,idx),focus_loc(i,2,idx),SERVE,P);
    t=0; for m=1:M, t=t+log2(1+SNR_t*abs(h(m,:)*ws(:,m))^2)/M; end
    best=max(best,t);
end
r=best;
end
