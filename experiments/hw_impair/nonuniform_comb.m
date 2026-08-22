% ========================================================================
% "PLACE THE PILOTS BETTER, DON'T USE MORE OF THEM" -- greedy offset selection
% ------------------------------------------------------------------------
% pilot_spacing_map.m confirmed the comb model 20/20, on the corner named in
% advance:
%       K =  8, Delta = 0.25 -> 100% of ideal
%       K = 16, Delta = 0.50 ->  58% of ideal
% Twice the pilots fails; half the pilots at the right spacing reaches ideal. So
% the variable is the pilot comb's RESOLUTION relative to the sub-array beamwidth
% 2/P, not the pilot count. With the baseline's hard-wired Delta = 2/K this forces
%       Delta = 2/K <= 2/P   =>   K >= P   =>   K * N_TTD = Nt
%
% A UNIFORM comb cannot beat K = P: covering a span of 2 at a resolution of 2/P
% needs P points wherever the resolution is required. But map_focus_shift showed
% the holes are NOT uniformly distributed -- they concentrate near boresight:
%       P=4, K=3 : biggest hole [-0.020, +0.137]
%       P=8, K=3 : biggest hole [-0.140, +0.252]
% Where sharing opens no hole, the fine resolution is wasted. So a NON-UNIFORM
% comb -- dense where holes appear, sparse elsewhere -- should reach full coverage
% with K < P, which a uniform comb provably cannot.
%
% The selection needs no channel knowledge and no extra measurement: each
% candidate offset's realized focus set is computed OFFLINE by exactly the same
% actual_focus table that the recalibrated lookup already requires. So this rides
% on the existing contribution at zero run-time cost -- greedy set cover over
% precomputed coverage sets.
%
% Falsifiable: if greedy K < P reaches ~100% while the uniform comb at the same K
% does not, pilot PLACEMENT is a second algorithm-side lever and the K >= P law
% applies only to uniform combs. If greedy never beats uniform, K >= P is
% fundamental to the beam family and the exchange law K * N_TTD = Nt stands
% unconditionally -- which is itself the cleaner paper.
%
% Runtime: the offline stage is 32 single-pilot focus tables (~ one K=32 run).
% ========================================================================
clear; clc;

ensure_baseline();
addpath(genpath(fileparts(mfilename('fullpath'))));

Nt=256; fc=30e9; B=5e9; M=1024; c=3e8; d=(c/fc)/2; Q=1;
f = fc + B/M*((1:M)-1-(M-1)/2);
Rmin=5; Rmax=30; th_min=-pi/3; th_max=pi/3;
TH1_0=349/11; TH2_0=-30; AL1=-427/800; AL2=1859/3200;
SNR_dB=15; SNR_t=10^(SNR_dB/10); N_iter=30; P=8; LIM=0.866;

rng(23); CH=cell(N_iter,1);
for it=1:N_iter
    CH{it}=near_field_channel(Nt,d,fc,B,M, Rmin+rand*(Rmax-Rmin), th_min+rand*(th_max-th_min));
end

% ---------- offline stage: coverage set of every candidate offset ----------
NC = 32; offs = (0:NC-1)*(2/NC);
edges = -LIM : (2/Nt) : LIM;  NB = numel(edges)-1;
COVSET = false(NC, NB);  FL = cell(NC,1);
fprintf('offline: %d candidate offsets, %d coverage bins of width 2/Nt\n', NC, NB);
for ci = 1:NC
    th1s = TH1_0 - offs(ci);
    wc = ddbs_beam_arch(Nt,B,fc,M,d,th1s,TH2_0,AL1,AL2,1,'shared',P,Inf,Inf,[]);
    flc = actual_focus(wc, cal_loc(fc,f,th1s,TH2_0,AL1,AL2,M,1), Nt,fc,B,M,d,1);
    FL{ci} = flc;
    b = discretize(flc(:,1,1), edges);
    COVSET(ci, b(~isnan(b))) = true;
end
fprintf('   per-offset coverage: min %.1f%%, median %.1f%%, max %.1f%% of bins\n\n', ...
        100*min(mean(COVSET,2)), 100*median(mean(COVSET,2)), 100*max(mean(COVSET,2)));

% ---------- greedy set cover ----------
order = zeros(1,NC); have = false(1,NB);
for step = 1:NC
    gainv = sum(COVSET & ~have, 2); gainv(order(order>0)) = -1;
    [~, j] = max(gainv); order(step) = j; have = have | COVSET(j,:);
end

r_ideal = ratefor(TH1_0 - (0:2)*(2/3), P,Nt,B,fc,M,d,TH1_0,TH2_0,AL1,AL2,f,CH,SNR_t,SNR_dB,Q,1);
fprintf('ideal (full TTD, baseline K=3) = %.3f\n', r_ideal);
fprintf('uniform comb needs K >= P = %d.  Can greedy do it with fewer?\n\n', P);
fprintf('%4s | %8s %6s %8s | %8s %6s %8s\n','K','uniform','%id','gap/BW','greedy','%id','gap/BW');
for K = 4:9
    ou = ratefor(TH1_0 - (0:K-1)*(2/K), P,Nt,B,fc,M,d,TH1_0,TH2_0,AL1,AL2,f,CH,SNR_t,SNR_dB,Q,0);
    og = ratefor(TH1_0 - offs(order(1:K)),  P,Nt,B,fc,M,d,TH1_0,TH2_0,AL1,AL2,f,CH,SNR_t,SNR_dB,Q,0);
    fprintf('%4d | %8.3f %5.0f%% %8.1f | %8.3f %5.0f%% %8.1f\n', ...
            K, ou(1), 100*ou(1)/r_ideal, ou(3), og(1), 100*og(1)/r_ideal, og(3));
end
fprintf('\nchosen offsets (greedy order): '); fprintf('%.3f ', offs(order(1:P))); fprintf('\n');
fprintf('uniform  offsets at K=%d      : ', P); fprintf('%.3f ', (0:P-1)*(2/P)); fprintf('\n');

% ------------------------------------------------------------------------
function out = ratefor(th1_list, P,Nt,B,fc,M,d,~,TH2,AL1,AL2,f,CH,SNR_t,SNR_dB,Q,ideal)
K = numel(th1_list);
arch='shared'; if ideal, arch='full'; end
w = zeros(Nt,K,M); focus_loc = zeros(M,2,K);
for s = 1:K
    w(:,s,:)         = ddbs_beam_arch(Nt,B,fc,M,d,th1_list(s),TH2,AL1,AL2,1,arch,P,Inf,Inf,[]);
    focus_loc(:,:,s) = cal_loc(fc,f,th1_list(s),TH2,AL1,AL2,M,1);
end
fl = actual_focus(w,focus_loc,Nt,fc,B,M,d,K);
r  = mean(cellfun(@(h) rate_ongrid(h,w,fl,Nt,B,fc,M,d,SNR_t,SNR_dB,Q,K), CH));
[gp, cov] = cov_gap(fl, Nt);
out = [r, cov, gp];
if ideal, out = r; end
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
