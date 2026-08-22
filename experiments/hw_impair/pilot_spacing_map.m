% ========================================================================
% THE DECISIVE TEST -- in the right variables this time: (K, Delta)
% ------------------------------------------------------------------------
% kspace_map.m refuted BOTH constraints. Read its table by ROWS and s barely
% moves anything; read it by COLUMNS and K decides everything:
%
%   K= 3 :  43%  46%  39%  42%      maxgap ~ 50 beamwidths
%   K= 6 :  72%  87%  86%  80%      maxgap ~ 9-10
%   K= 9 : 100% 102% 100%  96%      maxgap 0.3-0.6
%   K=12 : 101% 101%  99% 100%      maxgap 0.3-0.6
%
% K=9, s=1 has P*|theta_t| = 253.8, three times over the supposed limit, and it
% delivers 100%. So the P*|theta_t| <= 85 invariant was never real: on the
% redesign diagonal theta_t was a function of K, and the invariant was K in
% disguise. C1 (K*s >= K0) dies too -- K=9, s=0.25 has K*s = 2.25 and passes.
%
% What survives: K_min depends on P alone -- not on s (this table), not on B
% (frontier PART 2), not on theta_t. And K_min tracks P linearly:
%       P     2   4   8  16
%       K_min 3   5   9  18        ~ P + 1  (or 1.12 P -- not yet separable)
%
% MECHANISM THAT FITS ALL OF IT. The baseline places its K pilots by stepping the
% TTD intercept, theta_t,s = theta_t - 2(s-1)/K, so the pilots form a comb of
% SPACING Delta = 2/K spanning ~2. Sharing groups P antennas, whose sub-array
% pattern has beamwidth ~2/P. A hole opened by sharing is filled only if some
% pilot lands inside it, which needs the comb finer than the sub-array beam:
%
%       (R1) RESOLUTION   Delta <= 2/P          -> with Delta = 2/K:  K >= P
%       (R2) SPAN         K * Delta >= 2        -> the strips must still tile
%
% Neither contains B, theta_t or s -- which is exactly why K_min was measured to
% be independent of all three. Together they give K >= P, against a measured
% 1.12 P. No other hypothesis so far explains the null results as well as the
% positive one.
%
% This script breaks the baseline's hard-wired Delta = 2/K by building the pilots
% one at a time, so K and Delta scan INDEPENDENTLY. The pass region must be the
% corner {Delta <= 2/P} AND {K*Delta >= 2}:
%
%   P=8, K= 4, Delta=0.50 : span ok, TOO COARSE   -> predict FAIL
%   P=8, K= 4, Delta=0.25 : fine enough, TOO NARROW -> predict FAIL
%   P=8, K= 8, Delta=0.25 : both ok               -> predict PASS  (the K=P corner)
%   P=8, K=16, Delta=0.125: both ok               -> predict PASS
%
% If K=8/Delta=0.25 passes while K=16/Delta=0.5 fails, pilot COUNT is not the
% variable -- pilot RESOLUTION is, and the fix becomes "place the pilots better",
% which is cheaper than "use more of them".
% ========================================================================
clear; clc;

ensure_baseline();
addpath(genpath(fileparts(mfilename('fullpath'))));

Nt=256; fc=30e9; B=5e9; M=1024; c=3e8; d=(c/fc)/2; Q=1;
f = fc + B/M*((1:M)-1-(M-1)/2);
Rmin=5; Rmax=30; th_min=-pi/3; th_max=pi/3;
TH1_0=349/11; TH2_0=-30; AL1=-427/800; AL2=1859/3200;
SNR_dB=15; SNR_t=10^(SNR_dB/10); N_iter=30; P=8;

rng(23); CH=cell(N_iter,1);
for it=1:N_iter
    CH{it}=near_field_channel(Nt,d,fc,B,M, Rmin+rand*(Rmax-Rmin), th_min+rand*(th_max-th_min));
end
r_ideal = runcomb(TH1_0,TH2_0,1,3,2/3, Nt,B,fc,M,d,AL1,AL2,f,CH,SNR_t,SNR_dB,Q);
fprintf('P=%d (%d TTDs), s=1 throughout.  ideal = %.3f.  sub-array beamwidth 2/P = %.3f\n\n', ...
        P, Nt/P, r_ideal(1), 2/P);
fprintf('%4s %8s | %8s %8s | %8s %6s %9s | %s\n', ...
        'K','Delta','K*Delta','Delta*P/2','rate','%ideal','maxgap/BW','R1 R2 predict');

for K = [4 6 8 12 16]
  for Del = [0.5 0.25 0.125 2/K]
    out = runcomb(TH1_0,TH2_0,P,K,Del, Nt,B,fc,M,d,AL1,AL2,f,CH,SNR_t,SNR_dB,Q);
    r1 = Del <= 2/P + 1e-9;  r2 = K*Del >= 2 - 1e-9;
    pred='FAIL'; if r1 && r2, pred='PASS'; end
    tag=''; if abs(Del-2/K)<1e-9, tag=' <- baseline comb'; end
    fprintf('%4d %8.3f | %8.2f %8.2f | %8.3f %5.0f%% %9.1f | %3d %2d  %s%s\n', ...
            K, Del, K*Del, Del*P/2, out(1), 100*out(1)/r_ideal(1), out(3), r1, r2, pred, tag);
  end
end
fprintf(['\nRead: if the pass region is the (R1 and R2) corner, the design rule is\n' ...
         '"pilot comb finer than the sub-array beam, still spanning the sweep".\n' ...
         'If instead only K matters and Delta does not, the comb model is wrong too\n' ...
         'and K_min ~ 1.12 P stands as a purely empirical law.\n']);

% ------------------------------------------------------------------------
function out = runcomb(TH1,TH2,P,K,Del, Nt,B,fc,M,d,AL1,AL2,f,CH,SNR_t,SNR_dB,Q)
% Same as the baseline DDBS pilot set, except the intercept comb spacing Del is
% free instead of hard-wired to 2/K. Each pilot is generated as its own K=1 beam.
w = zeros(Nt,K,M); focus_loc = zeros(M,2,K);
arch='full'; if P>1, arch='shared'; end
for s = 1:K
    th1s = TH1 - Del*(s-1);
    w(:,s,:)      = ddbs_beam_arch(Nt,B,fc,M,d,th1s,TH2,AL1,AL2,1,arch,P,Inf,Inf,[]);
    focus_loc(:,:,s) = cal_loc(fc,f,th1s,TH2,AL1,AL2,M,1);
end
fl = actual_focus(w,focus_loc,Nt,fc,B,M,d,K);
r  = mean(cellfun(@(h) rate_ongrid(h,w,fl,Nt,B,fc,M,d,SNR_t,SNR_dB,Q,K), CH));
sn = sort(reshape(fl(:,1,:),[],1));
cov = (min(max(sn),0.866) - max(min(sn),-0.866))/1.732;
inner = sn(sn>=-0.866 & sn<=0.866);
out = [r, max(cov,0), max(diff(inner))/(2/Nt)];
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
