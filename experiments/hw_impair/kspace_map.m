% ========================================================================
% THE KEY FIGURE: the feasible region in (K, s) at a fixed TTD budget
% ------------------------------------------------------------------------
% map_focus_shift.m produced the decisive comparison, in a place I was not
% looking -- the PASSING configuration and the FAILING one have essentially
% IDENTICAL displacement statistics:
%
%    P=8, K=3  (42% of ideal) : stayed put 14.1%,  |k|>2 61.1%,  mean|dev| 0.617
%    P=8, K=9 (100% of ideal) : stayed put 15.7%,  |k|>2 56.7%,  mean|dev| 0.578
%
% So sharing scrambles the foci just as violently in the configuration that
% WORKS. "The beams are displaced" is therefore NOT the failure mode -- with the
% recalibrated lookup, displacement is harmless because it is known. What differs
% is only whether the displaced set still TILES: maxgap 0.3 vs 50.1 beamwidths.
%
% That resolves the law into TWO INDEPENDENT CONSTRAINTS:
%
%   (C1) COVERAGE  -- DD-BS's own design constraint. The sweep per pilot is
%        proportional to theta_p ~ s, so the total designed coverage is K*s, and
%        it must be held at the baseline value K0 = 3:          K * s >= K0
%
%   (C2) SHARING   -- the new one. Above P*|theta_t| ~ 85 the displaced foci stop
%        tiling and holes open:            P * |theta_t0| * s <= Theta ~ 85
%
%   Eliminating s:  K_min = K0 * P * |theta_t0| / Theta = 1.12 * P   -- the
%   measured frontier, now DERIVED from two constraints instead of fitted.
%
% redesign_for_shared.m only ever walked the diagonal K*s = K0, so it could not
% separate the two. decouple_theta.m broke the tie in the theta direction and
% failed (46% ceiling), which is (C1) biting. This script scans K and s as an
% INDEPENDENT 2-D grid, so the pass region should be exactly the intersection of
% the two half-planes -- a falsifiable prediction with an off-diagonal test:
%
%   K=9, s=1     : C1 satisfied (9 >= 3), C2 VIOLATED (253.8 > 85)  -> predict FAIL
%   K=3, s=1/3   : C1 VIOLATED (1 < 3),   C2 satisfied (84.6 <= 85) -> predict FAIL
%   K=9, s=1/3   : both satisfied                                   -> predict PASS
%
% Both single-constraint corners must fail for the decomposition to hold. If one
% of them passes, that constraint is not real and the law needs rebuilding.
%
% Runtime: 16 points, each an M*K focus search. K=12 rows are the slow ones.
% ========================================================================
clear; clc;

ensure_baseline();
addpath(genpath(fileparts(mfilename('fullpath'))));

Nt=256; fc=30e9; B=5e9; M=1024; c=3e8; d=(c/fc)/2; Q=1;
f = fc + B/M*((1:M)-1-(M-1)/2);
Rmin=5; Rmax=30; th_min=-pi/3; th_max=pi/3;
TH1_0=349/11; TH2_0=-30; AL1=-427/800; AL2=1859/3200; K0=3;
SNR_dB=15; SNR_t=10^(SNR_dB/10); N_iter=30; P=8; THETA=85;

rng(23); CH=cell(N_iter,1);
for it=1:N_iter
    CH{it}=near_field_channel(Nt,d,fc,B,M, Rmin+rand*(Rmax-Rmin), th_min+rand*(th_max-th_min));
end
r_ideal = runpt(TH1_0,TH2_0,1,3, Nt,B,fc,M,d,AL1,AL2,f,CH,SNR_t,SNR_dB,Q);
fprintf('P=%d (%d TTDs).  ideal (full TTD, baseline) = %.3f\n\n', P, Nt/P, r_ideal(1));
fprintf('%4s %7s | %6s %8s | %8s %6s %9s | %s\n', ...
        'K','s','K*s','P|th_t|','rate','%ideal','maxgap/BW','C1 C2 predict');

Ks = [3 6 9 12]; ss = [1 0.5 1/3 0.25];
R = nan(numel(Ks),numel(ss)); Gp = R;
for iK = 1:numel(Ks)
  for is = 1:numel(ss)
    K = Ks(iK); s = ss(is);
    TH1 = TH1_0*s; TH2 = TH2_0*s;
    out = runpt(TH1,TH2,P,K, Nt,B,fc,M,d,AL1,AL2,f,CH,SNR_t,SNR_dB,Q);
    R(iK,is) = out(1); Gp(iK,is) = out(3);
    c1 = K*s >= K0-1e-9;  c2 = P*abs(TH1) <= THETA;
    pred = 'FAIL'; if c1 && c2, pred = 'PASS'; end
    fprintf('%4d %7.3f | %6.2f %8.1f | %8.3f %5.0f%% %9.1f | %3d %2d  %s\n', ...
            K, s, K*s, P*abs(TH1), out(1), 100*out(1)/r_ideal(1), out(3), c1, c2, pred);
  end
end

save('kspace_map_results.mat','Ks','ss','R','Gp','r_ideal');
figure; imagesc(ss, Ks, 100*R/r_ideal(1)); set(gca,'YDir','normal'); colorbar;
xlabel('sweep scaling s'); ylabel('pilots K');
title(sprintf('%% of ideal rate, P=%d (%d TTDs) -- feasible region = C1 \\cap C2', P, Nt/P));
hold on;
sv = linspace(min(ss),max(ss),100);
plot(sv, K0./sv, 'w-', 'LineWidth', 2);                       % C1: K*s = K0
xline(THETA/(P*TH1_0), 'w--', 'LineWidth', 2);                % C2: P|theta_t| = Theta
legend('C1: K s = K_0','C2: P|\theta_t| = \Theta','Location','northwest');

% ------------------------------------------------------------------------
function out = runpt(TH1,TH2,P,K, Nt,B,fc,M,d,AL1,AL2,f,CH,SNR_t,SNR_dB,Q)
focus_loc = cal_loc(fc,f,TH1,TH2,AL1,AL2,M,K);
arch='full'; if P>1, arch='shared'; end
w  = ddbs_beam_arch(Nt,B,fc,M,d,TH1,TH2,AL1,AL2,K,arch,P,Inf,Inf,[]);
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
