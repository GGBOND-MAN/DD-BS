% ========================================================================
% Pin K_min at unit resolution, and test the master curve's PREDICTION
% ------------------------------------------------------------------------
% Every K_min so far came from a coarse ladder: 1.12*P (steps of 3), then P
% (steps of 3-6), then 7 at P=8 (from a spacing sweep). This scans K in UNIT
% steps around P so the constant stops moving.
%
% PRE-REGISTERED PREDICTION. With the baseline comb Delta = 2/K, the master
% curve's x-variable is Delta*P/2 = P/K, and it measured
%       P/K = 1.00 -> 100%   1.14 -> 98%   1.33 -> 73%
% so the 95% crossing sits just past P/K = 1.14, i.e. K_min = ceil(P/1.14):
%
%       P =  4  ->  K_min = 4      (K=3 gives P/K = 1.33, predicted ~73%)
%       P =  8  ->  K_min = 7      (already measured: K=7 -> 98%, K=6 -> 73%)
%       P = 16  ->  K_min = 14     (K=13 -> P/K = 1.23, borderline by design)
%
% P=8 is included as a CONTROL: K=6 and K=7 were already measured at 73% and
% 98% through a different script, so if they do not reproduce here the harness
% has drifted and nothing else in the table can be trusted.
%
% If all three land, K_min = ceil(P/1.14) = ceil(0.875 P) is the law and the
% exchange constant is 0.875, not 1 and not 1.12. If P=4 and P=16 disagree with
% P=8, the constant is P-dependent and the law must be stated as a curve rather
% than a single number -- which is worth knowing before it goes in a paper.
%
% The P/K column is the check that matters: rows from DIFFERENT P but equal P/K
% must give equal rate, since the master curve says P and K enter only through
% their ratio.
%
% Runtime: 8 points; the P=16 window (K = 13..15) dominates.
% ========================================================================
clear; clc;

ensure_baseline();
addpath(genpath(fileparts(mfilename('fullpath'))));

Nt=256; fc=30e9; B=5e9; M=1024; c=3e8; d=(c/fc)/2; Q=1;
f = fc + B/M*((1:M)-1-(M-1)/2);
Rmin=5; Rmax=30; th_min=-pi/3; th_max=pi/3;
TH1=349/11; TH2=-30; AL1=-427/800; AL2=1859/3200;
SNR_dB=15; SNR_t=10^(SNR_dB/10); N_iter=30; PASS=0.95;

rng(23); CH=cell(N_iter,1);
for it=1:N_iter
    CH{it}=near_field_channel(Nt,d,fc,B,M, Rmin+rand*(Rmax-Rmin), th_min+rand*(th_max-th_min));
end
r_ideal = evalK(1, 3, Nt,B,fc,M,d,TH1,TH2,AL1,AL2,f,CH,SNR_t,SNR_dB,Q);
fprintf('ideal (full TTD, K=3) = %.3f.  pass = %.0f%%\n\n', r_ideal(1), 100*PASS);
fprintf('%4s %4s %7s | %8s %6s %9s | %s\n','P','K','P/K','rate','%ideal','gap/BW','predicted K_min');

% Trimmed to the boundary neighbourhood at each P: cells far from the crossing
% cost actual_focus time without discriminating between the candidate laws.
windows = {[4, 3:5], [8, 6:7], [16, 13:15]};
for wi = 1:numel(windows)
    P = windows{wi}(1); Ks = windows{wi}(2:end);
    Kmin = NaN;
    for K = Ks
        out = evalK(P, K, Nt,B,fc,M,d,TH1,TH2,AL1,AL2,f,CH,SNR_t,SNR_dB,Q);
        pct = 100*out(1)/r_ideal(1);
        if isnan(Kmin) && out(1) >= PASS*r_ideal(1), Kmin = K; end
        fprintf('%4d %4d %7.3f | %8.3f %5.0f%% %9.1f |%s\n', ...
                P, K, P/K, out(1), pct, out(3), ...
                repmat(' <-- K_min',1, double((K==Kmin)&&~isnan(Kmin))));
    end
    fprintf('   P=%2d : measured K_min = %s , predicted %d , K_min/P = %.3f\n\n', ...
            P, num2str(Kmin), ceil(P/1.14), Kmin/P);
end
fprintf(['\nCross-check -- the master curve says P and K enter ONLY through their\n' ...
         'ratio, so these pairs, measured at DIFFERENT P, must give equal rate:\n' ...
         '   P/K = 1.333 : (P=4,K=3) vs (P=8,K=6)\n' ...
         '   P/K = 1.143 : (P=8,K=7) vs (P=16,K=14)\n' ...
         'If either pair disagrees, rate is not a function of P/K alone and the\n' ...
         'master curve is a coincidence of the P=8 sweep it was built from.\n']);

% ------------------------------------------------------------------------
function out = evalK(P, K, Nt,B,fc,M,d,TH1,TH2,AL1,AL2,f,CH,SNR_t,SNR_dB,Q)
focus_loc = cal_loc(fc,f,TH1,TH2,AL1,AL2,M,K);
arch='full'; if P>1, arch='shared'; end
w  = ddbs_beam_arch(Nt,B,fc,M,d,TH1,TH2,AL1,AL2,K,arch,P,Inf,Inf,[]);
fl = actual_focus(w,focus_loc,Nt,fc,B,M,d,K);
r  = mean(cellfun(@(h) rate_ongrid(h,w,fl,Nt,B,fc,M,d,SNR_t,SNR_dB,Q,K), CH));
[gp, cov] = cov_gap(fl, Nt);
out = [r, cov, gp];
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
