% ========================================================================
% Where do the focus points actually GO under shared TTD?
% ------------------------------------------------------------------------
% diag_shared_mechanism.m settled two things:
%   * gain truncation (M1) is REFUTED -- measured usable bandwidth is 5-25x
%     wider than the sub-array-factor prediction;
%   * coverage HOLES (M2) are the mechanism -- maxgap/(2/Nt) separates every
%     passing configuration (<=1) from every failing one (>=5), with no overlap.
% And decouple_theta.m refuted the theta_t/theta_p RATIO design equation: at
% P=8, K=3 no (theta_t, theta_p) pair exceeds 46% of ideal, and shrinking theta_t
% alone makes things WORSE (12% at s_t=1/16).
%
% So the holes are real but their cause is still open. This script stops guessing
% and MEASURES the focus map itself, testing one specific structural hypothesis:
%
%   GRATING-LOBE SELECTION. Sharing turns the TTD layer into an array of Nt/P
%   elements spaced P*d = P*lambda/2 apart, so its array factor has GRATING LOBES
%   every 1/P in sin(theta). The per-antenna PS forms the sub-array pattern that
%   normally selects the correct lobe; away from fc that pattern tilts, so a
%   DIFFERENT lobe wins and the realized focus JUMPS by an integer multiple of
%   1/P. Under this hypothesis the focus map is not blurred but PIECEWISE
%   DISPLACED, which is exactly how holes (rather than smearing) would appear.
%
% Falsifiable signature: (fl - designed) * P should cluster on INTEGERS.
%   * integers, |k| >= 1 for many m  -> grating-lobe selection confirmed
%   * continuous / non-integer       -> refuted; the map is smoothly warped
%
% Also reports which designed focus angles are the ones that go missing, i.e.
% what the biggest hole was supposed to contain.
% ========================================================================
clear; clc;

ensure_baseline();
addpath(genpath(fileparts(mfilename('fullpath'))));

Nt=256; fc=30e9; B=5e9; M=1024; c=3e8; d=(c/fc)/2;
f = fc + B/M*((1:M)-1-(M-1)/2);
TH1_0=349/11; TH2_0=-30; AL1=-427/800; AL2=1859/3200; K0=3;

for cfg = {[1 3],[4 3],[8 3],[8 6],[8 9]}
    P = cfg{1}(1); K = cfg{1}(2);
    s = K0/K; TH1 = TH1_0*s; TH2 = TH2_0*s;
    focus_loc = cal_loc(fc,f,TH1,TH2,AL1,AL2,M,K);
    arch='full'; if P>1, arch='shared'; end
    w  = ddbs_beam_arch(Nt,B,fc,M,d,TH1,TH2,AL1,AL2,K,arch,P,Inf,Inf,[]);
    fl = actual_focus(w,focus_loc,Nt,fc,B,M,d,K);

    des = squeeze(focus_loc(:,1,:));  act = squeeze(fl(:,1,:));
    dev = act - des;                          % displacement in sin(theta)
    kk  = dev * P;                            % in units of the grating-lobe step 1/P
    frac = abs(kk - round(kk));               % 0 => lands exactly on a lobe

    fprintf('\n=== P=%d (N_TTD=%d), K=%d, theta_t=%.2f, P|theta_t|=%.1f ===\n', ...
            P, Nt/P, K, TH1, P*abs(TH1));
    fprintf('displacement |dev|: mean %.4f, max %.4f   (grating step 1/P = %.4f)\n', ...
            mean(abs(dev(:))), max(abs(dev(:))), 1/P);
    fprintf('stayed put (|dev| < 1/(2P)) : %5.1f%%\n', 100*mean(abs(dev(:))<1/(2*P)));
    fprintf('lobe index k = dev*P  ->  k=0:%4.1f%%  |k|=1:%4.1f%%  |k|=2:%4.1f%%  |k|>2:%4.1f%%\n', ...
            100*mean(round(kk(:))==0), 100*mean(abs(round(kk(:)))==1), ...
            100*mean(abs(round(kk(:)))==2), 100*mean(abs(round(kk(:)))>2));
    fprintf('distance to nearest lobe |k - round(k)| : mean %.3f  (0 = on a lobe,\n', mean(frac(:)));
    fprintf('   0.25 = uniform/no lobe structure)\n');

    % where is the biggest hole, and what was supposed to fill it?
    sv = sort(act(:)); sv = sv(abs(sv)<=0.866);
    [gap, i] = max(diff(sv)); lo = sv(i); hi = sv(i+1);
    inhole = des >= lo & des <= hi;
    fprintf('biggest hole in the REALIZED map: [%+.3f, %+.3f], width %.4f = %.1f beamwidths\n', ...
            lo, hi, gap, gap/(2/Nt));
    fprintf('   designed foci that fell inside it: %d of %d (%.1f%%) -- these are the\n', ...
            sum(inhole(:)), numel(des), 100*mean(inhole(:)));
    fprintf('   users that lose their beam. They ended up at: mean %+.3f, spread %.3f\n', ...
            mean(act(inhole)), std(act(inhole)));
end
fprintf('\nRead: if |k - round(k)| is near 0 and |k|>=1 is common, sharing is\n');
fprintf('SELECTING THE WRONG GRATING LOBE and the fix is lobe disambiguation\n');
fprintf('(sub-array-level design), not more pilots. If it is near 0.25 the map is\n');
fprintf('smoothly warped and the fix is a redesigned sweep law.\n');
