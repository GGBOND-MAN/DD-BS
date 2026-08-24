% ========================================================================
% Rate vs SNR at the headline operating point -- the figure reviewers expect
% ------------------------------------------------------------------------
% Everything so far is a single SNR=20 dB slice. This sweeps SNR at L=8
% (32 TTDs, an 8x reduction, enforced in BOTH stages via serve_beam) and puts
% four curves on one axis:
%
%   MRT        per-subcarrier analog MRT with perfect CSI -- unconstrained bound
%   oracle     serving beam steered at the TRUE location on the SAME Nt/L delays
%              -- the ceiling this architecture can reach, so the gap to MRT is
%              the cost of grouping and the gap to it is the cost of training
%   ours       compensated PS + the (N_sec = 8, K_alpha = 1) split, 8 pilots
%   theirs     their (14) PS + their (N_sec = 4, K_alpha = 3) split, 12 pilots
%
% Two things this figure has to settle, and one trap it avoids:
%
% 1. DOES THE GAIN HOLD AT LOW SNR? The compensated arm has more array gain per
%    pilot, so it should degrade more gracefully -- but that is a prediction, and
%    the argmax decision could fail first for either arm. If theirs closes the gap
%    below some SNR, that is a limitation to report, not to omit.
% 2. WHERE DOES EACH ARM SATURATE? "ours reaches 98% of the ceiling" is a 20 dB
%    statement until the curve shows it holds across the range.
%
% THE TRAP: the pilot counts differ (8 vs 12), so this is NOT an equal-overhead
% comparison -- it is the best configuration each side has, which favours theirs
% on overhead and ours on rate. The equal-K comparison is sector_alloc's job and
% must not be conflated with this one. Both are reported.
%
% SNR enters twice, and both are swept together: the training decision (argmax
% under noise) and the served rate. Decoupling them would flatter whichever arm
% has the better decision rule.
% ========================================================================
clear; clc;

ensure_baseline();
addpath(genpath(fileparts(mfilename('fullpath'))));

Nt=256; fc=30e9; B=5e9; M=1024; c=3e8; d=(c/fc)/2; Q=1;
f = fc + B/M*((1:M)-1-(M-1)/2); km = 2*pi*f/c; kc = 2*pi*fc/c;
thlim=[-1 1]; alim=[0.0025 0.1]; gamma=0.9;
Rmin=5; Rmax=200; N_iter=200; L=8; SERVE='shared';
SNRs = 0:5:30;

rng(23); CH=cell(N_iter,1); TRU=zeros(N_iter,2);
for it=1:N_iter
    r = Rmin + rand*(Rmax-Rmin); th = asin(-0.9+1.8*rand);
    CH{it} = near_field_channel(Nt,d,fc,B,M,r,th);
    sn = sin(th); TRU(it,:) = [sn, (1-sn^2)/(2*r)];
end

% beams and focus tables are SNR-independent -- build once, reuse across the sweep
[secO, rhoO] = ojcoms_algorithm1(Nt,fc,B,M,d,L,thlim,alim,gamma,8,1);   % ours:   8 pilots
[secT, rhoT] = ojcoms_algorithm1(Nt,fc,B,M,d,L,thlim,alim,gamma,4,3);   % theirs: 12 pilots
[wO, flO] = buildarm(secO, rhoO, 'shared', L, Nt,B,fc,M,d,f,km,kc);
[wT, flT] = buildarm(secT, rhoT, 'ojcoms', L, Nt,B,fc,M,d,f,km,kc);
fprintf('L=%d (%d TTDs, both stages), N_iter=%d.  ours: %d pilots | theirs: %d pilots\n\n', ...
        L, Nt/L, N_iter, numel(secO), numel(secT));
fprintf('%6s | %9s %9s | %19s | %19s\n','SNR','MRT','oracle','ours (8 pilots)','theirs (12 pilots)');

R = zeros(numel(SNRs),4);
for is = 1:numel(SNRs)
    SNR_dB = SNRs(is);  SNR_t = 10^(SNR_dB/10);
    mrt = 0; orc = 0;
    for it = 1:N_iter
        h = CH{it};
        mrt = mrt + mean(arrayfun(@(m) ...
              log2(1+SNR_t*abs(h(m,:)*(exp(1j*angle(h(m,:)'))/sqrt(Nt)))^2), 1:M));
        ws = serve_beam(Nt,B,fc,M,d,TRU(it,1),TRU(it,2),SERVE,L);
        orc = orc + mean(arrayfun(@(m) log2(1+SNR_t*abs(h(m,:)*ws(:,m))^2), 1:M));
    end
    [rO,cO] = armrate(wO,flO,L,Nt,B,fc,M,d,CH,SNR_t,SNR_dB,Q,SERVE);
    [rT,cT] = armrate(wT,flT,L,Nt,B,fc,M,d,CH,SNR_t,SNR_dB,Q,SERVE);
    R(is,:) = [mrt/N_iter, orc/N_iter, rO, rT];
    fprintf('%4d dB | %9.3f %9.3f | %7.3f+-%-5.3f %3.0f%% | %7.3f+-%-5.3f %3.0f%%\n', ...
            SNR_dB, R(is,1), R(is,2), rO,cO,100*rO/R(is,2), rT,cT,100*rT/R(is,2));
end
save('rate_vs_snr_results.mat','SNRs','R','L','N_iter');

figure; hold on; box on; grid on;
plot(SNRs,R(:,1),'k--','LineWidth',1.4,'DisplayName','perfect-CSI analog MRT');
plot(SNRs,R(:,2),'k-.','LineWidth',1.4,'DisplayName',sprintf('oracle location, %d TTDs',Nt/L));
plot(SNRs,R(:,3),'-o','LineWidth',1.8,'DisplayName','proposed (8 pilots)');
plot(SNRs,R(:,4),'-s','LineWidth',1.8,'DisplayName','AoSA TTD-PS baseline (12 pilots)');
xlabel('SNR (dB)'); ylabel('average achievable rate (bit/s/Hz)');
title(sprintf('N_t=%d, L=%d (%d TTDs in both stages), B=%.0f GHz', Nt, L, Nt/L, B/1e9));
legend('Location','northwest');

fprintf(['\nRead: the oracle curve is the ceiling this hardware can reach, so the\n' ...
         'MRT-oracle gap is the cost of GROUPING and the oracle-ours gap is the cost\n' ...
         'of TRAINING. Pilot counts differ (8 vs 12) -- this is best-configuration,\n' ...
         'not equal-overhead; sector_alloc.m carries the equal-K comparison.\n']);

% ------------------------------------------------------------------------
function [w, fl] = buildarm(sec, rho, arch, L, Nt,B,fc,M,d,f,km,kc)
K = numel(sec); w = zeros(Nt,K,M); cf = zeros(M,2,K);
for s = 1:K
    w(:,s,:) = ddbs_beam_arch(Nt,B,fc,M,d, sec(s).theta_t, sec(s).theta_p, ...
                              sec(s).alpha_t, sec(s).alpha_p, 1, arch, L, Inf, Inf, []);
    cf(:,1,s) = (kc./km).'*sec(s).theta_p + rho.theta*sec(s).theta_t;
    cf(:,2,s) = (kc./km).'*sec(s).alpha_p + rho.alpha*sec(s).alpha_t;
end
fl = actual_focus(w, cf, Nt, fc, B, M, d, K);
end

function [r, ci] = armrate(w,fl,L,Nt,B,fc,M,d,CH,SNR_t,SNR_dB,Q,SERVE)
K = size(w,2); v = zeros(numel(CH),1);
for it = 1:numel(CH)
    h = CH{it}; best = -inf;
    for t = 1:K
        y = zeros(M,Q);
        for m=1:M, y(m,:)=repmat(h(m,:)*w(:,t,m),[1,Q]); end
        y = add_awgn(y, SNR_dB*2/sqrt(3));
        [~,i] = max(abs(sum(y,2)).^2);
        ws = serve_beam(Nt,B,fc,M,d,fl(i,1,t),fl(i,2,t),SERVE,L);
        tt = 0; for m=1:M, tt = tt + log2(1+SNR_t*abs(h(m,:)*ws(:,m))^2)/M; end
        best = max(best,tt);
    end
    v(it) = best;
end
r = mean(v); ci = 1.96*std(v)/sqrt(numel(v));
end
