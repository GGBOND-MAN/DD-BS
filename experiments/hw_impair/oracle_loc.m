% ========================================================================
% Is there ANY headroom left for a better estimator? -- the oracle bound
% ------------------------------------------------------------------------
% phase_refine.m produced a clean negative: a range estimate accurate to ~3 cm
% (gate: 5/5 within 0.32 m) makes the rate WORSE, not better --
%       ours  1.00x  0.91x  0.99x  0.80x  0.94x
% Before writing that off as a failed idea, the right question is whether the
% failure is in the refinement or in the premise. So: hand the pipeline the
% TRUE (theta, alpha) and see whether the rate moves at all.
%
%   ORACLE ~ COARSE  -> the recalibrated table is already at the beamforming
%                       optimum, NO estimator can help, and the whole
%                       estimator direction closes with a proof rather than a
%                       shrug. That also strengthens the lookup's story: it is
%                       not merely adequate, it is where extra geometric
%                       accuracy stops paying.
%   ORACLE >> COARSE -> there is real headroom, the phase refinement was simply
%                       designed wrong (alpha-only, from a theta it did not
%                       re-estimate), and one joint (theta, alpha) attempt is
%                       justified.
%
% WHY ALPHA-ONLY REFINEMENT CAN HURT, and why the oracle settles it. The serving
% beam uses the PAIR (theta_hat, alpha_hat) from the coarse table, and that pair
% is jointly consistent -- it is where the beam actually focuses, so the array
% can form it. Near-field focusing couples theta and range, so a theta error is
% partly absorbed by an alpha offset, and the coarse pair may already sit at that
% compensated optimum. Replacing alpha with a physically accurate value derived
% from an INACCURATE theta breaks the compensation. On top of that, the range
% matched filter de-embeds the beam term evaluated at the coarse point, so its
% own accuracy depends on the coarse estimate being good -- which is exactly what
% sharing destroys. It worked in the gate because the gate used an ideal L=1 beam.
% ========================================================================
clear; clc;

ensure_baseline();
addpath(genpath(fileparts(mfilename('fullpath'))));

Nt=256; fc=30e9; B=5e9; M=1024; c=3e8; d=(c/fc)/2; Q=1;
f = fc + B/M*((1:M)-1-(M-1)/2); km = 2*pi*f/c;
thlim=[-1 1]; alim=[0.0025 0.1]; gamma=0.9;
Rmin=5; Rmax=200; SNR_dB=20; SNR_t=10^(SNR_dB/10); N_iter=20;

rng(23); CH=cell(N_iter,1); TRU=zeros(N_iter,2);
for it=1:N_iter
    r = Rmin + rand*(Rmax-Rmin); th = asin(-0.9+1.8*rand);
    CH{it} = near_field_channel(Nt,d,fc,B,M,r,th);
    sn = sin(th);  TRU(it,:) = [sn, (1-sn^2)/(2*r)];
end

% absolute ceiling: per-subcarrier analog MRT with perfect CSI
mrt = mean(cellfun(@(h) mean(arrayfun(@(m) ...
        log2(1+SNR_t*abs(h(m,:)*(exp(1j*angle(h(m,:)'))/sqrt(Nt)))^2), 1:M)), CH));
fprintf('perfect-CSI analog MRT ceiling = %.3f bit/s/Hz (UNCONSTRAINED: per-antenna phase)\n', mrt);
fprintf('oracle = serving beam steered at the TRUE (theta, alpha); no search.\n\n');
fprintf('%3s %7s %6s %5s | %8s | %8s %8s %8s | %8s\n', ...
        'L','N_TTD','Nsec','K','arm','coarse','oracle','ratio','vs MRT');

cfg = { 1,1,12; 4,4,12; 8,4,12; 8,8,8; 16,12,12 };
for SERVEc = {'ideal','shared'}
SERVE = SERVEc{1};
fprintf('\n-- serving beam: %s --\n', SERVE);
for i = 1:size(cfg,1)
    L = cfg{i,1}; Nsec = cfg{i,2}; Ka = max(1,round(cfg{i,3}/Nsec));
    [sec, rho] = ojcoms_algorithm1(Nt,fc,B,M,d,L,thlim,alim,gamma,Nsec,Ka);
    for arch = {'ojcoms','shared'}
        a = arch{1};  if L==1, a2='full'; else, a2=a; end
        [rc, ro] = pairrate(sec,rho,a2,L,Nt,B,fc,M,d,f,km,CH,TRU,SNR_t,SNR_dB,Q,SERVE);
        nm = 'them'; if strcmp(a,'shared'), nm='ours'; end
        if L==1, nm='ungrouped'; end
        fprintf('%3d %7d %6d %5d | %8s | %8.3f %8.3f %7.2fx | %7.0f%%\n', ...
                L, Nt/L, Nsec, numel(sec), nm, rc, ro, ro/max(rc,eps), 100*ro/mrt);
        if L==1, break; end
    end
end
end
fprintf(['\nRead the ratio column. Near 1.00 means the coarse table already sits\n' ...
         'at the beamforming optimum and no estimator -- ours, theirs, MUSIC or\n' ...
         'otherwise -- can recover anything, which closes the estimator direction\n' ...
         'with a bound instead of a failed attempt. Well above 1.00 means the\n' ...
         'headroom is real and the phase refinement was simply mis-designed.\n']);

% ------------------------------------------------------------------------
function [rc, ro] = pairrate(sec,rho,arch,L,Nt,B,fc,M,d,f,km,CH,TRU,SNR_t,SNR_dB,Q,SERVE)
K = numel(sec); c=3e8; kc=2*pi*fc/c;
w = zeros(Nt,K,M); cf = zeros(M,2,K);
for s = 1:K
    w(:,s,:) = ddbs_beam_arch(Nt,B,fc,M,d, sec(s).theta_t, sec(s).theta_p, ...
                              sec(s).alpha_t, sec(s).alpha_p, 1, arch, L, Inf, Inf, []);
    cf(:,1,s) = (kc./km).'*sec(s).theta_p + rho.theta*sec(s).theta_t;
    cf(:,2,s) = (kc./km).'*sec(s).alpha_p + rho.alpha*sec(s).alpha_t;
end
fl = actual_focus(w, cf, Nt, fc, B, M, d, K);
acc = 0; acco = 0;
for it = 1:numel(CH)
    h = CH{it};
    acc  = acc  + coarse_rate(h,w,fl,Nt,B,fc,M,d,SNR_t,SNR_dB,Q,K,SERVE,L);
    % oracle: true location, served on the SAME hardware as the coarse arm, so
    % the ratio isolates estimation error and not the serving stage
    ws   = serve_beam(Nt,B,fc,M,d,TRU(it,1),TRU(it,2),SERVE,L);
    t=0; for m=1:M, t=t+log2(1+SNR_t*abs(h(m,:)*ws(:,m))^2)/M; end
    acco = acco + t;
end
rc = acc/numel(CH);  ro = acco/numel(CH);
end

function best = coarse_rate(h,w,fl,Nt,B,fc,M,d,SNR_t,SNR_dB,Q,K,SERVE,P)
best=-inf;
for t=1:K
    y=zeros(M,Q);
    for m=1:M, y(m,:)=repmat(h(m,:)*w(:,t,m),[1,Q]); end
    y = add_awgn(y, SNR_dB*2/sqrt(3));
    [~,i]=max(abs(sum(y,2)).^2);
    ws = serve_beam(Nt,B,fc,M,d,fl(i,1,t),fl(i,2,t),SERVE,P);
    tt=0; for m=1:M, tt=tt+log2(1+SNR_t*abs(h(m,:)*ws(:,m))^2)/M; end
    best=max(best,tt);
end
end
