% ========================================================================
% IS "L = 2 GIVES NO GAIN" A PROPERTY OF L, OR OF THE SPLIT? (THEORY sec 42.2)
% ------------------------------------------------------------------------
% ablation_ps.m sets  Nsec = min(4, max(1,L)),  Ka = round(12/Nsec),  so
%
%       L = 1  -> (Nsec, Ka) = (1, 12)
%       L = 2  -> (2, 6)          <-- the only row with this split
%       L = 4, 8, 16 -> (4, 3)
%
% L is therefore CONFOUNDED with the pilot split along that sweep. Limitation #2
% of the paper ("L = 2 gives no gain, so the claim is restricted to L >= 4") is
% read off the one row that has its own split. That is the sec 18 / sec 39
% pattern again -- an effect attributed to a variable that is tied to another
% along the swept locus -- and it has to be broken before the limitation is
% quoted as a property of L.
%
% sec 42.2 tried to break it with the offline geometry, which returned a flat
% 1.43 / 1.45 / 1.44 / 1.44 across four splits. THAT TEST IS NOT ADMISSIBLE: the
% same section shows the offline statistic mispredicts the L = 2 gain by a sign
% (it says 1.43, the measurement says 0.98). A tool that is wrong about the L = 2
% gain cannot be trusted about what moves the L = 2 gain. Only the end-to-end
% simulation can settle it, which is this file.
%
% PART 1  L = 2 only, five splits at K = 12 held fixed. If the gain stays ~1.0
%         across all of them, limitation #2 is real and is about L.
%         If it rises with Nsec, limitation #2 is an artefact of the split and
%         must be rewritten -- and the L >= 4 restriction may not be needed.
%
% PART 2  the full L sweep at a FIXED split (4,3), so the ablation becomes a
%         genuine single-variable experiment. NOTE this violates Nsec >= L
%         (sec 27) at L = 8 and 16 -- but it violates it IDENTICALLY on both
%         arms, so the RATIO is still a fair single-variable comparison. Read
%         PART 2 for the ratio only; for absolute rates use ablation_ps.m.
%
% PAIRED STATISTICS. Both arms see the same channel realisations, so the gain is
% a paired comparison and the unpaired CI that ablation_ps.m propagates is
% needlessly wide (it calls itself conservative). Here the per-channel difference
% d_i = r_comp,i - r_theirs,i is formed first, and its CI is what decides
% whether the two arms differ. At L = 2 that is exactly the question.
% ========================================================================
clear; clc;

ensure_baseline();
addpath(genpath(fileparts(mfilename('fullpath'))));

Nt=256; fc=30e9; B=5e9; M=1024; c=3e8; d=(c/fc)/2; Q=1;
f = fc + B/M*((1:M)-1-(M-1)/2);
thlim=[-1 1]; alim=[0.0025 0.1]; gamma=0.9;
Rmin=5; Rmax=200; SNR_dB=20; SNR_t=10^(SNR_dB/10); N_iter=200;
KTOT = 12;  SERVE = 'shared';        % serving beam on the same Nt/L delays

rng(23); CH=cell(N_iter,1);          % same seed/channels as ablation_ps.m
for it=1:N_iter
    CH{it} = near_field_channel(Nt,d,fc,B,M, Rmin+rand*(Rmax-Rmin), asin(-0.9+1.8*rand));
end

fprintf('Split control for the L=2 limitation.  Nt=%d, SNR=%d dB, K=%d, N_iter=%d\n', Nt,SNR_dB,KTOT,N_iter);
fprintf('serving beam on the SAME Nt/L delays as training.  Paired CI on the difference.\n\n');

fprintf('PART 1 -- L = 2, K = 12 held fixed, only the (N_sec, K_alpha) split changes\n');
fprintf('%12s | %15s %15s | %17s | %13s\n','split','their PS','compensated','paired diff','gain');
splits = [2 6; 3 4; 4 3; 6 2; 12 1];
for i = 1:size(splits,1)
    Ns = splits(i,1); Ka = splits(i,2);
    [sec, rho] = ojcoms_algorithm1(Nt,fc,B,M,d,2,thlim,alim,gamma,Ns,Ka);
    vT = armvec(sec,rho,'ojcoms',2,Nt,B,fc,M,d,f,CH,SNR_t,SNR_dB,Q,SERVE);
    vC = armvec(sec,rho,'shared',2,Nt,B,fc,M,d,f,CH,SNR_t,SNR_dB,Q,SERVE);
    report(sprintf('N_sec=%d,Ka=%d',Ns,Ka), vT, vC);
end
fprintf(['\nRead PART 1: if every row reads ~1.0 with a paired CI that CONTAINS zero\n' ...
         'difference, limitation #2 is real and is about L -- keep the L>=4 wording.\n' ...
         'If the gain climbs with N_sec, the limitation is the split and must be\n' ...
         'rewritten; the L>=4 restriction would then be unnecessary.\n\n']);

fprintf('PART 2 -- full L sweep at a FIXED split (4,3).  Single-variable in L.\n');
fprintf('(N_sec >= L is violated at L = 8, 16 -- identically on both arms, so the\n');
fprintf(' RATIO is fair; absolute rates there are not the operating point.)\n');
fprintf('%12s | %15s %15s | %17s | %13s\n','L','their PS','compensated','paired diff','gain');
for L = [1 2 4 8 16]
    [sec, rho] = ojcoms_algorithm1(Nt,fc,B,M,d,L,thlim,alim,gamma,4,3);
    aT = 'ojcoms'; aC = 'shared';
    if L == 1, aT = 'full'; aC = 'full'; end          % zero control: arms coincide
    vT = armvec(sec,rho,aT,L,Nt,B,fc,M,d,f,CH,SNR_t,SNR_dB,Q,SERVE);
    vC = armvec(sec,rho,aC,L,Nt,B,fc,M,d,f,CH,SNR_t,SNR_dB,Q,SERVE);
    report(sprintf('L=%d (%d TTD)',L,Nt/L), vT, vC);
end
fprintf(['\nRead PART 2: L=1 must give a paired difference of EXACTLY zero -- with no\n' ...
         'sharing the fc term is a global constant and both arms are the same beam.\n' ...
         'Any other value there is a bug, not a result.\n']);

% ------------------------------------------------------------------------
function report(label, vT, vC)
dv = vC - vT;  n = numel(dv);
mT = mean(vT); mC = mean(vC);
dm = mean(dv); dci = 1.96*std(dv)/sqrt(n);
g  = mC/max(mT,eps);
gci = g*sqrt((1.96*std(vC)/sqrt(n)/mC)^2 + (1.96*std(vT)/sqrt(n)/mT)^2);  % unpaired, for comparison
sig = 'no';  if abs(dm) > dci, sig = 'YES'; end
fprintf('%12s | %8.3f%8s %8.3f%8s | %+7.3f+-%-6.3f %3s | %5.2f+-%-5.2f\n', ...
        label, mT, '', mC, '', dm, dci, sig, g, gci);
end

function v = armvec(sec, rho, arch, L, Nt,B,fc,M,d,f,CH,SNR_t,SNR_dB,Q,SERVE)
K = numel(sec); c=3e8; kc=2*pi*fc/c; km=2*pi*f/c;
w = zeros(Nt,K,M); cf = zeros(M,2,K);
for s = 1:K
    w(:,s,:) = ddbs_beam_arch(Nt,B,fc,M,d, sec(s).theta_t, sec(s).theta_p, ...
                              sec(s).alpha_t, sec(s).alpha_p, 1, arch, L, Inf, Inf, []);
    cf(:,1,s) = (kc./km).'*sec(s).theta_p + rho.theta*sec(s).theta_t;
    cf(:,2,s) = (kc./km).'*sec(s).alpha_p + rho.alpha*sec(s).alpha_t;
end
fl = actual_focus(w, cf, Nt, fc, B, M, d, K);        % same lookup for both arms
v  = cellfun(@(h) rate_ongrid(h,w,fl,Nt,B,fc,M,d,SNR_t,SNR_dB,Q,K,SERVE,L), CH);
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
