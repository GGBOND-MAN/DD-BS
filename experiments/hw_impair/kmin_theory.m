% ========================================================================
% DIRECT TEST OF THE DERIVED COVERAGE LAW (THEORY.md sec 40)
% ------------------------------------------------------------------------
% sec 39 left  K_min = max(L, 3, C/s)  with  C ~ 4  FITTED on six cells and no
% mechanism. sec 40 derives it, and in doing so CORRECTS it: the coverage
% constraint binds on N_sec, not on K.
%
%   focus locus (exact):  theta_m = < u_m*theta'p + rho_theta*theta_t >_{2*u_m}
%                         alpha_m =   u_m*alpha'p + rho_alpha*alpha_t
%                         u_m = fc/f_m,  <x>_p = x mod p into [-p/2,p/2)
%
%   The modulus 2*u_m = 2*fc/f_m is the grating period AT f_m -- frequency
%   dependent. Eliminating u_m makes the locus a PENCIL OF LINES through the
%   pivot (rho_theta*theta_t, rho_alpha*alpha_t), one line per grating order.
%
%   alpha coverage per pilot = alpha'p * Lambda = H  (their eq (69)), for ANY s.
%   theta traversals per pilot:
%           n_tr = |q|*Lambda/2,   q = theta'p - 2*round((theta'p+theta_piv)/2)
%                ~ (fc/fH)*S*Lambda/2   ∝  s
%
%   COVERAGE LAW (derived, zero fitted parameters):
%
%           N_sec * n_tr  >=  H / (2*d_alpha)  =  H*d*Nt^2/4
%
%   with d_alpha = 2/(d*Nt^2) the near-field focusing depth (pi/2 edge phase).
%   In the sweep-rate form,  N_sec*s >= C  with
%
%           C = H*d*Nt^3*xi_H^2 / (3.52*gamma*M)  = 2.96  at these parameters
%
%   (sec 39's fitted C ~ 4 was the upper end of a factor-of-two bracket, AND it
%   was attributed to K because hw_budget held N_sec = 8 while sweeping K_alpha
%   -- the sec 18 mistake repeating: an invariant fitted on a locus where the
%   candidate variables are tied.)
%
% WHAT THIS SCRIPT TESTS
%   PART 1  sweeps N_sec at K_alpha = 1 and prints, per cell, the measured rate
%           beside the derived pass/fail call. The law is falsified if the
%           >=95% boundary does not sit where N_sec*n_tr crosses H*d*Nt^2/4.
%   PART 2  shows K_alpha is a DIFFERENT primitive, not a substitute: sec 40.6
%           finds the interleaved pilot works by pushing 40% of its sweep below
%           the realizable alpha range, where the recalibrated lookup pins it to
%           alpha_min and it becomes a flat angular scan. Offline (Python port)
%           this reaches the criterion at K = 16 where sectors alone need K ~ 42.
%
% PREDICTION TO CHECK IN PART 1 (offline coverage statistic, sec 40.5):
%   s = 1/8 : N_sec = 24 -> fail, 32 -> fail, 42 -> PASS, 48 -> PASS
%   s = 1/4 : N_sec = 8, 12 -> fail, 16 -> PASS
%   s = 1/2 : N_sec = 8 -> PASS (marginal, 9.68 vs threshold 8.21)
% ========================================================================
clear; clc;

ensure_baseline();
addpath(genpath(fileparts(mfilename('fullpath'))));

Nt=256; fc=30e9; B=5e9; M=1024; c=3e8; d=(c/fc)/2; Q=1; L=8;
f = fc + B/M*((1:M)-1-(M-1)/2);
thlim=[-1 1]; alim=[0.0025 0.1]; gamma=0.9;
Rmin=5; Rmax=200; SNR_dB=20; SNR_t=10^(SNR_dB/10); N_iter=100;
fL=fc-B/2; fH=fc+B/2; Lam=fc/fL-fc/fH; H=alim(2)-alim(1);
d_alpha = 2/(d*Nt^2);                 % focusing depth, pi/2 edge phase
T       = H/(2*d_alpha);              % = H*d*Nt^2/4
C_pred  = H*d*Nt^3*(fH/fc)^2/(3.52*gamma*M);

rng(23); CH=cell(N_iter,1);
for it=1:N_iter
    CH{it} = near_field_channel(Nt,d,fc,B,M, Rmin+rand*(Rmax-Rmin), asin(-0.9+1.8*rand));
end
[s0,r0] = ojcoms_algorithm1(Nt,fc,B,M,d,1,thlim,alim,gamma,1,12);
r_ref = armrate(s0,r0,'full',1,Nt,B,fc,M,d,f,CH,SNR_t,SNR_dB,Q,'ideal');

fprintf('Lambda = %.6f   H = %.4f   d_alpha = %.4e   threshold H/(2*d_alpha) = %.3f\n', Lam,H,d_alpha,T);
fprintf('derived C (N_sec*s form) = %.3f      reference (ungrouped, 12 pilots) = %.3f\n\n', C_pred, r_ref);

fprintf('PART 1 -- N_sec sweep at K_alpha = 1, L = %d (%d TTDs).  N_sec >= L is the sec 27 floor.\n', L, Nt/L);
fprintf('%6s %6s %6s | %8s %10s | %8s %5s | %s\n', ...
        's','N_sec','K','n_tr','N_sec*n_tr','rate','%ref','derived call');
for s_sw = [1 1/2 1/4 1/8]
  for Ns = [8 12 16 24 32 42 48]
      [sec, rho, dbg] = ojcoms_algorithm1(Nt,fc,B,M,d,L,thlim,alim,gamma,Ns,1,s_sw);
      tp  = sec(1).theta_p;  piv = rho.theta*sec(1).theta_t;
      q   = tp - 2*round((tp+piv)/2);
      ntr = abs(q)*Lam/2;
      r   = armrate(sec,rho,'shared',L,Nt,B,fc,M,d,f,CH,SNR_t,SNR_dB,Q,'shared');
      call = 'fail';  if Ns*ntr >= T, call = 'PASS'; end
      fprintf('%6.3f %6d %6d | %8.4f %10.2f | %8.3f %4.0f%% | %s\n', ...
              s_sw, Ns, numel(sec), ntr, Ns*ntr, r, 100*r/r_ref, call);
  end
  fprintf('\n');
end

fprintf('PART 2 -- K_alpha is a different primitive (sec 40.6), not more sectors.\n');
fprintf('%6s %6s %3s %6s | %10s | %8s %5s\n','s','N_sec','Ka','K','N_sec*n_tr','rate','%ref');
for s_sw = [1/4 1/8]
  for Ka = [1 2 4]
      [sec, rho, dbg] = ojcoms_algorithm1(Nt,fc,B,M,d,L,thlim,alim,gamma,8,Ka,s_sw);
      tp  = sec(1).theta_p;  piv = rho.theta*sec(1).theta_t;
      ntr = abs(tp - 2*round((tp+piv)/2))*Lam/2;
      r   = armrate(sec,rho,'shared',L,Nt,B,fc,M,d,f,CH,SNR_t,SNR_dB,Q,'shared');
      fprintf('%6.3f %6d %3d %6d | %10.2f | %8.3f %4.0f%%\n', s_sw,8,Ka,numel(sec),8*ntr,r,100*r/r_ref);
  end
end
fprintf(['\nRead PART 1: the >=95%% boundary should coincide with the PASS/fail column.\n' ...
         'Read PART 2: these cells reach >=95%% while N_sec*n_tr stays far below the\n' ...
         'threshold -- K_alpha does not buy sector coverage, it adds a flat alpha_min\n' ...
         'angular scan that only the RECALIBRATED lookup can exploit (sec 40.6).\n']);

% ------------------------------------------------------------------------
function r = armrate(sec, rho, arch, L, Nt,B,fc,M,d,f,CH,SNR_t,SNR_dB,Q,SERVE)
K = numel(sec); c=3e8; kc=2*pi*fc/c; km=2*pi*f/c;
w = zeros(Nt,K,M); cf = zeros(M,2,K);
for s = 1:K
    w(:,s,:) = ddbs_beam_arch(Nt,B,fc,M,d, sec(s).theta_t, sec(s).theta_p, ...
                              sec(s).alpha_t, sec(s).alpha_p, 1, arch, L, Inf, Inf, []);
    cf(:,1,s) = (kc./km).'*sec(s).theta_p + rho.theta*sec(s).theta_t;
    cf(:,2,s) = (kc./km).'*sec(s).alpha_p + rho.alpha*sec(s).alpha_t;
end
fl = actual_focus(w, cf, Nt, fc, B, M, d, K);
r  = mean(cellfun(@(h) rate1(h,w,fl,Nt,B,fc,M,d,SNR_t,SNR_dB,Q,K,SERVE,L), CH));
end

function best = rate1(h,w,fl,Nt,B,fc,M,d,SNR_t,SNR_dB,Q,K,SERVE,P)
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
