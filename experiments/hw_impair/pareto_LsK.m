% ========================================================================
% IS K = 32 THE FLOOR? -- the (L, s) plane, not the L = 8 slice
% ------------------------------------------------------------------------
% hw_budget swept the sweep rate s at a FIXED L = 8, and reported K = 32 at the
% smallest delay range. That is the minimum ON THAT SLICE. The design space has
% two axes and they buy different things:
%
%   L   -> TTD COUNT  (N_TTD = Nt/L)   and, through N_sec >= L (sec 27), the
%          PILOT FLOOR: at L = 8 no configuration can use fewer than 8 pilots
%   s   -> DELAY RANGE (range ~ |theta_t| ~ S ~ s) and, through K_alpha, the
%          pilot MULTIPLIER
%
% So K = N_sec * K_alpha ~ L * K_alpha(s), and the two costs are separable. The
% obvious untested move is to spend TTDs to buy pilots back at the SAME range:
%
%       L = 8,  s = 1/8 :  32 TTDs,  12.9 ns,  N_sec = 8,  K_alpha = 4  -> 32 pilots
%       L = 4,  s = 1/8 :  64 TTDs,  12.9 ns,  N_sec = 4,  K_alpha = 4  -> 16 pilots
%       L = 2,  s = 1/8 : 128 TTDs,  12.9 ns,  N_sec = 2,  K_alpha = 4  ->  8 pilots
%
% If that holds, K = 32 is not a floor at all -- it is the price of insisting on
% 32 TTDs while ALSO insisting on the smallest delay range. Halving L should halve
% the pilots at unchanged range.
%
% TWO HARD FLOORS that no configuration can go below, and they should be visible:
%   K >= N_sec >= ceil(theta_tot * L / 2) = L      sharing (sec 27, Delta <= 2/P)
%   K >= 3                                          the baseline's own DDBS floor
% so at 32 TTDs the absolute minimum is 8 pilots -- and s = 1 and s = 1/2 already
% reach it (98%). The 32 belongs to s = 1/8 alone.
%
% Reported per cell: TTD count, delay range, pilots, rate. The Pareto surface is
% three-dimensional, so no single row "wins" -- the output is the frontier, and
% the paper should print it as such rather than picking one operating point.
% ========================================================================
clear; clc;

ensure_baseline();
addpath(genpath(fileparts(mfilename('fullpath'))));

Nt=256; fc=30e9; B=5e9; M=1024; c=3e8; d=(c/fc)/2; Q=1;
f = fc + B/M*((1:M)-1-(M-1)/2);
thlim=[-1 1]; alim=[0.0025 0.1]; gamma=0.9;
Rmin=5; Rmax=200; SNR_dB=20; SNR_t=10^(SNR_dB/10); N_iter=100; DEV_PS=508;

rng(23); CH=cell(N_iter,1);
for it=1:N_iter
    CH{it} = near_field_channel(Nt,d,fc,B,M, Rmin+rand*(Rmax-Rmin), asin(-0.9+1.8*rand));
end
[s0,r0] = ojcoms_algorithm1(Nt,fc,B,M,d,1,thlim,alim,gamma,1,12);
r_ref = armrate(s0,r0,'full',1,Nt,B,fc,M,d,f,CH,SNR_t,SNR_dB,Q,'ideal');
fprintf('reference (ungrouped, 12 pilots) = %.3f\n', r_ref);
fprintf('hard floors: K >= N_sec >= L (sec 27) and K >= 3\n\n');
fprintf('%3s %7s | %6s %8s %9s | %5s %7s %6s | %8s %5s\n', ...
        'L','N_TTD','s','range','x over dev','Nsec','Kalpha','K_tot','rate','%ref');

nn=(-(Nt-1)/2:(Nt-1)/2)';
for s_sw = [1 1/4 1/8]
  for L = [2 4 8]
    for Ka = [1 2 4]
      Ns = L;                                   % the sec 27 floor
      [sec, rho, dbg] = ojcoms_algorithm1(Nt,fc,B,M,d,L,thlim,alim,gamma,Ns,Ka,s_sw);
      tt = sec(1).theta_t; at = sec(1).alpha_t;
      v  = (nn*d*tt - (nn*d).^2*at)/c;  rg = max(v)-min(v);
      r  = armrate(sec,rho,'shared',L,Nt,B,fc,M,d,f,CH,SNR_t,SNR_dB,Q,'shared');
      fprintf('%3d %7d | %6.3f %7.1fns %8.0fx | %5d %7d %6d | %8.3f %4.0f%%\n', ...
              L, Nt/L, s_sw, rg*1e9, rg*1e12/DEV_PS, Ns, Ka, numel(sec), r, 100*r/r_ref);
    end
  end
  fprintf('\n');
end
fprintf(['Read: compare rows at EQUAL delay range (same s) across L. If pilots\n' ...
         'halve when L halves, K=32 was never a floor -- it was the price of\n' ...
         'demanding 32 TTDs and the smallest range at the same time. Report the\n' ...
         'frontier, not one point: TTD count, delay range and pilots are three\n' ...
         'separate costs and the designer picks which to spend.\n']);

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
