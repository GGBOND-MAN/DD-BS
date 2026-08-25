% ========================================================================
% THE CLOSING FIGURE: the three-way budget (TTD count, delay RANGE, pilots)
% ------------------------------------------------------------------------
% Everything so far counted TTD ELEMENTS. That is only half the hardware bill,
% and the half this project has not yet paid:
%
%   training beam, theta_t ~ -31.8 (their sectors, and ours)  ->  135 ns
%   training beam, sweep slowed by s = 1/2, 1/4, 1/8          ->  68, 34, 17 ns
%   serving beam (focus only)                                 ->  4.2 ns
%   realizable TTD devices, sec 13 survey                     ->  1.47 - 508 ps
%
% 135 ns is **40 m of transmission line per element**, ~270x the best device in
% the survey. Sharing divides the COUNT by L; it does not touch the RANGE each
% surviving element must cover. So the sec 13 unrealizability problem is NOT
% solved by anything in sec 26-34, and saying otherwise would be false.
%
% What does reduce it is the sec 14 knob: slowing the angular sweep scales
% theta_t, and the range with it, at the cost of pilots. Sharing and sweep-rate
% are INDEPENDENT axes, so the real design space is three-dimensional:
%
%       N_TTD = Nt/L        (element count)   <- sharing, sec 27
%       range ~ |theta_t| ~ s (per-element)   <- sweep rate, sec 14
%       K = N_sec * K_alpha (overhead)        <- pays for both
%
% This script measures that space at the headline point L=8 and reports the full
% bill, so the paper can state a realizable operating point instead of an element
% count. N_sec is held at L (the sec 27 rule, Delta <= 2/P); K_alpha is swept,
% because slowing the sweep thins the strips and the distance dimension is what
% has to absorb it.
%
% The honest outcome may be that no point on this grid is realizable with today's
% devices. If so, that is the finding: quantify the gap and say which device
% technology would close it, rather than quietly reporting element counts alone.
% ========================================================================
clear; clc;

ensure_baseline();
addpath(genpath(fileparts(mfilename('fullpath'))));

Nt=256; fc=30e9; B=5e9; M=1024; c=3e8; d=(c/fc)/2; Q=1;
f = fc + B/M*((1:M)-1-(M-1)/2);
thlim=[-1 1]; alim=[0.0025 0.1]; gamma=0.9;
Rmin=5; Rmax=200; SNR_dB=20; SNR_t=10^(SNR_dB/10); N_iter=100;
L = 8;  DEV_PS = 508;                       % best realizable range in the survey [ps]

rng(23); CH=cell(N_iter,1);
for it=1:N_iter
    CH{it} = near_field_channel(Nt,d,fc,B,M, Rmin+rand*(Rmax-Rmin), asin(-0.9+1.8*rand));
end
[sec0, rho0] = ojcoms_algorithm1(Nt,fc,B,M,d,1,thlim,alim,gamma,1,12);
r_ref = armrate(sec0, rho0, 'full', 1, Nt,B,fc,M,d,f,CH,SNR_t,SNR_dB,Q,'ideal');
fprintf('reference (ungrouped, 12 pilots) = %.3f | L=%d -> N_TTD=%d in both stages\n', ...
        r_ref, L, Nt/L);
fprintf('best realizable TTD range in the survey: %g ps\n\n', DEV_PS);
fprintf('%6s %8s %8s %8s %8s | %5s %7s %6s | %8s %5s | %9s\n', ...
        's','S','theta_t','range','line/el','Nsec','Kalpha','K_tot','rate','%ref','x over dev');

nn=(-(Nt-1)/2:(Nt-1)/2)';
for s_sw = [1 1/2 1/4 1/8]
  for Ka = [1 2 4]
    % CORRECTED: slow the sweep by scaling S inside Algorithm 1, which re-solves
    % (63) so the sector centres stay put. The first version scaled theta_t and
    % theta_p here instead, which dragged every sector to boresight -- see the
    % s_sweep note in ojcoms_algorithm1.m.
    [sec, rho] = ojcoms_algorithm1(Nt,fc,B,M,d,L,thlim,alim,gamma,L,Ka,s_sw);
    tt = sec(1).theta_t; at = sec(1).alpha_t;
    v  = (nn*d*tt - (nn*d).^2*at)/c;  rg = (max(v)-min(v));
    r  = armrate(sec, rho, 'shared', L, Nt,B,fc,M,d,f,CH,SNR_t,SNR_dB,Q,'shared');
    [~,~,dbg] = ojcoms_algorithm1(Nt,fc,B,M,d,L,thlim,alim,gamma,L,Ka,s_sw);
    fprintf('%6.3f %8.2f %8.2f %7.1fns %7.2fm | %5d %7d %6d | %8.3f %4.0f%% | %8.0fx\n', ...
            s_sw, dbg.Sbound, tt, rg*1e9, rg*c, L, Ka, numel(sec), r, 100*r/r_ref, rg*1e12/DEV_PS);
  end
end
fprintf(['\nRead: the last column is how far each design point still is from the\n' ...
         'best TTD in the survey. Elements alone are not the bill -- a 32-element\n' ...
         'array whose elements each need 40 m of line is not buildable. Look for\n' ...
         'the row with >=95%% rate at the smallest "x over dev", and report the\n' ...
         'pilot cost that buys it.\n']);

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
