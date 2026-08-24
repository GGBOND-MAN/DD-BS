% ========================================================================
% BORROWED REFINEMENT: use the phase our pipeline throws away (Luo & Gao, TWC 2024)
% ------------------------------------------------------------------------
% Our decision rule is argmax_m |sum_q y|^2 followed by a table lookup. It uses
% MAGNITUDE ONLY. Luo & Gao (IEEE TWC 23(5), May 2024) refine distance from the
% multi-carrier PHASE DIFFERENCE between max-power subcarriers, reaching ~0.10 m
% RMSE at 15 dB. That information is free -- already received, currently discarded
% -- and distance is exactly where the grouped architecture is weakest (sec 27:
% L=16 tops out at 89%).
%
% THE OBSERVABLE. For a single-path LoS user at range r,
%       y_{m,t} = g * exp(-j 2 pi f_m r / c) * [ a(theta_u, alpha_u)^T w_{:,t,m} ]
% The bracket is KNOWN OFFLINE once a coarse (theta, alpha) is in hand -- it is the
% same beam table the recalibrated lookup already needs. De-embedding it leaves
%       z_m = y_m * conj(beam term)  ~  |.| * exp(-j 2 pi f_m r / c)
% a pure ranging phase ramp, and
%       r_hat = argmax_r | sum_{m in S} z_m exp(+j 2 pi f_m r / c) |
% is a matched filter in range over the strong subcarriers S. alpha follows as
% (1 - theta^2)/(2r).
%
% WRAPPING. The per-subcarrier phase step is 2 pi r (B/M) / c, which at r = 200 m
% is 2 pi * 3.25 -- ambiguous. The search is therefore seeded by the coarse
% estimate, and unambiguous only while the coarse range error stays below
% c/(2*(B/M)) = 30.7 m. That is a real precondition, not a formality.
%
% STAGE 0 IS A GATE, and it can fail. Everything above assumes the baseline
% channel model actually carries the common propagation phase exp(-j2*pi*f_m*r/c).
% Some near-field models normalise it away, in which case the ranging observable
% does not exist and this whole idea is void -- so stage 0 fits r from the phase
% of a noiseless measurement and compares it against the truth BEFORE any rate is
% reported. If the fit does not track the true range, stop: the fallback would be
% a joint (theta, alpha) phase refinement using the beam phase alone, which is
% weaker and needs its own design.
%
% NOVELTY NOTE. The refinement is BORROWED and must be cited (Luo & Gao TWC 2024;
% the coarse-to-fine structure is TVT 2026). It is applied to BOTH arms precisely
% so it cannot be mistaken for this work's contribution -- and so that a reviewer
% cannot claim the architecture gap would close under a stronger estimator.
% ========================================================================
clear; clc;

ensure_baseline();
addpath(genpath(fileparts(mfilename('fullpath'))));

Nt=256; fc=30e9; B=5e9; M=1024; c=3e8; d=(c/fc)/2; Q=1;
f = fc + B/M*((1:M)-1-(M-1)/2); km = 2*pi*f/c;
nn=(-(Nt-1)/2:(Nt-1)/2)'; N1=nn*d; N2=(nn*d).^2;
thlim=[-1 1]; alim=[0.0025 0.1]; gamma=0.9;
Rmin=5; Rmax=200; SNR_dB=20; SNR_t=10^(SNR_dB/10); N_iter=20;
rgrid = Rmin:0.05:Rmax;                       % ranging search grid

% ---------------- STAGE 0: does the ranging phase exist at all? ----------------
fprintf('=== STAGE 0 gate: is the propagation phase present in the model? ===\n');
fprintf('ambiguity period c*M/B = %.2f m; in-lobe resolution c/(2B) = %.3f m\n', c*M/B, c/(2*B));
fprintf('search is seeded and confined to one period -- unseeded search picks the\nwrong lobe, which is what the first run of this gate demonstrated.\n');
rng(7); ok = 0;
for t = 1:5
    r_true = Rmin + rand*(Rmax-Rmin);  th_true = asin(-0.9+1.8*rand);
    [h,~] = near_field_channel(Nt,d,fc,B,M,r_true,th_true);
    sn = sin(th_true); al = (1-sn^2)/(2*r_true);
    w  = ddbs_beam_arch(Nt,B,fc,M,d,349/11,-30,-427/800,1859/3200,1,'full',1,Inf,Inf,[]);
    y = zeros(M,1); bt = zeros(M,1);
    for m=1:M
        y(m)  = h(m,:)*w(:,1,m);
        bt(m) = exp(1j*km(m)*(N1.'*sn - N2.'*al))*w(:,1,m);
    end
    r_hat = range_mf(y, bt, f, rgrid, c, r_true);   % seeded: existence test only
    fprintf('  r_true = %6.2f m -> r_hat = %6.2f m   (err %+6.2f m)\n', r_true, r_hat, r_hat-r_true);
    ok = ok + (abs(r_hat-r_true) < 1);
end
if ok < 4
    fprintf(['\nGATE FAILED (%d/5 within 1 m). The model does not carry a usable\n' ...
             'ranging phase, so this refinement is void. STOP -- do not report rates.\n'], ok);
    return;
end
fprintf('GATE PASSED (%d/5 within 1 m). Proceeding.\n\n', ok);

% ---------------- STAGE 1: refinement on BOTH arms ----------------
rng(23); CH=cell(N_iter,1); TRUE=zeros(N_iter,2);
for it=1:N_iter
    r = Rmin + rand*(Rmax-Rmin); th = asin(-0.9+1.8*rand);
    CH{it} = near_field_channel(Nt,d,fc,B,M,r,th);  TRUE(it,:) = [r sin(th)];
end
[sec0, rho0] = ojcoms_algorithm1(Nt,fc,B,M,d,1,thlim,alim,gamma,1,12);
r_ref = armrate(sec0, rho0, 'full', 1, false, Nt,B,fc,M,d,f,km,N1,N2,CH,SNR_t,SNR_dB,Q,rgrid,c);
fprintf('ungrouped reference (L=1, 12 pilots) = %.3f\n\n', r_ref);
fprintf('%3s %7s %6s %5s | %8s %8s %6s | %8s %8s %6s\n', ...
        'L','N_TTD','Nsec','K', 'them','them+ph','gain','ours','ours+ph','gain');
cfg = { 4,4,3; 8,4,3; 8,8,1; 16,4,3; 16,12,1 };
for i = 1:size(cfg,1)
    L = cfg{i,1}; Nsec = cfg{i,2}; Ka = cfg{i,3};
    [sec, rho] = ojcoms_algorithm1(Nt,fc,B,M,d,L,thlim,alim,gamma,Nsec,Ka);
    a1 = armrate(sec,rho,'ojcoms',L,false, Nt,B,fc,M,d,f,km,N1,N2,CH,SNR_t,SNR_dB,Q,rgrid,c);
    a2 = armrate(sec,rho,'ojcoms',L,true , Nt,B,fc,M,d,f,km,N1,N2,CH,SNR_t,SNR_dB,Q,rgrid,c);
    b1 = armrate(sec,rho,'shared',L,false, Nt,B,fc,M,d,f,km,N1,N2,CH,SNR_t,SNR_dB,Q,rgrid,c);
    b2 = armrate(sec,rho,'shared',L,true , Nt,B,fc,M,d,f,km,N1,N2,CH,SNR_t,SNR_dB,Q,rgrid,c);
    fprintf('%3d %7d %6d %5d | %8.3f %8.3f %5.2fx | %8.3f %8.3f %5.2fx\n', ...
            L, Nt/L, Nsec, numel(sec), a1,a2,a2/max(a1,eps), b1,b2,b2/max(b1,eps));
end
fprintf(['\nRead: the refinement is applied identically to both arms. If the\n' ...
         '"ours" column still leads after both are refined, the architecture gap\n' ...
         'is not an artifact of a weak estimator -- which is the point of running\n' ...
         'it on both. A large "gain" on the ours column at L=16 would also close\n' ...
         'the distance-domain shortfall of sec 27.\n']);

% ------------------------------------------------------------------------
function [rhat, namb] = range_mf(y, bt, f, rgrid, c, r_seed)
% Matched filter in range over the strong subcarriers, after de-embedding the
% known beam term -- the phase our magnitude-only argmax discards.
%
% AMBIGUITY, and the bug this call signature exists to prevent. With
% f_m = fc + df*(m-m0) the matched filter is periodic in r with period c/df =
% c*M/B (61.44 m here), so searching the whole [Rmin, Rmax] picks a neighbouring
% lobe whenever the true range is not the strongest one. The first version of
% this file did exactly that, and its gate failed with errors of -61.44, +61.28
% and +122.72 m -- integer multiples of the period, which is what proved the
% observable was present rather than absent. The search is therefore confined to
% one period around a seed, as the file header always said it must be.
%
% Within a period the resolution is set by the bandwidth, c/(2B) = 3 cm here,
% which is why correctly-lobed estimates come back exact to 1 cm.
z = y .* conj(bt);
[~, ord] = sort(abs(bt), 'descend');
S = ord(1:max(8, round(numel(ord)/4)));          % strong-gain subcarriers only
half = c/(2*(f(2)-f(1)));                        % half the ambiguity period
win  = rgrid(rgrid >= r_seed-half & rgrid <= r_seed+half);
if isempty(win), rhat = r_seed; namb = 0; return; end
E = exp(1j*2*pi*(f(S).')*(win/c));
[~, i] = max(abs(z(S).' * E));
rhat = win(i);
namb = round((rhat - r_seed)/(2*half));          % lobes away from the seed (0 = in window)
end

function r = armrate(sec, rho, arch, L, refine, Nt,B,fc,M,d,f,km,N1,N2,CH,SNR_t,SNR_dB,Q,rgrid,c)
K = numel(sec); kc = 2*pi*fc/c;
w = zeros(Nt,K,M); cf = zeros(M,2,K);
for s = 1:K
    w(:,s,:) = ddbs_beam_arch(Nt,B,fc,M,d, sec(s).theta_t, sec(s).theta_p, ...
                              sec(s).alpha_t, sec(s).alpha_p, 1, arch, L, Inf, Inf, []);
    cf(:,1,s) = (kc./km).'*sec(s).theta_p + rho.theta*sec(s).theta_t;
    cf(:,2,s) = (kc./km).'*sec(s).alpha_p + rho.alpha*sec(s).alpha_t;
end
fl = actual_focus(w, cf, Nt, fc, B, M, d, K);
r = mean(cellfun(@(h) onechan(h,w,fl,refine,Nt,B,fc,M,d,f,km,N1,N2,SNR_t,SNR_dB,Q,K,rgrid,c), CH));
end

function best = onechan(h,w,fl,refine,Nt,B,fc,M,d,f,km,N1,N2,SNR_t,SNR_dB,Q,K,rgrid,c)
best=-inf;
for t=1:K
    y=zeros(M,Q);
    for m=1:M, y(m,:)=repmat(h(m,:)*w(:,t,m),[1,Q]); end
    y = add_awgn(y, SNR_dB*2/sqrt(3));
    ys = sum(y,2);
    [~,i]=max(abs(ys).^2);
    th = fl(i,1,t);  al = fl(i,2,t);
    if refine
        % de-embed the beam at the COARSE estimate, then range-match on phase
        bt = zeros(M,1);
        for m=1:M, bt(m) = exp(1j*km(m)*(N1.'*th - N2.'*al))*w(:,t,m); end
        r_coarse = (1-th^2)/(2*max(al,eps));      % seed from the coarse table entry
        rh = range_mf(ys, bt, f, rgrid, c, r_coarse);
        al = (1-th^2)/(2*rh);                     % alpha from the refined range
    end
    ws = TTD_beam(Nt,B,fc,M,d,th,al);
    tt=0; for m=1:M, tt=tt+log2(1+SNR_t*abs(h(m,:)*ws(:,m))^2)/M; end
    best=max(best,tt);
end
end
