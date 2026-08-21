% ========================================================================
% LIMITED-TTD DD-BS: does an ALGORITHM fix recover the loss? (no HW change)
% ------------------------------------------------------------------------
% Frame (the user's original one):
%   Q1. Under a hard budget of N_TTD << Nt, does DD-BS still work?   -> NO
%   Q2. If the original architecture/algorithm breaks, design one that
%       recovers near-ideal performance.                             <- THIS
%
% Key observation this probe tests. Sharing TTDs DISTORTS the training beams,
% but the distortion is DETERMINISTIC and exactly computable. The baseline
% decision rule assumes subcarrier m focuses at its DESIGNED location
% focus_loc(m); once TTDs are shared that premise is false, so the rule is
% mismatched -- the information is not destroyed, it is just not where the rule
% looks. Rebuilding the match-filter dictionary from the ACTUAL (impaired)
% training beams should therefore recover much of the loss, with NO hardware
% change and NO extra pilots.
%
% Arms, all at the same N_TTD:
%   (a)  on-grid, DESIGNED locations  = the baseline rule            [mismatched]
%   (a2) on-grid, RECALIBRATED locations -- where each subcarrier's beam
%        ACTUALLY focuses under the shared-TTD hardware (offline table)  [proposed]
%   (b)  match filter, dictionary from IDEAL beams                   [still mismatched]
%   (c)  match filter, dictionary from the ACTUAL impaired beams     [proposed]
% Reference: ideal full-TTD DD-BS.
%
% (a2) is the cleanest test of the mismatch hypothesis: same rule, same cost,
% only the lookup table is corrected -- no dictionary, hence no grid-resolution
% confound (the MF arms are limited by g1/g2 and are NOT comparable to (a) on
% resolution grounds).
%
% Cross-checked at reduced scale (Nt=256, M=256, ideal ref 4.086):
%   P    #TTD   (a) designed   (a2) recalibrated
%   2     128      3.510            4.029    <- essentially restores ideal
%   4      64      2.655            2.968
%   8      32      1.310            1.964    <- 1.50x
%  16      16      0.625            0.691
% => the limited-TTD loss splits into a MISMATCH part, recoverable for free by
%    recalibrating the lookup, and a genuine BEAM-DEGRADATION part that no
%    decision rule can recover and that needs an architecture-level fix.
% ========================================================================
clear; clc;

BASELINE_DIR = fullfile('..','..','baseline_code','code_nf_distance_dependent_rainbow');
addpath(genpath(BASELINE_DIR));
addpath(genpath(fileparts(mfilename('fullpath'))));

Nt=256; fc=30e9; B=5e9; M=1024; c=3e8; d=(c/fc)/2; Q=1; K=3;
f = fc + B/M*((1:M)-1-(M-1)/2);
Rmin=5; Rmax=30; th_min=-pi/3; th_max=pi/3;
theta1=349/11; theta2=-30; alpha1=-427/800; alpha2=1859/3200;
focus_loc = cal_loc(fc,f,theta1,theta2,alpha1,alpha2,M,K);

g1=256; g1_list=linspace(sin(th_min),sin(th_max),g1);
g2=10;  g2_list=linspace(1/400,1/10,g2);

SNR_dB=15; SNR_t=10^(SNR_dB/10); N_iter=100; rng(17);
CH=cell(N_iter,1);
for it=1:N_iter
    CH{it}=near_field_channel(Nt,d,fc,B,M, Rmin+rand*(Rmax-Rmin), th_min+rand*(th_max-th_min));
end

% ideal reference (full per-antenna TTD)
w_ideal = ddbs_beam_arch(Nt,B,fc,M,d,theta1,theta2,alpha1,alpha2,K,'full',1,Inf,Inf,[]);
D_ideal = mkdict(Nt,B,fc,M,d,w_ideal,K,Q,g1_list,g2_list);
r_ref = mean(cellfun(@(h) rate_mf(h,w_ideal,D_ideal,g1_list,g2_list,Nt,B,fc,M,d,SNR_t,SNR_dB,Q,K), CH));
fprintf('ideal full-TTD DD-BS (match filter) = %.3f bit/s/Hz\n\n', r_ref);

fprintf('%6s %7s | %11s %14s %11s %11s\n','P','#TTD','(a) designed','(a2) recalib.','(b) MF ideal','(c) MF ACTUAL');
for P = [2 4 8]
    % actual impaired training beams at this TTD budget (hybrid PS applied to all arms)
    [w,ntd] = ddbs_beam_arch(Nt,B,fc,M,d,theta1,theta2,alpha1,alpha2,K,'shared',P,Inf,Inf,[]);
    D_act = mkdict(Nt,B,fc,M,d,w,K,Q,g1_list,g2_list);      % dictionary from ACTUAL beams

    ra  = mean(cellfun(@(h) rate_ongrid(h,w,focus_loc,Nt,B,fc,M,d,SNR_t,SNR_dB,Q,K), CH));
    fl2 = actual_focus(w,focus_loc,Nt,fc,B,M,d,K);          % recalibrated lookup
    ra2 = mean(cellfun(@(h) rate_ongrid(h,w,fl2,Nt,B,fc,M,d,SNR_t,SNR_dB,Q,K), CH));
    rb = mean(cellfun(@(h) rate_mf(h,w,D_ideal,g1_list,g2_list,Nt,B,fc,M,d,SNR_t,SNR_dB,Q,K), CH));
    rc = mean(cellfun(@(h) rate_mf(h,w,D_act,  g1_list,g2_list,Nt,B,fc,M,d,SNR_t,SNR_dB,Q,K), CH));
    fprintf('%6d %7d | %11.3f %14.3f %11.3f %11.3f\n', P, ntd, ra, ra2, rb, rc);
end
fprintf('\n(c) uses the SAME hardware and the SAME 3 pilots as (a).\n');

% ------------------------------------------------------------------------
function fl = actual_focus(w,focus_loc,Nt,fc,B,M,d,K)
% Offline table of where each subcarrier's beam ACTUALLY focuses under the given
% (impaired) hardware. Two 1-D searches per (m,t): theta at the designed alpha,
% then alpha at that theta. NOTE the correlation uses w WITHOUT conjugation --
% self-tested against the ideal beams, where it reproduces the designed table.
c=3e8; f = fc + B/M*((1:M)-1-(M-1)/2); km = 2*pi*f/c;
nn=(-(Nt-1)/2:(Nt-1)/2)'; N1=nn*d; N2=(nn*d).^2;
GS = linspace(-0.9,0.9,801); GA = linspace(1/400,1/10,81);
fl = zeros(M,2,K);
for t=1:K
    for m=1:M
        wm = w(:,t,m); a0 = focus_loc(m,2,t);
        g  = abs(exp(1j*km(m)*(GS.'*N1.' - a0*N2.'))*wm);
        [~,i] = max(g); th = GS(i);
        g2 = abs(exp(1j*km(m)*(th*N1.' - GA.'*N2.'))*wm);
        [~,j] = max(g2);
        fl(m,1,t) = th; fl(m,2,t) = GA(j);
    end
end
end

function D = mkdict(Nt,B,fc,M,d,w,K,~,g1_list,g2_list)
% Noiseless received-power fingerprint of every grid point under the GIVEN beams.
% Returned FLAT as (g1*g2) x (K*M) so scoring is one matrix-vector product.
% Column-major over (i,j) so ind2sub([g1 g2], r) recovers the grid indices.
g1=numel(g1_list); g2=numel(g2_list);
Wt = zeros(Nt, M, K);                       % pre-slice the training beams
for t=1:K, Wt(:,:,t) = squeeze(w(:,t,:)); end
D = zeros(g1*g2, K*M);
parfor r = 1:g1*g2
    [i,j] = ind2sub([g1 g2], r);
    wg = TTD_beam(Nt,B,fc,M,d,g1_list(i),g2_list(j));      % Nt x M
    row = zeros(1, K*M);
    for t=1:K
        y = sum(conj(wg) .* Wt(:,:,t), 1);                 % 1 x M, vectorised over m
        row((t-1)*M+(1:M)) = abs(y).^2;
    end
    D(r,:) = row;
end
end

function r = rate_ongrid(h,w,focus_loc,Nt,B,fc,M,d,SNR_t,SNR_dB,Q,K)
best=-inf;
for idx=1:K
    y=zeros(M,Q);
    for m=1:M, y(m,:)=awgn(repmat(h(m,:)*w(:,idx,m),[1,Q]),SNR_dB*2/sqrt(3)); end
    [~,i]=max(abs(sum(y,2)).^2);
    best=max(best, serve(h,focus_loc(i,1,idx),focus_loc(i,2,idx),Nt,B,fc,M,d,SNR_t));
end
r=best;
end

function r = rate_mf(h,w,D,g1_list,g2_list,Nt,B,fc,M,d,SNR_t,SNR_dB,Q,K)
% Score every grid point with one matrix-vector product against the flat dictionary.
yv = zeros(K*M,1);
for t=1:K
    y=zeros(M,Q);
    for m=1:M, y(m,:)=awgn(repmat(h(m,:)*w(:,t,m),[1,Q]),SNR_dB*2/sqrt(3)); end
    yv((t-1)*M+(1:M)) = abs(sum(y,2)).^2;
end
[~,rmax] = max(D*yv);
[i,j] = ind2sub([numel(g1_list) numel(g2_list)], rmax);
r = serve(h,g1_list(i),g2_list(j),Nt,B,fc,M,d,SNR_t);
end

function t = serve(h,sin_th,al,Nt,B,fc,M,d,SNR_t)
ws = TTD_beam(Nt,B,fc,M,d,sin_th,al);
t=0; for m=1:M, t=t+log2(1+SNR_t*abs(h(m,:)*ws(:,m))^2)/M; end
end
