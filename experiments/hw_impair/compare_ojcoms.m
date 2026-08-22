% ========================================================================
% HEAD-TO-HEAD vs OJ-COMS AoSA -- at EQUAL TTD count and EQUAL pilot budget
% ------------------------------------------------------------------------
% The first thing a reviewer will ask. OJ-COMS 2026 (Qaid, Nasir, Al-Ahmadi, Liu)
% reduces the TTD COUNT with an AoSA/grouped-TTD array; this work characterises
% what that reduction costs in pilots and how to pay it. A comparison is only
% honest if both sides get the same hardware AND the same overhead, so every cell
% below fixes (N_TTD, K) and varies only the architecture and the decision rule.
%
% FIVE ARMS
%   A  full        one TTD per antenna, designed lookup        -- upper bound
%   B  shared      grouped TTD, DD-BS's ORIGINAL decision rule -- naive port
%   C  ojcoms      their eqs (13)-(14) form, designed lookup   -- their assumption
%   D  ojcoms      their architecture + OUR recalibrated lookup
%   E  shared      grouped TTD + OUR recalibrated lookup       -- ours
%
% Arm D matters as much as arm E. If recalibration lifts THEIR architecture too,
% the contribution is architecture-agnostic and must be reported that way -- that
% is a stronger and more citable claim than "our box beats their box", and hiding
% it would be the kind of thing a reviewer finds. If D stays low while E rises,
% the gain is specific to compensating the intra-group delay at fc.
%
% WHAT THIS IS NOT. 'ojcoms' reproduces their ARCHITECTURE (TTD carries the full
% 2*pi*f_m*tau on the sub-array-centre grid; the PS is the plain per-antenna DDBS
% term and does not compensate the intra-group difference). Their full scheme also
% redesigns the DDBS parameters and adds sector-shift / distance-interleaving /
% extra pilots. So arm C is NOT their published result and must never be labelled
% as such -- it is the architecture evaluated under this paper's protocol.
%
% Expected shape, from the master curve (rate depends on Delta*P/2 = P/K):
%   P/K > ~1.14  -> every shared arm fails, ours included: the pilots are unpaid
%   P/K <= ~1.14 -> E reaches the full-TTD bound; B stays down (wrong lookup)
%
% Runtime: 6 cells x 2 focus tables; the K=16 cell dominates.
% ========================================================================
clear; clc;

ensure_baseline();
addpath(genpath(fileparts(mfilename('fullpath'))));

Nt=256; fc=30e9; B=5e9; M=1024; c=3e8; d=(c/fc)/2; Q=1;
f = fc + B/M*((1:M)-1-(M-1)/2);
Rmin=5; Rmax=30; th_min=-pi/3; th_max=pi/3;
TH1=349/11; TH2=-30; AL1=-427/800; AL2=1859/3200;
SNR_dB=15; SNR_t=10^(SNR_dB/10); N_iter=30;

rng(23); CH=cell(N_iter,1);
for it=1:N_iter
    CH{it}=near_field_channel(Nt,d,fc,B,M, Rmin+rand*(Rmax-Rmin), th_min+rand*(th_max-th_min));
end

fprintf('Head-to-head at equal (N_TTD, K).  Nt=%d, SNR=%d dB, N_iter=%d\n', Nt, SNR_dB, N_iter);
fprintf('Rates in bit/s/Hz; the percentage is against arm A at the same K.\n\n');
fprintf('%6s %4s %6s | %8s | %13s %13s | %13s %13s\n', ...
        'N_TTD','K','P/K','A full','B shared','C ojcoms','D ojcoms','E shared');
fprintf('%6s %4s %6s | %8s | %13s %13s | %13s %13s\n', ...
        '','','','designed','designed','designed','RECALIBR.','RECALIBR.');

% Cells straddle the predicted boundary P/K ~ 1.14 at each P rather than filling a
% grid -- actual_focus cost grows with K, and off-boundary cells add runtime
% without adding evidence.
cells = [4 3; 4 4; 8 3; 8 8; 16 8; 16 16];
rA_k = containers.Map('KeyType','double','ValueType','double');

for ci = 1:size(cells,1)
    P = cells(ci,1); K = cells(ci,2);
    if ~isKey(rA_k, K)
        rA_k(K) = evalarm('full', 1, K, false, Nt,B,fc,M,d,TH1,TH2,AL1,AL2,f,CH,SNR_t,SNR_dB,Q);
        fprintf('%6d %4d %6s | %8.3f | %13s %13s | %13s %13s\n', ...
                Nt, K, '-', rA_k(K), '-','-','-','-');
    end
    rA = rA_k(K);
    rB = evalarm('shared', P, K, false, Nt,B,fc,M,d,TH1,TH2,AL1,AL2,f,CH,SNR_t,SNR_dB,Q);
    rC = evalarm('ojcoms', P, K, false, Nt,B,fc,M,d,TH1,TH2,AL1,AL2,f,CH,SNR_t,SNR_dB,Q);
    rD = evalarm('ojcoms', P, K, true,  Nt,B,fc,M,d,TH1,TH2,AL1,AL2,f,CH,SNR_t,SNR_dB,Q);
    rE = evalarm('shared', P, K, true,  Nt,B,fc,M,d,TH1,TH2,AL1,AL2,f,CH,SNR_t,SNR_dB,Q);
    fprintf('%6d %4d %6.2f |          | %7.3f %4.0f%% %7.3f %4.0f%% | %7.3f %4.0f%% %7.3f %4.0f%%\n', ...
            Nt/P, K, P/K, rB,100*rB/rA, rC,100*rC/rA, rD,100*rD/rA, rE,100*rE/rA);
end

fprintf(['\nRead ACROSS a row at fixed (N_TTD, K):\n' ...
         '  E - B   isolates the DECISION RULE   (same hardware, our lookup vs theirs)\n' ...
         '  E - D   isolates the ARCHITECTURE    (same lookup, our sharing vs theirs)\n' ...
         '  D - C   says whether our lookup is architecture-agnostic. If it lifts\n' ...
         '          their array too, report it that way -- it is the stronger claim.\n' ...
         'P/K <= ~1.14 is the predicted pass side. Arm C is the OJ-COMS ARCHITECTURE\n' ...
         'under this protocol, NOT their published result.\n']);

% ------------------------------------------------------------------------
function r = evalarm(arch, P, K, recal, Nt,B,fc,M,d,TH1,TH2,AL1,AL2,f,CH,SNR_t,SNR_dB,Q)
focus_loc = cal_loc(fc,f,TH1,TH2,AL1,AL2,M,K);
w  = ddbs_beam_arch(Nt,B,fc,M,d,TH1,TH2,AL1,AL2,K,arch,P,Inf,Inf,[]);
if recal, tbl = actual_focus(w,focus_loc,Nt,fc,B,M,d,K); else, tbl = focus_loc; end
r  = mean(cellfun(@(h) rate_ongrid(h,w,tbl,Nt,B,fc,M,d,SNR_t,SNR_dB,Q,K), CH));
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
