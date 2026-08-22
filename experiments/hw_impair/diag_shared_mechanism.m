% ========================================================================
% DIAGNOSTIC: what actually breaks under shared TTD?
% ------------------------------------------------------------------------
% Measures the two candidate mechanisms DIRECTLY, instead of inferring them from
% the rate, so the K_min = 1.12*P law gets a mechanism rather than a fit.
%
%   (M1) GAIN TRUNCATION -- the sub-array factor kills the band edges. Predicted
%        usable half-bandwidth (3 dB): |f-fc| <= 0.443 * 2*fc/(P*|theta_t|).
%        Measured here as the actual half-width of G(m) = realized gain at the
%        recalibrated focus.
%   (M2) COVERAGE LOSS -- after truncation the surviving subcarriers no longer
%        tile the user region, so users fall in gaps. Measured as the span and
%        the largest interior gap of the recalibrated focus angles, counting only
%        subcarriers that survive (M1).
%
% If measured usable bandwidth tracks the 2*fc/(P*|theta_t|) prediction and the
% surviving coverage collapses at exactly the (P,K) pairs that fail on rate, the
% coverage-truncation mechanism is confirmed and (*) in decouple_theta.m is the
% design equation.
% ========================================================================
clear; clc;

ensure_baseline();
addpath(genpath(fileparts(mfilename('fullpath'))));

Nt=256; fc=30e9; B=5e9; M=1024; c=3e8; d=(c/fc)/2;
f = fc + B/M*((1:M)-1-(M-1)/2); km=2*pi*f/c;
nn=(-(Nt-1)/2:(Nt-1)/2)'; N1=nn*d; N2=(nn*d).^2;
TH1_0=349/11; TH2_0=-30; AL1=-427/800; AL2=1859/3200; K0=3;

fprintf('%4s %4s %8s | %9s %9s | %7s %8s %8s\n', ...
        'P','K','P|th_t|','BW_meas','BW_pred','span','maxgap','BW/2/Nt');
for P = [1 2 4 8 16]
  for K = [3 6 9 12 18]
    if P>1 && K > max(3,ceil(1.5*P)), continue; end
    s = K0/K; TH1 = TH1_0*s; TH2 = TH2_0*s;
    focus_loc = cal_loc(fc,f,TH1,TH2,AL1,AL2,M,K);
    arch='full'; if P>1, arch='shared'; end
    w  = ddbs_beam_arch(Nt,B,fc,M,d,TH1,TH2,AL1,AL2,K,arch,P,Inf,Inf,[]);
    fl = actual_focus(w,focus_loc,Nt,fc,B,M,d,K);

    % (M1) realized gain at the recalibrated focus, per subcarrier
    G = zeros(M,K);
    for t=1:K
        for m=1:M
            a = exp(-1j*km(m)*(N1*fl(m,1,t) - N2*fl(m,2,t)))/sqrt(Nt);
            G(m,t) = abs(a'*w(:,t,m));
        end
    end
    g = mean(G,2); keep = g >= 0.707*max(g);
    BW_meas = (sum(keep)/M)*B;
    BW_pred = 2*0.443*2*fc/(P*abs(TH1));       % full usable width, 3 dB
    if P==1, BW_pred = B; end

    % (M2) coverage of the SURVIVING subcarriers
    sn = [];
    for t=1:K, sn = [sn; fl(keep,1,t)]; end %#ok<AGROW>
    sn = sort(sn(:)); sn = sn(abs(sn)<=0.9);
    span = (max(sn)-min(sn))/1.732;
    inner = sn(sn>=-0.866 & sn<=0.866);
    maxgap = max(diff(inner));

    fprintf('%4d %4d %8.1f | %6.2fGHz %6.2fGHz | %7.2f %8.4f %8.1f\n', ...
            P,K,P*abs(TH1), BW_meas/1e9, min(BW_pred,B)/1e9, span, maxgap, maxgap/(2/Nt));
  end
end
fprintf(['\nRead: BW_meas vs BW_pred tests (M1). maxgap/(2/Nt) is the largest\n' ...
         'coverage hole in beamwidths -- >1 means users fall between strips (M2).\n']);
