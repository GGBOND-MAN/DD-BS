% ========================================================================
% (B) RATE-DOMAIN ABLATION: the paper's main table/figure
% ------------------------------------------------------------------------
% The focusing-gain surface (hw_probe_surface) in the paper's own metric: average
% rate, with training AND serving beams built from ONE shared hardware model
% (same architecture, same absolute delay LSB) in every arm.
%
%   arch = 'ojcoms' : AoSA form (TTD carries full 2*pi*f_m*tau; plain DDBS PS)
%                     applied to the ORIGINAL DDBS parameters. NOT a reproduction
%                     of the OJ-COMS scheme (which redesigns parameters and uses
%                     ~12 pilots) -- never report this as beating them.
%   arch = 'shared' : same TTD sharing + hybrid TD-PS split (contribution A)
%                     -> the strongest FAIR sharing baseline.
%   arch = 'kron'   : + Kronecker two-stage TTD (contribution B).
%
% Expected (from the gain surface): at 32 TTDs and 12-bit the ablation should
% read roughly 0.06 -> 0.24 -> 0.94 of the ideal focusing gain, and the rate
% version should show the same ordering.
% ========================================================================
clear; clc;

ensure_baseline();          % locate the baseline code wherever it lives
addpath(genpath(fileparts(mfilename('fullpath'))));

Nt=256; fc=30e9; B=5e9; M=1024; c=3e8; d=(c/fc)/2; Q=1; K=3;
f = fc + B/M*((1:M)-1-(M-1)/2);
Rmin=5; Rmax=30; th_min=-pi/3; th_max=pi/3;
theta1=349/11; theta2=-30; alpha1=-427/800; alpha2=1859/3200;
focus_loc = cal_loc(fc,f,theta1,theta2,alpha1,alpha2,M,K);
nn=(-(Nt-1)/2:(Nt-1)/2)';
tau_rng = max((nn*d*theta1-(nn*d).^2*alpha1)/c) - min((nn*d*theta1-(nn*d).^2*alpha1)/c);

SNR_dB=15; SNR_t=10^(SNR_dB/10); N_iter=100; rng(7);
CH = cell(N_iter,1);
for it=1:N_iter
    CH{it} = near_field_channel(Nt,d,fc,B,M, Rmin+rand*(Rmax-Rmin), th_min+rand*(th_max-th_min));
end

% ideal reference: full per-antenna TTD, infinite resolution
w_id = ddbs_beam_arch(Nt,B,fc,M,d,theta1,theta2,alpha1,alpha2,K,'full',1,Inf,Inf,[]);
r_ideal = mean(cellfun(@(h) run1(h,w_id,'full',1,[],Nt,B,fc,f,M,d,focus_loc,SNR_t,SNR_dB,Q,K), CH));
rp = mean(cellfun(@(h) sum(arrayfun(@(m) ...
      log2(1+SNR_t*abs(h(m,:)*(exp(1j*angle(h(m,:)'))/sqrt(Nt)))^2)/M, 1:M)), CH));
fprintf('ideal full-TTD DDBS %.3f | Perfect CSI %.3f   (SNR %d dB, K=%d)\n\n', r_ideal, rp, SNR_dB, K);

Btd_list = [Inf 14 13 12 11 10];
cases = { {'ojcoms', 8}, {'shared', 8}, {'kron', 16} };   % all -> 32 TTDs
fprintf('All arms at 32 TTDs (ojcoms/shared: P=8 -> Nt/P; kron: P=16 -> Nt/P+P)\n');
fprintf('%-9s %6s |', 'arch', '#TTD'); fprintf('%9s', string(Btd_list)); fprintf('   <- B_td\n');
for ci = 1:numel(cases)
    arch = cases{ci}{1}; P = cases{ci}{2};
    [~,ntd] = ddbs_beam_arch(Nt,B,fc,M,d,theta1,theta2,alpha1,alpha2,K,arch,P,Inf,Inf,[]);
    fprintf('%-9s %6d |', arch, ntd);
    for Btd = Btd_list
        lsb = tau_rng/2^Btd;  if ~isfinite(Btd), lsb = []; end
        wtr = ddbs_beam_arch(Nt,B,fc,M,d,theta1,theta2,alpha1,alpha2,K,arch,P,Inf,Inf,lsb);
        r = mean(cellfun(@(h) run1(h,wtr,arch,P,lsb,Nt,B,fc,f,M,d,focus_loc,SNR_t,SNR_dB,Q,K), CH));
        fprintf('%9.3f', r);
    end
    fprintf('\n');
end
fprintf('\n(ideal = %.3f; divide by it for the fraction-of-ideal ablation)\n', r_ideal);

% ------------------------------------------------------------------------
function r = run1(h,w_train,arch,P,lsb,Nt,B,fc,f,M,d,focus_loc,SNR_t,SNR_dB,Q,K)
% DDBS training (argmax subcarrier -> its intended focus), then serve with a beam
% built from the SAME architecture and the SAME absolute delay LSB.
best = -inf;
for idx = 1:K
    y = zeros(M,Q);
    for m = 1:M
        y(m,:) = awgn( repmat( h(m,:)*w_train(:,idx,m), [1,Q] ), SNR_dB*2/sqrt(3) );
    end
    [~,i] = max(abs(sum(y,2)).^2);
    % serving beam = same generator with a single focus (theta2=alpha2=0, K=1)
    ws = ddbs_beam_arch(Nt,B,fc,M,d, focus_loc(i,1,idx), 0, focus_loc(i,2,idx), 0, 1, arch,P,Inf,Inf,lsb);
    t = 0;
    for m = 1:M
        t = t + log2(1 + SNR_t*abs(h(m,:)*ws(:,1,m))^2)/M;
    end
    best = max(best, t);
end
r = best;
end
