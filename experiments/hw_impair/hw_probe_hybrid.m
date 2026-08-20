% ========================================================================
% (B) PROPOSED FIX -- naive vs hybrid TD-PS under low-resolution TD
% ------------------------------------------------------------------------
% Compares the baseline-style DDBS beam (all delay in the TTD) against the
% proposed hybrid split (PS absorbs the frequency-independent 2*pi*fc*tau term;
% TTD carries only the frequency-dependent residual). Metric: mean realized
% focusing gain at each subcarrier's intended focus (ideal ~0.975).
%
% Expected (validated in Python at the same parameters):
%   B_td :   14      13      12      11      10
%   naive:  0.877   0.581   0.096   0.070   0.083
%   hybrid: 0.975   0.974   0.971   0.960   0.917
% => ~4.5 bits of TD resolution saved (~20x relaxation), PS stays cheap.
% ========================================================================
clear; clc;

BASELINE_DIR = fullfile('..','..','baseline_code','code_nf_distance_dependent_rainbow');
addpath(genpath(BASELINE_DIR));
addpath(genpath(fileparts(mfilename('fullpath'))));

Nt=256; fc=30e9; B=5e9; M=1024; c=3e8; d=(c/fc)/2; K=1;
theta1=-32.95; theta2=0.784; alpha1=-0.53; alpha2=0.58;
f = fc + B/M*((1:M)-1-(M-1)/2); km = 2*pi*f/c; nn = (-(Nt-1)/2:(Nt-1)/2)';

loc = cal_loc(fc, f, theta1, theta2, alpha1, alpha2, M, 1);
th = squeeze(loc(:,1,1)); al = squeeze(loc(:,2,1));
gainof = @(w) mean(arrayfun(@(m) ...
    abs( (exp(-1j*km(m)*(nn*d*th(m)-(nn*d).^2*al(m)))/sqrt(Nt))' * w(:,1,m) ), 1:M));

tau_rng = max((nn*d*theta1-(nn*d).^2*alpha1)/c) - min((nn*d*theta1-(nn*d).^2*alpha1)/c);
fprintf('delay range %.0f ns | 1/fH %.1f ps | 1/B %.0f ps\n', tau_rng*1e9, 1/(fc+B/2)*1e12, 1/B*1e12);
fprintf('ideal gain %.3f\n\n', gainof(ddbs_beam_impaired(Nt,B,fc,M,d,theta1,theta2,alpha1,alpha2,K,struct())));

fprintf('%5s %9s | %8s %8s\n','B_td','LSB[ps]','naive','HYBRID');
for Btd = [16 14 13 12 11 10 9 8]
    imp = struct('Btd',Btd);
    gn = gainof(ddbs_beam_impaired(Nt,B,fc,M,d,theta1,theta2,alpha1,alpha2,K,imp));
    gh = gainof(ddbs_beam_hybrid  (Nt,B,fc,M,d,theta1,theta2,alpha1,alpha2,K,imp));
    fprintf('%5d %9.1f | %8.3f %8.3f\n', Btd, tau_rng/2^Btd*1e12, gn, gh);
end

fprintf('\nHybrid, does the burden move to the PS?  (B_td = 10)\n');
for Bps = [Inf 6 5 4 3]
    g = gainof(ddbs_beam_hybrid(Nt,B,fc,M,d,theta1,theta2,alpha1,alpha2,K, ...
               struct('Btd',10,'Bps',Bps)));
    fprintf('  B_ps=%4s : %.3f\n', num2str(Bps), g);
end
