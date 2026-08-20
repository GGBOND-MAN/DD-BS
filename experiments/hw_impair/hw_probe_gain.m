% ========================================================================
% (B) Hardware-impairment robustness of DDBS -- PROBE 1: focusing gain
% ------------------------------------------------------------------------
% Improve-on-existing target: OJ-COMS 2026 (AoSA / grouped-TTD) reduced the TD
% COUNT but assumed IDEAL TD/PS values. The untouched axis is TD/PS RESOLUTION &
% error. This probe maps how the DDBS multi-strip focusing gain degrades under
% finite-bit TD/PS, TD fabrication jitter, and phase noise.
%
% Metric: mean realized array gain at each subcarrier's INTENDED focus (ideal ~0.97).
%
% Pre-validated finding (full Nt=256, M=1024) -- run to reproduce:
%   * TD bits (quantised over the 140 ns delay range): CATASTROPHIC. 8-bit -> 0.07.
%     DDBS needs ~13-bit adjustable TTD; this is WHY the baseline uses fixed lines.
%   * TD jitter (hits fixed lines too): tight -- <1 dB needs sigma_tau <~2 ps;
%     5 ps -> 0.63, 10 ps -> 0.18.
%   * PS bits: ROBUST -- 2-bit -> 0.88, >=3-bit -> >=0.95.
%   * Phase noise: moderate -- 20 deg -> 0.92, 30 deg -> 0.85.
% => The fragility is ASYMMETRIC and concentrated in TD precision. That is the
%    contribution: make DDBS tolerate low-resolution TD (delay-range reduction /
%    coarse-PS + fine-TTD / sub-array TTD), and quantify the bit/precision saving.
% ========================================================================
clear; clc;

BASELINE_DIR = fullfile('..','..','baseline_code','code_nf_distance_dependent_rainbow');
addpath(genpath(BASELINE_DIR));
addpath(genpath(fileparts(mfilename('fullpath'))));

Nt=256; fc=30e9; B=5e9; M=1024; c=3e8; d=(c/fc)/2;
% one representative DDBS pilot (Rate_snr beam-training section)
theta1=-32.95; theta2=0.784; alpha1=-0.53; alpha2=0.58; K=1;
f = fc + B/M*((1:M)-1-(M-1)/2); km = 2*pi*f/c;
nn = (-(Nt-1)/2:(Nt-1)/2)';
nrep = 8;   % Monte-Carlo repeats for random impairments

% intended per-subcarrier focus (baseline cal_loc)
loc = cal_loc(fc, f, theta1, theta2, alpha1, alpha2, M, 1);
th = squeeze(loc(:,1,1)); al = squeeze(loc(:,2,1));

meanfocusgain = @(w) mean(arrayfun(@(m) ...
    abs( (exp(-1j*km(m)*(nn*d*th(m)-(nn*d).^2*al(m)))/sqrt(Nt))' * w(:,1,m) ), 1:M));

sweep = @(imp) mean(arrayfun(@(~) meanfocusgain( ...
    ddbs_beam_impaired(Nt,B,fc,M,d,theta1,theta2,alpha1,alpha2,K,imp)), 1:nrep));

fprintf('delay range = %.0f ns, fH period = %.1f ps, ideal gain = %.3f\n\n', ...
        (max((nn*d*theta1-(nn*d).^2*alpha1)/c)-min((nn*d*theta1-(nn*d).^2*alpha1)/c))*1e9, ...
        1/(fc+B/2)*1e12, sweep(struct()));

fprintf('-- TD bits (quantise full delay range) --\n');
for Btd = [Inf 12 10 8 6 4]
    fprintf('  B_td=%4s : gain %.3f\n', num2str(Btd), sweep(struct('Btd',Btd)));
end
fprintf('-- TD jitter sigma_tau [ps] (fabrication) --\n');
for st = [0 0.5 1 2 5 10]
    fprintf('  sigma_tau=%4.1f ps : gain %.3f\n', st, sweep(struct('sigma_tau',st*1e-12)));
end
fprintf('-- PS bits --\n');
for Bps = [Inf 6 4 3 2]
    fprintf('  B_ps=%4s : gain %.3f\n', num2str(Bps), sweep(struct('Bps',Bps)));
end
fprintf('-- phase noise sigma [deg] --\n');
for sd = [0 5 10 20 30]
    fprintf('  sigma_pn=%2d deg : gain %.3f\n', sd, sweep(struct('sigma_pn',deg2rad(sd))));
end
