function p = derisk_params()
% DERISK_PARAMS  Central system parameters for the multipath de-risk study.
%
% Numbers copied verbatim from the baseline driver Rate_snr.m so that the
% single-path (L=1) behaviour of these scripts reproduces the paper exactly.
% Baseline: T. Zheng, M. Cui, Z. Wu, L. Dai, "Near-field wideband beam
% training based on distance-dependent beam split," IEEE TWC 24(2), 2025.

c        = 3e8;
p.c      = c;
p.Nt     = 256;            % # BS antennas (ULA)
p.fc     = 30e9;           % carrier frequency
p.B      = 5e9;            % bandwidth
p.M      = 1024;           % # subcarriers
p.Q      = 1;              % # repeated symbols per pilot (baseline uses 1)
p.lambda = c / p.fc;
p.d      = p.lambda / 2;   % antenna spacing

% subcarrier frequencies (same indexing as baseline)
p.f = zeros(1, p.M);
for m = 1:p.M
    p.f(m) = p.fc + p.B/p.M * (m - 1 - (p.M-1)/2);
end

% user region
p.Rmin      = 5;
p.Rmax      = 30;
p.theta_min = -pi/3;
p.theta_max =  pi/3;

% ---- DDBS training codebook parameters (exact values from Rate_snr.m) ----
p.k      = 3;              % number of pilots K
p.theta1 = 349/11;   p.theta2 = -30;
p.alpha1 = -427/800; p.alpha2 = 1859/3200;
end
