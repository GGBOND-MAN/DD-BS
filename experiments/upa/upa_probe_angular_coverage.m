% ========================================================================
% UPA / 3D-DDBS feasibility probe -- SUPPORT: 2D angular coverage locus
% ------------------------------------------------------------------------
% On a ULA, sweeping frequency sweeps ONE angle -> beam-split tiles the 1D angle
% line. On a UPA a separable design driven by the SAME frequency axis makes the
% focus trace a 1D CURVE (Omega_x(m), Omega_y(m)) through the 2D angle plane --
% NOT an automatic 2D fill. Covering the 2D plane therefore needs (i) different
% per-axis sweep rates so the curve is space-filling (Lissajous-like), and
% (ii) a few interleaved pilots (shifted intercepts). This probe visualizes that
% and gives a rough coverage fraction -- i.e. the UPA pilot-overhead intuition,
% which is the whole selling point of a beam-split scheme.
%
% Reuses the baseline cal_loc.m (per-axis focus locus Omega(m) = theta_t +
% (fc/f_m)(theta_p + 2p), wrapped to [-1,1]).
% ========================================================================
clear; clc;

BASELINE_DIR = fullfile('..','..','baseline_code','code_nf_distance_dependent_rainbow');
addpath(genpath(BASELINE_DIR));
addpath(genpath(fileparts(mfilename('fullpath'))));

fc = 30e9; B = 5e9; M = 1024; c = 3e8;
f = fc + B/M*((1:M) - 1 - (M-1)/2);

% per-axis DDBS angle params (theta_t, theta_p); different rates -> space-filling
xp = struct('tt',-30.0,'tp', 1.90);   % x-axis: fast sweep
yp = struct('tt',-11.0,'tp', 0.63);   % y-axis: ~1/3 the rate  (Lissajous)

Kx = 3; Ky = 3;                        % interleaved pilots per axis (intercept shifts)

% cal_loc(fc,f,theta1,theta2,alpha1,alpha2,M,k) -> loc(m,1,k) is the angle locus
locx = cal_loc(fc, f, xp.tt, xp.tp, 0, 0, M, Kx);   % M x 2 x Kx
locy = cal_loc(fc, f, yp.tt, yp.tp, 0, 0, M, Ky);

Wx = squeeze(locx(:,1,:));             % M x Kx
Wy = squeeze(locy(:,1,:));             % M x Ky

% pair every x-pilot with every y-pilot -> Kx*Ky pilots total
smax = sin(pi/3);
figure; hold on; box on;
allx = []; ally = [];
for a = 1:Kx
    for b = 1:Ky
        plot(Wx(:,a), Wy(:,b), '.', 'MarkerSize', 4);
        allx = [allx; Wx(:,a)]; ally = [ally; Wy(:,b)]; %#ok<AGROW>
    end
end
axis equal; xlim([-1 1]); ylim([-1 1]);
rectangle('Position',[-smax -smax 2*smax 2*smax],'EdgeColor',[.4 .4 .4],'LineStyle','--');
xlabel('\Omega_x'); ylabel('\Omega_y');
title(sprintf('UPA angular coverage locus  (K = %d pilots)', Kx*Ky));

% rough coverage: fraction of a 2D angle grid within a beamwidth of the locus set
Ng = 60; gx = linspace(-smax,smax,Ng); gy = linspace(-smax,smax,Ng);
Nt_lin = 64;                              % per-axis element count (match upa_probe_coupling)
bw = 0.886/Nt_lin;                        % ~3dB beamwidth in sine units (per axis)
covered = 0; total = 0;
for i = 1:Ng
    for j = 1:Ng
        if gx(i)^2+gy(j)^2 > smax^2, continue; end
        total = total + 1;
        dmin = min( max(abs(allx-gx(i)), abs(ally-gy(j))) );
        if dmin <= bw, covered = covered + 1; end
    end
end
fprintf('Approx 2D angular coverage with K=%d pilots: %.0f%%  (beamwidth ~%.3f)\n', ...
        Kx*Ky, 100*covered/total, bw);
fprintf('If coverage is low, raise per-axis pilots / retune sweep-rate ratio.\n');
