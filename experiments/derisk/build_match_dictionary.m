function [match_y, g1_list, g2_list] = build_match_dictionary(p, w_2D_near, g1, g2)
% BUILD_MATCH_DICTIONARY  Precompute the single-location match-filter dictionary.
%
% Identical construction to the dictionary in the baseline Rate_snr.m: for every
% sampled grid location (theta, alpha) it stores the *noiseless* received power
% fingerprint |a_mf(theta,alpha) * w_DDBS|^2 across all subcarriers and pilots,
% where a_mf is a single-location TTD focusing beam (TTD_beam).
%
% This is exactly the object whose SINGLE-PATH assumption we are stress-testing:
% under multipath the received fingerprint is a superposition of several such
% templates, so correlating against this dictionary is a template mismatch.
%
%   p         : struct from derisk_params()
%   w_2D_near : Nt x k x M DDBS training beams (from delay_polar_2d)
%   g1        : # angle grid points  (sine domain)
%   g2        : # distance-ring (alpha) grid points
%
% Outputs:
%   match_y : g1 x g2 x k x M x Q dictionary of received power fingerprints
%   g1_list : 1 x g1 sampled sine-angles
%   g2_list : 1 x g2 sampled alpha values

Nt = p.Nt; B = p.B; fc = p.fc; M = p.M; d = p.d; k = p.k; Q = p.Q;

g1_list = linspace(sin(p.theta_min), sin(p.theta_max), g1);
g2_list = linspace(1/400, 1/10, g2);     % same alpha range as baseline

match_y = zeros(g1, g2, k, M, Q);
parfor i = 1:g1
    row = zeros(g2, k, M, Q);
    for j = 1:g2
        w_near = TTD_beam(Nt, B, fc, M, d, g1_list(i), g2_list(j));  % Nt x M
        mf = w_near';                                               % M x Nt
        for t = 1:k
            y_near = zeros(M, Q);
            for m = 1:M
                y_near(m, :) = repmat(mf(m, :) * w_2D_near(:, t, m), [1, Q]);
            end
            y_near = abs(sum(y_near, 2)).^2;
            row(j, t, :, :) = y_near;
        end
    end
    match_y(i, :, :, :, :) = row;
end
end
