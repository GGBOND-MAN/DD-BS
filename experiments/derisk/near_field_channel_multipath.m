function [H, hc] = near_field_channel_multipath(Nt, d, fc, B, M, r_list, theta_phys_list, g_list)
% NEAR_FIELD_CHANNEL_MULTIPATH  Wideband near-field channel with L paths.
%
% Multipath generalisation of the baseline near_field_channel.m. Implements
% the Saleh-Valenzuela near-field model (paper Eq. (2)):
%
%     h_m = sum_{l=1}^L  g_l * exp(-j 2 pi f_m r_l / c) * a_m(theta_l, r_l)
%
% where a_m(.) is the spherical-wave array-response vector (near_field_manifold),
% which already carries the 1/sqrt(Nt) normalisation. The common phase term
% exp(-j 2 pi f_m r_l / c) uses tau_l = r_l/c, exactly matching how the baseline
% treats the LoS delay in near_field_channel.m.
%
% Inputs (per-path row vectors of length L):
%   r_list           : distance of each path to array centre [m]
%   theta_phys_list  : PHYSICAL direction of each path [rad], same convention as
%                      the baseline near_field_channel.m (manifold applies sin()).
%   g_list           : complex path gains. Convention: g_list(1)=1 for LoS.
%
% Outputs:
%   H  : M x Nt  channel over the M data subcarriers
%   hc : 1 x Nt  channel at the centre frequency fc (for the Perfect-CSI bound)
%
% Setting L=1 with g_list=1 reproduces near_field_channel.m exactly.

c = 3e8;
L = numel(r_list);
H = zeros(M, Nt);
hc = zeros(1, Nt);

for m = 1:M+1
    if m == M+1
        f = fc;                     % centre frequency -> hc
    else
        f = fc + B/M * (m - 1 - (M-1)/2);
    end

    hm = zeros(1, Nt);
    for l = 1:L
        at = near_field_manifold(Nt, d, f, r_list(l), theta_phys_list(l));  % 1 x Nt
        hm = hm + g_list(l) * exp(-1j*2*pi*f*r_list(l)/c) * at;
    end

    if m == M+1
        hc = hm;
    else
        H(m, :) = hm;
    end
end
end
