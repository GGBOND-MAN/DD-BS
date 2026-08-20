function s = gen_multipath_scenario(p, L, kappa_dB, min_sep_sin)
% GEN_MULTIPATH_SCENARIO  Draw one random L-path near-field user realisation.
%
%   p          : struct from derisk_params()
%   L          : number of paths (L=1 -> pure LoS, reproduces baseline)
%   kappa_dB   : per-NLoS-path power relative to LoS, in dB (e.g. -10).
%                Each NLoS gain has magnitude sqrt(10^(kappa_dB/10)) and a
%                uniformly random phase. LoS gain is fixed to 1.
%   min_sep_sin: (optional) minimum separation between paths in the sine-angle
%                domain, to keep paths resolvable for the Gate-B separability
%                test. Default 0 (fully random, used for the Gate-A honesty test).
%
% Output struct s with fields:
%   s.r         1xL distances [m]
%   s.theta     1xL PHYSICAL angles [rad]      (feed to near_field_channel_multipath)
%   s.theta_sin 1xL sine-domain angles         (compare against DDBS grid / dictionary)
%   s.g         1xL complex path gains (g(1)=1 for LoS)
%
% Path 1 is always the LoS path.

if nargin < 4 || isempty(min_sep_sin), min_sep_sin = 0; end

r      = zeros(1, L);
th     = zeros(1, L);
th_sin = zeros(1, L);
g      = zeros(1, L);

kappa = 10^(kappa_dB/10);          % linear NLoS/LoS power ratio

for l = 1:L
    % rejection-sample an angle that is at least min_sep_sin (in sine domain)
    % away from previously drawn paths
    ok = false;
    while ~ok
        r_l  = p.Rmin + rand*(p.Rmax - p.Rmin);
        th_l = p.theta_min + rand*(p.theta_max - p.theta_min);
        sin_l = sin(th_l);
        if l == 1
            ok = true;
        else
            ok = all(abs(sin_l - th_sin(1:l-1)) >= min_sep_sin);
        end
    end
    r(l)      = r_l;
    th(l)     = th_l;
    th_sin(l) = sin_l;

    if l == 1
        g(l) = 1;                                  % LoS reference gain
    else
        g(l) = sqrt(kappa) * exp(1j*2*pi*rand);    % NLoS: |g|=sqrt(kappa), random phase
    end
end

s.r         = r;
s.theta     = th;
s.theta_sin = th_sin;
s.g         = g;
end
