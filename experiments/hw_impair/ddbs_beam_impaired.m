function w = ddbs_beam_impaired(Nt, B, fc, M, d, theta1, theta2, alpha1, alpha2, K, imp)
% DDBS_BEAM_IMPAIRED  DDBS training beams with hardware impairments.
%
% Regenerates the baseline delay_polar_2d beams from the underlying TD delays and
% PS phases, then applies finite-resolution / error impairments. With imp all
% ideal (Inf bits, 0 error) this reproduces delay_polar_2d exactly.
%
% The DDBS per-antenna phase splits into:
%   TD (true-time-delay): 2*pi*f_m * tau_s(n),  tau_s(n) = (n d theta_t,s - n^2 d^2 alpha_t)/c
%   PS (phase-shift):      phi(n)             = kc (n d theta_p - n^2 d^2 alpha_p)
% Pilots differ only through the TD intercept theta_t,s = theta1 - 2(s-1)/K (baseline Alg.1).
%
% imp fields (any omitted -> ideal):
%   imp.Btd        : TD resolution in bits over the delay range (Inf = ideal)
%   imp.Bps        : PS resolution in bits over [0,2pi) (Inf = ideal)
%   imp.sigma_tau  : TD delay jitter std in SECONDS (fabrication error; 0 = none)
%   imp.sigma_pn   : per-element phase-noise std in RADIANS (0 = none)
%   imp.amp_std    : per-element amplitude error std (0 = none; models IL variation)
%
% Returns w (Nt x K x M).

if nargin < 11, imp = struct(); end
getf = @(s,def) ternary(isfield(imp,s) && ~isempty(imp.(s)), @() imp.(s), def);
Btd   = getf('Btd', Inf);
Bps   = getf('Bps', Inf);
s_tau = getf('sigma_tau', 0);
s_pn  = getf('sigma_pn', 0);
a_std = getf('amp_std', 0);

c = 3e8;
f = fc + B/M*((1:M) - 1 - (M-1)/2);
km = 2*pi*f/c; kc = 2*pi*fc/c;
nn = (-(Nt-1)/2:(Nt-1)/2)';

% ---- ideal TD delays (per pilot) and PS phase ----
tau = zeros(Nt, K);
for s = 1:K
    tt_s = theta1 - 2*(s-1)/K;
    tau(:,s) = (nn*d*tt_s - (nn*d).^2*alpha1) / c;
end
phi = kc*(nn*d*theta2 - (nn*d).^2*alpha2);       % Nt x 1

% ---- TD impairments ----
if isfinite(Btd)
    lo = min(tau(:)); hi = max(tau(:)); lsb = (hi-lo)/2^Btd;
    tau = lo + round((tau-lo)/lsb)*lsb;
end
if s_tau > 0
    tau = tau + s_tau*randn(Nt, K);
end

% ---- PS impairments ----
phi_q = phi;
if isfinite(Bps)
    lsb = 2*pi/2^Bps;
    phi_q = round(mod(phi,2*pi)/lsb)*lsb;
end

% ---- amplitude (IL variation) ----
if a_std > 0
    amp = max(0, 1 + a_std*randn(Nt, K));
else
    amp = ones(Nt, K);
end

% ---- compose beams ----
w = zeros(Nt, K, M);
for s = 1:K
    pn = (s_pn>0) * s_pn * randn(Nt, 1);          % phase noise per element/pilot
    for m = 1:M
        ph = 2*pi*f(m)*tau(:,s) + phi_q + pn;
        w(:,s,m) = amp(:,s) .* exp(-1j*ph) / sqrt(Nt);
    end
end
end

function y = ternary(cond, ftrue, def)
if cond, y = ftrue(); else, y = def; end
end
