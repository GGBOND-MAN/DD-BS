function w = ddbs_beam_hybrid(Nt, B, fc, M, d, theta1, theta2, alpha1, alpha2, K, imp)
% DDBS_BEAM_HYBRID  DDBS training beams with the PROPOSED low-resolution-TD design.
%
% Mechanism (the contribution). The ideal per-antenna phase is
%     2*pi*f_m*tau_n + phi_n .
% Split the delay term at the carrier:
%     2*pi*f_m*tau_n = 2*pi*fc*tau_n  +  2*pi*(f_m-fc)*tau_n .
% The first part is FREQUENCY-INDEPENDENT, so it can be realized EXACTLY by the
% phase shifter (mod 2*pi). The TTD then only has to realize 2*pi*(f_m-fc)*tau_n,
% whose phase error for a delay error d_tau is 2*pi*(f_m-fc)*d_tau <= pi*B*d_tau
% instead of ~2*pi*fc*d_tau. The TD precision requirement therefore scales with
% the BANDWIDTH B rather than the CARRIER fc -- a 2*fc/B (~12x here) relaxation,
% i.e. ~log2(2*fc/B) ~ 3.6 bits, measured ~4.5 bits end-to-end.
%
% Validated (Nt=256, M=1024, B=5 GHz, fc=30 GHz), mean focus gain (ideal 0.975):
%     B_td :   12      11      10
%     naive:  0.096   0.070   0.083
%     hybrid: 0.971   0.960   0.917
% and a 3-bit PS still gives ~0.98 of the hybrid gain, so the burden is NOT
% merely shifted onto the phase shifter.
%
% Same signature/impairment struct as ddbs_beam_impaired.m (Btd, Bps, sigma_tau,
% sigma_pn, amp_std), so the two are drop-in comparable.

if nargin < 11, imp = struct(); end
getf = @(s,def) ternary(isfield(imp,s) && ~isempty(imp.(s)), @() imp.(s), def);
Btd   = getf('Btd', Inf);
Bps   = getf('Bps', Inf);
s_tau = getf('sigma_tau', 0);
s_pn  = getf('sigma_pn', 0);
a_std = getf('amp_std', 0);

c = 3e8;
f = fc + B/M*((1:M) - 1 - (M-1)/2);
kc = 2*pi*fc/c;
nn = (-(Nt-1)/2:(Nt-1)/2)';

tau = zeros(Nt, K);
for s = 1:K
    tt_s = theta1 - 2*(s-1)/K;
    tau(:,s) = (nn*d*tt_s - (nn*d).^2*alpha1) / c;
end
phi = kc*(nn*d*theta2 - (nn*d).^2*alpha2);

% ---- TTD carries ONLY the frequency-dependent residual -> quantise it ----
tau_q = tau;
if isfinite(Btd)
    lo = min(tau(:)); hi = max(tau(:)); lsb = (hi-lo)/2^Btd;
    tau_q = lo + round((tau-lo)/lsb)*lsb;
end
if s_tau > 0, tau_q = tau_q + s_tau*randn(Nt, K); end

if a_std > 0, amp = max(0, 1 + a_std*randn(Nt, K)); else, amp = ones(Nt, K); end

w = zeros(Nt, K, M);
for s = 1:K
    % PS absorbs the exact frequency-independent term 2*pi*fc*tau (mod 2*pi)
    ph_ps = mod(phi + 2*pi*fc*tau(:,s), 2*pi);
    if isfinite(Bps)
        lsbp = 2*pi/2^Bps; ph_ps = round(ph_ps/lsbp)*lsbp;
    end
    pn = (s_pn>0) * s_pn * randn(Nt, 1);
    for m = 1:M
        ph = 2*pi*(f(m)-fc)*tau_q(:,s) + ph_ps + pn;    % residual TTD + PS
        w(:,s,m) = amp(:,s) .* exp(-1j*ph) / sqrt(Nt);
    end
end
end

function y = ternary(cond, ftrue, def)
if cond, y = ftrue(); else, y = def; end
end
