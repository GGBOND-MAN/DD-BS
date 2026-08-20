function w = ttd_beam_impaired(Nt, B, fc, M, d, sin_theta, alpha, lsb_sec, sigma_tau)
% TTD_BEAM_IMPAIRED  Serving (data) beam built with the SAME impaired hardware.
%
% Baseline TTD_beam assumes ideal delays. In reality the data beam uses the same
% TD hardware as training, so it inherits the same absolute time resolution.
% lsb_sec is the hardware delay LSB in SECONDS (not bits) -- the physically
% meaningful, architecture-independent knob; pass Inf/0 for ideal.
%
% NOTE the asymmetry this exposes: a serving beam focuses at sin_theta in [-1,1],
% so its delay range is only ~N*d/c (~2 ns here), whereas a DDBS TRAINING beam
% uses the large intercept theta_t (~-33) and spans ~140 ns. For the SAME hardware
% LSB the training beam is therefore far more fragile -- the difficulty is
% specific to the DDBS training waveform, not to near-field focusing in general.

c = 3e8;
f = fc + B/M*((1:M) - 1 - (M-1)/2);
nn = (-(Nt-1)/2:(Nt-1)/2)';
delay = nn*d*sin_theta - (nn*d).^2*alpha;      % metres of path
tau = delay / c;

if nargin >= 8 && isfinite(lsb_sec) && lsb_sec > 0
    tau = round(tau/lsb_sec)*lsb_sec;
end
if nargin >= 9 && sigma_tau > 0
    tau = tau + sigma_tau*randn(Nt,1);
end

w = zeros(Nt, M);
for m = 1:M
    w(:,m) = exp(-1j*2*pi*f(m)*tau) / sqrt(Nt);
end
end
