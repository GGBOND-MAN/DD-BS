function w = ttd_beam_impaired(Nt, B, fc, M, d, sin_theta, alpha, lsb_sec, sigma_tau, use_hybrid)
% TTD_BEAM_IMPAIRED  Serving (data) beam built with the SAME impaired hardware.
%
% Baseline TTD_beam assumes ideal delays. In reality the data beam uses the same
% TD hardware as training, so it inherits the same absolute time resolution.
% lsb_sec is the hardware delay LSB in SECONDS (architecture-independent knob);
% pass Inf for ideal.
%
% use_hybrid (default false): apply the SAME TD-PS split as ddbs_beam_hybrid --
%   2*pi*f_m*tau = 2*pi*fc*tau (exact, in the PS)  +  2*pi*(f_m-fc)*tau (in the TTD)
% The baseline's TTD_beam is pure-TD, but the array hardware is TD-PS (paper
% Fig. 1), so the phase shifter is available for the data beam as well. Applying
% the fix to ONLY the training beam leaves the serving beam as the bottleneck and
% erases the end-to-end benefit -- both stages must use the same architecture.
%
% Measured serving-beam focusing gain (range ~2.1 ns):
%   LSB [ps] :   16.5    32.9   131.7
%   pure-TD  :  0.650   0.044   0.043
%   hybrid   :  0.999   0.996   0.941
%
% Note the precision requirement is set by the CARRIER PERIOD (not the delay
% range): even this 2.1 ns beam is destroyed at a 32.9 ps LSB. The range only
% sets the BIT COUNT, bits = log2(range/LSB).

if nargin < 10 || isempty(use_hybrid), use_hybrid = false; end
if nargin < 9  || isempty(sigma_tau),  sigma_tau  = 0;     end
if nargin < 8  || isempty(lsb_sec),    lsb_sec    = Inf;   end

c = 3e8;
f = fc + B/M*((1:M) - 1 - (M-1)/2);
nn = (-(Nt-1)/2:(Nt-1)/2)';
tau = (nn*d*sin_theta - (nn*d).^2*alpha) / c;

tau_q = tau;
if isfinite(lsb_sec) && lsb_sec > 0
    tau_q = round(tau/lsb_sec)*lsb_sec;
end
if sigma_tau > 0
    tau_q = tau_q + sigma_tau*randn(Nt,1);
end

w = zeros(Nt, M);
if use_hybrid
    ph_ps = mod(2*pi*fc*tau, 2*pi);            % exact fc term via the phase shifter
    for m = 1:M
        w(:,m) = exp(-1j*(2*pi*(f(m)-fc)*tau_q + ph_ps)) / sqrt(Nt);
    end
else
    for m = 1:M
        w(:,m) = exp(-1j*2*pi*f(m)*tau_q) / sqrt(Nt);
    end
end
end
