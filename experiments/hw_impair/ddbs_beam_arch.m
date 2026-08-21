function [w, n_ttd] = ddbs_beam_arch(Nt,B,fc,M,d,theta1,theta2,alpha1,alpha2,K,arch,P,Btd,Bps)
% DDBS_BEAM_ARCH  DDBS training beams under a chosen TD ARCHITECTURE x RESOLUTION.
%
% Spans the 2-D hardware cost surface: TD COUNT (arch,P) x TD PRECISION (Btd).
%
% arch:
%   'full'   one TTD per antenna                       -> n_ttd = Nt
%   'shared' generic sub-array sharing, P per group    -> n_ttd = Nt/P
%            (one TTD per sub-array; PS compensates the intra-group difference
%             at fc). This is a GENERIC shared-TTD model -- it is NOT a
%             reproduction of the OJ-COMS AoSA scheme, which additionally
%             redesigns the training procedure.
%   'kron'   PROPOSED two-stage / Kronecker TTD        -> n_ttd = Nt/P + P
%
% Kronecker rationale. With n = n_g + n_q (n_g the sub-array offset, n_q the
% position inside it -- an exact split), the DDBS delay
%       tau_n = (n*d*theta_t - n^2*d^2*alpha_t)/c
% has a LINEAR term that is exactly separable, tau ~ A(n_g) + C(n_q), so
% Nt/P coarse TTDs plus P fine TTDs (reused by every sub-array) suffice. The only
% residual is the quadratic cross term
%       d_tau = 2*n_g*n_q*d^2*|alpha_t|/c   <=  Nt*P*d^2*|alpha_t|/(2c),
% which is small because the DDBS delay is linear-dominated (140 ns vs 0.72 ns
% here). Total TTDs Nt/P + P is minimised at P = sqrt(Nt) -> 2*sqrt(Nt).
%
% Btd/Bps: TD/PS resolution in bits (Inf = ideal). The PS always carries the exact
% frequency-independent term 2*pi*fc*tau_n (the hybrid split), so this function
% composes the count reduction with the precision relaxation.

if nargin<12||isempty(P),   P   = 1;   end
if nargin<13||isempty(Btd), Btd = Inf; end
if nargin<14||isempty(Bps), Bps = Inf; end

c = 3e8;
f  = fc + B/M*((1:M)-1-(M-1)/2);
kc = 2*pi*fc/c;
nn = (-(Nt-1)/2:(Nt-1)/2)';
G  = Nt/P;

w = zeros(Nt,K,M);
for s = 1:K
    tt_s = theta1 - 2*(s-1)/K;
    tau  = (nn*d*tt_s - (nn*d).^2*alpha1)/c;      % exact per-antenna delay

    switch arch
        case 'full'
            tg = tau;                    n_ttd = Nt;
        case 'shared'
            tg = repelem(mean(reshape(tau,P,G),1).',P);   n_ttd = G;
        case 'kron'
            idx = (0:Nt-1).';  g = floor(idx/P);  q = mod(idx,P);
            ng  = (g-(G-1)/2)*P;         nq = q-(P-1)/2;   % n = ng+nq exactly
            A   = (ng*d*tt_s - (ng*d).^2*alpha1)/c;        % G coarse values
            C   = (nq*d*tt_s - (nq*d).^2*alpha1)/c;        % P fine values
            tg  = A + C;                 n_ttd = G + P;
        otherwise
            error('arch must be full|shared|kron');
    end

    if isfinite(Btd)                                   % quantise the TD values
        lo = min(tau); R = max(tau)-min(tau); lsb = R/2^Btd;
        tg = lo + round((tg-lo)/lsb)*lsb;
    end

    % PS carries the exact per-antenna fc term + the DDBS phase-shift parameters
    ps = mod(kc*(nn*d*theta2 - (nn*d).^2*alpha2) + 2*pi*fc*tau, 2*pi);
    if arch=="shared"   % PS also compensates the intra-group delay difference
        ps = mod(ps + 2*pi*fc*(tau-tg), 2*pi);
    end
    if isfinite(Bps), lsbp = 2*pi/2^Bps; ps = round(ps/lsbp)*lsbp; end

    for m = 1:M
        w(:,s,m) = exp(-1j*(2*pi*(f(m)-fc)*tg + ps))/sqrt(Nt);
    end
end
end
