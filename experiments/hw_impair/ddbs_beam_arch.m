function [w, n_ttd, lsb_used] = ddbs_beam_arch(Nt,B,fc,M,d,theta1,theta2,alpha1,alpha2,K,arch,P,Btd,Bps,lsb_sec)
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
%   'ojcoms' literal OJ-COMS AoSA form (their eqs (13)-(14)): TTD carries the
%            full 2*pi*f_m*tau on the sub-array-centre grid; the PS is the plain
%            per-antenna DDBS term and does NOT compensate the intra-group delay.
%            Note their full scheme ALSO redesigns the DDBS parameters and adds
%            sector-shift / distance-interleaving / extra pilots, which this does
%            not reproduce -- so treat it as the architecture, not their result.
%   'kron'   Kronecker two-stage delay factorization    -> Nt/P + P DISTINCT
%            delay VALUES. NOTE (erratum): this is NOT a physical element-count
%            reduction. Sharing applies a delay before the fan-out so P antennas
%            truly share one element; the Kronecker fine term C(n_q) must instead
%            be realized per antenna, so the element count is Nt/P + Nt. Treat the
%            returned n_ttd for 'kron' as a distinct-value count only.
%
% A serving (data) beam focused at (sin_theta, alpha) is just this function with
% theta1=sin_theta, alpha1=alpha, theta2=alpha2=0, K=1 -- so the same hardware
% model covers both training and data beams.
%
% Kronecker rationale (distinct-value reduction only -- see the note above).
% With n = n_g + n_q (n_g the sub-array offset, n_q the position inside it -- an
% exact split), the DDBS delay
%       tau_n = (n*d*theta_t - n^2*d^2*alpha_t)/c
% has a LINEAR term that is exactly separable, tau ~ A(n_g) + C(n_q), so
% Nt/P coarse plus P fine delay VALUES suffice (the fine values recur across
% sub-arrays but need their own hardware). The only residual is the cross term
%       d_tau = 2*n_g*n_q*d^2*|alpha_t|/c   <=  Nt*P*d^2*|alpha_t|/(2c),
% which is small because the DDBS delay is linear-dominated (140 ns vs 0.72 ns
% here). The distinct-value count Nt/P + P is minimised at P = sqrt(Nt).
%
% Btd/Bps: TD/PS resolution in bits (Inf = ideal). The PS always carries the exact
% frequency-independent term 2*pi*fc*tau_n (the hybrid split), so this function
% composes the count reduction with the precision relaxation.

if nargin<12||isempty(P),   P   = 1;   end
if nargin<13||isempty(Btd), Btd = Inf; end
if nargin<14||isempty(Bps), Bps = Inf; end
if nargin<15, lsb_sec = []; end   % explicit hardware LSB [s]; overrides Btd

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
        case {'shared','ojcoms'}
            % one TTD per sub-array, evaluated at the sub-array centre
            tg = repelem(mean(reshape(tau,P,G),1).',P);   n_ttd = G;
        case 'kron'
            idx = (0:Nt-1).';  g = floor(idx/P);  q = mod(idx,P);
            ng  = (g-(G-1)/2)*P;         nq = q-(P-1)/2;   % n = ng+nq exactly
            A   = (ng*d*tt_s - (ng*d).^2*alpha1)/c;        % G coarse values
            C   = (nq*d*tt_s - (nq*d).^2*alpha1)/c;        % P fine values
            tg  = A + C;
            if P==1, n_ttd = Nt; else, n_ttd = G + P; end
        otherwise
            error('arch must be full|shared|ojcoms|kron');
    end

    % Quantise the TD values. lsb_sec (an absolute hardware LSB, in seconds) takes
    % precedence over Btd, which is defined relative to THIS beam's delay range --
    % use lsb_sec whenever training and serving beams must share one hardware.
    lo = min(tau); R = max(tau)-min(tau);
    if ~isempty(lsb_sec) && isfinite(lsb_sec) && lsb_sec > 0
        lsb = lsb_sec;
    elseif isfinite(Btd)
        lsb = R/2^Btd;
    else
        lsb = 0;
    end
    if lsb > 0, tg = lo + round((tg-lo)/lsb)*lsb; end
    lsb_used = lsb;

    % PS carries the exact per-antenna fc term + the DDBS phase-shift parameters
    if strcmp(arch,'ojcoms')
        % Literal OJ-COMS AoSA form (their eqs (13)-(14)): the TTD carries the
        % FULL 2*pi*f_m*tau_g on the sub-array-centre grid, and the PS is the
        % standard per-antenna DDBS phase term with its own design parameters --
        % it does NOT compensate the intra-group delay difference.
        ps = mod(kc*(nn*d*theta2 - (nn*d).^2*alpha2), 2*pi);
        if isfinite(Bps), lsbp = 2*pi/2^Bps; ps = round(ps/lsbp)*lsbp; end
        for m = 1:M
            w(:,s,m) = exp(-1j*(2*pi*f(m)*tg + ps))/sqrt(Nt);
        end
    else
        % Hybrid split: the PS carries the EXACT per-antenna frequency-independent
        % term 2*pi*fc*tau_n, so the residual error is 2*pi*(f_m-fc)*(tg - tau_n)
        % -- for 'shared' this already gives the PS full per-antenna resolution at
        % fc (the strongest fair sharing baseline; do NOT add another fc term).
        ps = mod(kc*(nn*d*theta2 - (nn*d).^2*alpha2) + 2*pi*fc*tau, 2*pi);
        if isfinite(Bps), lsbp = 2*pi/2^Bps; ps = round(ps/lsbp)*lsbp; end
        for m = 1:M
            w(:,s,m) = exp(-1j*(2*pi*(f(m)-fc)*tg + ps))/sqrt(Nt);
        end
    end
end
end
