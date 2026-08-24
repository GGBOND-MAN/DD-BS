function ws = serve_beam(Nt,B,fc,M,d,th,al,mode,P)
% SERVE_BEAM  Data beam at (th, al), built on the hardware the system actually has.
%
% THE GAP THIS CLOSES (THEORY sec 31.2). Every rate in this project was measured
% with TTD_beam, i.e. ONE TTD PER ANTENNA, no matter how few delays the training
% stage was restricted to. So the numbers described a system that trains on 16-32
% delays and then serves on 256. sec 5-11 had already found the serving beam
% dominates end-to-end rate; that lesson was not carried into the OJ-COMS work.
%
%   mode = 'ideal'   per-antenna TTD -- the old behaviour, kept as a reference row
%                    so the cost of grouping the serving stage is visible
%   mode = 'shared'  Nt/P delays on the sub-array common legs + per-antenna phase
%                    shifters carrying the exact fc term -- the same hardware the
%                    training stage is limited to
%
% WHY BOTH ARMS SERVE WITH THE COMPENSATED PROGRAMMING, including "theirs".
% A serving beam must focus AT fc, and with a grouped array only the PS can supply
% the intra-sub-array phase to do it. Programming the PS their way for the serving
% beam (their (14) with the focus parameters, or with zeros) leaves an
% uncompensated 2*pi*fc*(tau_g - tau_n) at EVERY subcarrier including fc, so the
% beam never focuses at all -- that is a strawman, not their architecture. Their
% (14) describes the TRAINING beam design. Holding the serving stage identical
% therefore isolates the training-beam difference, which is the actual claim, and
% it is the reading generous to them.
%
% Self-check: at P = 1 the shared path reduces to exp(-j 2 pi f_m tau), which is
% exactly TTD_beam. Verified once per session below.
persistent checked
if nargin < 9 || isempty(P), P = 1; end
if strcmp(mode,'ideal') || P == 1
    ws = TTD_beam(Nt,B,fc,M,d,th,al);
    return;
end
if isempty(checked)
    a = TTD_beam(Nt,B,fc,M,d,th,al);
    b = squeeze(ddbs_beam_arch(Nt,B,fc,M,d,th,0,al,0,1,'shared',1,Inf,Inf,[]));
    e = max(abs(a(:)-b(:)));
    fprintf('[serve_beam] self-check at P=1: max|TTD_beam - shared| = %.2e\n', e);
    if e > 1e-9, warning('serve_beam:selftest','P=1 path does NOT reduce to TTD_beam'); end
    checked = true;
end
ws = squeeze(ddbs_beam_arch(Nt,B,fc,M,d,th,0,al,0,1,'shared',P,Inf,Inf,[]));
end
