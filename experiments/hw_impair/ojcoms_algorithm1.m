function [sec, rho, dbg] = ojcoms_algorithm1(Nt, fc, B, M, d, L, thlim, alim, gamma, Nsec, Kalpha, s_sweep)
% OJCOMS_ALGORITHM1  Their sector-shifted / distance-interleaved pilot set.
%
% IEEE OJ-COMS v7 2026 (Qaid, Nasir, Al-Ahmadi, Liu), DOI 10.1109/OJCOMS.2026.3695965.
% Implements (31), (53)-(54), (61)-(65), (69)-(70) and, for L = 2, returns their
% published Table 3 verbatim so the baseline is THEIR configuration rather than a
% reconstruction of it.
%
% WHAT THE FIRST ATTEMPT GOT WRONG (recorded so it is not repeated):
%   * S was solved from boundary focusing, giving S = -11.92. Their text says
%     "S_m^(s) is chosen NEAR this upper bound", i.e. S ~ +S_bound ~ 35. Feeding
%     S = 35.53 and theta_M = 1 into their (63) reproduces theta't = -31.797
%     EXACTLY, which is what fixed the recipe.
%   * alpha't from (69)-(70) was already right: -0.5338 against their -0.533.
%
% s_sweep (optional, default 1) SLOWS THE ANGULAR SWEEP by scaling the composite
% sweep strength S, which is what shrinks |theta_t| and with it the per-element
% DELAY RANGE (sec 36). It must scale S ALONE: (63) is then re-solved so each
% sector's target direction theta_M stays where it was. Scaling theta_t and
% theta_p together instead -- the obvious move, and the one hw_budget did first --
% drags every sector centre toward boresight (theta_M goes 0.97, 0.49, 0.24, 0.12
% as s halves), so the array stops covering the search region and the rate
% collapses for a reason that has nothing to do with the sweep rate. Passing
% s_sweep ~= 1 forces the analytic branch, since Table 3 is an s = 1 design.
%
% STILL NOT DETERMINED BY THE TEXT, and flagged rather than guessed:
%   * the split of S into (theta'p, p_M) -- only their totals are recoverable;
%     Table 3 supplies theta'p directly for L = 2, so it is used there.
%   * what distinguishes the K_alpha = 3 distance pilots inside one sector.
%     Table 3 lists ONE (alpha't, alpha'p) per sector and the text says they are
%     "fixed across all sectors", which would make the three pilots identical.
%     Interpreted here as interleaving alpha't across the (70) interval; Kalpha=1
%     reproduces Table 3 literally. Report both.
if nargin == 0
    fprintf('ojcoms_algorithm1 generates their pilot set; call it from a script.\n');
    fprintf('Runnable scripts: rate_vs_snr, ablation_ps, sector_alloc, oracle_loc, ojcoms_baseline, phase_refine, kmin_fine, kspace_map, pilot_spacing_map, compare_ojcoms, hw_probe_*\n');
    return;
end
if nargin < 12 || isempty(s_sweep), s_sweep = 1; end
c = 3e8;  fL = fc - B/2;  fH = fc + B/2;  xiH = fH/fc;  G = Nt/L;

% ---- (31) LS coefficients ----
n   = (0:Nt-1).';
nh  = n - (Nt-1)/2;
nhs = floor(n/L) - (G-1)/2;
rho.theta = L*sum(nh.*nhs)/sum(nh.^2);
rho.alpha = L^2*sum(nh.^2 .* nhs.^2)/sum(nh.^4);

% ---- (69)-(70) distance parameters ----
U  = (alim(2)-alim(1)) / (fc/fL - fc/fH);
at = (alim(2) - (fc/fL)*U)/rho.alpha;              % (70) at equality
dbg.U = U; dbg.alpha_t_formula = at;

% ---- (59)-(60) sweep-strength bound, and (53)-(54) sector count ----
Sbound  = s_sweep * 1.76*gamma*fL*M/(Nt*B);      % s_sweep scales S ALONE
Wstrong = 5.56*Sbound / (L*pi*abs(31.73)*xiH^2);   % (53), nominal |theta't|
if isempty(Nsec), Nsec = max(1, ceil((thlim(2)-thlim(1))/abs(Wstrong))); end
dbg.Sbound = Sbound; dbg.Wstrong = Wstrong; dbg.Nsec_formula = max(1, ceil(2/abs(Wstrong)));

% ---- angle parameters ----
TAB3 = [-31.797 28.40; -33.840 31.93; -33.545 31.61; -32.102 28.51];  % their Table 3
if L == 2 && Nsec <= 4 && s_sweep == 1
    th_t = TAB3(1:Nsec,1);  th_p = TAB3(1:Nsec,2);
    dbg.source = 'Table 3 (verbatim)';
    % TABLE 3 CONTRADICTS THEIR OWN (69). With alpha'p = 0.158 and
    % alpha't = -0.533 the focus sweeps alpha in [-0.387, -0.361] -- the whole
    % pilot set points OUTSIDE the target interval [0.0025, 0.1], so no user can
    % ever be found. Their (69) gives alpha'p = U = 0.5809, which sweeps
    % [0.0033, 0.1008], i.e. exactly the target interval by construction; and
    % 0.5809 is also the DD-BS baseline's own alpha_p = 1859/3200. Three
    % independent confirmations that 0.158 is a typo. (69) is used, and the
    % table value is available via alpha_p_mode for comparison.
    if evalin('base','exist(''OJCOMS_ALPHA_P_TABLE'',''var'')') && ...
       evalin('base','OJCOMS_ALPHA_P_TABLE')
        alpha_p = 0.158;  dbg.source = 'Table 3 verbatim (alpha_p = 0.158, out of range)';
    else
        alpha_p = U;      dbg.source = 'Table 3 angles + eq (69) alpha_p';
    end
else
    % (61)-(63) with S near the upper bound; theta_M walks the sectors.
    th_t = zeros(Nsec,1); th_p = zeros(Nsec,1);
    for s = 1:Nsec
        thM = thlim(2) - (s-1)*(thlim(2)-thlim(1))/Nsec;
        S   = Sbound;
        th_t(s) = (thM - (fc/fH)*S)/rho.theta;      % (63)
        th_p(s) = S - (2/L)*round((L/2)*(S - 28.40));   % (62), p_M from their L=2 split
    end
    dbg.source = 'analytic (53)-(63)';
    alpha_p = U;
end
dbg.theta_t1 = th_t(1); dbg.theta_p1 = th_p(1); dbg.alpha_p = alpha_p;

% ---- assemble the pilot list ----
if isempty(Kalpha), Kalpha = 3; end
if Kalpha > 1
    % CENTRED, NESTED interleave: offsets 0, -1, +1, -2, +2, ... scaled to +/-spread.
    %
    % The previous version used linspace(-spread, spread, Kalpha), which for EVEN
    % Kalpha contains no zero -- so Kalpha = 1 sat exactly on the (70) optimum
    % while Kalpha = 2 sat on neither side of it, and adding a pilot made things
    % WORSE. hw_budget showed that in all four sweep rows (98->95, 98->95, 94->90,
    % 80->74), which reads as a non-monotonic K_alpha axis and is purely an
    % artifact of this placement, not of OJ-COMS's scheme (their text does not
    % specify the interleave -- see the note above).
    %
    % This ordering guarantees the set for Kalpha is a SUBSET of the set for
    % Kalpha+1, so an extra distance pilot can never remove the optimum and the
    % rate is monotone in Kalpha up to Monte-Carlo noise. Kalpha = 1 is unchanged.
    spread = abs(U)*0.15;
    k   = 0:(Kalpha-1);
    off = ceil(k/2) .* (-1).^k;                     % 0, -1, +1, -2, +2, ...
    ats = at + off * (spread/max(1, floor(Kalpha/2)));
else
    ats = at;
end
sec = struct('theta_t',{},'theta_p',{},'alpha_t',{},'alpha_p',{},'sector',{});
for s = 1:Nsec
    for ka = 1:Kalpha
        sec(end+1) = struct('theta_t',th_t(s),'theta_p',th_p(s), ...
                            'alpha_t',ats(ka),'alpha_p',alpha_p,'sector',s); %#ok<AGROW>
    end
end
dbg.Nsec = Nsec; dbg.Kalpha = Kalpha; dbg.Ktotal = numel(sec);
end
