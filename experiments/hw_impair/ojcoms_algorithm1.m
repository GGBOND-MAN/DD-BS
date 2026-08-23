function [sec, rho, dbg] = ojcoms_algorithm1(Nt, fc, B, M, d, L, thlim, alim, gamma, Nsec, Kalpha)
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
% STILL NOT DETERMINED BY THE TEXT, and flagged rather than guessed:
%   * the split of S into (theta'p, p_M) -- only their totals are recoverable;
%     Table 3 supplies theta'p directly for L = 2, so it is used there.
%   * what distinguishes the K_alpha = 3 distance pilots inside one sector.
%     Table 3 lists ONE (alpha't, alpha'p) per sector and the text says they are
%     "fixed across all sectors", which would make the three pilots identical.
%     Interpreted here as interleaving alpha't across the (70) interval; Kalpha=1
%     reproduces Table 3 literally. Report both.
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
Sbound  = 1.76*gamma*fL*M/(Nt*B);
Wstrong = 5.56*Sbound / (L*pi*abs(31.73)*xiH^2);   % (53), nominal |theta't|
if isempty(Nsec), Nsec = max(1, ceil((thlim(2)-thlim(1))/abs(Wstrong))); end
dbg.Sbound = Sbound; dbg.Wstrong = Wstrong; dbg.Nsec_formula = max(1, ceil(2/abs(Wstrong)));

% ---- angle parameters ----
TAB3 = [-31.797 28.40; -33.840 31.93; -33.545 31.61; -32.102 28.51];  % their Table 3
if L == 2 && Nsec <= 4
    th_t = TAB3(1:Nsec,1);  th_p = TAB3(1:Nsec,2);
    dbg.source = 'Table 3 (verbatim)';
    alpha_p = 0.158;                                % Table 3
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
    spread = abs(U)*0.15;                           % interleave alpha't -- see header
    ats = at + linspace(-spread, spread, Kalpha);
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
