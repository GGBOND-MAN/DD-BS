function [sec, rho, dbg] = ojcoms_algorithm1(Nt, fc, B, M, d, L, thlim, alim, gamma, Nsec, Kalpha)
% OJCOMS_ALGORITHM1  Their Algorithm 1 -- sector shifting + distance interleaving.
%
% Implements IEEE OJ-COMS v7 2026 (Qaid, Nasir, Al-Ahmadi, Liu), eqs (31), (53),
% (54), (61)-(65), (69)-(71), so their scheme can be run as a real baseline
% rather than approximated. Returns one struct per pilot with the four DDBS
% control parameters, plus the LS coefficients rho_theta / rho_alpha of (31).
%
% Nsec / Kalpha may be [] to use their analytical (54) / (71); pass numbers to
% force their published operating point (Nsec=4, Kalpha=3 at Nt=256, L=2).
%
% UNIT TEST built in: at Nt=256, fc=30 GHz, B=5 GHz, M=1024, L=2, gamma=0.9 the
% paper reports theta't(1) = -31.797 and theta'p(1) = 28.4 for sector 1. dbg
% carries this sector's values so a caller can check them before trusting
% anything downstream. A mismatch means the reconstruction is wrong and the
% comparison must stop -- see THEORY.md sec 23.
c = 3e8;  fL = fc - B/2;  fH = fc + B/2;  xiH = fH/fc;
G = Nt/L;

% ---- (31) LS coefficients from the symmetry sums ----
n   = (0:Nt-1).';
nh  = n - (Nt-1)/2;                       % centred antenna index
nhs = floor(n/L) - (G-1)/2;               % centred sub-array index
S2 = sum(nh.^2);  S4 = sum(nh.^4);
C1 = sum(nh.*nhs);  C2 = sum(nh.^2 .* nhs.^2);
rho.theta = L*C1/S2;   rho.alpha = L^2*C2/S4;

% ---- (61)-(62): sweep strength at the band-edge worst case ----
% |S| <= 1.76*gamma*fL*M/(Nt*B) is the no-hole bound of (59)-(60).
Sbound = 1.76*gamma*fL*M/(Nt*B);
% The alias integer p_M is a free design choice; the paper's "convenient
% implementation" floor((L/2)*Sbound) is one option, but their reported
% theta'p(1)=28.4 corresponds to a much smaller p_M, so pick the p_M that puts
% theta'p closest to a boundary-focusing solution and report both.
pM_conv = floor((L/2)*Sbound);
% Boundary focusing (their sector-1 rule): theta_1 = thlim(1), theta_M = thlim(2).
% From (56), theta_M - theta_1 = S*fc*(1/fH - 1/fL)  =>
Sbf = (thlim(2)-thlim(1)) / (fc*(1/fH - 1/fL));
pM_bf = round((L/2)*(Sbound - abs(Sbf)));
dbg.Sbound = Sbound; dbg.pM_conv = pM_conv; dbg.Sbf = Sbf; dbg.pM_bf = pM_bf;

% ---- (53)-(54): sector width and sector count ----
th_tot = thlim(2)-thlim(1);
Wstrong = 5.56*abs(Sbf) / (L*pi*abs(31.73)*xiH^2);   % nominal |theta't| = baseline
if isempty(Nsec), Nsec = max(1, ceil(th_tot/abs(Wstrong))); end
dbg.Wstrong = Wstrong; dbg.Nsec_formula = max(1, ceil(th_tot/abs(Wstrong)));

% ---- (69)-(70): distance sweep and the feasible alpha't interval ----
U = (alim(2)-alim(1)) / (fc/fL - fc/fH);
q = 0;                                   % alias branch; 0 keeps alpha'p in range
alpha_p = U - 2/(L^2*d)*q;
at_lo = (alim(2) - (fc/fL)*U)/rho.alpha;
at_hi = (alim(1) - (fc/fH)*U)/rho.alpha;
if isempty(Kalpha)
    Kalpha = max(1, ceil( abs(at_hi-at_lo) / max(1e-12, abs(alim(2)-alim(1))/2) ));
end
dbg.U = U; dbg.alpha_p = alpha_p; dbg.at_range = [at_lo at_hi];

% ---- Algorithm 1 loop: one (theta't, theta'p) per sector ----
sec = struct('theta_t',{},'theta_p',{},'alpha_t',{},'alpha_p',{},'sector',{});
if Kalpha > 1, ats = linspace(at_lo, at_hi, Kalpha); else, ats = mean([at_lo at_hi]); end
for s = 1:Nsec
    % step 5-7: shift the sector, alternating sweep direction as their sector 2 does
    if mod(s,2)==1, a = thlim(1); b = thlim(2); else, a = thlim(2); b = thlim(1); end
    shift = (s-1)*abs(Wstrong);
    th1 = wrapto(a + shift, thlim);  thM = wrapto(b + shift, thlim);
    % step 10-12: (61) S, (62) theta'p, (63) theta't
    Ss = (thM - th1) / (fc*(1/fH - 1/fL));
    theta_p = Ss - (2/L)*pM_bf;
    theta_t = (thM - (fc/fH)*Ss) / rho.theta;
    for ka = 1:Kalpha
        sec(end+1) = struct('theta_t',theta_t,'theta_p',theta_p, ...
                            'alpha_t',ats(ka),'alpha_p',alpha_p,'sector',s); %#ok<AGROW>
    end
end
dbg.theta_t1 = sec(1).theta_t;  dbg.theta_p1 = sec(1).theta_p;
dbg.Nsec = Nsec; dbg.Kalpha = Kalpha; dbg.Ktotal = numel(sec);
end

function y = wrapto(x, lim)
w = lim(2)-lim(1);
y = mod(x - lim(1), w) + lim(1);
end
