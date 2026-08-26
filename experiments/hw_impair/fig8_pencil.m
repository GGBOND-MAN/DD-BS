% ========================================================================
% FIGURE 8 -- the focus locus is a pencil of lines through one pivot
%             (THEORY sec 40.1-40.3).  COMPUTED, not hand-drawn.
% ------------------------------------------------------------------------
% This is the mechanism figure. It needs no channel model and no Monte Carlo --
% every curve is the closed form
%
%       theta_m = < u_m*theta'p + rho_theta*theta_t >_{2*u_m}     u_m = fc/f_m
%       alpha_m =   u_m*alpha'p + rho_alpha*alpha_t
%
% evaluated at the real design parameters. The modulus is 2*u_m = 2*fc/f_m, the
% grating period of a lambda_c/2 array AT f_m -- frequency dependent, and the
% term every earlier mechanism attempt missed (sec 40.1).
%
% Eliminating u_m between the two coordinates turns the locus into a PENCIL OF
% STRAIGHT LINES through the pivot (rho_theta*theta_t, rho_alpha*alpha_t), one
% line per grating order j, with slope alpha'p / (theta'p - 2j). By OJ-COMS's
% (63) the pivot sits at distance ~(fc/fH)*S from the sector target, so the
% sweep strength S -- what our s scales -- IS the pivot distance. That is the
% whole figure:
%
%   LEFT  (s = 1)    pivot far away  -> several steep branches, each crossing the
%                    whole angular window, sampling distance finely at each pass
%   RIGHT (s = 1/8)  pivot close     -> one shallow branch covering the whole
%                    distance range but only a slice of angle
%
% Top row: the pencil, drawn out to the pivot. Bottom row: the same locus after
% the mod-2u_m wrap, i.e. what actually lands in the search region.
%
% The invariant to point at in the caption: the distance coverage of one pilot
% is alpha'p * Lambda = H in BOTH columns -- (69) makes it so at every sweep
% rate. What s buys is angular traversals, n_tr ~ (fc/fH)*S*Lambda/2.
% ========================================================================
clear; clc;
addpath(genpath(fileparts(mfilename('fullpath'))));

Nt=256; fc=30e9; B=5e9; M=1024; c=3e8; d=(c/fc)/2;
f  = fc + B/M*((1:M)-1-(M-1)/2);
u  = fc./f;                                  % 1 x M, in [fc/fH, fc/fL]
thlim=[-1 1]; alim=[0.0025 0.1]; gamma=0.9;
L = 8; Nsec = 8;                             % 32 TTDs, the headline operating point
fL=fc-B/2; fH=fc+B/2; Lam = fc/fL - fc/fH; H = alim(2)-alim(1);

sweeps = [1 1/8];
figure('Position',[80 80 1000 640]);

for col = 1:2
    s_sw = sweeps(col);
    [sec, rho, dbg] = ojcoms_algorithm1(Nt,fc,B,M,d,L,thlim,alim,gamma,Nsec,1,s_sw);
    sc  = sec(1);
    thp = sc.theta_p;  alp = sc.alpha_p;
    pth = rho.theta*sc.theta_t;              % pivot, theta
    pal = rho.alpha*sc.alpha_t;              % pivot, alpha
    th_nom = u*thp + pth;
    al     = u*alp + pal;
    per    = 2*u;                            % frequency-dependent grating period
    jj     = round(th_nom./per);             % active grating order per subcarrier
    th_w   = th_nom - per.*jj;               % what actually focuses
    vis    = abs(th_w) <= 0.9;               % the search region actual_focus scans
    js     = unique(jj(vis));                % orders that land inside it
    cl     = lines(max(numel(js),3));
    bmin = inf; bmax = -inf;                 % extent of the branches, for xlim
    for i = 1:numel(js)
        bq = u*(thp - 2*js(i)) + pth;
        bmin = min(bmin, min(bq)); bmax = max(bmax, max(bq));
    end
    ntr    = abs(thp - 2*round((thp+pth)/2))*Lam/2;

    fprintf('s = %-6.3g  S = %7.3f  theta''p = %7.3f  pivot = (%8.3f, %7.4f)\n', ...
            s_sw, dbg.Sbound, thp, pth, pal);
    fprintf('            orders inside |theta|<=0.9: %s   n_tr = %.3f\n', mat2str(js), ntr);
    fprintf('            alpha swept by this pilot = %.4f   (H = %.4f)  <- invariant in s\n\n', ...
            max(al(vis))-min(al(vis)), H);

    % ---------- top row: the pencil ----------
    subplot(2,2,col); hold on; box on; grid on;
    uu = linspace(0, max(u)*1.02, 200);
    for i = 1:numel(js)
        q = thp - 2*js(i);
        plot(uu*q+pth, uu*alp+pal, ':', 'Color', cl(i,:), 'LineWidth', 1.0, ...
             'HandleVisibility','off');                       % extension to the pivot
        k = (jj == js(i)) & vis;
        plot(u(k)*q+pth, al(k), '-', 'Color', cl(i,:), 'LineWidth', 2.4, ...
             'DisplayName', sprintf('order j = %d', js(i)));   % the swept segment
    end
    plot(pth, pal, 'kp', 'MarkerSize', 14, 'MarkerFaceColor','k', ...
         'DisplayName','pivot (\rho_\theta\theta_t, \rho_\alpha\alpha_t)');
    yl = [min([al pal])-0.03, max(al)+0.03];
    plot([-0.9 -0.9], yl, 'k--','LineWidth',1.1,'HandleVisibility','off');
    plot([ 0.9  0.9], yl, 'k--','LineWidth',1.1,'DisplayName','search region |\theta| \leq 0.9');
    pad = 0.06*max(bmax-bmin, 1);
    ylim(yl); xlim([min(bmin,pth)-pad, bmax+pad]);
    xlabel('\theta = sin\vartheta   (unwrapped)'); ylabel('\alpha  (1/m)');
    title(sprintf('(%c)  s = %g:  pivot at \\theta = %.1f,  %d branch(es)', ...
                  'a'+col-1, s_sw, pth, numel(js)));
    legend('Location','northwest','FontSize',8);

    % ---------- bottom row: after the mod-2u_m wrap ----------
    subplot(2,2,2+col); hold on; box on; grid on;
    patch([-0.9 0.9 0.9 -0.9],[alim(1) alim(1) alim(2) alim(2)], ...
          [0.90 0.94 0.99],'EdgeColor',[0.4 0.5 0.7],'LineWidth',1.1, ...
          'DisplayName','search region');
    for i = 1:numel(js)
        k = (jj == js(i)) & vis;
        plot(th_w(k), al(k), '.', 'Color', cl(i,:), 'MarkerSize', 9, ...
             'DisplayName', sprintf('order j = %d', js(i)));
    end
    xlim([-1.02 1.02]); ylim([alim(1)-0.012, alim(2)+0.012]);
    xlabel('\theta = sin\vartheta   (after mod 2f_c/f_m)'); ylabel('\alpha  (1/m)');
    title(sprintf('(%c)  n_{tr} = %.2f angular traversals, \\alpha covers %.4f', ...
                  'c'+col-1, ntr, max(al(vis))-min(al(vis))));
    legend('Location','northeast','FontSize',8);
end

sgtitle(sprintf(['focus locus of ONE pilot,  N_t=%d, L=%d (%d TTDs), B=%.0f GHz  ' ...
                 '\\Lambda=%.4f'], Nt, L, Nt/L, B/1e9, Lam));

fprintf(['Caption to use: the focus locus of a single DDBS pilot is a pencil of\n' ...
         'straight lines through one pivot, one line per grating order, because the\n' ...
         'grating period of the focus map is 2*fc/f_m and therefore breathes across\n' ...
         'the band. OJ-COMS eq (63) places the pivot at distance ~(fc/fH)*S from the\n' ...
         'sector target, so the sweep strength IS the pivot distance. The distance\n' ...
         'coverage of one pilot is alpha''p*Lambda = H at every sweep rate, by their\n' ...
         '(69); what the sweep rate buys is angular traversals.\n']);
try, exportgraphics(gcf,'fig8_pencil.pdf','ContentType','vector');
catch, print(gcf,'-dpdf','fig8_pencil.pdf'); end
