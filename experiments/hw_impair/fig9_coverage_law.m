% ========================================================================
% FIGURE 9 -- the coverage law, end-to-end (THEORY sec 40-41)
% ------------------------------------------------------------------------
% Panel (a) is the paper figure: measured rate against the single derived
% quantity  N_sec * s,  with the derived threshold
%
%       C = H * d * Nt^3 * xi_H^2 / (3.52 * gamma * M) = 2.96
%
% drawn as a vertical line and the >=95% criterion as a horizontal one. 28 cells
% spanning four sweep rates and seven sector counts fall on the correct side of
% BOTH lines. Nothing in C is fitted.
%
% Panel (b) is the evidence for using the smooth form rather than the
% lattice-quantised one (sec 41.2). The same 28 cells against N_sec * n_tr with
% its own derived threshold H*d*Nt^2/4 = 7.99 misclassify three -- all
% conservative. Keep panel (b) for the appendix or the response letter; the main
% text needs only (a).
%
% DATA: kmin_theory.m, which now saves kmin_theory_results.mat. Run it first.
% Do NOT rebuild this figure from experiments/theory/ -- the offline coverage
% statistic there is systematically pessimistic at small s (sec 41.4) and would
% put three of these cells on the wrong side of the line.
% ========================================================================
clear; clc;

if exist('kmin_theory_results.mat','file') ~= 2
    error(['kmin_theory_results.mat not found. Run kmin_theory.m first ' ...
           '(it saves the 28 cells this figure plots).']);
end
S = load('kmin_theory_results.mat');
R = S.RES;                      % [s Nsec K n_tr Nsec*s Nsec*n_tr rate pct]
sv = R(:,1); Ns = R(:,2); xs = R(:,5); xn = R(:,6); pct = R(:,8);

svals = unique(sv,'stable');
mk = {'o','s','^','d','v','>','<'};
cl = lines(numel(svals));
CRIT = 95;

figure('Position',[100 100 980 400]);

for panel = 1:2
    subplot(1,2,panel); hold on; box on; grid on;
    if panel==1, x = xs; thr = S.C_pred;  xl = 'N_{sec} \cdot s';
    else,        x = xn; thr = S.T;       xl = 'N_{sec} \cdot n_{tr}';
    end
    yl = [min(pct)-3, max(pct)+2];
    % shade the region the derived law calls "fail"
    patch([min(x)*0.7 thr thr min(x)*0.7],[yl(1) yl(1) yl(2) yl(2)], ...
          [0.93 0.93 0.93],'EdgeColor','none','HandleVisibility','off');
    plot([thr thr], yl, 'k--','LineWidth',1.5,'DisplayName', ...
         sprintf('derived threshold = %.2f', thr));
    plot([min(x)*0.7 max(x)*1.4],[CRIT CRIT],'k:','LineWidth',1.2, ...
         'DisplayName','95% of reference');
    for i = 1:numel(svals)
        k = sv == svals(i);
        plot(x(k), pct(k), mk{min(i,numel(mk))}, 'Color', cl(i,:), ...
             'MarkerFaceColor', cl(i,:), 'MarkerSize', 7, 'LineStyle','none', ...
             'DisplayName', sprintf('s = %g', svals(i)));
    end
    wrong = (x >= thr) ~= (pct >= CRIT);
    if any(wrong)
        plot(x(wrong), pct(wrong), 'ro','MarkerSize',13,'LineWidth',1.6, ...
             'DisplayName', sprintf('misclassified (%d)', sum(wrong)));
    end
    set(gca,'XScale','log'); xlim([min(x)*0.7 max(x)*1.4]); ylim(yl);
    xlabel(xl); ylabel('rate, % of ungrouped reference');
    title(sprintf('(%c)  %d/%d correct', 'a'+panel-1, numel(pct)-sum(wrong), numel(pct)));
    legend('Location','southeast','FontSize',8);
end

sgtitle(sprintf('N_t=%d, L=%d (%d TTDs both stages), B=%.0f GHz, N_{iter}=%d', ...
        S.Nt, S.L, S.Nt/S.L, S.B/1e9, S.N_iter));

w1 = sum((xs >= S.C_pred) ~= (pct >= CRIT));
w2 = sum((xn >= S.T)      ~= (pct >= CRIT));
fprintf('N_sec*s   vs C = %.3f : %d/%d correct\n', S.C_pred, numel(pct)-w1, numel(pct));
fprintf('N_sec*n_tr vs T = %.3f : %d/%d correct\n', S.T,      numel(pct)-w2, numel(pct));
fprintf('highest failing cell: %s = %.2f | lowest passing: %s = %.2f\n', ...
        'N_sec*s', max(xs(pct<CRIT)), 'N_sec*s', min(xs(pct>=CRIT)));
fprintf(['\nCaption to use: 28 configurations, four sweep rates x seven sector\n' ...
         'counts, plotted against the single derived quantity N_sec*s. The dashed\n' ...
         'line is the closed-form threshold C = H d Nt^3 xi_H^2 / (3.52 gamma M),\n' ...
         'which contains no fitted parameter. Report the margin: the measured\n' ...
         'bracket is (2, 3], so C is resolved only to a factor 1.5.\n']);
try, exportgraphics(gcf,'fig9_coverage_law.pdf','ContentType','vector');
catch, print(gcf,'-dpdf','fig9_coverage_law.pdf'); end
