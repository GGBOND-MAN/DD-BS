% ========================================================================
% FIGURE 2 -- the primary ablation: one term in the phase shifter (THEORY sec 44)
% ------------------------------------------------------------------------
% DATA: ablation_split.m, NOT ablation_ps.m. sec 44.6 found three defects in the
% latter's L = 2 row -- it used the split (2,6) where OJ-COMS publish (4,3), it
% therefore ran on Table 3 while every other L ran on the analytic branch, and
% its paired test was diluted by independent noise draws. ablation_split.m holds
% the split at (4,3), forces the analytic branch at every L, and re-seeds the
% noise so the two arms are compared on identical realisations.
%
% (a) rate vs sub-array size, both arms, with their Fig. 8 read off for
%     reference. HOW FAITHFUL THE REPRODUCTION IS, stated plainly rather than
%     claimed: within 0.1 at L = 1 and 4, about 0.45 LOW at L = 2, and 0.4-0.7
%     HIGH at L = 8 and 16. Do NOT repeat the earlier "tracks their Fig. 8 at
%     five points": that was read off the (2,6) split, which is not the one they
%     publish, and even there L = 8 and 16 were 0.4-0.8 high.
%
% (b) the PAIRED difference r_compensated - r_theirs, in bit/s/Hz, with its own
%     CI. This is the statistically exact quantity: both arms saw the same
%     channels AND the same noise, so the difference is a matched pair. The gain
%     RATIO is printed on each bar because that is what the text quotes, but the
%     error bar belongs to the difference. L = 1 must sit at exactly zero.
%
% (c) the L = 2 control that settled sec 44: gain against split, for both
%     parameter sources. The split moves the gain by 0.03; switching Table 3 for
%     the analytic reconstruction moves it from 1.43 to 0.83. Table 3 is tuned
%     for the UNCOMPENSATED phase shifter -- which is what its authors designed
%     it for -- so the compensated architecture wants its own parameters. That
%     is limitation #2 in the paper, and this panel is its evidence.
%
% Run ablation_split.m first; it saves ablation_split_results.mat.
% ========================================================================
clear; clc;

if exist('ablation_split_results.mat','file') ~= 2
    error('ablation_split_results.mat not found. Run ablation_split.m first.');
end
S  = load('ablation_split_results.mat');
P2 = S.P2;                 % [L N_TTD mT cT mC cC dm dci g gci]
P1 = S.P1;                 % [analytic? Nsec Ka mT cT mC cC dm dci g gci table3?]
L  = P2(:,1); mT = P2(:,3); cT = P2(:,4); mC = P2(:,5); cC = P2(:,6);
dm = P2(:,7); dci = P2(:,8); gg = P2(:,9);

fig8L = [1 2 4 8 16];  fig8R = [6.5 5.0 2.2 0.6 0.4];   % read off their Fig. 8

figure('Position',[60 60 1180 380]);

% ---------- (a) rate vs L ----------
subplot(1,3,1); hold on; box on; grid on;
errorbar(L, mT, cT, '-s','LineWidth',1.8,'MarkerSize',7,'MarkerFaceColor','auto', ...
         'DisplayName','AoSA TTD-PS, their (14)');
errorbar(L, mC, cC, '-o','LineWidth',1.8,'MarkerSize',7,'MarkerFaceColor','auto', ...
         'DisplayName','+ f_c compensation (proposed)');
plot(fig8L, fig8R, 'k^','MarkerSize',8,'LineWidth',1.2, ...
     'DisplayName','their Fig. 8 (read off)');
set(gca,'XScale','log','XTick',fig8L,'XTickLabel',{'1','2','4','8','16'});
xlim([0.85 19]); xlabel('sub-array size L   (N_{TTD} = N_t/L)');
ylabel('average rate (bit/s/Hz)');
title('(a)  rate vs sharing');
legend('Location','southwest','FontSize',8);

% ---------- (b) paired difference ----------
subplot(1,3,2); hold on; box on; grid on;
bar(1:numel(L), dm, 0.55, 'FaceColor',[0.30 0.55 0.80],'EdgeColor','none');
errorbar(1:numel(L), dm, dci, 'k.','LineWidth',1.3,'CapSize',9);
plot([0.4 numel(L)+0.6],[0 0],'k-','LineWidth',1.0);
for i = 1:numel(L)
    yoff = 0.22; if dm(i) < 0, yoff = -0.35; end
    text(i, dm(i)+dci(i)+yoff, sprintf('%.2f\\times', gg(i)), ...
         'HorizontalAlignment','center','FontSize',9,'FontWeight','bold');
end
set(gca,'XTick',1:numel(L),'XTickLabel',arrayfun(@(v)sprintf('%d',v),L,'UniformOutput',false));
xlim([0.4 numel(L)+0.6]); ylim([min(0,min(dm)-0.8) max(dm)+1.1]);
xlabel('sub-array size L'); ylabel('paired difference  r_{comp} - r_{theirs}  (bit/s/Hz)');
title('(b)  paired gain, identical channels and noise');

% ---------- (c) the L = 2 control ----------
subplot(1,3,3); hold on; box on; grid on;
A = P1(P1(:,1)==1, :);                       % analytic everywhere
T = P1(P1(:,1)==0 & P1(:,12)==1, :);         % rows that actually used Table 3
xa = 1:size(A,1);
bar(xa, A(:,9), 0.55, 'FaceColor',[0.30 0.55 0.80],'EdgeColor','none', ...
    'DisplayName','analytic (61)-(63)');
[tf, loc] = ismember(T(:,2), A(:,2));        % match Table-3 rows to their split
xt = xa(loc(tf));
bar(xt, T(tf,9), 0.28, 'FaceColor',[0.85 0.45 0.25],'EdgeColor','none', ...
    'DisplayName','their Table 3');
plot([0.4 size(A,1)+0.6],[1 1],'k--','LineWidth',1.2,'DisplayName','no gain');
set(gca,'XTick',xa,'XTickLabel',arrayfun(@(a,b)sprintf('(%d,%d)',a,b),A(:,2),A(:,3),'UniformOutput',false));
xlim([0.4 size(A,1)+0.6]); ylim([0 max(A(:,9))*1.25]);
xlabel('(N_{sec}, K_\alpha) split at L = 2,  K = 12'); ylabel('gain,  r_{comp} / r_{theirs}');
title('(c)  L = 2: the source decides, not the split');
legend('Location','northwest','FontSize',8);

sgtitle(sprintf(['N_t=%d, K=%d pilots at every L, split (4,3), analytic parameters, ' ...
                 'SNR=%d dB, N_{iter}=%d'], S.Nt, S.KTOT, S.SNR_dB, S.N_iter));

fprintf('reproduction of their Fig. 8 (their arm vs read-off):\n');
for i = 1:numel(L)
    k = find(fig8L == L(i), 1);
    if ~isempty(k)
        fprintf('   L=%-2d  ours %.3f  vs  theirs ~%.1f   (%+.2f)\n', L(i), mT(i), fig8R(k), mT(i)-fig8R(k));
    end
end
fprintf(['\nCaption: identical pilots, identical lookup, identical channels and noise;\n' ...
         'only the phase-shifter term differs. The L=1 bar in (b) is exactly zero --\n' ...
         'without sharing the f_c term is a global constant and the two arms are the\n' ...
         'same beam, so that bar is a hard check on the experiment, not a result.\n' ...
         'State the reproduction honestly: within 0.1 at L=1 and 4, ~0.45 low at\n' ...
         'L=2, and 0.4-0.7 high at L=8 and 16. Quote the gains as an ablation ratio\n' ...
         'at a FIXED parameter set,\n' ...
         'never as a win over their published system (sec 44.5).\n']);
try, exportgraphics(gcf,'fig2_ablation.pdf','ContentType','vector');
catch, print(gcf,'-dpdf','fig2_ablation.pdf'); end
