% ========================================================================
% (A) Reopen multipath: does relaxing the pilot budget K lift the K=3 ceiling?
% ------------------------------------------------------------------------
% Gate C found complex-OMP multipath separation plateaus at K=3 (flat in SNR).
% Hypothesis (improve-on-existing framing): spend MORE DDBS pilots -> lift the
% identifiability ceiling -> beat DFT/polar-sweeping multipath methods while
% keeping K far below the sweeping overhead. This sweeps K and measures whether
% detection rises toward full recovery.
%
% Pilots are added the baseline way: delay_polar_2d generates K interleaved
% strips by shifting the intercept theta_t by -2/K per pilot (baseline Alg. 1).
%
% FINDING (pre-validated in Python at full Nt=256, M=1024): detection does NOT
% rise with K -- it is flat (~0.49) then declines (K=10 -> 0.37). The extra DDBS
% pilots are strongly CORRELATED (the same multi-strip pattern, slightly shifted),
% so they add little sensing diversity; near-collinear path signatures stay
% unresolved, and denser strips even raise dictionary mutual coherence. Hence
% relaxing K does NOT rescue multipath -- the limit is the DDBS measurement
% operator's structure, not the pilot count.  => hypothesis (A) NEGATIVE.
% ========================================================================
clear; clc;

BASELINE_DIR = fullfile('..','..','baseline_code','code_nf_distance_dependent_rainbow');
addpath(genpath(BASELINE_DIR));
addpath(genpath(fileparts(mfilename('fullpath'))));

p = derisk_params();

Gtheta = 81;  Gr = 12;              % dictionary grid
K_list = [3 6 10];
L_true = 3; kappa_dB = -3; minsep = 0.08;
SNR_list = [5 15 20];               % also checks the flat-in-SNR signature
N_iter = 50;  tol_sin = 0.03; tol_rrel = 0.25;
rng(1);

th_grid = linspace(p.theta_min, p.theta_max, Gtheta);   % physical angle
r_grid  = linspace(p.Rmin, p.Rmax, Gr);
[TT, RR] = ndgrid(th_grid, r_grid);
gsin = sin(TT(:)).'; gr = RR(:).'; G = numel(gsin);

fprintf('K-sweep | L=%d, kappa=%d dB, grid=%d atoms, N_iter=%d\n', ...
        L_true, kappa_dB, G, N_iter);
det = zeros(numel(K_list), numel(SNR_list));

for iK = 1:numel(K_list)
    K = K_list(iK);
    w = delay_polar_2d(p.Nt,p.B,p.fc,p.M,p.d, p.theta1,p.theta2,p.alpha1,p.alpha2, K);

    % complex dictionary for THIS K  (atom length K*M)
    D = zeros(K*p.M, G);
    for gi = 1:G
        H = near_field_channel(p.Nt,p.d,p.fc,p.B,p.M, gr(gi), asin(gsin(gi)));  % M x Nt
        atom = zeros(K*p.M,1);
        for t = 1:K
            wt = squeeze(w(:,t,:)).';
            atom((t-1)*p.M+(1:p.M)) = sum(H.*wt, 2);
        end
        D(:,gi) = atom;
    end
    D = D ./ vecnorm(D);

    for iS = 1:numel(SNR_list)
        SNR = 10^(SNR_list(iS)/10); ndet = 0; ntot = 0;
        for it = 1:N_iter
            s = gen_multipath_scenario(p, L_true, kappa_dB, minsep);
            [h,~] = near_field_channel_multipath(p.Nt,p.d,p.fc,p.B,p.M,s.r,s.theta,s.g);
            y = zeros(K*p.M,1);
            for t = 1:K
                wt = squeeze(w(:,t,:)).';
                y((t-1)*p.M+(1:p.M)) = sum(h.*wt, 2);
            end
            Ps = mean(abs(y).^2);
            y = y + sqrt(Ps/SNR/2)*(randn(size(y))+1j*randn(size(y)));
            idx = omp_complex(D, y, L_true);
            es = gsin(idx); er = gr(idx); used = false(1,numel(idx));
            for l = 1:L_true
                ntot = ntot+1; st = s.theta_sin(l); rt = s.r(l);
                best = inf; bi = 0;
                for e = 1:numel(idx)
                    if used(e), continue; end
                    dd = abs(es(e)-st) + abs(er(e)-rt)/rt;
                    if dd < best, best = dd; bi = e; end
                end
                if bi>0
                    used(bi)=true;
                    if abs(es(bi)-st)<=tol_sin && abs(er(bi)-rt)/rt<=tol_rrel, ndet=ndet+1; end
                end
            end
        end
        det(iK,iS) = ndet/ntot;
        fprintf('K=%2d, SNR=%2d dB : detection %.2f  (pilots %d vs sweeping ~%d)\n', ...
                K, SNR_list(iS), det(iK,iS), K, p.Nt*10);
    end
end

save('derisk_Ksweep_results.mat','K_list','SNR_list','det');

figure; hold on; box on; grid on;
for iS = 1:numel(SNR_list)
    plot(K_list, det(:,iS), '-o','LineWidth',1.6,'DisplayName',sprintf('SNR=%d dB',SNR_list(iS)));
end
yline(1,'--','full recovery (goal)'); ylim([0 1]);
xlabel('Number of DDBS pilots K'); ylabel('Detection probability (of L=3 paths)');
title('(A) Does relaxing K lift the multipath ceiling?  -- it does not');
legend('Location','southwest');

% ------------------------------------------------------------------------
function idx = omp_complex(D, y, K)
r = y; idx = [];
for kk = 1:K
    c = abs(D'*r); c(idx) = -inf; [~,j] = max(c); idx = [idx, j]; %#ok<AGROW>
    Ds = D(:,idx); x = Ds \ y; r = y - Ds*x;
end
end
