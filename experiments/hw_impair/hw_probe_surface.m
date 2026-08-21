% ========================================================================
% (B) THE 2-D HARDWARE COST SURFACE:  TD count  x  TD precision
% ------------------------------------------------------------------------
% Two axes:
%   * PRECISION (bits)          -- relaxed by the hybrid TD-PS split (Contribution A)
%   * DISTINCT DELAY VALUES     -- reduced by the Kronecker factorization (B)
% ERRATUM: the second axis was originally labelled "TD element count". That is
% wrong -- Kronecker reduces distinct delay VALUES, not physical elements (sharing
% applies the delay before the fan-out and genuinely shares hardware; the
% Kronecker fine term does not). The precision axis is the load-bearing result. Generic sub-array sharing is included as the reference way to
% cut the count (NOT a reproduction of the OJ-COMS AoSA scheme, which also
% redesigns the training procedure -- see ddbs_beam_arch.m).
%
% Arms: 'ojcoms' (literal AoSA form of the OJ-COMS architecture -- their TTD on
% the sub-array-centre grid with a plain per-antenna DDBS phase shifter),
% 'shared' (the same TTD sharing but with the hybrid PS giving exact per-antenna
% fc phase -- the STRONGEST fair sharing baseline), and 'kron' (proposed).
%
% Pre-validated (Nt=256, ideal-TD, hybrid PS; ideal full-TTD gain = 0.975):
%   #TTD :   128     64      32
%   shared: 0.668  0.402   0.238
%   kron  : 0.975  0.973   0.941      <- 4x better at an equal DISTINCT-VALUE count
% and Kronecker at 32 distinct values still gives 0.938 at 12-bit TD (vs 14+ bits
% needed by a pure-TD design).
% ========================================================================
clear; clc;

BASELINE_DIR = fullfile('..','..','baseline_code','code_nf_distance_dependent_rainbow');
addpath(genpath(BASELINE_DIR));
addpath(genpath(fileparts(mfilename('fullpath'))));

Nt=256; fc=30e9; B=5e9; M=1024; c=3e8; d=(c/fc)/2; K=1;
theta1=-32.95; theta2=0.784; alpha1=-0.53; alpha2=0.58;
f = fc + B/M*((1:M)-1-(M-1)/2); km = 2*pi*f/c; nn=(-(Nt-1)/2:(Nt-1)/2)';
loc = cal_loc(fc,f,theta1,theta2,alpha1,alpha2,M,1);
th = squeeze(loc(:,1,1)); al = squeeze(loc(:,2,1));
gainof = @(w) mean(arrayfun(@(m) ...
    abs( (exp(-1j*km(m)*(nn*d*th(m)-(nn*d).^2*al(m)))/sqrt(Nt))' * w(:,1,m) ), 1:M));

P_list   = [1 2 4 8 16 32];
Btd_list = [Inf 14 13 12 11 10 9];

fprintf('ideal full-TTD, ideal bits : %.3f\n\n', ...
  gainof(ddbs_beam_arch(Nt,B,fc,M,d,theta1,theta2,alpha1,alpha2,K,'full',1,Inf,Inf)));

for arch = ["ojcoms","shared","kron"]
    fprintf('=== arch = %s ===\n', arch);
    fprintf('%6s %7s |', 'P', '#TTD');
    fprintf('%8s', string(Btd_list)); fprintf('   <- B_td\n');
    for P = P_list
        if (arch=="shared"||arch=="ojcoms") && P==1, continue; end
        [w,ntd] = ddbs_beam_arch(Nt,B,fc,M,d,theta1,theta2,alpha1,alpha2,K,char(arch),P,Inf,Inf);
        fprintf('%6d %7d |', P, ntd);
        for Btd = Btd_list
            w = ddbs_beam_arch(Nt,B,fc,M,d,theta1,theta2,alpha1,alpha2,K,char(arch),P,Btd,Inf);
            fprintf('%8.3f', gainof(w));
        end
        fprintf('\n');
    end
    fprintf('\n');
end

% ---- contour of the proposed (Kronecker + hybrid) surface ----
Pg = [1 2 4 8 16 32]; Bg = [16 14 13 12 11 10 9];
Zk = zeros(numel(Pg),numel(Bg)); Ntd = zeros(size(Pg));
for i=1:numel(Pg)
    for j=1:numel(Bg)
        [w,ntd] = ddbs_beam_arch(Nt,B,fc,M,d,theta1,theta2,alpha1,alpha2,K,'kron',Pg(i),Bg(j),Inf);
        Zk(i,j) = gainof(w); Ntd(i) = ntd;
    end
end
save('hw_probe_surface_results.mat','Pg','Bg','Zk','Ntd');
figure; imagesc(Bg, 1:numel(Pg), Zk); set(gca,'YTick',1:numel(Pg),'YTickLabel',compose('P=%d (%d TTD)',Pg(:),Ntd(:)));
set(gca,'XTick',fliplr(Bg)); xlabel('TD resolution (bits)'); colorbar; caxis([0 1]);
title('Proposed Kronecker + hybrid TD-PS: focusing gain over the hardware cost surface');
