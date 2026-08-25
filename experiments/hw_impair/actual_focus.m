function fl = actual_focus(w,focus_loc,Nt,fc,B,M,d,K)
% Offline table of where each subcarrier's beam ACTUALLY focuses under the given
% (impaired) hardware. Two 1-D searches per (m,t): theta at the designed alpha,
% then alpha at that theta. NOTE the correlation uses w WITHOUT conjugation --
% self-tested against the ideal beams, where it reproduces the designed table.
if nargin == 0
    fprintf('actual_focus builds a focus table; call it from a script.\n');
    fprintf('Runnable scripts: rate_vs_snr, ablation_ps, sector_alloc, oracle_loc, ojcoms_baseline, phase_refine, kmin_fine, kspace_map, pilot_spacing_map, compare_ojcoms, hw_probe_*\n');
    return;
end
c=3e8; f = fc + B/M*((1:M)-1-(M-1)/2); km = 2*pi*f/c;
nn=(-(Nt-1)/2:(Nt-1)/2)'; N1=nn*d; N2=(nn*d).^2;
GS = linspace(-0.9,0.9,801); GA = linspace(1/400,1/10,81);
fl = zeros(M,2,K);
for t=1:K
    for m=1:M
        wm = w(:,t,m); a0 = focus_loc(m,2,t);
        g  = abs(exp(1j*km(m)*(GS.'*N1.' - a0*N2.'))*wm);
        [~,i] = max(g); th = GS(i);
        g2 = abs(exp(1j*km(m)*(th*N1.' - GA.'*N2.'))*wm);
        [~,j] = max(g2);
        fl(m,1,t) = th; fl(m,2,t) = GA(j);
    end
end
