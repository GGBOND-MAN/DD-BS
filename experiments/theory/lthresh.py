exec(open('collapse.py').read().split('def cell(')[0])
import numpy as np
# ablation_ps.m config: Nsec = min(4,max(1,L)),  Ka = round(12/Nsec)
def cell2(L, arch):
    Ns = min(4, max(1, L)); Ka = max(1, round(12/Ns))
    sec, rh, dbg = D.alg1(Nt, fc, B, M, d, L=L, Nsec=Ns, Kalpha=Ka, s_sweep=1.0)
    K = len(sec)
    W = np.zeros((Nt, K, M), complex)
    for t, s in enumerate(sec):
        W[:, t, :] = D.beam(Nt, B, fc, M, d, s['theta_t'], s['theta_p'],
                            s['alpha_t'], s['alpha_p'], arch, L)
    a0 = np.stack([u*s['alpha_p'] + rh['alpha']*s['alpha_t'] for s in sec], axis=1)
    fl = D.focus_table(W, a0, Nt, fc, B, M, d, NF=2048, na=201)
    g  = gmax(fl, thu, alu)
    return K, g

meas = {1:(6.524,6.538), 2:(5.481,5.383), 4:(2.313,6.490), 8:(0.982,4.613), 16:(1.212,3.128)}
print('%3s %4s | %-24s | %-24s | %8s %8s'
      %('L','K','their PS (ojcoms)','compensated (shared)','pred','meas'))
print('%3s %4s | %8s %8s %6s | %8s %8s %6s | %8s %8s'
      %('','','mean g','P(g>.9)','rate','mean g','P(g>.9)','rate','ratio','ratio'))
for L in (1,2,4,8,16):
    arch_t = 'full' if L==1 else 'ojcoms'
    arch_c = 'full' if L==1 else 'shared'
    K,gT = cell2(L, arch_t)
    _,gC = cell2(L, arch_c)
    rT,rC = meas[L]
    # crude rate proxy at the calibrated operating point (SNR*Nt ~ 92)
    pr = lambda g: np.mean(np.log2(1+92*g**2))
    print('%3d %4d | %8.4f %8.3f %6.3f | %8.4f %8.3f %6.3f | %8.2f %8.2f'
          %(L,K,gT.mean(),(gT>0.9).mean(),pr(gT),gC.mean(),(gC>0.9).mean(),pr(gC),
            pr(gC)/pr(gT), rC/rT), flush=True)
