import numpy as np

c = 3e8
def params(Nt=256, fc=30e9, B=5e9, M=1024):
    d = (c/fc)/2
    f = fc + B/M*(np.arange(1,M+1)-1-(M-1)/2)
    return dict(Nt=Nt, fc=fc, B=B, M=M, d=d, f=f)

def rho(Nt, L):
    G = Nt//L
    n  = np.arange(Nt)
    nh = n - (Nt-1)/2
    nhs = np.floor(n/L) - (G-1)/2
    rt = L*np.sum(nh*nhs)/np.sum(nh**2)
    ra = L**2*np.sum(nh**2*nhs**2)/np.sum(nh**4)
    return rt, ra

def alg1(Nt, fc, B, M, d, L, thlim=(-1,1), alim=(0.0025,0.1), gamma=0.9,
         Nsec=None, Kalpha=1, s_sweep=1.0):
    fL, fH = fc-B/2, fc+B/2
    rt, ra = rho(Nt, L)
    U  = (alim[1]-alim[0])/(fc/fL - fc/fH)
    at = (alim[1] - (fc/fL)*U)/ra
    S  = s_sweep*1.76*gamma*fL*M/(Nt*B)
    th_t = np.zeros(Nsec); th_p = np.zeros(Nsec)
    for s in range(Nsec):
        thM = thlim[1] - s*(thlim[1]-thlim[0])/Nsec
        th_t[s] = (thM - (fc/fH)*S)/rt
        th_p[s] = S - (2/L)*np.round((L/2)*(S-28.40))
    alpha_p = U
    if Kalpha > 1:
        spread = abs(U)*0.15
        k = np.arange(Kalpha)
        off = np.ceil(k/2)*(-1.0)**k
        ats = at + off*(spread/max(1, Kalpha//2))
    else:
        ats = np.array([at])
    sec = [dict(theta_t=th_t[s], theta_p=th_p[s], alpha_t=ats[ka], alpha_p=alpha_p, sector=s)
           for s in range(Nsec) for ka in range(Kalpha)]
    return sec, dict(theta=rt, alpha=ra), dict(S=S, U=U, at=at, theta_p=th_p[0])

def beam(Nt,B,fc,M,d,theta1,theta2,alpha1,alpha2,arch,P):
    """one pilot (K=1), returns w (Nt,M).  Mirrors ddbs_beam_arch."""
    f  = fc + B/M*(np.arange(1,M+1)-1-(M-1)/2)
    kc = 2*np.pi*fc/c
    nn = (np.arange(Nt) - (Nt-1)/2)
    G  = Nt//P
    tau = (nn*d*theta1 - (nn*d)**2*alpha1)/c
    if arch == 'full':
        tg = tau
    elif arch in ('shared','ojcoms'):
        tg = np.repeat(tau.reshape(G,P).mean(axis=1), P)
    else:
        raise ValueError(arch)
    if arch == 'ojcoms':
        ps = np.mod(kc*(nn*d*theta2 - (nn*d)**2*alpha2), 2*np.pi)
        ph = 2*np.pi*f[None,:]*tg[:,None] + ps[:,None]
    else:
        ps = np.mod(kc*(nn*d*theta2 - (nn*d)**2*alpha2) + 2*np.pi*fc*tau, 2*np.pi)
        ph = 2*np.pi*(f[None,:]-fc)*tg[:,None] + ps[:,None]
    return np.exp(-1j*ph)/np.sqrt(Nt)

def focus_table(W, a0, Nt, fc, B, M, d, NF=8192, alim=(1/400,1/10), na=401):
    """W: (Nt,K,M).  a0: (M,K) designed alpha.  Returns fl (M,2,K)."""
    f  = fc + B/M*(np.arange(1,M+1)-1-(M-1)/2)
    km = 2*np.pi*f/c
    nn = (np.arange(Nt) - (Nt-1)/2)
    N1 = nn*d; N2 = (nn*d)**2
    K  = W.shape[1]
    GA = np.linspace(alim[0], alim[1], na)
    fl = np.zeros((M,2,K))
    psi = -2*np.pi*np.arange(NF)/NF
    psi = np.where(psi < -np.pi, psi+2*np.pi, psi)          # wrap to (-pi,pi]
    for m in range(M):
        Wm = W[:,:,m]                                        # (Nt,K)
        # --- theta search at designed alpha ---
        v = np.exp(-1j*km[m]*a0[m,:][None,:]*N2[:,None])*Wm  # (Nt,K)
        v = v*np.exp(-1j*psi_shift(km[m],d,nn))[:,None] if False else v
        F = np.fft.fft(v, n=NF, axis=0)                      # sum_idx v e^{-2pi i k idx/NF}
        th_grid = psi/(km[m]*d)
        ok = np.abs(th_grid) <= 0.9
        Fm = np.where(ok[:,None], np.abs(F), -1.0)
        i  = np.argmax(Fm, axis=0)
        th = th_grid[i]
        fl[m,0,:] = th
        # --- alpha search at that theta ---
        E2 = np.exp(-1j*km[m]*GA[:,None]*N2[None,:])         # (na,Nt)
        V  = Wm*np.exp(1j*km[m]*th[None,:]*N1[:,None])       # (Nt,K)
        g2 = np.abs(E2@V)                                    # (na,K)
        fl[m,1,:] = GA[np.argmax(g2, axis=0)]
    return fl

def psi_shift(km,d,nn):  # unused
    return np.zeros_like(nn)
