import numpy as np, ddbs as D
c=3e8
p=D.params(); Nt,fc,B,M,d,f=p['Nt'],p['fc'],p['B'],p['M'],p['d'],p['f']; u=fc/f
nn=np.arange(Nt)-(Nt-1)/2; kc=2*np.pi*fc/c

# ---------- gain LUT  g(dtheta, dalpha) ----------
TT, TA, NTH, NAL = 0.025, 0.05, 2049, 1201
dal = np.linspace(-TA,TA,NAL)
NF  = 8192
base = np.exp(1j*kc*(-1)*np.outer(dal,(nn*d)**2))          # (NAL,Nt)
Fq   = np.fft.fft(base, n=NF, axis=1)/Nt                    # freq = -2pi k/NF vs index
psi  = -2*np.pi*np.arange(NF)/NF; psi=np.where(psi<-np.pi,psi+2*np.pi,psi)
dth_g= psi/(kc*d)                                           # theta grid (unsorted)
o    = np.argsort(dth_g); dth_g=dth_g[o]; Fq=np.abs(Fq[:,o])
sel  = np.abs(dth_g)<=TT; dth_g=dth_g[sel]; LUT=Fq[:,sel]   # (NAL, NTHsel)
print('LUT %s  dtheta step %.2e  dalpha step %.2e  peak %.4f'%(LUT.shape,dth_g[1]-dth_g[0],dal[1]-dal[0],LUT.max()))

def gmax(fl, thu, alu, chunk=64):
    K=fl.shape[2]; TH=fl[:,0,:].ravel(); AL=fl[:,1,:].ravel()
    out=np.zeros(len(thu))
    for i in range(0,len(thu),chunk):
        sl=slice(i,min(i+chunk,len(thu)))
        dt=TH[None,:]-thu[sl][:,None]; da=AL[None,:]-alu[sl][:,None]
        ok=(np.abs(dt)<TT)&(np.abs(da)<TA)
        ia=np.clip(np.searchsorted(dth_g,dt),0,len(dth_g)-1)
        ib=np.clip(((da+TA)/(2*TA)*(NAL-1)).round().astype(int),0,NAL-1)
        g=np.where(ok,LUT[ib,ia],0.0)
        out[sl]=g.max(axis=1)
    return out

# ---------- users ----------
rng=np.random.default_rng(23); Nu=3000
thu=-0.9+1.8*rng.random(Nu); r=5+195*rng.random(Nu); alu=(1-thu**2)/(2*r)

def cell(L,Ns,Ka,sw,arch=None):
    sec,rh,dbg=D.alg1(Nt,fc,B,M,d,L=L,Nsec=Ns,Kalpha=Ka,s_sweep=sw); K=len(sec)
    arch = arch or ('full' if L==1 else 'shared')
    W=np.zeros((Nt,K,M),complex)
    for t,s in enumerate(sec):
        W[:,t,:]=D.beam(Nt,B,fc,M,d,s['theta_t'],s['theta_p'],s['alpha_t'],s['alpha_p'],arch,L)
    a0=np.stack([u*s['alpha_p']+rh['alpha']*s['alpha_t'] for s in sec],axis=1)
    fl=D.focus_table(W,a0,Nt,fc,B,M,d,NF=4096,na=401)
    g=gmax(fl,thu,alu)
    return K,g,fl,dbg

SNR=100.0
def prate(g): return np.mean(np.log2(1+SNR*Nt*g**2))

Kr,gr,_,_=cell(1,12,1,1.0)                      # ungrouped reference, 12 pilots
ref=prate(gr); print('reference L=1 K=%d : mean g %.4f  P(g>0.9) %.3f  proxy rate %.3f\n'%(Kr,gr.mean(),(gr>0.9).mean(),ref))
meas={(1.0,1):98,(1.0,2):100,(1.0,4):100,(0.5,1):98,(0.5,2):100,(0.5,4):100,
      (0.25,1):94,(0.25,2):98,(0.25,4):100,(0.125,1):80,(0.125,2):92,(0.125,4):97}
print('%6s %3s %4s | %7s %8s %9s | %8s'%('s','Ka','K','mean g','P(g>.9)','pred %ref','meas %ref'))
for sw in (1.0,0.5,0.25,0.125):
    for Ka in (1,2,4):
        K,g,fl,dbg=cell(8,8,Ka,sw)
        print('%6.3f %3d %4d | %7.4f %8.3f %8.0f%% | %7d%%'%(sw,Ka,K,g.mean(),(g>0.9).mean(),100*prate(g)/ref,meas[(sw,Ka)]))
