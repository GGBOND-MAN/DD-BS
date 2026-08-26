exec(open('collapse.py').read().split('rows=[]')[0])
import numpy as np
for Ka in (1,2):
    sec,rh,dbg=D.alg1(Nt,fc,B,M,d,L=8,Nsec=8,Kalpha=Ka,s_sweep=0.125); K=len(sec)
    W=np.zeros((Nt,K,M),complex)
    for t,s in enumerate(sec):
        W[:,t,:]=D.beam(Nt,B,fc,M,d,s['theta_t'],s['theta_p'],s['alpha_t'],s['alpha_p'],'shared',8)
    a0=np.stack([u*s['alpha_p']+rh['alpha']*s['alpha_t'] for s in sec],axis=1)
    fl=D.focus_table(W,a0,Nt,fc,B,M,d,NF=2048,na=201)
    for t in range(min(2,K)):
        des=a0[:,t]; got=fl[:,1,t]
        print('Ka=%d pilot %d: designed alpha %.4f..%.4f | reported %.4f..%.4f | at grid floor %.0f%% | theta span of floor entries %.3f'
              %(Ka,t,des.min(),des.max(),got.min(),got.max(),100*np.mean(got<=0.0026),
                np.ptp(fl[got<=0.0026,0,t]) if (got<=0.0026).any() else 0))
