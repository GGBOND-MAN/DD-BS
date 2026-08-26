exec(open('collapse.py').read().split('rows=[]')[0])
import numpy as np
print('%6s %3s | %14s %14s'%('s','Ka','recalibrated','nominal (no recal)'))
for sw in (1.0,0.125):
  for Ka in (1,2,4):
    sec,rh,dbg=D.alg1(Nt,fc,B,M,d,L=8,Nsec=8,Kalpha=Ka,s_sweep=sw); K=len(sec)
    W=np.zeros((Nt,K,M),complex)
    for t,s in enumerate(sec):
        W[:,t,:]=D.beam(Nt,B,fc,M,d,s['theta_t'],s['theta_p'],s['alpha_t'],s['alpha_p'],'shared',8)
    a0=np.stack([u*s['alpha_p']+rh['alpha']*s['alpha_t'] for s in sec],axis=1)
    fl=D.focus_table(W,a0,Nt,fc,B,M,d,NF=2048,na=201)
    nom=np.zeros_like(fl); per=(2*u)[:,None]
    thn=np.stack([u*s['theta_p']+rh['theta']*s['theta_t'] for s in sec],axis=1)
    nom[:,0,:]=thn-per*np.round(thn/per); nom[:,1,:]=a0
    print('%6.3f %3d | %14.4f %14.4f'%(sw,Ka,(gmax(fl,thu,alu)>0.9).mean(),(gmax(nom,thu,alu)>0.9).mean()),flush=True)
