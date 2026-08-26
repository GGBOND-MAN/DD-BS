import numpy as np, ddbs as D
c=3e8
p=D.params(); Nt,fc,B,M,d,f=p['Nt'],p['fc'],p['B'],p['M'],p['d'],p['f']; u=fc/f
for L,Ns,sw in [(1,8,1.0),(8,8,1.0),(8,8,0.25),(8,8,0.125)]:
    sec,rh,dbg=D.alg1(Nt,fc,B,M,d,L=L,Nsec=Ns,Kalpha=1,s_sweep=sw); K=len(sec)
    arch='full' if L==1 else 'shared'
    W=np.zeros((Nt,K,M),complex)
    for t,s in enumerate(sec):
        W[:,t,:]=D.beam(Nt,B,fc,M,d,s['theta_t'],s['theta_p'],s['alpha_t'],s['alpha_p'],arch,L)
    a0=np.stack([u*s['alpha_p']+rh['alpha']*s['alpha_t'] for s in sec],axis=1)
    fl=D.focus_table(W,a0,Nt,fc,B,M,d,NF=4096,na=401)
    per=(2*u)[:,None]
    thn=np.stack([u*s['theta_p']+rh['theta']*s['theta_t'] for s in sec],axis=1)
    near=thn-per*np.round((thn-fl[:,0,:])/per)
    e=np.abs(near-fl[:,0,:]); ea=np.abs(fl[:,1,:]-a0)
    print('L=%d s=%-5.3f arch=%-6s |dtheta| med %.2e p95 %.2e | |dalpha| med %.2e p95 %.2e | alias order j: %s'
          %(L,sw,arch,np.median(e),np.percentile(e,95),np.median(ea),np.percentile(ea,95),
            np.unique(np.round((thn[:,0]-fl[:,0,0])/per[:,0]).astype(int))[:8]))
