exec(open('cover.py').read().split('# ---------- users ----------')[0])
import numpy as np, pickle
rng=np.random.default_rng(23); Nu=1500
thu=-0.9+1.8*rng.random(Nu); r=5+195*rng.random(Nu); alu=(1-thu**2)/(2*r)
fL,fH=fc-B/2,fc+B/2; Lam=fc/fL-fc/fH; H=0.0975

def cell(L,Ns,Ka,sw):
    sec,rh,dbg=D.alg1(Nt,fc,B,M,d,L=L,Nsec=Ns,Kalpha=Ka,s_sweep=sw); K=len(sec)
    W=np.zeros((Nt,K,M),complex)
    for t,s in enumerate(sec):
        W[:,t,:]=D.beam(Nt,B,fc,M,d,s['theta_t'],s['theta_p'],s['alpha_t'],s['alpha_p'],'shared',L)
    a0=np.stack([u*s['alpha_p']+rh['alpha']*s['alpha_t'] for s in sec],axis=1)
    fl=D.focus_table(W,a0,Nt,fc,B,M,d,NF=2048,na=201)
    g=gmax(fl,thu,alu)
    # exact n_tr: lattice-quantised traversal count of sector 1
    tp=dbg['theta_p']; piv=rh['theta']*sec[0]['theta_t']
    q=tp-2*np.round((tp+piv)/2); ntr=abs(q)*Lam/2
    return K,g,ntr

rows=[]
for sw in (1.0,0.5,0.25,0.125):
    for Ns in (8,12,16):
        for Ka in (1,2,3,4,6):
            K,g,ntr=cell(8,Ns,Ka,sw)
            rows.append((sw,Ns,Ka,K,ntr,float(g.mean()),float((g>0.9).mean())))
            print('%6.3f Ns=%2d Ka=%d K=%3d  n_tr=%.4f  K*s=%6.2f  K*n_tr=%6.2f  meang %.4f  P(g>.9) %.4f'
                  %(sw,Ns,Ka,K,ntr,K*sw,K*ntr,g.mean(),(g>0.9).mean()),flush=True)
pickle.dump(rows,open('collapse.pkl','wb'))
