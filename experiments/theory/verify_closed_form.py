import numpy as np, ddbs as D
p = D.params(); Nt,fc,B,M,d,f = p['Nt'],p['fc'],p['B'],p['M'],p['d'],p['f']
sec, rh, dbg = D.alg1(Nt,fc,B,M,d,L=1,Nsec=1,Kalpha=1,s_sweep=1.0)
s = sec[0]; u = fc/f
W = D.beam(Nt,B,fc,M,d,s['theta_t'],s['theta_p'],s['alpha_t'],s['alpha_p'],'full',1)[:,None,:]
a0 = (u*s['alpha_p'] + rh['alpha']*s['alpha_t'])[:,None]
th_nom = u*s['theta_p'] + rh['theta']*s['theta_t']
fl = D.focus_table(W, a0, Nt,fc,B,M,d)
# alias period at subcarrier m is 2*u(m)
per = 2*u
res = th_nom - per*np.round((th_nom-fl[:,0,0])/per)     # nearest alias to measurement
sel = np.abs(fl[:,0,0])<0.85
print('theta: max|err| over %d in-range subcarriers = %.2e'%(sel.sum(), np.abs(res[sel]-fl[sel,0,0]).max()))
print('alpha: max|err| = %.2e   (grid step %.2e)'%(np.abs(fl[sel,1,0]-a0[sel,0]).max(), (0.1-1/400)/400))
