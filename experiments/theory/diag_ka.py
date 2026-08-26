exec(open('collapse.py').read().split('rows=[]')[0])
import numpy as np
res={}
for Ka in (1,2,4):
    K,g,ntr=cell(8,8,Ka,0.125); res[Ka]=g
q=np.quantile(alu,np.linspace(0,1,9))
print('user alpha deciles (r uniform 5-200 m):')
print('%12s'%'alpha bin', ''.join('%9s'%('Ka=%d'%k) for k in (1,2,4)), '  median r')
for i in range(8):
    m=(alu>=q[i])&(alu<q[i+1]); 
    print('%6.4f-%6.4f'%(q[i],q[i+1]), ''.join('%9.3f'%res[k][m].mean() for k in (1,2,4)),
          '  %6.1f m  (n=%d)'%(np.median(r[m]),m.sum()))
print('\noverall mean g  ', ''.join('%9.3f'%res[k].mean() for k in (1,2,4)))
print('P(g>0.9)        ', ''.join('%9.3f'%(res[k]>0.9).mean() for k in (1,2,4)))
# where do the alpha samples sit at a fixed theta?
for Ka in (1,2):
    sec,rh,dbg=D.alg1(Nt,fc,B,M,d,L=8,Nsec=8,Kalpha=Ka,s_sweep=0.125)
    print('Ka=%d alpha_t values:'%Ka, ['%.4f'%s['alpha_t'] for s in sec[:4]],
          ' alpha locus min per pilot:', ['%.4f'%(min(u)*s['alpha_p']+rh['alpha']*s['alpha_t']) for s in sec[:4]])
