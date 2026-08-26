exec(open('collapse.py').read().split('rows=[]')[0])
import numpy as np
H=0.0975; dal=5.939e-3; thr=H/(2*dal)
print('threshold N_sec*n_tr >= %.2f ; at s=1/8 n_tr=0.1977 -> N_sec >= %.1f\n'%(thr,thr/0.1977))
print('%6s %5s %3s %4s %10s %9s   %s'%('s','Nsec','Ka','K','Nsec*n_tr','P(g>.9)','predicted'))
for sw,Ns,Ka in [(0.125,24,1),(0.125,32,1),(0.125,42,1),(0.125,48,1),(0.125,8,2),(0.125,8,4),(0.0625,8,1),(0.0625,16,1)]:
    K,g,ntr=cell(8,Ns,Ka,sw)
    print('%6.4f %5d %3d %4d %10.2f %9.4f   %s'%(sw,Ns,Ka,K,Ns*ntr,(g>0.9).mean(),'PASS' if Ns*ntr>=thr else 'fail'),flush=True)
