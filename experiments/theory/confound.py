exec(open('lthresh.py').read().split('meas = {')[0].replace(
  "Ns = min(4, max(1, L)); Ka = max(1, round(12/Ns))","Ns, Ka = split"))
import numpy as np
pr = lambda g: np.mean(np.log2(1+92*g**2))
print('K=12 held fixed; only the (N_sec, K_alpha) SPLIT changes\n')
print('%3s %-10s | %14s %14s | %6s'%('L','split','their PS','compensated','ratio'))
for L,splits in [(2,[(2,6),(4,3),(3,4),(6,2)]), (4,[(4,3),(2,6)]), (8,[(4,3),(8,1)])]:
    for split in splits:
        globals()['split']=split
        _,gT = cell2(L,'ojcoms'); _,gC = cell2(L,'shared')
        print('%3d N_sec=%d,Ka=%-2d | %7.4f %6.3f | %7.4f %6.3f | %6.2f'
              %(L,split[0],split[1],gT.mean(),pr(gT),gC.mean(),pr(gC),pr(gC)/pr(gT)),flush=True)
