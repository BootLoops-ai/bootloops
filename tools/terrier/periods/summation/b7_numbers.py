import math
def li_neg(d,x):
    eul={0:None,1:[1],2:[1,1],3:[1,4,1],4:[1,11,11,1],5:[1,26,66,26,1]}
    if d==0: return 1/(1-x)
    return x*sum(c*x**j for j,c in enumerate(eul[d]))/(1-x)**(d+1)
def tail(d,M,x):  # exact-ish direct sum
    s,m=0.0,M+1
    while True:
        t=(m**d)*x**m; s+=t; m+=1
        if t<1e-22*max(s,1e-300): break
    return s
for lbl,S in (("box center",0.746),("box worst corner",0.83)):
    x=S*S
    print(f"{lbl}: S={S} x={x:.4f} Li_-4={li_neg(4,x):.1f} tail_4(M=50)={tail(4,50,x):.3e} tail_4(M=100)={tail(4,100,x):.3e} tail_0(M=50)={tail(0,50,x):.3e}")
from math import comb
print("support pts C(M+6,6): M=50:",comb(56,6),"M=100:",comb(106,6))
