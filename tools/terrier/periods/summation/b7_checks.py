import math, random, itertools
random.seed(7)
# 1) multinomial Cauchy-Schwarz majorant with jet weights, n=6
def multinom(m,k):
    r=math.factorial(m)
    for ki in k: r//=math.factorial(ki)
    return r
def check_cs(m,y,alpha):
    d=sum(alpha); lhs=0.0
    for k in itertools.product(range(m+1),repeat=6):
        if sum(k)!=m: continue
        w=1.0
        for ki,ai in zip(k,alpha): w*=ki**ai
        lhs+=w*multinom(m,list(k))**2*math.prod(yi**ki for yi,ki in zip(y,k))
    rhs=(m**d)*(sum(math.sqrt(yi) for yi in y))**(2*m)
    return lhs<=rhs*(1+1e-12), lhs, rhs
ok=True
for trial in range(20):
    m=random.randint(1,5); y=[random.uniform(0.001,0.05) for _ in range(6)]
    alpha=[random.randint(0,2) for _ in range(6)]
    while sum(alpha)>4: alpha=[random.randint(0,2) for _ in range(6)]
    good,l,r=check_cs(m,y,alpha); ok&=good
print("CS+jet majorant 20 trials:", "PASS" if ok else "FAIL")
# 2) Eulerian closed forms Li_{-d}(x)=sum m^d x^m, d=0..5
eul={0:[1],1:[1],2:[1,1],3:[1,4,1],4:[1,11,11,1],5:[1,26,66,26,1]}
def li_neg(d,x):
    if d==0: return x/(1-x)  # careful: sum_{m>=1} x^m; m=0 term is 0 except d=0 where m^0=1 at m=0 -> add 1
    num=x*sum(c*x**j for j,c in enumerate(eul[d]))
    return num/(1-x)**(d+1)
ok2=True
for d in range(6):
    for x in (0.3,0.69,0.83):
        s=sum((m**d)*x**m for m in range(0,4000))
        cf=(1/(1-x)) if d==0 else li_neg(d,x)
        ok2 &= abs(s-cf)<1e-9*max(1,abs(cf))
print("Eulerian closed forms d=0..5:", "PASS" if ok2 else "FAIL")
# 3) tail corollary: sum_{m>M} m^d x^m <= (M+1)^d x^(M+1)/(1-x*e^(d/(M+1))) when x*e^(d/(M+1))<1
ok3=True
for d in range(6):
    for x in (0.3,0.69,0.83):
        for M in (20,50):
            q=x*math.exp(d/(M+1))
            if q>=1: continue
            tail=sum((m**d)*x**m for m in range(M+1,6000))
            bound=((M+1)**d)*x**(M+1)/(1-q)
            ok3 &= tail<=bound*(1+1e-12)
print("tail corollary bound:", "PASS" if ok3 else "FAIL")
# 4) |log(1+w)| <= -log(1-|w|), complex w
import cmath
ok4=True
for trial in range(2000):
    r=random.uniform(0,0.95); th=random.uniform(0,2*math.pi); w=r*cmath.exp(1j*th)
    ok4 &= abs(cmath.log(1+w))<=-math.log(1-r)+1e-14
print("log modulus inequality:", "PASS" if ok4 else "FAIL")
