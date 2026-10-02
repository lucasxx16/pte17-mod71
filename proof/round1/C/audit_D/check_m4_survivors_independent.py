"""Independently verify every D m=4 screen survivor via Hasse orders."""
from math import comb
from pathlib import Path
import json

HERE=Path(__file__).parent
P,N,M=71,17,4

def mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            c[i+j]=(c[i+j]+x*y)%P
    return c

def basis():
    t=[1]
    for _ in range(M-1):
        t=mul(t,[P-1,0,1])
    J=[]
    for k in range(6):
        product=mul(t,[0]*(2*k)+[1])
        J.append([0]+[N*x*pow(i+1,-1,P)%P for i,x in enumerate(product)])
    return J

def poly_for(R,J):
    f=[0]*18
    for k,r in enumerate(R):
        for i,x in enumerate(J[k]):
            f[i]=(f[i]+r*x)%P
    c=sum(f)%P
    f[0]=-c%P
    return f,c

def evalpoly(f,a):
    out=0
    for x in reversed(f):
        out=(out*a+x)%P
    return out

def root_order(f,a):
    for k in range(len(f)):
        h=sum(comb(i,k)*f[i]*pow(a,i-k,P) for i in range(k,len(f)))%P
        if h:
            return k
    raise AssertionError('zero polynomial')

def main():
    src=HERE.parent.parent/'D'/'multiplicity4_result.json'
    data=json.loads(src.read_text(encoding='utf-8'))
    J=basis()
    assert len(data['survivors'])==6534
    mismatches=[]
    total_hist={}
    distinct_hist={}
    for index,rec in enumerate(data['survivors']):
        f,c=poly_for(rec['R_ascending'],J)
        roots=[]
        for a in range(P):
            if evalpoly(f,a)==0:
                roots.extend([a]*root_order(f,a))
        if c!=rec['c'] or roots!=rec['roots_with_multiplicity'] or \
           len(set(roots))!=rec['distinct_roots'] or \
           len(roots)!=rec['total_rational_multiplicity']:
            mismatches.append(index)
        total_hist[len(roots)]=total_hist.get(len(roots),0)+1
        distinct_hist[len(set(roots))]=distinct_hist.get(len(set(roots)),0)+1
    result={'survivors_checked':len(data['survivors']),
            'mismatches':len(mismatches),'first_mismatch_indices':mismatches[:10],
            'rational_root_multiplicity_histogram':total_hist,
            'distinct_root_histogram':distinct_hist,
            'fully_split':sum(k==N for k in total_hist for _ in range(total_hist[k]))}
    (HERE/'check_m4_survivors_independent_result.json').write_text(
        json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    main()
