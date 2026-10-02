"""Complete exact check for a root of multiplicity >=4, p=71, degree 17.

Parametrization: f'=17*(x^2-1)^3*(x^10+r4*x^8+...+r0),
P=f-f(1), f odd and monic. Enumerate all 71^5 parameter choices through
71^4 histograms. NumPy integer arithmetic only; no floating point or search.
"""
import json, time
import numpy as np
P=71
started=time.monotonic()

def value(q,x):
    y=0
    for a in reversed(q): y=(y*x+a)%P
    return y

def part(k):
    q=[0]*18
    for j,a in enumerate([-1,3,-3,1]):
        d=2*(j+k)+1
        q[d]=17*a*pow(d,-1,P)%P
    return q

parts=[part(k) for k in range(6)]
xs=list(range(2,P))
E=np.array([[(value(q,x)-sum(q))%P for x in xs] for q in parts],dtype=np.int64)
regular=np.flatnonzero(E[0])
exceptional=np.flatnonzero(E[0]==0)
inverse=np.array([pow(int(E[0,i]),-1,P) for i in regular],dtype=np.int64)
survivors=[]
solutions=[]
histogram_survivors=0

def split(q):
    q=q[:]
    roots=[]
    for x in range(P):
        while len(q)>1 and value(q,x)==0:
            z=[0]*(len(q)-1)
            z[-1]=q[-1]
            for j in range(len(z)-2,-1,-1): z[j]=(q[j+1]+x*z[j+1])%P
            assert (q[0]+x*z[0])%P==0
            roots.append(x)
            q=z
    return roots,q

for offset in range(0,P**4,8192):
    codes=np.arange(offset,min(offset+8192,P**4),dtype=np.int64)
    digits=np.column_stack([(codes//(P**j))%P for j in range(4)])
    base=(digits@E[1:5]+E[5])%P
    desired=(-base[:,regular]*inverse)%P
    rowoffset=np.arange(len(codes),dtype=np.int64)*P
    hist=np.bincount((desired+rowoffset[:,None]).ravel(),minlength=len(codes)*P).reshape(-1,P)
    if len(exceptional): hist+=(base[:,exceptional]==0).sum(axis=1)[:,None]
    ii,jj=np.nonzero(hist>=8)
    histogram_survivors+=len(ii)
    for i,r0 in zip(ii.tolist(),jj.tolist()):
        parameters=[r0]+digits[i].tolist()+[1]
        f=[sum(parameters[k]*parts[k][d] for k in range(6))%P for d in range(18)]
        c=sum(f)%P
        if not c: continue
        q=f[:]
        q[0]=-c%P
        roots,remainder=split(q)
        rec={'R_ascending':parameters,'c':c,'distinct_roots':int(hist[i,r0])+1,'total_rational_multiplicity':len(roots),'roots_with_multiplicity':roots,'remaining_factor_ascending':remainder}
        survivors.append(rec)
        if len(roots)==17:
            assert all(sum(pow(a,j,P) for a in roots)%P==0 for j in range(1,16,2))
            solutions.append(rec)
    if offset and offset//8192%512==0:
        print(json.dumps({'outer_done':int(offset),'outer_total':P**4,'survivors_nonzero_c':len(survivors),'elapsed_seconds':time.monotonic()-started}),flush=True)

result={'p':P,'normalization':'root 1 has multiplicity >=4','parameter_polynomials':P**5,'histogram_outer_quadruples':P**4,'last_parameter_direction_extra_zeros':[xs[i] for i in exceptional],'histogram_survivors_including_c_zero':histogram_survivors,'survivor_count':len(survivors),'survivors':survivors,'solutions':solutions,'elapsed_seconds':time.monotonic()-started}
with open('multiplicity4_result.json','w',encoding='utf8') as file: json.dump(result,file,indent=2)
print(json.dumps({k:v for k,v in result.items() if k!='survivors'},ensure_ascii=False),flush=True)
