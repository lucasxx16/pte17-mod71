"""E independent full m>=4 scan: eliminate r4, not D's r0; derivative multiplicities."""
from pathlib import Path
import json,time,math,collections
import numpy as np
HERE=Path(__file__).resolve().parent
started=time.monotonic(); p=71
reference=json.loads((HERE/'multiplicity4_result.json').read_text())
expected={tuple(x['R_ascending'][:5]):x for x in reference['survivors']}
assert len(expected)==6534

def mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b): c[i+j]=(c[i+j]+x*y)%p
    return c
def eval_poly(a,x):
    value=0
    for c in a[::-1]: value=(value*x+c)%p
    return value
def integral_of_derivative(r):
    even=[0]*(2*len(r)-1)
    for i,c in enumerate(r): even[2*i]=c
    derivative=mul(mul(mul([-1,0,1],[-1,0,1]),[-1,0,1]),even)
    out=[0]+[17*c*pow(i+1,-1,p)%p for i,c in enumerate(derivative)]
    return out+[0]*(18-len(out))
basis=[]
for k in range(6):
    r=[0]*(k+1);r[k]=1
    basis.append(integral_of_derivative(r))
at1=np.array([eval_poly(a,1) for a in basis],dtype=np.int64)
evaluation=np.array([[eval_poly(a,x) for x in range(2,71)] for a in basis],dtype=np.int64)
differences=(evaluation-at1[:,None])%p
direction=differences[4]; nonzero=np.where(direction!=0)[0]; zeros=np.where(direction==0)[0]
inverses=np.array([pow(int(direction[j]),-1,p) for j in nonzero],dtype=np.int64)
seen=set(); all_pass=0; nonzero_constant=0; root_mass_hist=collections.Counter(); distinct_hist=collections.Counter()
for block in range(0,p**4,4096):
    ids=np.arange(block,min(block+4096,p**4),dtype=np.int64)
    digits=np.empty((len(ids),4),dtype=np.int64)
    remainder=ids.copy()
    for k in range(4): digits[:,k]=remainder%p; remainder//=p
    residual=(digits@differences[:4]+differences[5])%p
    needed=(-residual[:,nonzero]*inverses)%p
    labels=needed+p*np.arange(len(ids),dtype=np.int64)[:,None]
    hist=np.bincount(labels.flatten(),minlength=p*len(ids)).reshape(len(ids),p)
    if len(zeros): hist+=(residual[:,zeros]==0).sum(axis=1)[:,None]
    row,col=np.where(hist>=8)
    all_pass+=len(row)
    for i,r4 in zip(row.tolist(),col.tolist()):
        r=digits[i].tolist()+[r4]
        constant=int((sum(r[k]*at1[k] for k in range(5))+at1[5])%p)
        if constant==0: continue
        nonzero_constant+=1; key=tuple(r)
        assert key in expected, ('new candidate',key)
        assert key not in seen
        seen.add(key)
        # Rebuild by polynomial products and exact formal integration.
        f=integral_of_derivative(r+[1]); f[0]=-eval_poly(f,1)%p
        assert f[-1]==1 and f[0]==-constant%p
        # Multiplicities from independent formal derivatives, with no factor division.
        derivatives=[f]
        while len(derivatives[-1])>1:
            previous=derivatives[-1]
            derivatives.append([(j*previous[j])%p for j in range(1,len(previous))])
        roots=[]
        for a in range(p):
            if eval_poly(f,a): continue
            m=next(k for k in range(1,18) if eval_poly(derivatives[k],a))
            roots.extend([a]*m)
        assert roots==expected[key]['roots_with_multiplicity']
        assert len(set(roots))==int(hist[i,r4])+1
        assert len(roots)<17
        root_mass_hist[len(roots)]+=1;distinct_hist[len(set(roots))]+=1
    if block and block//4096%1024==0:
        print(json.dumps({'outer_processed':block,'survivors':len(seen),'elapsed':time.monotonic()-started}),flush=True)
assert seen==set(expected)
assert all_pass==52910 and nonzero_constant==6534
result={'status':'PASS','eliminated_parameter':'r4','outer_quadruples':p**4,'covered_polynomials':p**5,
        'all_histogram_survivors':all_pass,'nonzero_constant_survivors':nonzero_constant,
        'direction_extra_zero_points':[int(j)+2 for j in zeros],
        'mass_histogram':dict(root_mass_hist),'distinct_root_histogram':dict(distinct_hist),
        'complete_survivor_set_matches_D':True,'solutions':0,'elapsed_seconds':time.monotonic()-started}
(HERE/'independent_m4_audit_result.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2),flush=True)
