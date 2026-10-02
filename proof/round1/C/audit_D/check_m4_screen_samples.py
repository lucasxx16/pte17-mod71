"""Brute all 71 r0 values and all 69 other field points on sampled m4 outer tuples."""
import json
import random
from pathlib import Path
import numpy as np

HERE=Path(__file__).parent
P=71

def multiply(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            c[i+j]=(c[i+j]+x*y)%P
    return c

def evalpoly(a,x):
    t=0
    for c in reversed(a):
        t=(t*x+c)%P
    return t

def basis():
    factor=[1]
    for _ in range(3):
        factor=multiply(factor,[P-1,0,1])
    out=[]
    for k in range(6):
        g=multiply(factor,[0]*(2*k)+[1])
        out.append([0]+[17*c*pow(i+1,-1,P)%P for i,c in enumerate(g)])
    return out

def main():
    data=json.loads((HERE.parent.parent/'D'/'multiplicity4_result.json').read_text(encoding='utf-8'))
    recorded={tuple(rec['R_ascending'][:5]) for rec in data['survivors']}
    survivor_outer={tuple(rec['R_ascending'][1:5]) for rec in data['survivors']}
    rng=random.Random(20261002)
    outers={(0,0,0,0),(70,70,70,70)}
    outers.update(tuple(rng.randrange(P) for _ in range(4)) for _ in range(500))
    outers.update(sorted(survivor_outer)[:100])
    Js=basis()
    xs=list(range(2,P))
    E=np.array([[(evalpoly(J,x)-evalpoly(J,1))%P for x in xs] for J in Js],dtype=np.int64)
    assert [xs[i] for i in np.flatnonzero(E[0]==0)]==[39]
    other=np.arange(P,dtype=np.int64)
    mismatches=[]
    total_direct_candidates=0
    for outer in sorted(outers):
        base=(E[5]+sum(outer[j-1]*E[j] for j in range(1,5)))%P
        vals=(base[None,:]+other[:,None]*E[0][None,:])%P
        counts=(vals==0).sum(axis=1)
        assert np.all(counts<=69)
        j1=(evalpoly(Js[5],1)+sum(outer[j-1]*evalpoly(Js[j],1) for j in range(1,5))+
            other*evalpoly(Js[0],1))%P
        direct={tuple([int(r0)]+list(outer)) for r0 in range(P)
                if counts[r0]>=8 and j1[r0]!=0}
        total_direct_candidates+=len(direct)
        listed={x for x in recorded if x[1:]==outer}
        if direct!=listed:
            mismatches.append({"outer":outer,"missing":list(direct-listed),
                               "extra":list(listed-direct)})
    result={"outer_tuples_checked":len(outers),"r0_points_each":71,
            "field_points_each":69,"exceptional_E0_zero_points":[39],
            "direct_candidates_in_sample":total_direct_candidates,
            "record_mismatches":len(mismatches),"first_mismatches":mismatches[:5]}
    (HERE/'check_m4_screen_samples_result.json').write_text(
        json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    main()
