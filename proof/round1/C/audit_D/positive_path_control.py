"""Trace a known short lattice vector through every integer-pruning interval."""
import json
from math import isqrt
from pathlib import Path

HERE=Path(__file__).parent
P,N,S,R=71,35,1024,41
cert=json.loads((HERE/'lattice_certificate_input.json').read_text(encoding='utf-8'))
record=json.loads((HERE.parent.parent/'D'/'lattice_records.json').read_text(encoding='utf-8'))[0]
B,A,M,D,bounds=(cert[k] for k in
                ('B','inverse_times_71','M','D','coordinate_bounds'))
RECORD_V=record['vector']

def ceildiv(a,b):
    return -((-a)//b)

def main():
    v=RECORD_V[:]
    assert len(v)==N
    assert sum(map(abs,v))==18 and sum(x*x for x in v)==22
    assert max(map(abs,v))==2
    assert all(sum(v[a-1]*pow(a,j,P) for a in range(1,N+1))%P==0
               for j in range(1,16,2))
    numer=[sum(v[j]*A[j][i] for j in range(N)) for i in range(N)]
    assert all(x%P==0 for x in numer)
    z=[x//P for x in numer]
    top=next(i for i in reversed(range(N)) if z[i])
    if z[top]<0:
        z=[-x for x in z]
        v=[-x for x in v]
    assert all(sum(z[i]*B[i][j] for i in range(N))==v[j] for j in range(N))
    assert all(abs(z[i])<=bounds[i] for i in range(N))
    remaining=R*S**3
    tail_abs=0
    tail_nonzero=False
    trace=[]
    for k in reversed(range(N)):
        center=sum(M[j][k]*z[j] for j in range(k+1,N))
        reach=isqrt(remaining//D[k])+tail_abs
        lo=max(-bounds[k],ceildiv(-reach-center,S))
        hi=min(bounds[k],(reach-center)//S)
        if not tail_nonzero:
            lo=max(0,lo)
        assert lo<=z[k]<=hi,(k,z[k],lo,hi)
        if k>0:
            a=max(0,abs(S*z[k]+center)-tail_abs)
            remaining-=D[k]*a*a
            assert remaining>=0
            tail_abs+=abs(z[k])
            tail_nonzero|=z[k]!=0
        trace.append({"k":k,"chosen_z":z[k],"lo":lo,"hi":hi,
                      "remaining_before":remaining if k==0 else None})
    result={"status":"PASS","l1":sum(map(abs,v)),"l2_squared":sum(x*x for x in v),
            "H_zero":True,"z_top_index":top,"z_top_positive":z[top]>0,
            "all_35_intervals_contain_known_vector":True,
            "original_parity_filter_would_reject":sum(v)%2==0,
            "leaf_mass_18_positive_control":True,
            "min_margin_to_bound":min(min(x['chosen_z']-x['lo'],x['hi']-x['chosen_z']) for x in trace),
            "trace":trace}
    (HERE/'positive_path_control_result.json').write_text(
        json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='trace'},indent=2))

if __name__=='__main__':
    main()
