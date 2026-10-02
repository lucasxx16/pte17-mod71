"""Independent exact verification of D's saved lattice certificate."""
from fractions import Fraction
from math import isqrt, prod
from pathlib import Path
import json
import random

HERE=Path(__file__).parent
data=json.loads((HERE/'lattice_certificate_input.json').read_text(encoding='utf-8'))
P,N,S,R=71,35,1024,41
B=data['B'];A=data['inverse_times_71'];M=data['M'];D=data['D'];bounds=data['coordinate_bounds']

def determinant_mod(matrix,p):
    a=[[x%p for x in row] for row in matrix]
    r=0;det=1
    for j in range(len(a[0])):
        pivot=next((i for i in range(r,len(a)) if a[i][j]),None)
        if pivot is None:
            continue
        if pivot!=r:
            a[r],a[pivot]=a[pivot],a[r]
            det=-det%p
        x=a[r][j]
        det=det*x%p
        xi=pow(x,-1,p)
        for i in range(r+1,len(a)):
            factor=a[i][j]*xi%p
            for k in range(j,len(a[0])):
                a[i][k]=(a[i][k]-factor*a[r][k])%p
        r+=1
        if r==len(a):
            break
    return r,det

def main():
    assert (data['p'],data['S'],data['radius_squared'])==(P,S,R)
    assert len(B)==len(A)==len(M)==len(D)==len(bounds)==N
    H=[[pow(a,2*j+1,P) for a in range(1,N+1)] for j in range(8)]
    assert all(sum(x*y for x,y in zip(row,h))%P==0 for row in B for h in H)
    rank,_=determinant_mod([h[:8] for h in H],P)
    assert rank==8
    assert all(sum(B[i][k]*A[k][j] for k in range(N))==P*(i==j)
               for i in range(N) for j in range(N))
    norms=[Fraction(*x) for x in data['gso_norms_exact']]
    mus=[[Fraction(*x) for x in row] for row in data['mu_exact']]
    assert all(x>0 for x in norms)
    assert prod(norms)==P**16
    assert all((S*norms[i]).numerator//(S*norms[i]).denominator==D[i]
               for i in range(N))
    assert all((S*mus[i][j]).numerator//(S*mus[i][j]).denominator==M[i][j]
               for i in range(N) for j in range(i))
    # Independent Gram = L diag(norms) L^T identity, L_ii=1.
    def l(i,k):
        return Fraction(1) if i==k else mus[i][k]
    assert all(sum(B[i][t]*B[j][t] for t in range(N))==
               sum(l(i,k)*l(j,k)*norms[k] for k in range(min(i,j)+1))
               for i in range(N) for j in range(N))
    assert all(bounds[j]==isqrt(R*sum(A[k][j]**2 for k in range(N)))//P
               for j in range(N))
    expected=[N,S,R]+[x for row in B for x in row]+[x for row in M for x in row]+D+bounds
    actual=list(map(int,(HERE/'lattice_certificate_input.txt').read_text(encoding='ascii').split()))
    assert expected==actual
    rng=random.Random(20261002)
    lower_bound_checks=0
    for _ in range(1000):
        z=[rng.randint(-bounds[i],bounds[i]) for i in range(N)]
        v=[sum(z[i]*B[i][j] for i in range(N)) for j in range(N)]
        exact=sum(x*x for x in v)
        lb_num=0
        for i in range(N):
            center=S*z[i]+sum(M[j][i]*z[j] for j in range(i+1,N))
            err=sum(abs(z[j]) for j in range(i+1,N))
            a=max(0,abs(center)-err)
            lb_num+=D[i]*a*a
        assert lb_num<=S**3*exact
        lower_bound_checks+=1
    result={"status":"PASS","n":N,"rank_H":rank,"det_abs":P**8,
            "GS_product":str(prod(norms)),"min_D":min(D),"max_D":max(D),
            "input_text_matches_json":True,"independent_Gram_identity":True,
            "random_signed_lower_bound_checks":lower_bound_checks,
            "coordinate_bounds":bounds}
    (HERE/'check_certificate_independent_result.json').write_text(
        json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    main()
