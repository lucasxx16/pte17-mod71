"""E audit: independent direct GS and Bareiss determinant; exact safety bounds."""
from pathlib import Path
from fractions import Fraction
import json, math, hashlib

HERE=Path(__file__).resolve().parent
data=json.loads((HERE/'lattice_certificate_input.json').read_text())
B=data['B']; A=data['inverse_times_71']; M=data['M']; D=data['D']; bounds=data['coordinate_bounds']
N=35; scale=1024; radius=41; budget=radius*scale**3
assert len(B)==N and all(len(row)==N for row in B)
assert all(sum(b[a-1]*pow(a,j,71) for a in range(1,36))%71==0 for b in B for j in range(1,16,2))
assert all(sum(B[i][k]*A[k][j] for k in range(N))==71*(i==j) for i in range(N) for j in range(N))

# Bareiss elimination directly on B, independent of GS determinant product.
T=[row[:] for row in B]; previous=1; sign=1
for k in range(N-1):
    pivot=next(i for i in range(k,N) if T[i][k])
    if pivot!=k: T[k],T[pivot]=T[pivot],T[k]; sign=-sign
    p=T[k][k]
    for i in range(k+1,N):
        for j in range(k+1,N):
            numerator=p*T[i][j]-T[i][k]*T[k][j]
            assert numerator%previous==0
            T[i][j]=numerator//previous
        T[i][k]=0
    previous=p
determinant=sign*T[-1][-1]
assert abs(determinant)==71**8

# Explicit rational orthogonal vectors, rather than D's LDL implementation.
orthogonal=[]; norms=[]; mu=[]
for i,b in enumerate(B):
    w=list(map(Fraction,b)); row=[Fraction(0)]*N
    for j in range(i):
        row[j]=sum(Fraction(b[k])*orthogonal[j][k] for k in range(N))/norms[j]
        w=[x-row[j]*y for x,y in zip(w,orthogonal[j])]
        assert M[i][j]==(scale*row[j]).numerator//(scale*row[j]).denominator
    norm=sum(x*x for x in w)
    assert norm>0 and D[i]==(scale*norm).numerator//(scale*norm).denominator
    assert data['gso_norms_exact'][i]==[norm.numerator,norm.denominator]
    assert data['mu_exact'][i]==[[x.numerator,x.denominator] for x in row]
    orthogonal.append(w); norms.append(norm); mu.append(row)
assert math.prod(norms)==71**16
assert bounds==[math.isqrt(radius*sum(A[j][i]**2 for j in range(N)))//71 for i in range(N)]
expected=[N,scale,radius]+sum(B,[])+sum(M,[])+D+bounds
assert list(map(int,(HERE/'lattice_certificate_input.txt').read_text().split()))==expected

center=[sum(abs(M[j][i])*bounds[j] for j in range(i+1,N)) for i in range(N)]
tail_abs=[sum(bounds[i+1:]) for i in range(N)]
sqrt_argument=[budget//x for x in D]
root_bound=[math.isqrt(x) for x in sqrt_argument]
reach=[a+b for a,b in zip(root_bound,tail_abs)]
linear_form=[scale*b+c for b,c in zip(bounds,center)]
coordinate_abs=[sum(bounds[j]*abs(B[j][i]) for j in range(N)) for i in range(N)]
rowsums=list(map(sum,B))
safety={
    'budget':budget,'min_D':min(D),'max_D':max(D),
    'max_basis_coordinate_bound':max(bounds),'sum_basis_coordinate_bounds':sum(bounds),
    'max_abs_center':max(center),'max_tail_abs':max(tail_abs),
    'max_isqrt_argument':max(sqrt_argument),'max_exact_isqrt':max(root_bound),
    'max_isqrt_correction_square_if_initial_off_by_one':(max(root_bound)+2)**2,
    'max_reach':max(reach),'max_abs_floor_div_numerator':max(a+b for a,b in zip(reach,center)),
    'max_abs_scaled_center':max(linear_form),
    'max_recursive_D_times_a_squared':budget,
    'max_abs_leaf_coordinate':max(coordinate_abs),
    'max_leaf_coordinate_square':max(coordinate_abs)**2,
    'max_leaf_total_mass':sum(coordinate_abs),
    'max_leaf_total_squared_norm':sum(x*x for x in coordinate_abs),
    'max_abs_parity_update':max(1+b*abs(t) for b,t in zip(bounds,rowsums)),
    'center_by_level':center,'coordinate_abs_bounds':coordinate_abs,
}
assert all(v<2**63 for v in safety.values() if isinstance(v,int))
result={'status':'PASS','basis_determinant':determinant,'safety_bounds':safety,
        'sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in HERE.iterdir() if p.is_file() and p.name in ['lattice_certificate_input.json','lattice_certificate_input.txt','IntegerLatticeEnumerator.cs','check_multiplicity4.py']}}
(HERE/'independent_certificate_audit_result.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
