"""Read-only independent exact validation of the saved enumeration constants."""
from fractions import Fraction as Q
from pathlib import Path
import json, math

data=json.loads(Path('lattice_certificate_input.json').read_text(encoding='utf8'))
B=data['B']; M=data['M']; D=data['D']; A=data['inverse_times_71']
N=35; S=1024; R=41; p=71
assert data['p']==p and data['S']==S and data['radius_squared']==R
assert len(B)==N and all(len(r)==N for r in B)
for b in B:
    assert all(sum(b[a-1]*pow(a,j,p) for a in range(1,N+1))%p==0 for j in range(1,16,2))
for i in range(N):
    for j in range(N):
        assert sum(B[i][k]*A[k][j] for k in range(N))==p*int(i==j)

G=[[Q(sum(x*y for x,y in zip(a,b))) for b in B] for a in B]
mu=[[Q(0) for j in range(N)] for i in range(N)]
d=[]
for i in range(N):
    for j in range(i):
        mu[i][j]=(G[i][j]-sum(mu[i][k]*mu[j][k]*d[k] for k in range(j)))/d[j]
        assert M[i][j]==(S*mu[i][j]).__floor__()
    d.append(G[i][i]-sum(mu[i][k]**2*d[k] for k in range(i)))
    assert d[i]>0 and D[i]==(S*d[i]).__floor__()
assert math.prod(d)==p**16
for j in range(N):
    bound=math.isqrt(R*sum(A[k][j]**2 for k in range(N)))//p
    assert bound==data['coordinate_bounds'][j]

expected=[N,S,R]+[x for row in B for x in row]+[x for row in M for x in row]+D+data['coordinate_bounds']
actual=list(map(int,Path('lattice_certificate_input.txt').read_text(encoding='ascii').split()))
assert actual==expected
print('PASS: exact moment kernel, exact inverse, exact determinant, rational Gram-Schmidt floors, coefficient bounds, text/JSON equality.')
print('The 8-by-8 minor at columns 1..8 is nonzero: nonzero column factors times a Vandermonde in distinct squares modulo 71.')
