"""Exact local odd-size symmetric PTE witness checker. Standard library only.

This is a witness checker plus small control tests, NOT an exhaustive search
for n=17,p=71. Repeated residues are deliberately supported.
Input JSON: {"p":71,"n":17,"A":[17 integer entries]}.
"""
from collections import Counter
from itertools import combinations_with_replacement
from math import isqrt, gcd, comb
from pathlib import Path
import argparse
import json

def prime(p):
    return p>=2 and all(p%d for d in range(2,isqrt(p)+1))

def root_polynomial(A,p):
    # Ascending coefficients; multiply linear factors exactly.
    f=[1]
    for a in A:
        out=[0]*(len(f)+1)
        for j,c in enumerate(f):
            out[j]=(out[j]-a*c)%p
            out[j+1]=(out[j+1]+c)%p
        f=out
    return f

def inspect(A,p,n):
    if type(p) is not int or type(n) is not int or not prime(p) or not 3<=n<p or n%2!=1:
        raise ValueError('Require an odd integer n>=3 and prime p>n.')
    if not isinstance(A,(list,tuple)) or any(type(a) is not int for a in A):
        raise ValueError('A must be a list of integers, retaining multiplicities.')
    if len(A)!=n:raise ValueError('Incorrect multiset cardinality.')
    A=[a%p for a in A];counts=Counter(A)
    moments=[sum(pow(a,j,p) for a in A)%p for j in range(1,n,2)]
    moments_ok=not any(moments)
    f=root_polynomial(A,p)
    shape_ok=all(f[j]==0 for j in range(2,n,2))
    # Second formulation: all nonconstant even coefficients vanish.
    assert moments_ok==shape_ok
    symmetric=counts==Counter((-a)%p for a in A)
    disjoint=all(counts.get((-a)%p,0)==0 for a in counts)
    if moments_ok:
        assert symmetric==(f[0]==0)
        assert disjoint==(f[0]!=0)
    valid=moments_ok and not symmetric
    result={'p':p,'n':n,'A':sorted(A),'odd_exponents':list(range(1,n,2)),
            'odd_moments':moments,'root_polynomial_ascending':f,
            'A_equals_negative_A':symmetric,'disjoint_from_negative_A':disjoint,
            'nontrivial_local_witness':valid,
            'scope':'Local finite-field witness check, not an integer PTE lifting proof.'}
    if valid:
        z=[counts.get(a,0)-counts.get(p-a,0) for a in range(1,(p+1)//2)]
        assert sum(abs(x) for x in z)==n
        assert all(sum(z[a-1]*pow(a,j,p) for a in range(1,(p+1)//2))%p==0
                   for j in range(1,n,2))
        assert sum(x!=0 for x in z)>=(n+1)//2
        result['signed_multiplicities']=z
        result['constant_difference']=2*f[0]%p
        if gcd(n,p-1)==1:
            u=pow(pow(f[0],-1,p),pow(n,-1,p-1),p)
            scaled=[u*a%p for a in A]
            g=root_polynomial(scaled,p)
            assert g[0]==1
            assert not any(sum(pow(a,j,p) for a in scaled)%p for j in range(1,n,2))
            result['constant_one_normalization']={'scale':u,'A':sorted(scaled)}
    return result

def rank_mod(matrix,p):
    B=[row[:] for row in matrix];r=0
    for j in range(len(B[0])):
        pivot=next((i for i in range(r,len(B)) if B[i][j]%p),None)
        if pivot is None:continue
        B[r],B[pivot]=B[pivot],B[r]
        t=pow(B[r][j],-1,p);B[r]=[t*v%p for v in B[r]]
        for i in range(len(B)):
            if i!=r:
                c=B[i][j];B[i]=[(a-c*b)%p for a,b in zip(B[i],B[r])]
        r+=1
        if r==len(B):break
    return r

def selftest():
    records=[]
    for p,n in [(5,3),(7,3),(11,5)]:
        total=witnesses=0
        for A in combinations_with_replacement(range(p),n):
            result=inspect(A,p,n);total+=1
            witnesses+=result['nontrivial_local_witness']
        assert total==comb(p+n-1,n)
        records.append({'p':p,'n':n,'all_multisets_checked':total,'nontrivial':witnesses})
    # A valid repeated-root control: a set-only checker would miss it.
    assert inspect([1,1,3],5,3)['nontrivial_local_witness']
    # Independently generated multiplicative subgroups, not target witnesses.
    controls=[]
    for p,n in [(23,11),(31,15)]:
        A=[a for a in range(1,p) if pow(a,n,p)==1]
        assert len(A)==n
        result=inspect(A,p,n)
        assert result['nontrivial_local_witness']
        controls.append({'p':p,'n':n,'A':A})
    trivial=[0]+[x for a in range(1,9) for x in (a,71-a)]
    assert not inspect(trivial,71,17)['nontrivial_local_witness']
    assert not inspect([1]*17,71,17)['nontrivial_local_witness']
    H=[[pow(a,j,71) for a in range(1,36)] for j in range(1,17,2)]
    assert rank_mod(H,71)==8
    # The complete constant-normalization scalar identity.
    assert 17*33%70==1
    for c in range(1,71):
        u=pow(pow(c,-1,71),33,71)
        assert pow(u,17,71)*c%71==1
    order=next(d for d in range(1,71) if pow(2,d,71)==1)
    assert order==35
    return {'status':'all exact assertions passed','small_exhaustive_controls':records,
            'positive_subgroup_controls':controls,'target_matrix_rank':8,
            'order_of_2_mod_71':order,'target_exhaustive_search_performed':False,
            'target_status':'not resolved by this program',
            'size_estimates':{'all_nonzero_multisets':comb(86,17),
              'raw_eight_half':comb(77,8),'raw_nine_half':comb(78,9),
              'constant_one_coefficient_space':71**8}}

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path)
    parser.add_argument('--out',type=Path)
    args=parser.parse_args()
    if args.input:
        source=json.loads(args.input.read_text(encoding='utf-8-sig'))
        result=inspect(source['A'],source.get('p',71),source.get('n',17))
    else:result=selftest()
    text=json.dumps(result,ensure_ascii=False,indent=2)
    if args.out:args.out.write_text(text,encoding='utf-8')
    print(text)
    return 0 if not args.input or result['nontrivial_local_witness'] else 1

if __name__=='__main__':
    if not __debug__:raise RuntimeError('Do not disable assertions with -O.')
    raise SystemExit(main())
