"""Independent arithmetic and small direct orbit checks for E8."""
import json
from itertools import combinations
from math import comb
from pathlib import Path

BITS=35
MASK=(1<<BITS)-1

def rotate_right(v,k):
    k%=BITS
    return ((v>>k)|(v<<(BITS-k)))&MASK if k else v

def canonical(v):
    return min(rotate_right(v,k) for k in range(BITS))

def burnside(o):
    terms={"identity":comb(35,o),"order5":0,"order7":0,"order35":0}
    if o%5==0:
        terms["order5"]=4*comb(7,o//5)
    if o%7==0:
        terms["order7"]=6*comb(5,o//7)
    if o%35==0:
        terms["order35"]=24*comb(1,o//35)
    numerator=sum(terms.values())
    assert numerator%35==0
    return numerator//35,terms

def brute_orbits(o):
    reps=set()
    all_min_have_bit0=True
    for positions in combinations(range(35),o):
        mask=sum(1<<i for i in positions)
        c=canonical(mask)
        all_min_have_bit0 &= bool(c&1)
        reps.add(c)
    return len(reps),all_min_have_bit0

def main():
    rows=[]
    for o in range(1,18,2):
        count,terms=burnside(o)
        item={"odd_size":o,"orbit_count":count,"fixed_subset_terms":terms}
        if o<=5:
            brute,all_bit0=brute_orbits(o)
            item["direct_brute_orbit_count"]=brute
            item["all_canonical_reps_have_bit0"]=all_bit0
            assert brute==count and all_bit0
        rows.append(item)
    result={"generator_square_root_12":pow(12,2,71),
            "generator_order":next(k for k in range(1,71) if pow(2,k,71)==1),
            "orbit_counts":rows,"sum_o_le9":sum(x["orbit_count"] for x in rows if x["odd_size"]<=9),
            "sum_all":sum(x["orbit_count"] for x in rows),
            "periodic_examples":{"o5":[0,7,14,21,28],"o7":[0,5,10,15,20,25,30]}}
    Path(__file__).with_name("check_E8_orbits_result.json").write_text(
        json.dumps(result,indent=2),encoding="utf-8")
    print(json.dumps({"sum_o_le9":result["sum_o_le9"],
                      "sum_all":result["sum_all"],
                      "small_direct_counts":[(x["odd_size"],x.get("direct_brute_orbit_count"))
                                             for x in rows if x["odd_size"]<=5]},indent=2))

if __name__=="__main__":
    main()
