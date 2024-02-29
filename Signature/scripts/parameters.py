#!/usr/bin/env python3
from sage.all import *
proof.all(False)  # faster

import re
for l in open('sqisign_parameters.txt'):
    for k in ('lvl', 'p', 'B'):
        m = re.search(rf'^\s*{k}\s*=\s*([x0-9a-f]+)', l)
        if m:
            v = ZZ(m.groups()[0], 0)
            globals()[k] = v

L = {l for l,_ in (p**2 - 1).factor(limit=B+5) if l <= B}
assert 2 in L
L.remove(2)
f = (p+1).valuation(2)
if (p-1).valuation(2) > f:
    raise NotImplementedError('2-power torsion is on twist')
exp3 = (p+1).valuation(3)
if (p-1).valuation(3) > exp3:
    raise NotImplementedError('3-power torsion is on twist')
Lpls = {l for l in L if (p+1).valuation(l) >= (p-1).valuation(l)}
Lmin = L - Lpls
Lpls, Lmin = map(sorted, (Lpls, Lmin))
Epls = [(p+1).valuation(l) for l in Lpls]
Emin = [(p-1).valuation(l) for l in Lmin]
Tpls = prod(l**e for l,e in zip(Lpls,Epls))
Tmin = prod(l**e for l,e in zip(Lmin,Emin))

Dcom = (Tpls*Tmin).prime_to_m_part(2*3)
Dchall = prod(l**(p+1).valuation(l) for l in (2,3))

__all__ = ['lvl', 'p', 'B', 'f', 'exp3', 'Tpls', 'Tmin', 'Dcom', 'Dchall']

