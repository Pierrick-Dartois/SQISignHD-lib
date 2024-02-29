"""
This code is based on a copy of:
https://github.com/ThetaIsogenies/two-isogenies

MIT License

Copyright (c) 2023 Pierrick Dartois, Luciano Maino, Giacomo Pope and Damien Robert

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

from sage.all import *
from theta_structures.Tuple_point import TuplePoint
from theta_structures.Theta_dim2 import ThetaPointDim2
from isogenies_dim2.gluing_isogeny_dim2 import GluingThetaIsogenyDim2
from isogenies_dim2.isogeny_dim2 import ThetaIsogenyDim2
from utilities.strategy import optimised_strategy


class IsogenyChainDim2:
    r"""
    Given (P1, P2), (Q1, Q2) in (E1 x E2)[2^(n+2)] as the generators of a kernel
    of a (2^n, 2^n)-isogeny

    ker(Phi) = <(P1, P2), (Q1, Q2)>

    Input:

    - kernel = TuplePoint(P1, P2), TuplePoint(Q1, Q2):
      where points are on the elliptic curves E1, E2 of order 2^(n+2)
    - n: the length of the chain
    - strategy: the optimises strategy to compute a walk through the graph of
      images and doublings with a quasli-linear number of steps
    """

    def __init__(self, kernel, Theta12, M, n, strategy=None):
        self.n = n
        self.E1, self.E2 = kernel[0].parent_curves()
        assert kernel[1].parent_curves() == [self.E1, self.E2]

        self._domain = (self.E1, self.E2)

        if strategy is None:
            strategy = self.get_strategy()
        self.strategy = strategy

        self._phis = self.isogeny_chain(kernel, Theta12, M)

        self._codomain=self._phis[-1]._codomain

    def get_strategy(self):
        return optimised_strategy(self.n)

    def isogeny_chain(self, kernel, Theta12, M):
        """
        Compute the isogeny chain and store intermediate isogenies for evaluation
        """
        # Extract the CouplePoints from the Kernel
        Tp1, Tp2 = kernel

        # Store chain of (2,2)-isogenies
        isogeny_chain = []

        # Bookkeeping for optimal strategy
        strat_idx = 0
        level = [0]
        ker = (Tp1, Tp2)
        kernel_elements = [ker]

        for k in range(self.n):
            prev = sum(level)
            ker = kernel_elements[-1]

            while prev != (self.n - 1 - k):
                level.append(self.strategy[strat_idx])

                # Perform the doublings
                Tp1 = ker[0].double_iter(self.strategy[strat_idx])
                Tp2 = ker[1].double_iter(self.strategy[strat_idx])

                ker = (Tp1, Tp2)

                # Update kernel elements and bookkeeping variables
                kernel_elements.append(ker)
                prev += self.strategy[strat_idx]
                strat_idx += 1

            # Compute the codomain from the 8-torsion
            Tp1, Tp2 = ker
            if k == 0:
                phi = GluingThetaIsogenyDim2(Tp1, Tp2, Theta12, M)
            else:
                phi = ThetaIsogenyDim2(Th, Tp1, Tp2)

            # Update the chain of isogenies
            Th = phi._codomain
            isogeny_chain.append(phi)

            # Remove elements from list
            kernel_elements.pop()
            level.pop()

            # Push through points for the next step
            kernel_elements = [(phi(T1), phi(T2)) for T1, T2 in kernel_elements]

        return isogeny_chain

    def evaluate_isogeny(self, P):
        """
        Given a point P, of type TuplePoint on the domain E1 x E2, computes the
        ThetaPointDim2 on the codomain ThetaStructureDim2.
        """
        if not isinstance(P, TuplePoint):
            raise TypeError(
                "IsogenyChainDim2 isogeny expects as input a TuplePoint on the domain product E1 x E2"
            )
        n=len(self._phis)
        for i in range(n):
            P = self._phis[i](P)
        return P

    def __call__(self, P):
        """
        Evaluate a TuplePoint under the action of this isogeny.
        """
        return self.evaluate_isogeny(P)

    def dual(self):
        domain = self._codomain
        codomain = self._domain
        n=len(self._phis)
        isogenies=[]
        for i in range(n):
            isogenies.append(self._phis[n-1-i].dual())
        return DualIsogenyChainDim2(domain, codomain, isogenies)


class DualIsogenyChainDim2:
    def __init__(self, domain, codomain, isogenies):
        self._domain = domain
        self._codomain = codomain
        self._phis = isogenies

    def evaluate_isogeny(self, P):
        """
        Given a ThetaPointDim2 point P on the codomain ThetaStructureDim2, 
        computes the image TuplePoint on the codomain E1 x E2.
        """
        if not isinstance(P, ThetaPointDim2):
            raise TypeError(
                "DualIsogenyChainDim2 isogeny expects as input a ThetaPointDim2."
            )
        n=len(self._phis)
        for i in range(n):
            P = self._phis[i](P)
        return P

    def __call__(self, P):
        """
        Evaluate a ThetaPointDim2 under the action of this isogeny.
        """
        return self.evaluate_isogeny(P)
