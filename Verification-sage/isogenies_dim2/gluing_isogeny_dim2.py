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
from theta_structures.Theta_dim2 import ThetaStructureDim2, ThetaPointDim2
from theta_structures.theta_helpers_dim2 import batch_inversion
from basis_change.base_change_dim2 import montgomery_to_theta_matrix_dim2, apply_base_change_theta_dim2
from theta_structures.montgomery_theta import lift_kummer_montgomery_point

class GluingThetaIsogenyDim2:
    """
    Compute the gluing isogeny from E1 x E2 (Elliptic Product) -> A (Theta Model)

    Expected input:

    - (K1_8, K2_8) The 8-torsion above the kernel generating the isogeny
    - M (Optional) a base change matrix, if this is not including, it can
      be derived from [2](K1_8, K2_8)
    """

    def __init__(self, K1_8, K2_8, Theta12, N):
        # Double points to get four-torsion, we always need one of these, used
        # for the image computations but we'll need both if we wish to derived
        # the base change matrix as well
        K1_4 = 2*K1_8

        # Initalise self
        # This is the base change matrix for product Theta coordinates (not used, except in the dual)
        self._base_change_matrix_theta = N
        # Here, base change matrix directly applied to the Montgomery coordinates. null_point_bc is the 
        # theta null point obtained after applying the base change to the product Theta-structure.
        self._base_change_matrix, null_point_bc = montgomery_to_theta_matrix_dim2(Theta12.zero().coords(),N, return_null_point = True)
        self._domain_bc = ThetaStructureDim2(null_point_bc)
        self.T_shift = K1_4
        self._precomputation = None
        self._zero_idx = 0
        self._domain_product = Theta12
        self._domain=(K1_8[0].curve(), K1_8[1].curve())

        # Map points from elliptic product onto the product theta structure
        # using the base change matrix
        T1_8 = self.base_change(K1_8)
        T2_8 = self.base_change(K2_8)

        # Compute the codomain of the gluing isogeny
        self._codomain = self._special_compute_codomain(T1_8, T2_8)

    def apply_base_change(self, coords):
        """
        Apply the basis change by acting with matrix multiplication, treating
        the coordinates as a vector
        """
        N = self._base_change_matrix
        x, y, z, t = coords
        X = N[0, 0] * x + N[0, 1] * y + N[0, 2] * z + N[0, 3] * t
        Y = N[1, 0] * x + N[1, 1] * y + N[1, 2] * z + N[1, 3] * t
        Z = N[2, 0] * x + N[2, 1] * y + N[2, 2] * z + N[2, 3] * t
        T = N[3, 0] * x + N[3, 1] * y + N[3, 2] * z + N[3, 3] * t

        return (X, Y, Z, T)

    def base_change(self, P):
        """
        Compute the basis change on a TuplePoint to recover a ThetaPointDim2 of
        compatible form
        """
        if not isinstance(P, TuplePoint):
            raise TypeError("Function assumes that the input is of type `TuplePoint`")

        # extract X,Z coordinates on pairs of points
        P1, P2 = P.points()
        X1, Z1 = P1[0], P1[2]
        X2, Z2 = P2[0], P2[2]

        # Correct in the case of (0 : 0)
        if X1 == 0 and Z1 == 0:
            X1 = 1
            Z1 = 0
        if X2 == 0 and Z2 == 0:
            X2 = 1
            Z2 = 0

        # Apply the basis transformation on the product
        coords = self.apply_base_change([X1 * X2, X1 * Z2, Z1 * X2, Z1 * Z2])
        return coords

    def _special_compute_codomain(self, T1, T2):
        """
        Given twzero_matro isotropic points of 8-torsion T1 and T2, compatible with
        the theta null point, compute the level two theta null point A/K_2
        """
        xAxByCyD = ThetaPointDim2.to_squared_theta(*T1)
        zAtBzYtD = ThetaPointDim2.to_squared_theta(*T2)

        # Find the value of the non-zero index
        zero_idx = next((i for i, x in enumerate(xAxByCyD) if x == 0), None)
        self._zero_idx = zero_idx

        # Dumb check to make sure everything is OK
        assert xAxByCyD[self._zero_idx] == zAtBzYtD[self._zero_idx] == 0

        # Initialize lists
        # The zero index described the permutation
        ABCD = [0 for _ in range(4)]
        precomp = [0 for _ in range(4)]

        # Compute non-trivial numerators (Others are either 1 or 0)
        num_1 = zAtBzYtD[1 ^ self._zero_idx]
        num_2 = xAxByCyD[2 ^ self._zero_idx]
        num_3 = zAtBzYtD[3 ^ self._zero_idx]
        num_4 = xAxByCyD[3 ^ self._zero_idx]

        # Compute and invert non-trivial denominators
        den_1, den_2, den_3, den_4 = batch_inversion([num_1, num_2, num_3, num_4])

        # Compute A, B, C, D
        ABCD[0 ^ self._zero_idx] = 0
        ABCD[1 ^ self._zero_idx] = num_1 * den_3
        ABCD[2 ^ self._zero_idx] = num_2 * den_4
        ABCD[3 ^ self._zero_idx] = 1

        # Compute precomputation for isogeny images
        precomp[0 ^ self._zero_idx] = 0
        precomp[1 ^ self._zero_idx] = den_1 * num_3
        precomp[2 ^ self._zero_idx] = den_2 * num_4
        precomp[3 ^ self._zero_idx] = 1
        self._precomputation = precomp

        # Final Hadamard of the above coordinates
        a, b, c, d = ThetaPointDim2.to_hadamard(*ABCD)

        return ThetaStructureDim2([a, b, c, d])

    def special_image(self, P, translate):
        """
        When the domain is a non product theta structure on a product of
        elliptic curves, we will have one of A,B,C,D=0, so the image is more
        difficult. We need to give the coordinates of P but also of
        P+Ti, Ti one of the point of 4-torsion used in the isogeny
        normalisation
        """
        AxByCzDt = ThetaPointDim2.to_squared_theta(*P)

        # We are in the case where at most one of A, B, C, D is
        # zero, so we need to account for this
        #
        # To recover values, we use the translated point to get
        AyBxCtDz = ThetaPointDim2.to_squared_theta(*translate)

        # Directly compute y,z,t
        y = AxByCzDt[1 ^ self._zero_idx] * self._precomputation[1 ^ self._zero_idx]
        z = AxByCzDt[2 ^ self._zero_idx] * self._precomputation[2 ^ self._zero_idx]
        t = AxByCzDt[3 ^ self._zero_idx]

        # We can compute x from the translation
        # First we need a normalisation
        if z != 0:
            zb = AyBxCtDz[3 ^ self._zero_idx]
            lam = z / zb
        else:
            tb = AyBxCtDz[2 ^ self._zero_idx] * self._precomputation[2 ^ self._zero_idx]
            lam = t / tb

        # Finally we recover x
        xb = AyBxCtDz[1 ^ self._zero_idx] * self._precomputation[1 ^ self._zero_idx]
        x = xb * lam

        xyzt = [0 for _ in range(4)]
        xyzt[0 ^ self._zero_idx] = x
        xyzt[1 ^ self._zero_idx] = y
        xyzt[2 ^ self._zero_idx] = z
        xyzt[3 ^ self._zero_idx] = t

        image = ThetaPointDim2.to_hadamard(*xyzt)
        return self._codomain(image)

    def __call__(self, P):
        """
        Take into input the theta null point of A/K_2, and return the image
        of the point by the isogeny
        """
        if not isinstance(P, TuplePoint):
            raise TypeError(
                "Isogeny image for the gluing isogeny is defined to act on TuplePoints"
            )

        # Compute sum of points on elliptic curve
        P_sum_T = P + self.T_shift

        # Push both the point and the translation through the
        # completion
        iso_P = self.base_change(P)
        iso_P_sum_T = self.base_change(P_sum_T)

        return self.special_image(iso_P, iso_P_sum_T)

    def dual(self):
        domain = self._codomain.hadamard()
        codomain_bc = self._domain_bc.hadamard()
        codomain = self._domain

        precomputation = batch_inversion(codomain_bc.null_point_dual())

        N_split = self._base_change_matrix.inverse()

        return DualGluingThetaIsogenyDim2(domain, codomain_bc, codomain, N_split, precomputation)


class DualGluingThetaIsogenyDim2:
    def __init__(self, domain, codomain_bc, codomain, N_split, precomputation):
        self._domain = domain
        self._codomain_bc = codomain_bc # Theta structure
        self._codomain = codomain # Elliptic curves E1 and E2
        self._precomputation = precomputation
        self._splitting_matrix = N_split

    def __call__(self,P):
        # Returns a TuplePoint.
        if not isinstance(P, ThetaPointDim2):
            raise TypeError("Isogeny evaluation expects a ThetaPointDim2 as input")

        xx, yy, zz, tt = P.squared_theta()

        Ai, Bi, Ci, Di = self._precomputation

        xx = xx * Ai
        yy = yy * Bi
        zz = zz * Ci
        tt = tt * Di

        image_coords = (xx, yy, zz, tt)

        X1X2, X1Z2, Z1X2, Z1Z2 = apply_base_change_theta_dim2(self._splitting_matrix, image_coords)

        E1, E2 = self._codomain

        if Z1Z2!=0:
            #Z1=1, Z2=Z1Z2

            Z2_inv=1/Z1Z2
            X2=Z1X2*Z2_inv# Normalize (X2:Z2)=(X2/Z2:1)

            X1=X1Z2*Z2_inv

            assert X1*Z1X2==X1X2
            P1 = lift_kummer_montgomery_point(E1, X1)
            P2 = lift_kummer_montgomery_point(E2, X2)
            return TuplePoint(P1,P2)
        elif Z1X2==0 and X1Z2!=0:
            # Case (X1:Z1)=0, X1!=0 and (X2:Z2)!=0

            X2=X1X2/X1Z2
            P2 = lift_kummer_montgomery_point(E2, X2)
            return TuplePoint(E1(0),P2)
        elif Z1X2!=0 and X1Z2==0:
            # Case (X1:Z1)!=0 and (X2:Z2)=0, X2!=0

            X1=X1X2/Z1X2
            P1 = lift_kummer_montgomery_point(E1, X1)
            return TuplePoint(P1,E2(0))
        else:
            return TuplePoint(E1(0),E2(0))





