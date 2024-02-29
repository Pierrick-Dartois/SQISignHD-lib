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

from sage.all import ZZ
from theta_structures.Theta_dim2 import ThetaStructureDim2, ThetaPointDim2
from theta_structures.Tuple_point import TuplePoint
from utilities.theta_helpers import batch_inversion


class ThetaIsogeny:
    def __init__(self, domain, T1_8, T2_8, hadamard=(False, True)):
        """
        Compute a (2,2)-isogeny in the theta model. Expects as input:

        - domain: the ThetaStructureDim2 from which we compute the isogeny
        - (T1_8, T2_8): points of 8-torsion above the kernel generating the isogeny

        When the 8-torsion is not available (for example at the end of a long
        (2,2)-isogeny chain), the the helper functions in isogeny_sqrt.py
        must be used.

        NOTE: on the hadamard bools:

        The optional parameter 'hadamard' controls if we are in standard or dual
        coordinates, and if the codomain is in standard or dual coordinates. By
        default this is (False, True), meaning we use standard coordinates on
        the domain A and the codomain B.

        The kernel is then the kernel K_2 where the action is by sign. Other
        possibilities: - (False, False): standard coordinates on A, dual
        coordinates on B - (True, True): start in dual coordinates on A
        (alternatively: standard coordinates on A but quotient by K_1 whose
        action is by permutation), and standard coordinates on B. - (True,
        False): dual coordinates on A and B

        These can be composed as follows for A -> B -> C:

        - (False, True) -> (False, True) (False, False) -> (True, True):
          - standard coordinates on A and C,
          - standard/resp dual coordinates on B
        - (False, True) -> (False, False) (False, False) -> (True, False):
          - standard coordinates on A,
          - dual coordinates on C,
          - standard/resp dual coordinates on B
        - (True, True) -> (False, True) (True, False) -> (True, True):
          - dual coordinates on A,
          - standard coordinates on C,
          - standard/resp dual coordiantes on B
        - (True, True) -> (False, False) (True, False) -> (True, False):
          - dual coordinates on A and C
          - standard/resp dual coordinates on B

        On the other hand, these gives the multiplication by [2] on A:

        - (False, False) -> (False, True) (False, True) -> (True, True):
          - doubling in standard coordinates on A
          - going through dual/standard coordinates on B=A/K_2
        - (True, False) -> (False, False) (True, True) -> (True, False):
          - doubling in dual coordinates on A
          - going through dual/standard coordinates on B=A/K_2
            (alternatively: doubling in standard coordinates on A going
            through B'=A/K_1)
        - (False, False) -> (False, False) (False, True) -> (True, False):
          - doubling from standard to dual coordinates on A
        - (True, False) -> (False, True) (True, True) -> (True, True):
          - doubling from dual to standard coordinates on A
        """
        if not isinstance(domain, ThetaStructureDim2):
            raise ValueError
        self._domain = domain

        self._hadamard = hadamard
        self._precomputation = None
        self._codomain = self._compute_codomain(T1_8, T2_8)

    def _compute_codomain(self, T1, T2):
        """
        Given two isotropic points of 8-torsion T1 and T2, compatible with
        the theta null point, compute the level two theta null point A/K_2
        """
        if self._hadamard[0]:
            xA, xB, _, _ = ThetaPointDim2.to_squared_theta(
                *ThetaPointDim2.to_hadamard(*T1.coords())
            )
            zA, tB, zC, tD = ThetaPointDim2.to_squared_theta(
                *ThetaPointDim2.to_hadamard(*T2.coords())
            )
        else:
            xA, xB, _, _ = T1.squared_theta()
            zA, tB, zC, tD = T2.squared_theta()

        if not self._hadamard[0] and self._domain._precomputation:
            # Batch invert denominators
            xA_inv, zA_inv, tB_inv = batch_inversion(xA, zA, tB)

            # Compute A, B, C, D
            A = ZZ(1)
            B = xB * xA_inv
            C = zC * zA_inv
            D = tD * tB_inv * B

            _, _, _, BBinv, CCinv, DDinv = self._domain._arithmetic_precomputation()
            B_inv = BBinv * B
            C_inv = CCinv * C
            D_inv = DDinv * D
        else:
            # Batch invert denominators
            xA_inv, zA_inv, tB_inv, xB_inv, zC_inv, tD_inv = batch_inversion(
                xA, zA, tB, xB, zC, tD
            )

            # Compute A, B, C, D
            A = ZZ(1)
            B = xB * xA_inv
            C = zC * zA_inv
            D = tD * tB_inv * B
            B_inv = xB_inv * xA
            C_inv = zC_inv * zA
            D_inv = tD_inv * tB * B_inv

            # NOTE: some of the computations we did here could be reused for the
            # arithmetic precomputations of the codomain However, we are always
            # in the mode (False, True) except the very last isogeny, so we do
            # not lose much by not doing this optimisation Just in case we need
            # it later:
            # - for hadamard=(False, True): we can reuse the arithmetic
            #   precomputation; we do this already above
            # - for hadamard=(False, False): we can reuse the arithmetic
            #   precomputation as above, and furthermore we could reuse B_inv,
            #   C_inv, D_inv for the precomputation of the codomain
            # - for hadamard=(True, False): we could reuse B_inv, C_inv, D_inv
            #   for the precomputation of the codomain
            # - for hadamard=(True, True): nothing to reuse!

        self._precomputation = (B_inv, C_inv, D_inv)
        if self._hadamard[1]:
            a, b, c, d = ThetaPointDim2.to_hadamard(A, B, C, D)
            return ThetaStructureDim2([a, b, c, d])
        else:
            return ThetaStructureDim2([A, B, C, D])

    def __call__(self, P):
        """
        Take into input the theta null point of A/K_2, and return the image
        of the point by the isogeny
        """
        if not isinstance(P, ThetaPointDim2):
            raise TypeError("Isogeny evaluation expects a TuplePoint as input")

        if self._hadamard[0]:
            xx, yy, zz, tt = ThetaPointDim2.to_squared_theta(
                *ThetaPointDim2.to_hadamard(*P.coords())
            )
        else:
            xx, yy, zz, tt = P.squared_theta()

        Bi, Ci, Di = self._precomputation

        yy = yy * Bi
        zz = zz * Ci
        tt = tt * Di

        image_coords = (xx, yy, zz, tt)
        if self._hadamard[1]:
            image_coords = ThetaPointDim2.to_hadamard(*image_coords)
        return self._codomain(image_coords)

    def domain(self):
        """Return the domain of the morphism"""
        return self._domain

    def codomain(self):
        """Return the codomain of the morphism"""
        return self._codomain

class GluingThetaIsogeny(ThetaIsogeny):
    """
    Compute the gluing isogeny from E1 x E2 (Elliptic Product) -> A (Theta Model)

    Expected input:

    - (K1_8, K2_8) The 8-torsion above the kernel generating the isogeny
    - M (Optional) a base change matrix, if this is not including, it can
      be derived from [2](K1_8, K2_8)
    """

    def __init__(self, K1_8, K2_8, M=None):
        # Double points to get four-torsion, we always need one of these, used
        # for the image computations but we'll need both if we wish to derived
        # the base change matrix as well
        K1_4 = K1_8.double()

        # If M is not included, compute the matrix on the fly from the four
        # torsion.
        if M is None:
            K2_4 = K2_8.double()
            M = self.get_base_change_matrix(K1_4, K2_4)

        # Initalise self
        self._base_change_matrix = M
        self.T_shift = K1_4
        self._precomputation = None
        self._zero_idx = 0

        # Map points from elliptic product onto the product theta structure
        # using the base change matrix
        T1_8 = self.base_change(K1_8)
        T2_8 = self.base_change(K2_8)

        # Compute the codomain of the gluing isogeny
        self._codomain = self._special_compute_codomain(T1_8, T2_8)

    @staticmethod
    def get_base_change_matrix(T1, T2):
        """
        Given the four torsion above the kernel generating the gluing isogeny,
        compute the matrix M which allows us to map points on an elliptic
        product to the compatible theta structure.
        """

        def get_matrix(T):
            """
            Compute the matrix [a, b]
                               [c, d]
            From the (X : Z) coordinates of a point
            T and its double TT = (T + T)

            NOTE: Assumes that Z = 1 or 0
            """
            TT = T + T
            x = T[0]
            u = TT[0]

            # Our Z-coordinates are always 1 or 0
            # as Sage uses affine addition
            assert T[2] == 1
            assert TT[2] == 1

            # Compute the matrix coefficents
            det = x - u
            inv_det = 1 / det
            a = -u * inv_det
            b = -inv_det
            c = x * (u - det) * inv_det
            d = -a
            return Matrix(2, 2, [a, b, c, d])

        # Extract elliptic curve points from TuplePoints
        P1, P2 = T1.points()
        Q1, Q2 = T2.points()

        # Compute matrices from points
        # TODO: totally overkill, but if these were all computed together
        #       you could batch all 4 inversions into one.
        g1 = get_matrix(P1)
        g2 = get_matrix(P2)
        h1 = get_matrix(Q1)
        h2 = get_matrix(Q2)

        # the matrices gi, hi does not commute, but g1 \tens g2 should commute with h1 \tens h2
        gh1 = g1 * h1
        gh2 = g2 * h2

        # Access Coefficients once
        # Notice some coeffs are never used
        g00_1, g01_1, g10_1, g11_1 = g1.list()
        g00_2, _, g10_2, _ = g2.list()
        h00_1, _, h10_1, _ = h1.list()
        h00_2, h01_2, h10_2, h11_2 = h2.list()

        # First row from the product of gi * hi
        gh00_1, _, gh10_1, _ = gh1.list()
        gh00_2, _, gh10_2, _ = gh2.list()

        # start the trace with id
        a = 1
        b = 0
        c = 0
        d = 0

        # T1
        a += g00_1 * g00_2
        b += g00_1 * g10_2
        c += g10_1 * g00_2
        d += g10_1 * g10_2

        # T2
        a += h00_1 * h00_2
        b += h00_1 * h10_2
        c += h10_1 * h00_2
        d += h10_1 * h10_2

        # T1+T2
        a += gh00_1 * gh00_2
        b += gh00_1 * gh10_2
        c += gh10_1 * gh00_2
        d += gh10_1 * gh10_2

        # Now we act by (0, Q2)
        a1 = h00_2 * a + h01_2 * b
        b1 = h10_2 * a + h11_2 * b
        c1 = h00_2 * c + h01_2 * d
        d1 = h10_2 * c + h11_2 * d

        # Now we act by (P1, 0)
        a2 = g00_1 * a + g01_1 * c
        b2 = g00_1 * b + g01_1 * d
        c2 = g10_1 * a + g11_1 * c
        d2 = g10_1 * b + g11_1 * d

        # Now we act by (P1, Q2)
        a3 = g00_1 * a1 + g01_1 * c1
        b3 = g00_1 * b1 + g01_1 * d1
        c3 = g10_1 * a1 + g11_1 * c1
        d3 = g10_1 * b1 + g11_1 * d1

        return Matrix(
            [[a, b, c, d], [a1, b1, c1, d1], [a2, b2, c2, d2], [a3, b3, c3, d3]]
        )

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
        Compute the basis change on a TuplePoint to recover a ThetaPoint of
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
        Given two isotropic points of 8-torsion T1 and T2, compatible with
        the theta null point, compute the level two theta null point A/K_2
        """
        xAxByCyD = ThetaPoint.to_squared_theta(*T1)
        zAtBzYtD = ThetaPoint.to_squared_theta(*T2)

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
        den_1, den_2, den_3, den_4 = batched_inversion(num_1, num_2, num_3, num_4)

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
        a, b, c, d = ThetaPoint.to_hadamard(*ABCD)

        return ThetaStructure([a, b, c, d])

    def special_image(self, P, translate):
        """
        When the domain is a non product theta structure on a product of
        elliptic curves, we will have one of A,B,C,D=0, so the image is more
        difficult. We need to give the coordinates of P but also of
        P+Ti, Ti one of the point of 4-torsion used in the isogeny
        normalisation
        """
        AxByCzDt = ThetaPoint.to_squared_theta(*P)

        # We are in the case where at most one of A, B, C, D is
        # zero, so we need to account for this
        #
        # To recover values, we use the translated point to get
        AyBxCtDz = ThetaPoint.to_squared_theta(*translate)

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

        image = ThetaPoint.to_hadamard(*xyzt)
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

