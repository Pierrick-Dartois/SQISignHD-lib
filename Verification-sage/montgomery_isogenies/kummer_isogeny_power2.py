from collections import deque
from montgomery_isogenies.kummer_line import KummerLine
from utilities.strategy import optimised_strategy
#from utilities.supersingular import montgomery_coefficient
#from dim2_wrapper import fix_minus_sign

# =========================================== #
#   Extract coefficent from Montgomery curve  #
# =========================================== #
def montgomery_coefficient(E):
    a_inv = E.a_invariants()
    A = a_inv[1]
    if a_inv != (0, A, 0, 1, 0):
        raise ValueError(
            "Parent function assumes the curve E is in the Montgomery model"
        )
    return A


def fix_minus_sign(P, Q, PQ):
    """
    TODO: this is redundant, we have similar code in isogenies x only
    clean this up

    We have ±P, ±Q, ±(P - Q).
    This functions outputs either (-P, -Q)
    or (P, Q)
    """
    A = P.curve().a_invariants()[1]

    xP, xQ, xPQ = P[0], Q[0], PQ[0]
    yP, yQ = P[1], Q[1]

    lhs = A + xP + xQ + xPQ
    rhs = ((yQ - yP) / (xQ - xP)) ** 2

    if lhs == rhs:
        Q = -Q

    # Make sure everything now works
    assert A + P[0] + Q[0] + PQ[0] != ((Q[1] - P[1]) / (Q[0] - P[0])) ** 2

    return P, Q

def xDBL(X, Z, A24, C24):
    #input: projective coordinates xP=X/Z
    #       curve constant A24/C24 = (A/C+2)/4
    #output: projective coordinates x(2P)=X2/Z2
    t0 = X - Z
    t1 = X + Z
    t0 = t0*t0
    t1 = t1*t1
    Z2 = C24 * t0
    X2 = Z2 * t1
    t1 = t1 - t0
    t0 = A24 * t1
    Z2 = Z2 + t0
    Z2 = Z2 * t1
    return X2, Z2   #cost: 4M+2S+4a

#function for computing [2^e](X:Z) via repeated doublings
def xDBLe(XP,ZP,A24,C24,e):
    #input: projective coordinates xP=XP/ZP
    #       curve constant A24:C24
    #output: projective coordinates of x(2^e*P)=XeP/ZeP
    XeP = XP
    ZeP = ZP
    for i in range(e):
        XeP, ZeP = xDBL(XeP, ZeP, A24, C24)
    return XeP, ZeP                            #cost:4eM+2eS+(4e+3)a

def degree_2_isog(X2,Z2):
    #input: point of order 2 X2:Z2
    #output: projective coordinates of image curve E/<X2:Z2>
    A24 = X2*X2
    C24 = Z2*Z2
    A24 = C24 - A24
    return A24, C24

def eval_2_isog(X2, Z2, QX, QZ):
    #input: X2:Z2 a point of order 2 and the point QX:QZ to be pushed through the isogeny
    #output: image of QX:QZ
    t0 = X2 + Z2
    t1 = X2 - Z2
    t2 = QX + QZ
    t3 = QX - QZ
    t0 = t0*t3
    t1 = t1*t2
    t2 = t0 + t1
    t3 = t0 - t1
    QX = QX*t2
    QZ = QZ*t3
    return QX, QZ

def degree_4_isog(X4, Z4):
    #input: projective point of order four (X4:Z4)
    #output: 4-isog curve with projective coefficient A/C and 3 coefficients for evaluating
    K2 = X4 - Z4
    K3 = X4 + Z4
    K1 = Z4*Z4
    K1 = K1 + K1
    C24 = K1 * K1
    K1 = K1+K1
    A24 = X4*X4
    A24 = A24+ A24
    A24 = A24* A24
    return A24, C24, [K1,K2,K3]   #cost:5S+7a

#evaluate 4-isogenies: given coefficients from get_4_isog, evaluate at point in domain
def eval_4_isog(coeff, X, Z):
    #input: coefficients from get_4_isog
    #       projective point P=(X:Z)
    #output: projective point phi(P)=(X:Z)
    K1,K2,K3 = coeff
    t0 = X + Z
    t1 = X - Z
    X = t0*K2
    Z = t1*K3
    t0 = t0*t1
    t0 = t0*K1
    t1 = X + Z
    Z = X - Z
    t1 = t1*t1
    Z = Z*Z
    X = t0 + t1
    t0 = Z - t0
    X = X*t1
    Z = Z*t0
    return X, Z              #cost:9M+1S+6a

def isogeny_2e(A24, C24, RX, RZ, e2, points=[], strategy=[]):
    #inputs: curve constants (A24: C24) where (A24 : C24) = (A + 2C : 4C) for projective curve (A:C) and a torsion point RX, RZ of order e2
    #output: image curve, and optional image of torsion points on that curve
    r=[]
    #print(e2, strategy)

    # if e is odd, compute the first isogeny - (which must be a 2 isogeny)
    # ToOptimize: we lose all our doublings in this case (e odd)
    if (e2 % 2) == 1:
        TX, TZ = xDBLe(RX,RZ, A24, C24, e2-1)
        A24, C24 = degree_2_isog(TX, TZ)
        RX, RZ = eval_2_isog(TX, TZ, RX, RZ)
        points = [eval_2_isog(TX, TZ, PX, PZ) for (PX, PZ) in points]
        e2=e2-1

    iso_queue = deque()
    iso_queue.append((e2/2, RX, RZ))
    i = 0
    while iso_queue:
        h, X, Z = iso_queue.pop()
        #print(h, i, strategy[i])
        assert h == 1 or strategy[i] < h, "Dim 1 strategy is invalid."
        if h == 1:
            A24, C24, consts = degree_4_isog(X, Z)
            iso_queue_2 = deque()
            while iso_queue:
                h, X, Z = iso_queue.popleft()
                X, Z = eval_4_isog(consts, X, Z)
                iso_queue_2.append((h-1, X, Z))
            iso_queue = iso_queue_2
            points = [eval_4_isog(consts, PX, PZ) for (PX, PZ) in points]
        else:
           iso_queue.append((h, X, Z))
           X, Z = xDBLe(X, Z, A24, C24, 2*strategy[i])
           iso_queue.append((h-strategy[i], X, Z))
           i = i+1
    return A24, C24, points

def compute_isogeny_2e(K, P, e, strategy=None, points=()):
    if strategy is None:
        strategy = optimised_strategy(e//2, mul_c=2)

    #A24=A+2; C24=4*C
    #if P is EllipticCurvePoint_field:
    #    RX=P[0]; RZ=P[2]
    #else:
    #    RX, RZ=P
    #points = ((P[0],P[2]) for P in points) #warning: this gives (0:0) for 0_E...
    A24, C24 = K.AC24()
    RX, RZ = P

    A24new, C24new, points = isogeny_2e(A24, C24, RX, RZ, e, points=points, strategy=strategy)
    Anew=4*A24new-2*C24new; Cnew=C24new
    F=K.base_ring()
    K = KummerLine(F, (Anew, Cnew))
    points = (K(P) for P in points)
    return K, points

def isogeny_2e_basis(K, P, e, P0, Q0, PQ0=None, **kwds):
    if PQ0 is None:
        #assume P0, Q0 are elliptic points
        PQ0=P0-Q0
        points=(K(P) for P in (P0, Q0, P0-Q0))
    codom, images = compute_isogeny_2e(K, P, e, points=(P0, Q0, PQ0), **kwds)
    imE=codom.curve()
    # Now we need to lift the images back to E
    # ToOptimize: we only need imP0 on E, so we could save 2 sqrts
    images = (P.curve_point() for P in images)
    imP0, imQ0, imPQ0 = images
    imP0, imQ0 = fix_minus_sign(imP0, imQ0, imPQ0)
    return imE, imP0, imQ0

def get_isogeny_2e(E, e, c, P1, P2, alpha_beta = None, **kwds):
    F=E.base_ring(); p=F.characteristic()
    A, C = montgomery_coefficient(E), 1
    K = KummerLine(F, (A,C))
    P1, P2, P12 = K(P1), K(P2), K(P1-P2)
    A24, C24 = K.AC24()
    P=K(xDBLe(P1.X(), P1.Z(), A24, C24, e-c))

    if alpha_beta is not None: #we need to apply an isomorphism
        alpha, beta = alpha_beta
        #A=3*alpha; C=beta #No this is for when the domain has a short Weierstrass equation
        A=alpha*alpha+beta*beta; C=alpha*beta #this is what we want when the domain is in Montgomery
        K=KummerLine(F, (A,C))
        def iso(R):
            X,Z = R
            return K((X-alpha*Z, beta*Z))
            # x=X/Z
            # return K(((x-alpha)/beta, 1))
        P, P1, P2, P12=(iso(R) for R in (P, P1, P2, P12))
        ## if __debug__:
        ##     A24, C24 = K.AC24()
        ##     P_=K(xDBLe(P1.X(), P1.Z(), A24, C24, e-c))
        ##     assert P==P_

    assert (P1.X()**((p**2-1)/2) != 1) #when 1, we are above (0:1) so have a problem; this should have been fixed by the above isomorphism
    return isogeny_2e_basis(K, P, c, P1, P2, P12, **kwds)
