
# ============================================ #
#     Fast square root and quadratic roots     #
# ============================================ #

"""
Most of this code has been taken from:
https://github.com/FESTA-PKE/FESTA-SageMath

Copyright (c) 2023 Andrea Basso, Luciano Maino and Giacomo Pope.

Functions with another Copyright mention are not from the above authors.
"""

def sqrt_Fp2(a):
    """
    Efficiently computes the sqrt
    of an element in Fp2 using that
    we always have a prime p such that
    p ≡ 3 mod 4.
    """
    Fp2 = a.parent()
    p = Fp2.characteristic()
    i = Fp2.gen() # i = √-1

    a1 = a ** ((p - 3) // 4)
    x0 = a1 * a
    alpha = a1 * x0

    if alpha == -1:
        x = i * x0
    else:
        b = (1 + alpha) ** ((p - 1) // 2)
        x = b * x0

    return x

def n_sqrt(a, n):
    for _ in range(n):
        a = sqrt_Fp2(a)
    return a

def sqrt_Fp(a):
    """
    Efficiently computes the sqrt
    of an element in Fp using that
    we always have a prime p such that
    p ≡ 3 mod 4.

    Copyright (c) Pierrick Dartois 2025.
    """
    Fp = a.parent()
    p = Fp.characteristic()

    return a**((p+1)//4)
