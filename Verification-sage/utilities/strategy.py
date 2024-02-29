# ============================================================================ #
#     Compute optimised strategy for 2-isogeny chains (in dimensions 2 and 4)  #
# ============================================================================ #

"""
The function optimised_strategy has been taken from:
https://github.com/FESTA-PKE/FESTA-SageMath

Copyright (c) 2023 Andrea Basso, Luciano Maino and Giacomo Pope.

Other functions are original work.
"""

def optimised_strategy(n, mul_c=1):
    """
    Algorithm 60: https://sike.org/files/SIDH-spec.pdf
    Shown to be appropriate for (l,l)-chains in 
    https://ia.cr/2023/508
    
    Note: the costs we consider are:
       eval_c: the cost of one isogeny evaluation
       mul_c:  the cost of one element doubling
    """

    eval_c = 1.000
    mul_c  = mul_c

    S = {1:[]}
    C = {1:0 }
    for i in range(2, n+1):
        b, cost = min(((b, C[i-b] + C[b] + b*mul_c + (i-b)*eval_c) for b in range(1,i)), key=lambda t: t[1])
        S[i] = [b] + S[i-b] + S[b]
        C[i] = cost

    return S[n]

def optimised_strategy_with_first_eval(n,mul_c=1,first_eval_c=1):
    r"""
    Adapted from optimised_strategy when the fist isogeny evaluation is more costly.
    This is well suited to gluing comptations.
    NB: This takes into account the fact that doublings on the codomain of the first isogeny
    are impossible (because of zero dual theta constants).

    CAUTION: When splittings are involved, do not use this function. Use 
    optimised_strategy_with_first_eval_and_splitting instead.

    INPUT:
    - n: number of leaves of the strategy (length of the isogeny).
    - mul_c: relative cost of one doubling compared to one generic 2-isogeny evaluation.
    - first_eval_c: relative cost of an evaluation of the first 2-isogeny (gluing) 
    compared to one generic 2-isogeny evaluation.

    OUTPUT:
    - S_left[n]: an optimal strategy represented as a sequence [s_0,...,s_{t-2}], where
    there is an index i for every internal node of the strategy, where indices are 
    ordered depth-first left-first (as the way we move on the strategy) and s_i is 
    the number of leaves to the left of internal node i 
    (see https://sike.org/files/SIDH-spec.pdf, pp. 16-17).
    """


    eval_c = 1.000
    first_eval_c = first_eval_c
    mul_c  = mul_c

    S_left = {1:[], 2:[1]} # Optimal strategies "on the left" i.e. meeting the first right edge (more costly) 
    S_right = {1:[]} # Optimal strategies "on the right" i.e. not meeting the fisrt right edge
    C_left = {1:0, 2:mul_c+first_eval_c } # Cost of strategies on the left
    C_right = {1:0 } # Cost of strategies on the right
    for i in range(2, n+1):
        # Optimisation on the right
        b, cost = min(((b, C_right[i-b] + C_right[b] + b*mul_c + (i-b)*eval_c) for b in range(1,i)), key=lambda t: t[1])
        S_right[i] = [b] + S_right[i-b] + S_right[b]
        C_right[i] = cost

    for i in range(3,n+1):
        # Optimisation on the left (b<i-1 to avoid doublings on the codomain of the first isogeny)
        b, cost = min(((b, C_left[i-b] + C_right[b] + b*mul_c + (i-1-b)*eval_c+first_eval_c) for b in range(1,i-1)), key=lambda t: t[1])
        S_left[i] = [b] + S_left[i-b] + S_right[b]
        C_left[i] = cost

    return S_left[n]

def optimised_strategy_with_first_eval_and_splitting(n,m,mul_c=1,first_eval_c=1):

    eval_c = 1.000
    first_eval_c = first_eval_c
    mul_c  = mul_c

    S_left = {1:[], 2:[1]} # Optimal strategies "on the left" i.e. meeting the first right edge (more costly) 
    S_right = {1:[]} # Optimal strategies "on the right" i.e. not meeting the fisrt right edge
    S_right_split = {1:[]} # Optimal strategies "on the right" avoiding the (n-2*m-2)-th right edge (splitting isogeny)
    C_left = {1:0, 2:mul_c+first_eval_c } # Cost of strategies on the left
    C_right = {1:0 } # Cost of strategies on the right
    C_right_split = {1:0 }
    for i in range(2, n+1):
        # Optimisation on the right
        b, cost = min(((b, C_right[i-b] + C_right[b] + b*mul_c + (i-b)*eval_c) for b in range(1,i)), key=lambda t: t[1])
        S_right[i] = [b] + S_right[i-b] + S_right[b]
        C_right[i] = cost

        if i<=m+1:
            S_right_split[i] = S_right[i]
            C_right_split[i] = cost
        else:
            b, cost = min(((b, C_right[i-b] + C_right_split[b] + b*mul_c + (i-b)*eval_c) for b in list(range(1,m+1))+list(range(m+2,i))), key=lambda t: t[1])
            S_right_split[i] = [b] + S_right[i-b] + S_right_split[b]
            C_right_split[i] = cost


    for i in range(3,n):
        # Optimisation on the left (b<i-1 to avoid doublings on the codomain of the first isogeny)
        b, cost = min(((b, C_left[i-b] + C_right[b] + b*mul_c + (i-1-b)*eval_c+first_eval_c) for b in range(1,i-1)), key=lambda t: t[1])
        S_left[i] = [b] + S_left[i-b] + S_right[b]
        C_left[i] = cost

    # Optimisation on the left (b<n-1 to avoid doublings on the codomain of the first isogeny)
    b, cost = min(((b, C_left[n-b] + C_right_split[b] + b*mul_c + (i-1-b)*eval_c+first_eval_c) for b in list(range(1,m+1))+list(range(m+2,n))), key=lambda t: t[1])
    S_left[n] = [b] + S_left[n-b] + S_right_split[b]
    C_left[n] = cost

    return S_left[n]


