# ============================================================================ #
#     Compute optimised strategy for 2-isogeny chains (in dimensions 2 and 4)  #
# ============================================================================ #

"""
The function optimised_strategy has been taken from:
https://github.com/FESTA-PKE/FESTA-SageMath

Copyright (c) 2023 Andrea Basso, Luciano Maino and Giacomo Pope.

Other functions are original work.
"""

from sage.all import log, cached_function
import logging
logger = logging.getLogger(__name__)
#logger.setLevel("DEBUG")

@cached_function
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

@cached_function
def optimised_strategy_with_first_eval(n,mul_c=1,first_eval_c=1):
    r"""
    Adapted from optimised_strategy when the fist isogeny evaluation is more costly.
    This is well suited to gluing comptations. Computes optimal strategies with constraint
    at the beginning. This takes into account the fact that doublings on the codomain of 
    the first isogeny are impossible (because of zero dual theta constants).

    CAUTION: When splittings are involved, do not use this function. Use 
    optimised_strategy_with_first_eval_and_splitting instead.

    INPUT:
    - n: number of leaves of the strategy (length of the isogeny).
    - mul_c: relative cost of one doubling compared to one generic 2-isogeny evaluation.
    - first_eval_c: relative cost of an evaluation of the first 2-isogeny (gluing) 
    compared to one generic 2-isogeny evaluation.

    OUTPUT:
    - S_left[n]: an optimal strategy of depth n with constraint at the beginning
    represented as a sequence [s_0,...,s_{t-2}], where there is an index i for every 
    internal node of the strategy, where indices are ordered depth-first left-first 
    (as the way we move on the strategy) and s_i is the number of leaves to the right 
    of internal node i (see https://sike.org/files/SIDH-spec.pdf, pp. 16-17).
    """

    S_left, _, _, _ = optimised_strategies_with_first_eval_new(n,mul_c,first_eval_c)

    return S_left[n]

@cached_function
def optimised_strategies_with_first_eval_new(n,mul_c=1,first_eval_c=1):
    r"""
    Adapted from optimised_strategy when the fist isogeny evaluation is more costly.
    This is well suited to gluing comptations.

    INPUT:
    - n: number of leaves of the strategy (length of the isogeny).
    - mul_c: relative cost of one doubling compared to one generic 2-isogeny evaluation.
    - first_eval_c: relative cost of an evaluation of the first 2-isogeny (gluing) 
    compared to one generic 2-isogeny evaluation.

    OUTPUT:
    - S_left: Optimal strategies "on the left"/with higher cost at the beginning.
    - S_right: Optimal strategies "on the right" i.e. not meeting the fisrt left edge (uniform cost).
    """

    # print(f"Strategy eval: n={n}, mul_c={mul_c}, first_eval_c={first_eval_c}")

    eval_c = 1.000
    first_eval_c = first_eval_c
    mul_c  = mul_c

    S_left = {1:[], 2:[1]} # Optimal strategies "on the left" i.e. meeting the first left edge 
    S_right = {1:[]} # Optimal strategies "on the right" i.e. not meeting the first left edge
    C_left = {1:0, 2:mul_c+first_eval_c } # Cost of strategies on the left
    C_right = {1:0 } # Cost of strategies on the right
    for i in range(2, n+1):
        # Optimisation on the right
        b, cost = min(((b, C_right[i-b] + C_right[b] + b*mul_c + (i-b)*eval_c) for b in range(1,i)), key=lambda t: t[1])
        S_right[i] = [b] + S_right[i-b] + S_right[b]
        C_right[i] = cost

    for i in range(3,n+1):
        # Optimisation on the left
        b, cost = min(((b, C_left[i-b] + C_right[b] + b*mul_c + (i-1-b)*eval_c+first_eval_c) for b in range(1,i)), key=lambda t: t[1])
        S_left[i] = [b] + S_left[i-b] + S_right[b]
        C_left[i] = cost

    return S_left, S_right, C_left, C_right

@cached_function
def optimised_strategies_with_first_eval(n,mul_c=1,first_eval_c=1):
    r"""
    Deprecated: forbidding doublings on the codomain of the first isogeny is unnecessary.

    Adapted from optimised_strategy when the fist isogeny evaluation is more costly.
    This is well suited to gluing comptations. Computes optimal strategies with constraint
    at the beginning. This takes into account the fact that doublings on the codomain of 
    the first isogeny are impossible (because of zero dual theta constants).

    CAUTION: When splittings are involved, do not use this function. Use 
    optimised_strategy_with_first_eval_and_splitting instead.

    INPUT:
    - n: number of leaves of the strategy (length of the isogeny).
    - mul_c: relative cost of one doubling compared to one generic 2-isogeny evaluation.
    - first_eval_c: relative cost of an evaluation of the first 2-isogeny (gluing) 
    compared to one generic 2-isogeny evaluation.

    OUTPUT:
    - S_left: Optimal strategies "on the left"/with constraint at the beginning i.e. meeting the 
    first left edge that do not contain any left edge on the line y=sqrt(3)*(x-1).
    - S_right: Optimal strategies "on the right" i.e. not meeting the fisrt left edge (no constraint).
    """

    # print(f"Strategy eval: n={n}, mul_c={mul_c}, first_eval_c={first_eval_c}")

    eval_c = 1.000
    first_eval_c = first_eval_c
    mul_c  = mul_c

    S_left = {1:[], 2:[1]} # Optimal strategies "on the left" i.e. meeting the first left edge 
    S_right = {1:[]} # Optimal strategies "on the right" i.e. not meeting the first left edge
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

    return S_left, S_right, C_left, C_right

@cached_function
def optimised_strategy_with_first_eval_and_splitting(n,m,mul_c=1,first_eval_c=1):
    r""" Deprecated
    """
    eval_c = 1.000
    first_eval_c = first_eval_c
    mul_c  = mul_c

    S_left, S_middle, C_left, C_middle = optimised_strategies_with_first_eval(n-m,mul_c,first_eval_c)

    ## Optimization of the right part (with constraint at the end only)

    # Dictionnary of dictionnaries of translated strategies "on the right".
    # trans_S_right[d][i] is an optimal strategy of depth i 
    # without left edge on the line y=sqrt(3)*(x-(i-1-d))
    trans_S_right={}
    trans_C_right={}

    for d in range(1,m+1):
        trans_S_right[d]={1:[]}
        trans_C_right[d]={1:0}
        if d==1:
            for i in range(3,n-m+d):
                b, cost = min(((b, C_middle[i-b] + trans_C_right[1][b] + b*mul_c + (i-b)*eval_c) for b in [1]+list(range(3,i))), key=lambda t: t[1])
                trans_S_right[1][i] = [b] + S_middle[i-b] + trans_S_right[1][b]
                trans_C_right[1][i] = cost
        else:
            for i in range(2,n-m+d):
                if i!=d+1:
                    b = 1
                    cost = trans_C_right[d-b][i-b] + C_middle[b] + b*mul_c + (i-b)*eval_c
                    for k in range(2,min(i,d)):
                        cost_k = trans_C_right[d-k][i-k] + C_middle[k] + k*mul_c + (i-k)*eval_c
                        if cost_k<cost:
                            b = k
                            cost = cost_k
                    # k=d
                    if i>d:
                        cost_k = C_middle[i-d] + C_middle[d] + d*mul_c + (i-d)*eval_c
                        if cost_k<cost:
                            b = d
                            cost = cost_k
                    for k in range(d+2,i):
                        #print(d,i,k)
                        cost_k = C_middle[i-k] + trans_C_right[d][k] + k*mul_c + (i-k)*eval_c
                        if cost_k<cost:
                            b = k
                            cost = cost_k
                    if b<d:
                        trans_S_right[d][i] = [b] + trans_S_right[d-b][i-b] + S_middle[b]
                        trans_C_right[d][i] = cost
                    else:
                        trans_S_right[d][i] = [b] + S_middle[i-b] + trans_S_right[d][b]
                        trans_C_right[d][i] = cost

    ## Optimization on the left (last part) taking into account the constraints at the beginning and at the end
    for i in range(n-m+1,n+1):
        d = i-(n-m)
        b = 1
        cost = C_left[i-b] + trans_C_right[d][b] + b*mul_c + (i-1-b)*eval_c + first_eval_c
        for k in range(2,i):
            if k!=d+1 and k!=n-1:
                cost_k = C_left[i-k] + trans_C_right[d][k] + k*mul_c + (i-1-k)*eval_c + first_eval_c
                if cost_k<cost:
                    b = k
                    cost = cost_k

        S_left[i] = [b] + S_left[i-b] + trans_S_right[d][b]
        C_left[i] = cost

    return S_left[n]

@cached_function
def precompute_strategy_with_first_eval(e,m,M=1,S=0.8,I=100):
    r"""
    INPUT:
    - e: isogeny chain length.
    - m: length of the chain in dimension 2 before gluing in dimension 4.
    - M: multiplication cost.
    - S: squaring cost.
    - I: inversion cost.

    OUTPUT: Optimal strategy to compute an isogeny chain without splitting of
    length e with m steps in dimension 2 before gluing in dimension 4.
    """
    n = e - m
    eval_c = 4*(16*M+16*S)
    mul_c = 4*(32*M+32*S)+(48*M+I)/log(n*1.0)
    first_eval_c = 4*(2*I+(244+6*m)*M+(56+8*m)*S)

    return optimised_strategy_with_first_eval(n, mul_c = mul_c/eval_c, first_eval_c = first_eval_c/eval_c)

@cached_function
def precompute_strategy_with_first_eval_and_splitting(e,m,M=1,S=0.8,I=100):
    r"""
    Deprecated: forbidding doublings on the codomain of the first isogeny and on 
    the domain of the first splitting isogeny (m steps before the end) is unnecessary.
    
    INPUT:
    - e: isogeny chain length.
    - m: length of the chain in dimension 2 before gluing in dimension 4.
    - M: multiplication cost.
    - S: squaring cost.
    - I: inversion cost.

    OUTPUT: Optimal strategy to compute an isogeny chain of length e 
    with m steps in dimension 2 before gluing in dimension 4 and
    with splitting m steps before the end.
    """
    logger.debug(f"Strategy eval split: e={e}, m={m}")
    n = e - m
    eval_c = 4*(16*M+16*S)
    mul_c = 4*(32*M+32*S)+(48*M+I)/log(n*1.0)
    first_eval_c = 4*(2*I+(244+6*m)*M+(56+8*m)*S)

    return optimised_strategy_with_first_eval_and_splitting(n, m, mul_c = mul_c/eval_c, first_eval_c = first_eval_c/eval_c)
