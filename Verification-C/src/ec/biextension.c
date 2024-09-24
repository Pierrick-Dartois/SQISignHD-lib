#include <biextension.h>
#include <mp.h>
#include <assert.h>
#include <stdlib.h>

/*
 * We implement the biextension arithmetic by using the cubical torsor representation. For now only
 * implement the 2^e-ladder.
 *
 * Warning: both cubicalDBL and cubicalADD are off by a factor x4 with respect
 * to the cubical arithmetic.
 * Since this factor is the same, this means that the biextension
 * arithmetic is correct, so the pairings are ok (they only rely on the
 * biextension arithmetic).
 * In the special case where Q=P (self pairings), we use the cubical ladder
 * rather than the biextension ladder because this is faster. In that case,
 * when we do a ladder we are off by a factor 4^m, m the number of bits.
 * This factor thus disappear in the Weil pairing since we take a quotient,
 * and also in the Tate pairing due to the final exponentiation; so
 * everything is ok too.
 * (Note that when the curves are supersingular as in our case, the Tate
 * self pairing is always trivial anyway because the Galois structure of the
 * isogeneous curves are all the same, so the étale torsor representing the
 * Tate pairing has to be trivial).
 */

// this is exactly like xDBL_A24, except we use the fact that P is normalised
// to gain a multiplication
// Cost: 3M + 2S + 2a + 2s
static void
cubicalDBL(ec_point_t *Q, const ec_point_t *P, const ec_point_t *A24)
{
    // A24 = (A+2C/4C: 1)
    assert(fp2_is_one(&A24->z));

    fp2_t t0, t1, t2;

    fp2_add(&t0, &P->x, &P->z);
    fp2_sqr(&t0, &t0);
    fp2_sub(&t1, &P->x, &P->z);
    fp2_sqr(&t1, &t1);
    fp2_sub(&t2, &t0, &t1);
    // fp2_mul(&t1, &t1, &A24->z);
    fp2_mul(&Q->x, &t0, &t1);
    fp2_mul(&t0, &t2, &A24->x);
    fp2_add(&t0, &t0, &t1);
    fp2_mul(&Q->z, &t0, &t2);
}

// this would be exactly like xADD if PQ was 'antinormalised' as (1,z)
// Cost: 3M + 2S + 3a + 3s
static void
cubicalADD(ec_point_t *R, const ec_point_t *P, const ec_point_t *Q, const fp2_t *ixPQ)
{
    fp2_t t0, t1, t2, t3;

    fp2_add(&t0, &P->x, &P->z);
    fp2_sub(&t1, &P->x, &P->z);
    fp2_add(&t2, &Q->x, &Q->z);
    fp2_sub(&t3, &Q->x, &Q->z);
    fp2_mul(&t0, &t0, &t3);
    fp2_mul(&t1, &t1, &t2);
    fp2_add(&t2, &t0, &t1);
    fp2_sub(&t3, &t0, &t1);
    fp2_sqr(&R->z, &t3);
    fp2_sqr(&t2, &t2);
    fp2_mul(&R->x, ixPQ, &t2);
}

// Given cubical reps of P, Q and x(P - Q) = (1 : ixPQ)
// compute P + Q, [2]Q
// Cost: 6M + 4S + 4a + 4s
static void
cubicalDBLADD(ec_point_t *PpQ,
              ec_point_t *QQ,
              const ec_point_t *P,
              const ec_point_t *Q,
              const fp2_t *ixPQ,
              const ec_point_t *A24)
{
    // A24 = (A+2C/4C: 1)
    assert(fp2_is_one(&A24->z));

    fp2_t t0, t1, t2, t3;

    fp2_add(&t0, &P->x, &P->z);
    fp2_sub(&t1, &P->x, &P->z);
    fp2_add(&PpQ->x, &Q->x, &Q->z);
    fp2_sub(&t3, &Q->x, &Q->z);
    fp2_sqr(&t2, &PpQ->x);
    fp2_sqr(&QQ->z, &t3);
    fp2_mul(&t0, &t0, &t3);
    fp2_mul(&t1, &t1, &PpQ->x);
    fp2_add(&PpQ->x, &t0, &t1);
    fp2_sub(&t3, &t0, &t1);
    fp2_sqr(&PpQ->z, &t3);
    fp2_sqr(&PpQ->x, &PpQ->x);
    fp2_mul(&PpQ->x, ixPQ, &PpQ->x);
    fp2_sub(&t3, &t2, &QQ->z);
    fp2_mul(&QQ->x, &t2, &QQ->z);
    fp2_mul(&t0, &t3, &A24->x);
    fp2_add(&t0, &t0, &QQ->z);
    fp2_mul(&QQ->z, &t0, &t3);
}

// iterative biextension doubling
static void
biext_ladder_2e(uint64_t e,
                ec_point_t *PnQ,
                ec_point_t *nQ,
                const ec_point_t *PQ,
                const ec_point_t *Q,
                const fp2_t *ixP,
                const ec_point_t *A24)
{
    copy_point(PnQ, PQ);
    copy_point(nQ, Q);
    for (uint64_t i = 0; i < e; i++) {
        cubicalDBLADD(PnQ, nQ, PnQ, nQ, ixP, A24);
    }
}

// iterative biextension double and add
static void
biext_ladder(uint64_t* n, const unsigned int nwords,
    ec_point_t *PnQ,
    ec_point_t *nQ,
    const ec_point_t *PQ,
    const ec_point_t *Q,
    const fp2_t *ixP,
    const fp2_t *ixQ,
    const ec_point_t *A24, const ec_point_t *P)
{
    // non-cubical points (to test)
    ec_point_t ncnm1Q, ncnQ, ncPnQ;

    unsigned int bj,j,k;
    uint64_t mask;
    uint64_t one[nwords], nm1[nwords];
    ec_point_t nm1Q;

    mp_set_small(one,1,nwords);
    mp_sub(nm1,n,one,nwords);
    uint16_t nbits=mp_nbits(nm1, nwords);

    copy_point(PnQ, PQ);
    copy_point(nQ, Q);
    copy_point(&nm1Q, Q);
    cubicalDBLADD(PnQ, nQ, PnQ, nQ, ixP, A24);
    for(int i=1;i<nbits;i++){
        j=nbits-i-1;
        k=j%RADIX;
        mask=1ULL<<k;
        bj=(mask&nm1[j/RADIX])>>k;

        if(bj){
            cubicalADD(&nm1Q, nQ, &nm1Q, ixQ);
            cubicalDBLADD(PnQ, nQ, PnQ, nQ, ixP, A24);
        }
        else{
            cubicalADD(nQ, nQ, &nm1Q, ixQ);
            cubicalDBLADD(PnQ, &nm1Q, PnQ, &nm1Q, ixP, A24);
        }
    }
    assert(ec_is_equal(&nm1Q,Q));

    // Operations on non-cubical points
    copy_point(&ncPnQ, PQ);
    copy_point(&ncnQ, Q);
    copy_point(&ncnm1Q, Q);
    //cubicalDBLADD(PnQ, nQ, PnQ, nQ, ixP, A24);
    xDBLADD(&ncnQ, &ncPnQ, &ncnQ, &ncPnQ, P, A24);
    for(int i=1;i<nbits;i++){
        j=nbits-i-1;
        k=j%RADIX;
        mask=1ULL<<k;
        bj=(mask&nm1[j/RADIX])>>k;

        if(bj){
            //cubicalADD(&nm1Q, nQ, &nm1Q, ixQ);
            //cubicalDBLADD(PnQ, nQ, PnQ, nQ, ixP, A24);
            xADD(&ncnm1Q, &ncnQ, &ncnm1Q, Q);
            xDBLADD(&ncnQ, &ncPnQ, &ncnQ, &ncPnQ, P, A24);
        }
        else{
            //cubicalADD(nQ, nQ, &nm1Q, ixQ);
            //cubicalDBLADD(PnQ, &nm1Q, PnQ, &nm1Q, ixP, A24);
            xADD(&ncnQ, &ncnQ, &ncnm1Q, Q);
            xDBLADD(&ncnm1Q, &ncPnQ, &ncnm1Q, &ncPnQ, P, A24);
        }
    }
    printf("%u\n",ec_is_zero(&ncnQ));
    printf("%u\n",ec_is_equal(&ncPnQ,P));
}

// Compute the ratio X/Z above as a (X:Z) point to avoid a division
static void
point_ratio(ec_point_t *R, const ec_point_t *PnQ, const ec_point_t *nQ, const ec_point_t *P)
{
    // Sanity tests
    assert(ec_is_zero(nQ));
    assert(ec_is_equal(PnQ, P));

    fp2_mul(&R->x, &nQ->x, &P->x);
    fp2_copy(&R->z, &PnQ->x);
}

// Compute the cubical translation of P by a point of 2-torsion T
static void
translate(ec_point_t *P, const ec_point_t *T)
{
    // When we translate, the following three things can happen:
    // T = (A : 0) then the translation of P should be P
    // T = (0 : B) then the translation of P = (X : Z) should be (Z : X)
    // Otherwise T = (A : B) and P translates to (AX - BZ : BX - AZ)
    // We compute this in constant time by computing the generic case
    // and then using constant time swaps.
    fp2_t PX_new, PZ_new;

    {
        fp2_t t0, t1;

        // PX_new = AX - BZ
        fp2_mul(&t0, &T->x, &P->x);
        fp2_mul(&t1, &T->z, &P->z);
        fp2_sub(&PX_new, &t0, &t1);

        // PZ_new = BX - AZ
        fp2_mul(&t0, &T->z, &P->x);
        fp2_mul(&t1, &T->x, &P->z);
        fp2_sub(&PZ_new, &t0, &t1);
    }

    // When we have A zero we should return (Z : X)
    uint32_t TA_is_zero = fp2_is_zero(&T->x);
    fp2_select(&PX_new, &PX_new, &P->z, TA_is_zero);
    fp2_select(&PZ_new, &PZ_new, &P->x, TA_is_zero);

    // When we have B zero we should return (X : Z)
    uint32_t TB_is_zero = fp2_is_zero(&T->z);
    fp2_select(&PX_new, &PX_new, &P->x, TB_is_zero);
    fp2_select(&PZ_new, &PZ_new, &P->z, TB_is_zero);

    // Set the point to the desired result
    fp2_copy(&P->x, &PX_new);
    fp2_copy(&P->z, &PZ_new);
}

// Compute the monodromy P+2^e Q (in level 1)
// The suffix _i means that we are given 1/x(P) as parameter.
// Warning: to get meaningful result when using the monodromy to compute
// pairings, we need P, Q, PQ, A24 to be normalised
// (this is not strictly necessary, but care need to be taken when they are not normalised. Only
// handle the normalised case for now)
static void
monodromy_i(ec_point_t *R, const weil_params_t *weil_data, bool swap_PQ)
{
    fp2_t ixP;
    ec_point_t P, Q, PnQ, nQ;

    // When we compute the Weil pairing we need both P + [2^e]Q and
    // Q + [2^e]P which we can do easily with biext_ladder_2e() below
    // we use a bool to decide wether to use Q, ixP or P, ixQ in the
    // ladder and P or Q in translation.
    if (!swap_PQ) {
        copy_point(&P, &weil_data->P);
        copy_point(&Q, &weil_data->Q);
        fp2_copy(&ixP, &weil_data->ixP);
    } else {
        copy_point(&P, &weil_data->Q);
        copy_point(&Q, &weil_data->P);
        fp2_copy(&ixP, &weil_data->ixQ);
    }

    // Compute the biextension ladder P + [2^e]Q
    biext_ladder_2e(weil_data->e - 1, &PnQ, &nQ, &weil_data->PQ, &Q, &ixP, &weil_data->A24);
    translate(&PnQ, &nQ);
    translate(&nQ, &nQ);
    point_ratio(R, &PnQ, &nQ, &P);
}

// Compute the monodromy of P+[n]Q with n odd
static void
monodromy_odd_i(ec_point_t *R, const weil_params_t *weil_data, bool swap_PQ)
{
    fp2_t ixP, ixQ;
    ec_point_t P, Q, PnQ, nQ;
    const unsigned int nwords=weil_data->nwords;
    uint64_t m[nwords];
    mp_copy(m,weil_data->n,nwords);

    // When we compute the Weil pairing we need both P + [n]Q and
    // Q + [n]P which we can do easily with biext_ladder() below
    // we use a bool to decide wether to use Q, ixP or P, ixQ in the
    // ladder and P or Q in translation.
    if (!swap_PQ) {
        copy_point(&P, &weil_data->P);
        copy_point(&Q, &weil_data->Q);
        fp2_copy(&ixP, &weil_data->ixP);
        fp2_copy(&ixQ, &weil_data->ixQ);
    } else {
        copy_point(&P, &weil_data->Q);
        copy_point(&Q, &weil_data->P);
        fp2_copy(&ixP, &weil_data->ixQ);
        fp2_copy(&ixQ, &weil_data->ixP);
    }

    // Compute the biextension ladder P + [n]Q
    biext_ladder(m, nwords, &PnQ, &nQ, &weil_data->PQ, &Q, &ixP, &ixQ, &weil_data->A24, &P);
    point_ratio(R, &PnQ, &nQ, &P);
}

// Normalize the points and also store 1/x(P), 1/x(Q)
static void
cubical_normalization(weil_params_t *weil_data, const ec_point_t *P, const ec_point_t *Q)
{
    fp2_t t[4];
    fp2_copy(&t[0], &P->x);
    fp2_copy(&t[1], &P->z);
    fp2_copy(&t[2], &Q->x);
    fp2_copy(&t[3], &Q->z);
    fp2_batched_inv(t, 4);

    // Store PZ / PX and QZ / QX
    fp2_mul(&weil_data->ixP, &P->z, &t[0]);
    fp2_mul(&weil_data->ixQ, &Q->z, &t[2]);

    // Store x(P), x(Q) normalised to (X/Z : 1)
    fp2_mul(&weil_data->P.x, &P->x, &t[1]);
    fp2_mul(&weil_data->Q.x, &Q->x, &t[3]);
    fp2_set_one(&weil_data->P.z);
    fp2_set_one(&weil_data->Q.z);
}

// Weil pairing, PQ should be P+Q in (X:Z) coordinates
// We assume the points are normalised correctly
static void
weil_2e_n(fp2_t *r, weil_params_t *weil_data)
{
    ec_point_t R0, R1;
    monodromy_i(&R0, weil_data, true);
    monodromy_i(&R1, weil_data, false);

    // TODO: check if that's the Weil pairing or its inverse
    // Seems OK (see paper http://www.normalesup.org/~robert/pro/publications/articles/biextensions.pdf p. 57)
    fp2_mul(r, &R0.x, &R1.z);
    fp2_inv(r);
    fp2_mul(r, r, &R0.z);
    fp2_mul(r, r, &R1.x);
}

// Weil pairing, PQ should be P+Q in (X:Z) coordinates
// Normalise the points and call the code above
// The code will crash (division by 0) if either P or Q is (0:1)
void
weil_2e(fp2_t *r, uint64_t e, ec_point_t *P, ec_point_t *Q, ec_point_t *PQ, ec_curve_t *E)
{
    weil_params_t weil_data;
    // Construct the structure for the Weil pairing
    // Set (PX/PZ : 1), (QX : QZ : 1), PZ/PX and QZ/QX
    weil_data.e = e;
    cubical_normalization(&weil_data, P, Q);
    copy_point(&weil_data.PQ, PQ);

    // Ensure the input curve has A24 normalised and store
    // in a struct
    ec_curve_normalize_A24(E);
    copy_point(&weil_data.A24, &E->A24);

    // Compute the Weil pairing e_(2^n)(P, Q)
    weil_2e_n(r, &weil_data);
}


// For odd values of n, this is the Weil pairing for 2(O_E), we need 
// to power it to the exponent (n+1)/2 (i.e. take a square root in mu_n) 
// to get the Weil pairing for (O_E). In practice it does not matter (for dlogs). 
static void
weil_n(fp2_t *r, weil_params_t *weil_data,  bool odd)
{
    ec_point_t R0, R1;
    if(odd){
        monodromy_odd_i(&R0, weil_data, true);
        monodromy_odd_i(&R1, weil_data, false);
    }
    else{
        printf("Not implemented for n values of n.\n");
        exit(EXIT_FAILURE);
    }

    // TODO: check if that's the Weil pairing or its inverse
    // Seems OK (see paper http://www.normalesup.org/~robert/pro/publications/articles/biextensions.pdf p. 57)
    fp2_mul(r, &R0.x, &R1.z);
    fp2_inv(r);
    fp2_mul(r, r, &R0.z);
    fp2_mul(r, r, &R1.x);
}

void
weil(fp2_t *r, uint64_t *n, unsigned int nwords, ec_point_t *P, ec_point_t *Q, ec_point_t *PQ, ec_curve_t *E)
{
    weil_params_t weil_data;
    // Construct the structure for the Weil pairing
    // Set (PX/PZ : 1), (QX : QZ : 1), PZ/PX and QZ/QX
    weil_data.n=n;
    weil_data.nwords=nwords;
    cubical_normalization(&weil_data, P, Q);
    copy_point(&weil_data.PQ, PQ);

    // Ensure the input curve has A24 normalised and store
    // in a struct
    ec_curve_normalize_A24(E);
    copy_point(&weil_data.A24, &E->A24);

    // Compute the Weil pairing e_n(P, Q)
    weil_n(r, &weil_data, n[0]&1);
}

// Functions to compute discrete logs by computing the Weil pairing of points
// followed by computing the dlog in Fp^2
// TODO:
// Is it ever cheaper to work with full order points and then solve the dlog with
// the Tate pairing?

// recursive dlog function
static bool
fp2_dlog_2e_rec(digit_t *a, long len, fp2_t *pows_f, fp2_t *pows_g, long stacklen)
{
    if (len == 0) {
        // *a = 0;
        for (int i = 0; i < NWORDS_ORDER; i++) {
            a[i] = 0;
        }
        return true;
    } else if (len == 1) {
        if (fp2_is_one(&pows_f[stacklen - 1])) {
            // a = 0;
            for (int i = 0; i < NWORDS_ORDER; i++) {
                a[i] = 0;
            }
            for (int i = 0; i < stacklen - 1; ++i) {
                fp2_sqr(&pows_g[i], &pows_g[i]); // new_g = g^2
            }
            return true;
        } else if (fp2_is_equal(&pows_f[stacklen - 1], &pows_g[stacklen - 1])) {
            // a = 1;
            a[0] = 1;
            for (int i = 1; i < NWORDS_ORDER; i++) {
                a[i] = 0;
            }
            fp2_t tmp;
            for (int i = 0; i < stacklen - 1; ++i) {
                fp2_mul(&pows_f[i], &pows_f[i], &pows_g[i]); // new_f = f*g
                fp2_sqr(&pows_g[i], &pows_g[i]);             // new_g = g^2
            }
            return true;
        } else {
            return false;
        }
    } else {
        long right = (double)len * 0.5;
        long left = len - right;
        pows_f[stacklen] = pows_f[stacklen - 1];
        pows_g[stacklen] = pows_g[stacklen - 1];
        for (int i = 0; i < left; i++) {
            fp2_sqr(&pows_f[stacklen], &pows_f[stacklen]);
            fp2_sqr(&pows_g[stacklen], &pows_g[stacklen]);
        }
        // uint64_t dlp1 = 0, dlp2 = 0;
        digit_t dlp1[NWORDS_ORDER], dlp2[NWORDS_ORDER];
        bool ok;
        ok = fp2_dlog_2e_rec(dlp1, right, pows_f, pows_g, stacklen + 1);
        if (!ok)
            return false;
        ok = fp2_dlog_2e_rec(dlp2, left, pows_f, pows_g, stacklen);
        if (!ok)
            return false;
        // a = dlp1 + 2^right * dlp2
        multiple_mp_shiftl(dlp2, right, NWORDS_ORDER);
        mp_add(a, dlp2, dlp1, NWORDS_ORDER);

        return true;
    }
}

// compute DLP: compute scal such that f = g^scal with f, 1/g as input
static bool
fp2_dlog_2e(digit_t *scal, const fp2_t *f, const fp2_t *g_inverse, int e)
{
    long log, len = e;
    for (log = 0; len > 1; len >>= 1)
        log++;
    log += 1;

    fp2_t pows_f[log], pows_g[log];
    pows_f[0] = *f;
    pows_g[0] = *g_inverse;

    for (int i = 0; i < NWORDS_ORDER; i++) {
        scal[i] = 0;
    }

    bool ok = fp2_dlog_2e_rec(scal, e, pows_f, pows_g, 1);
    assert(ok);

    return ok;
}

// Normalize the bases (P, Q), (R, S) and store their inverse
// and additionally normalise the curve to (A/C : 1)
static void
cubical_normalization_dlog(weil_dlog_params_t *weil_dlog_data, ec_curve_t *curve)
{
    fp2_t t[9];
    ec_basis_t *PQ = &weil_dlog_data->PQ;
    ec_basis_t *RS = &weil_dlog_data->RS;
    fp2_copy(&t[0], &PQ->P.x);
    fp2_copy(&t[1], &PQ->P.z);
    fp2_copy(&t[2], &PQ->Q.x);
    fp2_copy(&t[3], &PQ->Q.z);
    fp2_copy(&t[4], &RS->P.x);
    fp2_copy(&t[5], &RS->P.z);
    fp2_copy(&t[6], &RS->Q.x);
    fp2_copy(&t[7], &RS->Q.z);
    fp2_copy(&t[8], &curve->C);

    fp2_batched_inv(t, 9);

    fp2_mul(&weil_dlog_data->ixP, &PQ->P.z, &t[0]);
    fp2_mul(&PQ->P.x, &PQ->P.x, &t[1]);
    fp2_set_one(&PQ->P.z);

    fp2_mul(&weil_dlog_data->ixQ, &PQ->Q.z, &t[2]);
    fp2_mul(&PQ->Q.x, &PQ->Q.x, &t[3]);
    fp2_set_one(&PQ->Q.z);

    fp2_mul(&weil_dlog_data->ixR, &RS->P.z, &t[4]);
    fp2_mul(&RS->P.x, &RS->P.x, &t[5]);
    fp2_set_one(&RS->P.z);

    fp2_mul(&weil_dlog_data->ixS, &RS->Q.z, &t[6]);
    fp2_mul(&RS->Q.x, &RS->Q.x, &t[7]);
    fp2_set_one(&RS->Q.z);

    fp2_mul(&curve->A, &curve->A, &t[8]);
    fp2_set_one(&curve->C);
}

// Given two bases <P, Q> and basis = <R, S> compute
// x(P - R), x(P - S), x(R - Q), x(S - Q)
static void
compute_difference_points(weil_dlog_params_t *weil_dlog_data, ec_curve_t *curve)
{
    jac_point_t xyP, xyQ, xyR, xyS, temp;

    // lifting the two basis points, assumes that x(P) and x(R)
    // and the curve itself are normalised to (X : 1)
    lift_basis_normalized(&xyP, &xyQ, &weil_dlog_data->PQ, curve);
    lift_basis_normalized(&xyR, &xyS, &weil_dlog_data->RS, curve);

    // computation of the differences
    // x(P - R)
    jac_neg(&temp, &xyR);
    ADD(&temp, &temp, &xyP, curve);
    jac_to_xz(&weil_dlog_data->diff.PmR, &temp);

    // x(P - S)
    jac_neg(&temp, &xyS);
    ADD(&temp, &temp, &xyP, curve);
    jac_to_xz(&weil_dlog_data->diff.PmS, &temp);

    // x(R - Q)
    jac_neg(&temp, &xyQ);
    ADD(&temp, &temp, &xyR, curve);
    jac_to_xz(&weil_dlog_data->diff.RmQ, &temp);

    // x(S - Q)
    jac_neg(&temp, &xyQ);
    ADD(&temp, &temp, &xyS, curve);
    jac_to_xz(&weil_dlog_data->diff.SmQ, &temp);
}

// Inline all the Weil pairing computations done in ec_dlog_2_weil
static void
weil_dlog(digit_t *r1, digit_t *r2, digit_t *s1, digit_t *s2, weil_dlog_params_t *weil_dlog_data)
{

    ec_point_t nP, nQ, nR, nS, nPQ, PnQ, nPR, PnR, nPS, PnS, nRQ, RnQ, nSQ, SnQ;

    copy_point(&nP, &weil_dlog_data->PQ.P);
    copy_point(&nQ, &weil_dlog_data->PQ.Q);
    copy_point(&nR, &weil_dlog_data->RS.P);
    copy_point(&nS, &weil_dlog_data->RS.Q);
    copy_point(&nPQ, &weil_dlog_data->PQ.PmQ);
    copy_point(&PnQ, &weil_dlog_data->PQ.PmQ);
    copy_point(&nPR, &weil_dlog_data->diff.PmR);
    copy_point(&nPS, &weil_dlog_data->diff.PmS);
    copy_point(&PnR, &weil_dlog_data->diff.PmR);
    copy_point(&PnS, &weil_dlog_data->diff.PmS);
    copy_point(&nRQ, &weil_dlog_data->diff.RmQ);
    copy_point(&nSQ, &weil_dlog_data->diff.SmQ);
    copy_point(&RnQ, &weil_dlog_data->diff.RmQ);
    copy_point(&SnQ, &weil_dlog_data->diff.SmQ);

    for (uint64_t i = 0; i < weil_dlog_data->e - 1; i++) {
        cubicalADD(&nPQ, &nPQ, &nP, &weil_dlog_data->ixQ);
        cubicalADD(&nPR, &nPR, &nP, &weil_dlog_data->ixR);
        cubicalDBLADD(&nPS, &nP, &nPS, &nP, &weil_dlog_data->ixS, &weil_dlog_data->A24);

        cubicalADD(&PnQ, &PnQ, &nQ, &weil_dlog_data->ixP);
        cubicalADD(&RnQ, &RnQ, &nQ, &weil_dlog_data->ixR);
        cubicalDBLADD(&SnQ, &nQ, &SnQ, &nQ, &weil_dlog_data->ixS, &weil_dlog_data->A24);

        cubicalADD(&PnR, &PnR, &nR, &weil_dlog_data->ixP);
        cubicalDBLADD(&nRQ, &nR, &nRQ, &nR, &weil_dlog_data->ixQ, &weil_dlog_data->A24);

        cubicalADD(&PnS, &PnS, &nS, &weil_dlog_data->ixP);
        cubicalDBLADD(&nSQ, &nS, &nSQ, &nS, &weil_dlog_data->ixQ, &weil_dlog_data->A24);
    }

    // weil(&w0,e,&PQ->P,&PQ->Q,&PQ->PmQ,&A24);
    translate(&nPQ, &nP);
    translate(&nPR, &nP);
    translate(&nPS, &nP);
    translate(&PnQ, &nQ);
    translate(&RnQ, &nQ);
    translate(&SnQ, &nQ);
    translate(&PnR, &nR);
    translate(&nRQ, &nR);
    translate(&PnS, &nS);
    translate(&nSQ, &nS);

    translate(&nP, &nP);
    translate(&nQ, &nQ);
    translate(&nR, &nR);
    translate(&nS, &nS);

    // computation of the reference weil pairing
    ec_point_t T0, T1;
    fp2_t w1[5], w2[5];

    // e(P, Q) = w0
    point_ratio(&T0, &nPQ, &nP, &weil_dlog_data->PQ.Q);
    point_ratio(&T1, &PnQ, &nQ, &weil_dlog_data->PQ.P);
    // For the first element we need it's inverse for
    // fp2_dlog_2e so we swap w1 and w2 here to save inversions
    fp2_mul(&w2[0], &T0.x, &T1.z);
    fp2_mul(&w1[0], &T1.x, &T0.z);

    // e(P,R) = w0^r2
    point_ratio(&T0, &nPR, &nP, &weil_dlog_data->RS.P);
    point_ratio(&T1, &PnR, &nR, &weil_dlog_data->PQ.P);
    fp2_mul(&w1[1], &T0.x, &T1.z);
    fp2_mul(&w2[1], &T1.x, &T0.z);

    // e(R,Q) = w0^r1
    point_ratio(&T0, &nRQ, &nR, &weil_dlog_data->PQ.Q);
    point_ratio(&T1, &RnQ, &nQ, &weil_dlog_data->RS.P);
    fp2_mul(&w1[2], &T0.x, &T1.z);
    fp2_mul(&w2[2], &T1.x, &T0.z);

    // e(P,S) = w0^s2
    point_ratio(&T0, &nPS, &nP, &weil_dlog_data->RS.Q);
    point_ratio(&T1, &PnS, &nS, &weil_dlog_data->PQ.P);
    fp2_mul(&w1[3], &T0.x, &T1.z);
    fp2_mul(&w2[3], &T1.x, &T0.z);

    // e(S,Q) = w0^s1
    point_ratio(&T0, &nSQ, &nS, &weil_dlog_data->PQ.Q);
    point_ratio(&T1, &SnQ, &nQ, &weil_dlog_data->RS.Q);
    fp2_mul(&w1[4], &T0.x, &T1.z);
    fp2_mul(&w2[4], &T1.x, &T0.z);

    fp2_batched_inv(w1, 5);
    for (int i = 0; i < 5; i++) {
        fp2_mul(&w1[i], &w1[i], &w2[i]);
    }

    fp2_dlog_2e(r2, &w1[1], &w1[0], weil_dlog_data->e);
    fp2_dlog_2e(r1, &w1[2], &w1[0], weil_dlog_data->e);
    fp2_dlog_2e(s2, &w1[3], &w1[0], weil_dlog_data->e);
    fp2_dlog_2e(s1, &w1[4], &w1[0], weil_dlog_data->e);
}

void
ec_dlog_2_weil(digit_t *r1,
               digit_t *r2,
               digit_t *s1,
               digit_t *s2,
               ec_basis_t *PQ,
               ec_basis_t *RS,
               ec_curve_t *curve,
               int e)
{
    assert(test_point_order_twof(&PQ->Q, curve, e));

    // precomputing the correct curve data
    ec_curve_normalize_A24(curve);

    weil_dlog_params_t weil_dlog_data;
    weil_dlog_data.e = e;
    weil_dlog_data.PQ = *PQ;
    weil_dlog_data.RS = *RS;
    weil_dlog_data.A24 = curve->A24;

    cubical_normalization_dlog(&weil_dlog_data, curve);
    compute_difference_points(&weil_dlog_data, curve);

    weil_dlog(r1, r2, s1, s2, &weil_dlog_data);

#ifndef NDEBUG
    ec_point_t test;
    ec_biscalar_mul(&test, r1, r2, e, PQ, curve);
    // R = [r1]P + [r2]Q
    assert(ec_is_equal(&test, &RS->P));
    ec_biscalar_mul(&test, s1, s2, e, PQ, curve);
    // S = [s1]P + [s2]Q
    assert(ec_is_equal(&test, &RS->Q));
#endif
}
