#include "isog.h"
#include <assert.h>

/*
// since we use degree 4 isogeny steps, we need to handle the odd case with care
static void
ec_eval_even_strategy(ec_curve_t *image,
                      ec_point_t *points,
                      unsigned short points_len,
                      ec_point_t *A24,
                      const ec_point_t *kernel,
                      const int isog_len)
{

    ec_kps4_t kps;

    uint8_t log2_of_e, tmp;
    digit_t e_half = (isog_len) >> 1;
    for (tmp = e_half, log2_of_e = 0; tmp > 0; tmp >>= 1, ++log2_of_e)
        ;
    log2_of_e *= 2; // In order to ensure each splits is at most size log2_of_e

    ec_point_t SPLITTING_POINTS[log2_of_e], K2;
    copy_point(&SPLITTING_POINTS[0], kernel);

    int strategy = 0, // Current element of the strategy to be used
        i, j;

    int BLOCK = 0,        // Keeps track of point order
        current = 0;      // Number of points being carried
    int XDBLs[log2_of_e]; // Number of doubles performed

    // If walk length is odd, we start with a 2-isogeny
    int is_odd = isog_len % 2;

    // if the length is long enough we normalize the first step
    // TODO: should we normalise throughout the chain too? This only
    // helps for the first step
    if (isog_len > 50) {
        ec_normalize_point(A24);
    }

    // Chain of 4-isogenies
    for (j = 0; j < (e_half - 1); j++) {
        // Get the next point of order 4
        while (BLOCK != (e_half - 1 - j)) {
            // A new split will be added
            current += 1;
            // We set the seed of the new split to be computed and saved
            copy_point(&SPLITTING_POINTS[current], &SPLITTING_POINTS[current - 1]);
            // if we copied from the very first element, then we perform one additional doubling
            if (is_odd && current == 1) {
                if (j == 0) {
                    assert(fp2_is_one(&A24->z));
                    xDBL_A24_normalized(
                        &SPLITTING_POINTS[current], &SPLITTING_POINTS[current], A24);
                } else {
                    xDBL_A24(&SPLITTING_POINTS[current], &SPLITTING_POINTS[current], A24);
                }
            }
            for (i = 0; i < 2 * STRATEGY4[TORSION_PLUS_EVEN_POWER - isog_len][strategy]; i++)
                if (j == 0) {
                    assert(fp2_is_one(&A24->z));
                    xDBL_A24_normalized(
                        &SPLITTING_POINTS[current], &SPLITTING_POINTS[current], A24);
                } else {
                    xDBL_A24(&SPLITTING_POINTS[current], &SPLITTING_POINTS[current], A24);
                }
            XDBLs[current] = STRATEGY4[TORSION_PLUS_EVEN_POWER - isog_len]
                                      [strategy]; // The number of doublings performed is saved
            BLOCK += STRATEGY4[TORSION_PLUS_EVEN_POWER - isog_len]
                              [strategy]; // BLOCK is increased by the number of doublings performed
            strategy += 1;                // Next, we move to the next element of the strategy
        }
        if (j == 0) {
            assert(current > 0);
            ec_point_t T;
            assert(fp2_is_one(&A24->z));
            xDBL_A24_normalized(&T, &SPLITTING_POINTS[current], A24);
            if (fp2_is_zero(&T.x)) {
                xisog_4_singular(&kps, A24, SPLITTING_POINTS[current], *A24);
                xeval_4_singular(
                    SPLITTING_POINTS, SPLITTING_POINTS, current, SPLITTING_POINTS[current], &kps);

                // Evaluate points
                if (points_len)
                    xeval_4_singular(points, points, points_len, SPLITTING_POINTS[current], &kps);
            } else {
                xisog_4(&kps, A24, SPLITTING_POINTS[current]);
                xeval_4(SPLITTING_POINTS, SPLITTING_POINTS, current, &kps);

                // Evaluate points
                if (points_len)
                    xeval_4(points, points, points_len, &kps);
            }
        } else {
            if (is_odd && current == 0) {
                xDBL_A24(&SPLITTING_POINTS[current], &SPLITTING_POINTS[current], A24);
            }
#ifndef NDEBUG
            // printf("%d \n",current);
            assert(!fp2_is_zero(&SPLITTING_POINTS[current].z));
            ec_point_t test;
            copy_point(&test, &SPLITTING_POINTS[current]);
            xDBL_A24(&test, &test, A24);
            assert(!fp2_is_zero(&test.z));
            xDBL_A24(&test, &test, A24);
            assert(fp2_is_zero(&test.z));
#endif
            // Evaluate 4-isogeny
            xisog_4(&kps, A24, SPLITTING_POINTS[current]);
            xeval_4(SPLITTING_POINTS, SPLITTING_POINTS, current, &kps);
            if (points_len)
                xeval_4(points, points, points_len, &kps);
        }

        BLOCK -= XDBLs[current];
        XDBLs[current] = 0;
        current -= 1;
    }
    // Final 4-isogeny
    if (is_odd) {
        current = 1;
        copy_point(&SPLITTING_POINTS[1], &SPLITTING_POINTS[0]);
        xDBL_A24(&SPLITTING_POINTS[current], &SPLITTING_POINTS[current], A24);
    }
    xisog_4(&kps, A24, SPLITTING_POINTS[current]);
    if (points_len)
        xeval_4(points, points, points_len, &kps);

    // current-=1;
    // final 2-isogeny
    if (is_odd) {
        xeval_4(SPLITTING_POINTS, SPLITTING_POINTS, 1, &kps);

#ifndef NDEBUG
        assert(!fp2_is_zero(&SPLITTING_POINTS[0].z));
        ec_point_t test;
        copy_point(&test, &SPLITTING_POINTS[0]);
        xDBL_A24(&test, &test, A24);
        assert(fp2_is_zero(&test.z));

#endif

        ec_kps2_t kps;
        xisog_2(&kps, A24, SPLITTING_POINTS[0]);
        if (points_len)
            xeval_2(points, points, points_len, &kps);
    }

    // Output curve in the form (A:C)
    A24_to_AC(image, A24);

    // TODO:
    // The curve does not have A24 normalised though
    // should we normalise it here, or do it later?
    image->is_A24_computed_and_normalized = false;
}

void
ec_eval_even(ec_curve_t *image, ec_isog_even_t *phi, ec_point_t *points, unsigned short length)
{
    ec_curve_normalize_A24(&phi->curve);
    ec_eval_even_strategy(image, points, length, &phi->curve.A24, &phi->kernel, phi->length);
}
*/

// naive implementation
void
ec_eval_small_chain(ec_curve_t *image,
                    const ec_point_t *kernel,
                    int len,
                    ec_point_t *points,
                    int len_points)
{

    ec_point_t A24;
    AC_to_A24(&A24, image);

    ec_kps2_t kps;

    ec_point_t small_K, big_K;
    copy_point(&big_K, kernel);

    for (int i = 0; i < len; i++) {
        copy_point(&small_K, &big_K);
        // small_K = big_K;
        for (int j = 0; j < len - i - 1; j++) {
            xDBL_A24(&small_K, &small_K, &A24);
        }
        if (fp2_is_zero(&small_K.x)) {
            ec_point_t B24;
            xisog_2_singular(&kps, &B24, A24);
            xeval_2_singular(&big_K, &big_K, 1, &kps);
            xeval_2_singular(points, points, len_points, &kps);
            copy_point(&A24, &B24);
        } else {
            xisog_2(&kps, &A24, small_K);
            xeval_2(&big_K, &big_K, 1, &kps);
            xeval_2(points, points, len_points, &kps);
        }
    }
    A24_to_AC(image, &A24);

    // TODO:
    // The curve does not have A24 normalised though
    // should we normalise it here, or do it later?
    image->is_A24_computed_and_normalized = false;
}

void
ec_isomorphism(ec_isom_t *isom, const ec_curve_t *from, const ec_curve_t *to)
{
    fp2_t t0, t1, t2, t3, t4;
    fp2_mul(&t0, &from->A, &to->C);
    fp2_sqr(&t0, &t0); // fromA^2toC^2
    fp2_mul(&t1, &to->A, &from->C);
    fp2_sqr(&t1, &t1); // toA^2fromC^2
    fp2_mul(&t2, &to->C, &from->C);
    fp2_sqr(&t2, &t2); // toC^2fromC^2
    fp2_add(&t3, &t2, &t2);
    fp2_add(&t2, &t3, &t2); // 3toC^2fromC^2
    fp2_sub(&t3, &t2, &t0); // 3toC^2fromC^2-fromA^2toC^2
    fp2_sub(&t4, &t2, &t1); // 3toC^2fromC^2-toA^2fromC^2
    fp2_inv(&t3);
    fp2_mul(&t4, &t4, &t3);
    fp2_sqrt(&t4); // lambda^2 constant for SW isomorphism
    fp2_sqr(&t3, &t4);
    fp2_mul(&t3, &t3, &t4); // lambda^6

    // Check sign of lambda^2, such that lambda^6 has the right sign
    fp2_sqr(&t0, &from->C);
    fp2_add(&t1, &t0, &t0);
    fp2_add(&t1, &t1, &t1);
    fp2_add(&t1, &t1, &t1);
    fp2_add(&t0, &t0, &t1); // 9fromC^2
    fp2_sqr(&t2, &from->A);
    fp2_add(&t2, &t2, &t2); // 2fromA^2
    fp2_sub(&t2, &t2, &t0);
    fp2_mul(&t2, &t2, &from->A); // -9fromC^2fromA+2fromA^3
    fp2_sqr(&t0, &to->C);
    fp2_mul(&t0, &t0, &to->C);
    fp2_mul(&t2, &t2, &t0); // toC^3* [-9fromC^2fromA+2fromA^3]
    fp2_mul(&t3, &t3, &t2); // lambda^6*(-9fromA+2fromA^3)*toC^3
    fp2_sqr(&t0, &to->C);
    fp2_add(&t1, &t0, &t0);
    fp2_add(&t1, &t1, &t1);
    fp2_add(&t1, &t1, &t1);
    fp2_add(&t0, &t0, &t1); // 9toC^2
    fp2_sqr(&t2, &to->A);
    fp2_add(&t2, &t2, &t2); // 2toA^2
    fp2_sub(&t2, &t2, &t0);
    fp2_mul(&t2, &t2, &to->A); // -9toC^2toA+2toA^3
    fp2_sqr(&t0, &from->C);
    fp2_mul(&t0, &t0, &from->C);
    fp2_mul(&t2, &t2, &t0); // fromC^3* [-9toC^2toA+2toA^3]
    if (!fp2_is_equal(&t2, &t3))
        fp2_neg(&t4, &t4);

    // Mont -> SW -> SW -> Mont
    fp2_set_one(&t0);
    fp2_add(&isom->D, &t0, &t0);
    fp2_add(&isom->D, &isom->D, &t0);
    fp2_mul(&isom->D, &isom->D, &from->C);
    fp2_mul(&isom->D, &isom->D, &to->C);
    fp2_mul(&isom->Nx, &isom->D, &t4);
    fp2_mul(&t4, &t4, &from->A);
    fp2_mul(&t4, &t4, &to->C);
    fp2_mul(&t0, &to->A, &from->C);
    fp2_sub(&isom->Nz, &t0, &t4);
}

void
ec_iso_eval(ec_point_t *P, ec_isom_t *isom)
{
    fp2_t tmp;
    fp2_mul(&P->x, &P->x, &isom->Nx);
    fp2_mul(&tmp, &P->z, &isom->Nz);
    fp2_sub(&P->x, &P->x, &tmp);
    fp2_mul(&P->z, &P->z, &isom->D);
}

void
ec_2_isog_chain(ec_2_isog_chain_t *chain,const ec_point_t *kernel, const ec_curve_t *domain, unsigned int len,
    const unsigned int *strategy)
{
    ec_kps2_t kps[len];
    ec_point_t A24[len+1];
    bool is_singular[len];

    unsigned int tmp,
        log2_of_e, // Height of the strategy tree topology
        strat_idx = 0, // Current element of the strategy to be used
        block = 0, // Keeps track of point order
        current = 0; // Number of points being carried

    for (tmp = len, log2_of_e = 0; tmp > 0; tmp >>= 1, ++log2_of_e)
        ;
    log2_of_e *= 2;

    ec_point_t kernel_elements[log2_of_e];
    unsigned int XDBLs[log2_of_e]; // Number of doubles performed

    AC_to_A24(&A24[0],domain);
    copy_point(&kernel_elements[0], kernel);

    for(int k=0;k<len;k++){
        while(block!=len-1-k){
            current += 1;

            // Append the last kernel element and performs the doublings
            copy_point(&kernel_elements[current],&kernel_elements[current-1]);
            for(int j=0;j<strategy[strat_idx];j++){
                xDBL_A24(&kernel_elements[current],&kernel_elements[current],&A24[k]);
            }

            // Update bookkeeping variables
            XDBLs[current]=strategy[strat_idx];
            block+=strategy[strat_idx];
            strat_idx+=1;
        }

        if(fp2_is_zero(&kernel_elements[current].x)){
            xisog_2_singular(&kps[k], &A24[k+1], A24[k]);
            xeval_2_singular(kernel_elements, kernel_elements, current, &kps[k]);
            is_singular[k]=true;
        }
        else{
            xisog_2(&kps[k], &A24[k+1], kernel_elements[current]);
            xeval_2(kernel_elements, kernel_elements, current, &kps[k]);
            ec_point_t test;
            xeval_2(&test, &kernel_elements[current], 1, &kps[k]);
            printf("%i %u\n",k,ec_is_zero(&test));
            is_singular[k]=false;
        }

        block -= XDBLs[current];
        XDBLs[current] = 0;
        current -= 1;
    }


    chain->kps=kps;
    chain->A24=A24;
    chain->is_singular=is_singular;
    chain->len=len;
    copy_curve(&chain->domain,domain);
    A24_to_AC(&chain->codomain,&A24[len]);
}

void 
ec_eval_2_isog_chain(ec_point_t *Q, const ec_point_t *P, const ec_2_isog_chain_t *chain)
{
    copy_point(Q,P);
    for(int i=0;i<chain->len;i++){
        if(chain->is_singular[i]){
            xeval_2_singular(Q, Q, 1, &chain->kps[i]);
        }
        else{
            xeval_2(Q, Q, 1, &chain->kps[i]);
        }
    }
}
