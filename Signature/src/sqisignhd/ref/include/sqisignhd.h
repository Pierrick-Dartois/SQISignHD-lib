/** @file
 *
 * @authors Antonin Leroux
 *
 * @brief The protocols
 */

#ifndef SQISIGNHD_H
#define SQISIGNHD_H

#include <intbig.h>
#include <quaternion.h>
#include <klpt_constants.h>
#include <quaternion_data.h>
#include <torsion_constants.h>
#include <rng.h>
#include <endomorphism_action.h>
#include <encoded_sizes.h>
#include <fp_constants.h>

#include <stdio.h>
#include <klpt.h>
#include <id2iso.h>
#include <hd.h>
#include <dim2id2iso.h>
#include <biextension.h>

/** @defgroup sqisignhd_sqisignhd SQIsigndim2_heuristic protocols
 * @{
 */
/** @defgroup sqisignhd_t Types for SQIsignHD protocols
 * @{
 */

/** @brief Type for the signature
 *
 * @typedef signature_t
 *
 * @struct signature
 *
 */

typedef struct signature
{
    ec_curve_t E_com; /// the montgomery A coefficient for commitment curve

    // Hints to recover the basis
    int *hint_com;
    int *hint_chal;

    // Image point coefficients
    ibz_t a;
    ibz_t b;
    ibz_t c_or_d;

    // Degree of the response
    ibz_t q;
    
    /// REMOVE: this old stuff
    ///int two_resp_length;
    ///ibz_t x;
    ///int hint_b;
    ///ibz_t b0;
    ///ibz_t d0;
    ///ibz_t b1;
    ///ibz_t d1;
    ///ibz_t c0_adjust;
    ///ibz_t e0_adjust;

} signature_t;

/** @brief Type for the public keys
 *
 * @typedef public_key_t
 *
 * @struct public_key
 *
 */
typedef struct public_key
{
    ec_curve_t curve; /// the normalized A coefficient of the Montgomery curve
    int *hint_pk;
} public_key_t;

/** @brief Type for the secret keys
 *
 * @typedef secret_key_t
 *
 * @struct secret_key
 *
 */
typedef struct secret_key
{
    ec_basis_t canonical_basis; // the canonical basis of the public key curve
    ec_curve_t curve;           /// the public curve, but with little precomputations
    quat_left_ideal_t secret_ideal;
    ibz_mat_2x2_t mat_BAcan_to_BA0_two; /// mat_BA0_to_BAcan*BA0 = BAcan, where BAcan is the
                                        /// canonical basis of EA[2^e], and BA0 the image of the
                                        /// basis of E0[2^e] through the secret isogeny
} secret_key_t;

/** @}
 */

/*************************** Functions *****************************/

void protocols_keygen(public_key_t *pk, secret_key_t *sk);
int protocols_sign(signature_t *sig,
                   const public_key_t *pk,
                   secret_key_t *sk,
                   const unsigned char *m,
                   size_t l,
                   int verbose);

void public_key_init(public_key_t *pk);
void public_key_finalize(public_key_t *pk);

void secret_key_init(secret_key_t *sk);
void secret_key_finalize(secret_key_t *sk);
void secret_sig_init(signature_t *sig);
void secret_sig_finalize(signature_t *sig);

//void print_signature(const signature_t *sig);
//void print_public_key(const public_key_t *pk);
void fprint_signature(FILE *p_file, const signature_t *sig);
void fprint_public_key(FILE *p_file, const public_key_t *pk);


/** @defgroup signature The signature protocol
 * @{
 */

/** @}
 */

/** @}
 */

#endif
