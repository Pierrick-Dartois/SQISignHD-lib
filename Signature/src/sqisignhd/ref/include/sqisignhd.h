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

/** @defgroup sqisignhd_sqisignhd SQIsignHD protocols
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


typedef struct signature {
    ec_curve_t E_com; /// commitment curve
    ibz_mat_2x2_t mat_sigma_phichall; /// the matrix of sigma o phi_chall from canonical basis EA to canonical basis of E_com
} signature_t;

/** @brief Type for the public keys
 * 
 * @typedef public_key_t
 * 
 * @struct public_key
 * 
*/
typedef struct public_key {
	ec_curve_t curve; /// the normalized A coefficient of the Montgomery curve
} public_key_t;

/** @brief Type for the secret keys
 * 
 * @typedef secret_key_t
 * 
 * @struct secret_key
 * 
*/
typedef struct secret_key {
	ec_curve_t curve; /// the public curve
    quat_left_ideal_t secret_ideal_two;
    quat_alg_elem_t two_to_three_transporter;
    ibz_mat_2x2_t mat_BAcan_to_BA0_two; /// mat_BA0_to_BAcan*BA0 = BAcan, where BAcan is the canonical basis of EA[2^e], and BA0 the image of the basis of E0[2^e] through the secret odd isogeny
    ibz_mat_2x2_t mat_BAcan_to_BA0_three; /// mat_BA0_to_BAcan*BA0 = BAcan, where BAcan is the canonical basis of EA[2^e], and BA0 the image of the basis of E0[2^e] through the secret odd isogeny
} secret_key_t;


/** @}
*/


/*************************** Functions *****************************/

// basis_two is set to dual_three(basis of E0)
void doublepath(quat_alg_elem_t *gamma, quat_left_ideal_t *lideal_even, quat_left_ideal_t *lideal_odd, 
    ec_basis_t *basis_three,
    ec_basis_t *basis_two,
    ec_curve_t *E_target,
    int verbose);

void protocols_keygen(public_key_t *pk, secret_key_t *sk);
int protocols_sign(signature_t *sig, const public_key_t *pk, const secret_key_t *sk, const unsigned char* m, size_t l, int verbose);

void secret_key_init(secret_key_t *sk);
void secret_key_finalize(secret_key_t *sk);
void secret_sig_init(signature_t *sig);
void secret_sig_finalize(signature_t *sig);

void print_signature(const signature_t *sig);
void print_public_key(const public_key_t *pk);

/** @defgroup signature The signature protocol
 * @{
*/



/** @}
*/


/** @}
*/

#endif
