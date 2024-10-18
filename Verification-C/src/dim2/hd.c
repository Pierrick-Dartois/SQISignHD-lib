#include <hd.h>
#include <assert.h>

/**
 * @brief Compute the double of the theta couple point in on the elliptic product E1E2
 *
 * @param out Output: the theta_couple_point
 * @param in the theta couple point in the elliptic product
 * @param E1E2 an elliptic product
 * in  = (P1,P2)
 * out = [2] (P1,P2)
 *
 **/
void
double_couple_point(theta_couple_point_t *out,
                    const theta_couple_point_t *in,
                    const theta_couple_curve_t *E1E2)
{
    ec_dbl(&out->P1, &in->P1, &E1E2->E1);
    ec_dbl(&out->P2, &in->P2, &E1E2->E2);
}

/**
 * @brief Compute the iterated double of the theta couple point in on the elliptic product E1E2
 *
 * @param out Output: the theta_couple_point
 * @param n : the number of iteration
 * @param in the theta couple point in the elliptic product
 * @param E1E2 an elliptic product
 * in  = (P1,P2)
 * out = [2^n] (P1,P2)
 *
 */
void
double_couple_point_iter(theta_couple_point_t *out,
                         unsigned n,
                         const theta_couple_point_t *in,
                         const theta_couple_curve_t *E1E2)
{
    if (n == 0) {
        *out = *in;
    } else {
        double_couple_point(out, in, E1E2);
        for (int i = 0; i < n - 1; i++) {
            double_couple_point(out, out, E1E2);
        }
    }
}

/**
 * @brief Compute the addition of two points in (X : Y : Z) coordinates on the elliptic product E12
 *
 * @param out Output: the theta_couple_jac_point
 * @param T1 the theta couple jac point in the elliptic product
 * @param T2 the theta couple jac point in the elliptic product
 * @param E1E2 an elliptic product
 * in  = (P1, P2), (Q1, Q2)
 * out = (P1 + Q1, P2 + Q2)
 *
 **/
void
add_couple_jac_points(theta_couple_jac_point_t *out,
                      const theta_couple_jac_point_t *T1,
                      const theta_couple_jac_point_t *T2,
                      const theta_couple_curve_t *E1E2)
{
    ADD(&out->P1, &T1->P1, &T2->P1, &E1E2->E1);
    ADD(&out->P2, &T1->P2, &T2->P2, &E1E2->E2);
}

/**
 * @brief Compute the double of the theta couple jac point in on the elliptic product E12
 *
 * @param out Output: the theta_couple_jac_point
 * @param E12 an elliptic product
 * @param in the theta couple jac point in the elliptic product
 * in = (P1,P2)
 * out = [2] (P1,P2)
 *
 */
void
double_couple_jac_point(theta_couple_jac_point_t *out,
                        const theta_couple_jac_point_t *in,
                        const theta_couple_curve_t *E1E2)
{
    DBL(&out->P1, &in->P1, &E1E2->E1);
    DBL(&out->P2, &in->P2, &E1E2->E2);
}

/**
 * @brief Compute the iterated double of the theta couple jac point in on the elliptic product E12
 *
 * @param out Output: the theta_couple_jac_point
 * @param n : the number of iteration
 * @param in the theta couple jac point in the elliptic product
 * @param E1E2 an elliptic product
 * in = (P1,P2)
 * out = [2^n] (P1,P2)
 *
 */
void
double_couple_jac_point_iter(theta_couple_jac_point_t *out,
                             unsigned n,
                             const theta_couple_jac_point_t *in,
                             const theta_couple_curve_t *E1E2)
{
    if (n == 0) {
        *out = *in;
    } else {
        double_couple_jac_point(out, in, E1E2);
        for (int i = 0; i < n - 1; i++) {
            double_couple_jac_point(out, out, E1E2);
        }
    }
}

/**
 * @brief A forgetful function which returns (X : Z) points given a pair of (X : Y : Z) points
 *
 * @param P Output: the theta_couple_point
 * @param xyP : the theta_couple_jac_point
 **/
void
couple_jac_to_xz(theta_couple_point_t *P, const theta_couple_jac_point_t *xyP)
{
    jac_to_xz(&P->P1, &xyP->P1);
    jac_to_xz(&P->P2, &xyP->P2);
}
