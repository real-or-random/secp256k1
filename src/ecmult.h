/***********************************************************************
 * Copyright (c) 2013, 2014, 2017 Pieter Wuille, Andrew Poelstra       *
 * Distributed under the MIT software license, see the accompanying    *
 * file COPYING or https://www.opensource.org/licenses/mit-license.php.*
 ***********************************************************************/

#ifndef SECP256K1_ECMULT_H
#define SECP256K1_ECMULT_H

#include "group.h"
#include "scalar.h"
#include "scratch.h"

#ifndef ECMULT_WINDOW_SIZE
#  define ECMULT_WINDOW_SIZE 15
#  ifdef DEBUG_CONFIG
#     pragma message DEBUG_CONFIG_MSG("ECMULT_WINDOW_SIZE undefined, assuming default value")
#  endif
#endif

#ifdef DEBUG_CONFIG
#  pragma message DEBUG_CONFIG_DEF(ECMULT_WINDOW_SIZE)
#endif

/** Larger values for ECMULT_WINDOW_SIZE result in possibly better
 *  performance at the cost of an exponentially larger precomputed
 *  table. The exact table size is
 *      (1 << (WINDOW_G - 2)) * sizeof(secp256k1_ge_storage)  bytes,
 *  where sizeof(secp256k1_ge_storage) is typically 64 bytes but can
 *  be larger due to platform-specific padding and alignment.
 *  Two tables of this size are used (due to the endomorphism
 *  optimization).
 *
 * Noone will ever need more than a window size of 24. The code might
 * be correct for larger values of ECMULT_WINDOW_SIZE but this is not
 * tested.
 *
 * The following limitations are known, and there are probably more:
 * If WINDOW_G > 27 and size_t has 32 bits, then the code is incorrect
 * because the size of the memory object that we allocate (in bytes)
 * will not fit in a size_t.
 * If WINDOW_G > 31 and int has 32 bits, then the code is incorrect
 * because certain expressions will overflow.
 */
#if ECMULT_WINDOW_SIZE < 2 || ECMULT_WINDOW_SIZE > 24
#  error Set ECMULT_WINDOW_SIZE to an integer in range [2..24].
#endif

#if defined(EXHAUSTIVE_TEST_ORDER)
/* We need to lower these values for exhaustive tests because
 * the tables cannot have infinities in them (this breaks the
 * affine-isomorphism stuff which tracks z-ratios) */
#  if EXHAUSTIVE_TEST_ORDER > 128
#    define WINDOW_A 5
#    define WINDOW_G 8
#  elif EXHAUSTIVE_TEST_ORDER > 8
#    define WINDOW_A 4
#    define WINDOW_G 4
#  else
#    define WINDOW_A 2
#    define WINDOW_G 2
#  endif
#else
#  define WINDOW_A 5 /* optimal for 128-bit and 256-bit exponents. */
#  define WINDOW_G ECMULT_WINDOW_SIZE
#endif

/** The number of entries a table with precomputed multiples needs to have. */
#define ECMULT_TABLE_SIZE(w) (1L << ((w)-2))

/** Double multiply: R = na*A + ng*G */
static void secp256k1_ecmult(secp256k1_gej *r, const secp256k1_gej *a, const secp256k1_scalar *na, const secp256k1_scalar *ng);

typedef int (secp256k1_ecmult_multi_callback)(secp256k1_scalar *sc, secp256k1_ge *pt, size_t idx, void *data);

/**
 * Multi-multiply: R = inp_g_sc * G + sum_i ni * Ai.
 * Chooses the right algorithm for a given number of points and scratch space
 * size. Resets and overwrites the given scratch space. If the points do not
 * fit in the scratch space the algorithm is repeatedly run with batches of
 * points. If no scratch space is given then a simple algorithm is used that
 * simply multiplies the points with the corresponding scalars and adds them up.
 * Returns: 1 on success (including when inp_g_sc is NULL and n is 0)
 *          0 if there is not enough scratch space for a single point or
 *          callback returns 0
 */
static int secp256k1_ecmult_multi_var(const secp256k1_callback* error_callback, secp256k1_scratch *scratch, secp256k1_gej *r, const secp256k1_scalar *inp_g_sc, secp256k1_ecmult_multi_callback cb, void *cbdata, size_t n);

#endif /* SECP256K1_ECMULT_H */
