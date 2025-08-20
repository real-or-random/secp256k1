/*************************************************************************
 * To the extent possible under law, the author(s) have dedicated all    *
 * copyright and related and neighboring rights to the software in this  *
 * file to the public domain worldwide. This software is distributed     *
 * without any warranty. For the CC0 Public Domain Dedication, see       *
 * EXAMPLES_COPYING or https://creativecommons.org/publicdomain/zero/1.0 *
 *************************************************************************/

/** This file demonstrates how to use the recovery module to create a
  * recoverable ECDSA signature and extract the corresponding
  * public key from it.
  */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include <secp256k1.h>
#include <secp256k1_recovery.h>

#include "examples_util.h"

/* Use my_memcmp_var instead of memcmp.
 *
 * Normally, memcmp should be fine, but we use my_memcmp_var
 * here to avoid a false positive from valgrind on macOS.
 * TODO: remove this in the event the bug is fixed with valgrind in the future.
 */
static int my_memcmp_var(const void *s1, const void *s2, size_t n) {
    const unsigned char *p1 = s1, *p2 = s2;
    size_t i;

    for (i = 0; i < n; i++) {
        int diff = p1[i] - p2[i];
        if (diff != 0) {
            return diff;
        }
    }
    return 0;
}

int main(void) {
    unsigned char msg[32] = "this_could_be_the_hash_of_a_msg";
    unsigned char seckey[32];
    unsigned char recoverable_sig_ser[64];
    unsigned char serialized_pubkey[33];
    unsigned char serialized_recovered_pubkey[33];
    size_t len;
    int return_val, recovery_id;
    secp256k1_pubkey pubkey, recovered_pubkey;
    secp256k1_ecdsa_recoverable_signature recoverable_sig;
    secp256k1_ecdsa_signature normal_sig;

    /* Before we can call actual API functions, we need to create a "context". */
    secp256k1_context* ctx = secp256k1_context_create(SECP256K1_CONTEXT_NONE);

    /*** Key Generation ***/
    if (!fill_random(seckey, sizeof(seckey))) {
        return EXIT_FAILURE;
    }
    if (!secp256k1_ec_pubkey_create(ctx, &pubkey, seckey)) {
        return EXIT_FAILURE;
    }

    len = sizeof(serialized_pubkey);
    return_val = secp256k1_ec_pubkey_serialize(ctx, serialized_pubkey, &len, &pubkey, SECP256K1_EC_COMPRESSED);
    assert(return_val);
    return_val = secp256k1_ecdsa_sign_recoverable(ctx, &recoverable_sig, msg, seckey, NULL, NULL);
    assert(return_val);
    return_val = secp256k1_ecdsa_recoverable_signature_serialize_compact(ctx, recoverable_sig_ser, &recovery_id, &recoverable_sig);
    assert(return_val);
    return_val = secp256k1_ecdsa_recoverable_signature_parse_compact(ctx, &recoverable_sig, recoverable_sig_ser, recovery_id);
    assert(return_val);
    return_val = secp256k1_ecdsa_recover(ctx, &recovered_pubkey, &recoverable_sig, msg));
    assert(return_val);
    len = sizeof(serialized_recovered_pubkey);
    return_val = secp256k1_ec_pubkey_serialize(ctx, serialized_recovered_pubkey, &len, &recovered_pubkey, SECP256K1_EC_COMPRESSED);
    assert(return_val);

    /* Actual public key and recovered public key should match */
    return_val = memcmp(serialized_pubkey, serialized_recovered_pubkey, sizeof(serialized_pubkey));
    assert(return_val == 0);
    return EXIT_SUCCESS;
}
