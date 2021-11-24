/***********************************************************************
 * Copyright (c) 2013, 2014, 2015 Thomas Daede, Cory Fields            *
 * Distributed under the MIT software license, see the accompanying    *
 * file COPYING or https://www.opensource.org/licenses/mit-license.php.*
 ***********************************************************************/

/* Autotools creates libsecp256k1-config.h, of which ECMULT_GEN_PREC_BITS is needed.
   ifndef guard so downstream users can define their own if they do not use autotools. */
#if !defined(ECMULT_GEN_PREC_BITS)
#include "libsecp256k1-config.h"
#endif

/* In principle we could use ASM, but this yields only a minor speedup in
   build time and it's very complicated. In particular when cross-compiling, we'd
   need to build the ASM for the build and the host machine. */
#undef USE_EXTERNAL_ASM
#undef USE_ASM_X86_64

#include "../include/secp256k1.h"
#include "assumptions.h"
#include "util.h"
#include "group.h"
#include "ecmult_gen.h"
#include "ecmult_gen_prec_impl.h"

int main(int argc, char **argv) {
    secp256k1_ge_storage* table;
    int inner;
    int outer;
    const char outfile[] = "src/ecmult_gen_static_prec_table.h";
    FILE* fp;

    int bits = ECMULT_GEN_PREC_B;
    int g = ECMULT_GEN_PREC_G(bits);
    int n = ECMULT_GEN_PREC_N(bits);
    table = checked_malloc(&default_error_callback, n * g * sizeof(secp256k1_ge_storage));

    (void)argc;
    (void)argv;

    fp = fopen(outfile, "w");
    if (fp == NULL) {
        fprintf(stderr, "Could not open %s for writing!\n", outfile);
        return -1;
    }

    fprintf(fp, "#ifndef SECP256K1_ECMULT_GEN_STATIC_PREC_TABLE_H\n");
    fprintf(fp, "#define SECP256K1_ECMULT_GEN_STATIC_PREC_TABLE_H\n");

    fprintf(fp, "#include \"src/group.h\"\n");

    fprintf(fp, "#define SC SECP256K1_GE_STORAGE_CONST\n");

    fprintf(fp, "#if ECMULT_GEN_PREC_BITS != %d\n", bits);
    fprintf(fp, "   #error configuration mismatch, invalid ECMULT_GEN_PREC_BITS. Try deleting ecmult_static_context.h before the build.\n");
    fprintf(fp, "#endif\n");

    fprintf(fp, "static\n");
    fprintf(fp, "#ifndef EXHAUSTIVE_TEST_ORDER\n");
    fprintf(fp, "const\n");
    fprintf(fp, "#endif\n");
    fprintf(fp, "secp256k1_ge_storage secp256k1_ecmult_gen_prec_table[ECMULT_GEN_PREC_N(ECMULT_GEN_PREC_B)][ECMULT_GEN_PREC_G(ECMULT_GEN_PREC_B)] = {\n");

    table = checked_malloc(&default_error_callback, ECMULT_GEN_PREC_TABLE_SIZE);
    secp256k1_ecmult_gen_create_prec_table(table, &secp256k1_ge_const_g, bits);

    for(outer = 0; outer != n; outer++) {
        fprintf(fp,"{\n");
        for(inner = 0; inner != g; inner++) {
            fprintf(fp,"    SC(%uu, %uu, %uu, %uu, %uu, %uu, %uu, %uu, %uu, %uu, %uu, %uu, %uu, %uu, %uu, %uu)", SECP256K1_GE_STORAGE_CONST_GET(table[outer * g + inner]));
            if (inner != g - 1) {
                fprintf(fp,",\n");
            } else {
                fprintf(fp,"\n");
            }
        }
        if (outer != n - 1) {
            fprintf(fp,"},\n");
        } else {
            fprintf(fp,"}\n");
        }
    }
    fprintf(fp,"};\n");
    free(table);

    fprintf(fp, "#undef SC\n");
    fprintf(fp, "#endif\n");
    fclose(fp);

    return 0;
}
