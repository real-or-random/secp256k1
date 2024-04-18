/***********************************************************************
 * Copyright (c) 2021 Russell O'Connor, Jonas Nick                     *
 * Distributed under the MIT software license, see the accompanying    *
 * file COPYING or https://www.opensource.org/licenses/mit-license.php.*
 ***********************************************************************/

#ifndef SECP256K1_HSORT_IMPL_H
#define SECP256K1_HSORT_IMPL_H

#include "hsort.h"

/* An array is a heap when, for all non-zero indexes i, the element at index i
 * compares as less than or equal to the element at index parent(i) = (i-1)/2.
 */

static SECP256K1_INLINE size_t secp256k1_heap_child1(size_t i) {
    VERIFY_CHECK(i <= (SIZE_MAX - 1)/2);
    return 2*i + 1;
}

static SECP256K1_INLINE size_t secp256k1_heap_child2(size_t i) {
    VERIFY_CHECK(i <= SIZE_MAX/2 - 1);
    return secp256k1_heap_child1(i)+1;
}

static SECP256K1_INLINE void secp256k1_heap_swap64(unsigned char *a, size_t i, size_t j, size_t stride, size_t swap_size) {
    unsigned char tmp[64];
    VERIFY_CHECK(swap_size <= 64);
    memcpy(tmp, a + i*stride, swap_size);
    memmove(a + i*stride, a + j*stride, swap_size);
    memcpy(a + j*stride, tmp, swap_size);
}

/* Swap the elements of a at indices i and j, assuming that the size of each element is stride. */
static SECP256K1_INLINE void secp256k1_heap_swap(unsigned char *a, size_t i, size_t j, size_t stride) {
    size_t remaining = stride;
    while (64 < remaining) {
        secp256k1_heap_swap64(a + (remaining - 64), i, j, stride, 64);
        remaining -= 64;
    }
    secp256k1_heap_swap64(a, i, j, stride, remaining);
}

static SECP256K1_INLINE void secp256k1_heap_down(unsigned char *a, size_t i, size_t heap_size, size_t stride,
                            int (*cmp)(const void *, const void *, void *), void *cmp_data) {
    while (i < heap_size/2) {
        VERIFY_CHECK(i <= SIZE_MAX/2 - 1);
        /* Proof:
         * i < heap_size/2
         * i + 1 <= heap_size/2
         * 2*i + 2 <= heap_size <= SIZE_MAX
         * 2*i <= SIZE_MAX - 2
         */

        VERIFY_CHECK(secp256k1_heap_child1(i) < heap_size);
        /* Proof:
         * i < heap_size/2
         * i + 1 <= heap_size/2
         * 2*i + 2 <= heap_size
         * 2*i + 1 < heap_size
         * child1(i) < heap_size
         */

        /* Let [x] be notation for the contents at a[x*stride].
         *
         * If [child1(i)] > [i] and [child2(i)] > [i],
         *   swap [i] with the larger child to ensure the new parent is larger
         *   than both children. When [child1(i)] == [child2(i)], swap [i] with
         *   [child2(i)].
         * Else if [child1(i)] > [i], swap [i] with [child1(i)].
         * Else if [child2(i)] > [i], swap [i] with [child2(i)].
         */
        if (secp256k1_heap_child2(i) < heap_size
                && 0 <= cmp(a + secp256k1_heap_child2(i)*stride, a + secp256k1_heap_child1(i)*stride, cmp_data)) {
            if (0 < cmp(a + secp256k1_heap_child2(i)*stride, a +         i*stride, cmp_data)) {
                secp256k1_heap_swap(a, i, secp256k1_heap_child2(i), stride);
                i = secp256k1_heap_child2(i);
            } else {
                /* At this point we have [child2(i)] >= [child1(i)] and we have
                 * [child2(i)] <= [i], and thus [child1(i)] <= [i] which means
                 * that the next comparison can be skipped. */
                return;
            }
        } else if (0 < cmp(a + secp256k1_heap_child1(i)*stride, a +         i*stride, cmp_data)) {
            secp256k1_heap_swap(a, i, secp256k1_heap_child1(i), stride);
            i = secp256k1_heap_child1(i);
        } else {
            return;
        }
    }
    /* heap_size/2 <= i
     * heap_size/2 < i + 1
     * heap_size < 2*i + 2
     * heap_size <= 2*i + 1
     * heap_size <= child1(i)
     * Thus child1(i) and child2(i) are now out of bounds and we are at a leaf.
     */
}

/* In-place heap sort. */
static void secp256k1_hsort(void *ptr, size_t count, size_t size,
                            int (*cmp)(const void *, const void *, void *),
                            void *cmp_data ) {
    size_t i;

    for(i = count/2; 0 < i; --i) {
        secp256k1_heap_down(ptr, i-1, count, size, cmp, cmp_data);
    }
    for(i = count; 1 < i; --i) {
        /* Extract the largest value from the heap */
        secp256k1_heap_swap(ptr, 0, i-1, size);

        /* Repair the heap condition */
        secp256k1_heap_down(ptr, 0, i-1, size, cmp, cmp_data);
    }
}

#endif
