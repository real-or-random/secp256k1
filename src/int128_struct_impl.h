#ifndef SECP256K1_INT128_64_IMPL_H
#define SECP256K1_INT128_64_IMPL_H

#include "int128.h"

#if defined(_M_X64) | defined(_M_ARM64) | defined(_WIN64) /* MSVC */
  #include <intrin.h>
  #define secp256k1_umulh __umulh
  #define secp256k1_mulh __mulh
#else
static SECP256K1_INLINE uint64_t secp256k1_umulh(uint64_t a, uint64_t b) {
   uint64_t t1 = (uint64_t)(uint32_t)a * (uint32_t)b;
   uint64_t t2 = (a >> 32) * (uint32_t)b;
   uint64_t t3 = (uint32_t)a*(b >> 32) + (t1 >> 32) + (uint32_t)t2;
   return (a >> 32)*(b >> 32) + (t2 >> 32) + (t3 >> 32);
}

static SECP256K1_INLINE int64_t secp256k1_mulh(int64_t a, int64_t b) {
   uint64_t t1 = (uint64_t)(uint32_t)a * (uint32_t)b;
   int64_t t2 = (a >> 32) * (uint32_t)b;
   int64_t t3 = (uint32_t)a * (b >> 32);
   uint64_t t4 = (t1 >> 32) + (uint32_t)t2 + (uint32_t)t3;
   return (a >> 32) * (b >> 32) + (t2 >> 32) + (t3 >> 32) + (int64_t)(t4 >> 32);
}
#endif

static SECP256K1_INLINE void secp256k1_u128_mul(secp256k1_uint128 *r, uint64_t a, uint64_t b) {
   r->hi = secp256k1_umulh(a, b);
   r->lo = a * b;
}

static SECP256K1_INLINE void secp256k1_u128_accum_mul(secp256k1_uint128 *r, uint64_t a, uint64_t b) {
   uint64_t lo = a * b;
   r->hi += secp256k1_umulh(a, b) + (~lo < r->lo);
   r->lo += lo;
}

static SECP256K1_INLINE void secp256k1_u128_accum_u64(secp256k1_uint128 *r, uint64_t a) {
   r->hi += (r->lo > ~a);
   r->lo += a;
}

/* Unsigned (logical) right shift.
 * Non-constant time in n.
 */
static SECP256K1_INLINE void secp256k1_u128_rshift(secp256k1_uint128 *r, unsigned int n) {
   if (n >= 64) {
     r->lo = (r->hi) >> (n-64);
     r->hi = 0;
   } else if (n > 0) {
     r->lo = ((r->hi & (((uint64_t)1<<n)-1)) << (64-n)) | r->lo >> n;
     r->hi >>= n;
   }
}

static SECP256K1_INLINE uint64_t secp256k1_u128_to_u64(const secp256k1_uint128 *a) {
   return a->lo;
}

static SECP256K1_INLINE uint64_t secp256k1_u128_hi_u64(const secp256k1_uint128 *a) {
   return a->hi;
}

static SECP256K1_INLINE void secp256k1_u128_from_u64(secp256k1_uint128 *r, uint64_t a) {
   r->hi = 0;
   r->lo = a;
}

static SECP256K1_INLINE int secp256k1_u128_check_bits(secp256k1_uint128 *r, unsigned int n) {
   return n >= 64 ? r->hi >> (n - 64) == 0
                  : r->hi == 0 && r->lo >> n == 0;
}

static SECP256K1_INLINE void secp256k1_i128_mul(secp256k1_int128 *r, int64_t a, int64_t b) {
   r->hi = (uint64_t)secp256k1_mulh(a, b);
   r->lo = (uint64_t)a * (uint64_t)b;
}

static SECP256K1_INLINE void secp256k1_i128_accum_mul(secp256k1_int128 *r, int64_t a, int64_t b) {
   uint64_t lo = (uint64_t)a * (uint64_t)b;
   uint64_t hi = (uint64_t)secp256k1_mulh(a, b) + (~lo < r->lo);
   VERIFY_CHECK((r->hi <= 0x7fffffffffffffffu && hi <= 0x7fffffffffffffffu) <= (r->hi + hi <= 0x7fffffffffffffffu));
   VERIFY_CHECK((r->hi > 0x7fffffffffffffffu && hi > 0x7fffffffffffffffu) <= (r->hi + hi > 0x7fffffffffffffffu));
   r->hi += hi;
   r->lo += lo;
}

static SECP256K1_INLINE void secp256k1_i128_dissip_mul(secp256k1_int128 *r, int64_t a, int64_t b) {
   uint64_t lo = (uint64_t)a * (uint64_t)b;
   uint64_t hi = (uint64_t)secp256k1_mulh(a, b) + (r->lo < lo);
   VERIFY_CHECK((r->hi <= 0x7fffffffffffffffu && hi > 0x7fffffffffffffffu) <= (r->hi - hi <= 0x7fffffffffffffffu));
   VERIFY_CHECK((r->hi > 0x7fffffffffffffffu && hi <= 0x7fffffffffffffffu) <= (r->hi - hi > 0x7fffffffffffffffu));
   r->hi -= hi;
   r->lo -= lo;
}

static SECP256K1_INLINE void secp256k1_i128_det(secp256k1_int128 *r, int64_t a, int64_t b, int64_t c, int64_t d) {
   secp256k1_i128_mul(r, a, d);
   secp256k1_i128_dissip_mul(r, b, c);
}

/* Signed (arithmetic) right shift.
 * Non-constant time in b.
 */
static SECP256K1_INLINE void secp256k1_i128_rshift(secp256k1_int128 *r, unsigned int b) {
   if (b >= 64) {
     r->lo = (uint64_t)((int64_t)(r->hi) >> (b-64));
     r->hi = (uint64_t)((int64_t)(r->hi) >> 63);
   } else if (b > 0) {
     r->lo = ((r->hi & (((uint64_t)1<<b)-1)) << (64-b)) | r->lo >> b;
     r->hi = (uint64_t)((int64_t)(r->hi) >> b);
   }
}

static SECP256K1_INLINE int64_t secp256k1_i128_to_i64(const secp256k1_int128 *a) {
   return (int64_t)a->lo;
}

static SECP256K1_INLINE void secp256k1_i128_from_i64(secp256k1_int128 *r, int64_t a) {
   r->hi = (uint64_t)(a >> 63);
   r->lo = (uint64_t)a;
}

static SECP256K1_INLINE int secp256k1_i128_eq(const secp256k1_int128 *a, const secp256k1_int128 *b) {
   return a->hi == b->hi && a->lo == b->lo;
}

static SECP256K1_INLINE int secp256k1_i128_check_bit(secp256k1_int128 *r, unsigned int n) {
   return n >= 64 ? r->hi == (uint64_t)1 << (n - 64) && r->lo == 0
                  : r->hi == 0 && r->lo == (uint64_t)1 << n;
}

#endif
