load("secp256k1_params.sage")

MAX_ORDER = 1000

# Set of (curve) orders we have encountered so far.
orders_done = set()

# Map from (subgroup) orders to [b, int(gen.x), int(gen.y), gen, lambda] for those subgroups.
gens = {}

# Iterate over curves of the form y^2 = x^3 + B.
for b in range(1, P):
    # There are only 6 curves (up to isomorphism) of the form y^2 = x^3 + B. Stop once we have tried all.
    if len(orders_done) == 6:
        break

    E = EllipticCurve(F, [0, b])
    print("Analyzing curve y^2 = x^3 + %i" % b)
    n = E.order()

    # Skip curves with an order we've already tried
    if n in orders_done:
        print("- Isomorphic to earlier curve")
        print()
        continue
    orders_done.add(n)

    # Skip curves isomorphic to the real secp256k1
    if n.is_pseudoprime():
        assert E.is_isomorphic(C)
        print(" - Isomorphic to secp256k1")
        print()
        continue

    print("- Finding subgroups")

    # Iterate over the generators of this curve.
    for g in E.gens():
        # Find what prime subgroups of group generated by g exist.
        g_order = g.order()
        for f, _ in g.order().factor():
            print("- Analyzing subgroup of order %i" % f)
            # Skip subgroups that are too large.
            if f < 4 or f > MAX_ORDER:
                print("  - Bad size")
                continue

            # Construct a generator for that subgroup.
            gen = g * (g_order // f)
            assert(gen.order() == f)

            # Find lambda for endomorphism. Skip if none can be found.
            lam = None
            for l in Integers(f)(1).nth_root(3, all=True):
                if int(l)*gen == E(BETA*gen[0], gen[1]):
                    lam = int(l)
                    break
            if lam is None:
                print("  - No endomorphism for this subgroup")
                break

            # Add to gens all the generators for that subgroup.
            gens.setdefault(f, [])
            gens[f].extend((b, int(P[0]), int(P[1]), P, lam) for P in [n*gen for n in range(1, f)])
            print("  - Found solutions")

    print()

def output_generator(g, name):
    print(f"#define {name} SECP256K1_GE_CONST(\\")
    print("    0x%08x, 0x%08x, 0x%08x, 0x%08x,\\" % tuple((int(g[0]) >> (32 * (7 - i))) & 0xffffffff for i in range(4)))
    print("    0x%08x, 0x%08x, 0x%08x, 0x%08x,\\" % tuple((int(g[0]) >> (32 * (7 - i))) & 0xffffffff for i in range(4, 8)))
    print("    0x%08x, 0x%08x, 0x%08x, 0x%08x,\\" % tuple((int(g[1]) >> (32 * (7 - i))) & 0xffffffff for i in range(4)))
    print("    0x%08x, 0x%08x, 0x%08x, 0x%08x\\" % tuple((int(g[1]) >> (32 * (7 - i))) & 0xffffffff for i in range(4, 8)))
    print(")")

def output_b(b):
    print(f"#define SECP256K1_B {int(b)}")

print()
print("To be put in src/group_impl.h:")
print()
print("/* Begin of section generated by sage/gen_exhaustive_groups.sage. */")
for f in sorted(gens.keys()):
    # Use as generator/2 the one with lowest b, and lowest (x, y) generator (interpreted as non-negative integers).
    b, _, _, HALF_G, lam = min(gens[f])
    output_generator(2 * HALF_G, f"SECP256K1_G_ORDER_{f}")
print("/** Generator for secp256k1, value 'g' defined in")
print(" *  \"Standards for Efficient Cryptography\" (SEC2) 2.7.1.")
print(" */")
output_generator(G, "SECP256K1_G")
print("/* These exhaustive group test orders and generators are chosen such that:")
print(" * - The field size is equal to that of secp256k1, so field code is the same.")
print(" * - The curve equation is of the form y^2=x^3+B for some small constant B.")
print(" * - The subgroup has a generator 2*P, where P.x is as small as possible.")
print(f" * - The subgroup has size less than {MAX_ORDER} to permit exhaustive testing.")
print(" * - The subgroup admits an endomorphism of the form lambda*(x,y) == (beta*x,y).")
print(" */")
print("#if defined(EXHAUSTIVE_TEST_ORDER)")
first = True
for f in sorted(gens.keys()):
    b, _, _, _, lam = min(gens[f])
    print(f"#  {'if' if first else 'elif'} EXHAUSTIVE_TEST_ORDER == {f}")
    first = False
    print()
    print(f"static const secp256k1_ge secp256k1_ge_const_g = SECP256K1_G_ORDER_{f};")
    output_b(b)
    print()
print("#  else")
print("#    error No known generator for the specified exhaustive test group order.")
print("#  endif")
print("#else")
print()
print("static const secp256k1_ge secp256k1_ge_const_g = SECP256K1_G;")
output_b(7)
print()
print("#endif")
print("/* End of section generated by sage/gen_exhaustive_groups.sage. */")


print()
print()
print("To be put in src/scalar_impl.h:")
print()
print("/* Begin of section generated by sage/gen_exhaustive_groups.sage. */")
first = True
for f in sorted(gens.keys()):
    _, _, _, _, lam = min(gens[f])
    print("#  %s EXHAUSTIVE_TEST_ORDER == %i" % ("if" if first else "elif", f))
    first = False
    print("#    define EXHAUSTIVE_TEST_LAMBDA %i" % lam)
print("#  else")
print("#    error No known lambda for the specified exhaustive test group order.")
print("#  endif")
print("/* End of section generated by sage/gen_exhaustive_groups.sage. */")
