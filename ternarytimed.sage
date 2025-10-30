# multi_h_arora_ge_final_timed_ternary.sage
# Run inside Sage (sage -python or from within sage cell)

from time import perf_counter
import random
from sage.rings.integer import Integer

# -------------------------
# Parameters (tune for your machine)
n = 13       # ring dimension (degree)
q = 31       # modulus
B = 1        # bound on coefficients (we'll draw from [-B, B] for ternary)
E = [-1, 0, 1]  # allowed small convolution outputs (ternary)
m = 10        # number of public h's to generate
# -------------------------

random.seed(42)  # reproducible runs

# Setup Rq = Z_q[x]/(x^n - 1)
F = Integers(q)
R.<x> = PolynomialRing(F)
Rq = R.quotient(x^n - 1, 'xbar')
xbar = Rq.gen()

# Helper: random small polynomial (coeffs in [-B..B])
def random_small_poly():
    coeffs = [random.randint(-B, B) for _ in range(n)]
    return sum(Integer(coeffs[i]) * xbar**i for i in range(n))

# Generate invertible g and inverse
def generate_g():
    while True:
        g = random_small_poly()
        try:
            # In Sage, .inverse_of_unit() available on ring elements that are units
            g_inv = g.inverse_of_unit()
            return g, g_inv
        except (ZeroDivisionError, ValueError, Exception):
            # keep trying until unit found
            continue

# Generate f
def generate_f():
    return random_small_poly()

# Build multivariate polynomial ring for x0..x_{n-1} over GF(q)
var_names = [f"x{i}" for i in range(n)]
Rpoly = PolynomialRing(GF(q), var_names)
x_vars = Rpoly.gens()

# Generate secret g once
print("Generating secret g ...")
g, g_inv = generate_g()
print("Secret g (ring element):", g)
print("g coefficients:", [int(g[i]) for i in range(n)])

# Pre-generate public h's and their Arora–Ge polynomials (we'll reuse for both runs)
all_h_coeffs = []
pre_arora_polys = []
for idx in range(m):
    f_i = generate_f()
    h_i = f_i * g_inv       # public key
    h_coeffs = [Integer(h_i[i]) for i in range(n)]
    all_h_coeffs.append(h_coeffs)

    # build t_l polynomials for this h_i
    for l in range(n):
        L = Rpoly(0)
        for j in range(n):
            k = (l - j) % n
            L += Integer(h_coeffs[j]) * x_vars[k]
        t_l = Rpoly(1)
        for b in E:
            t_l *= (L - Integer(b))
        pre_arora_polys.append(t_l)

# Ternary constraints: xi*(xi-1)*(xi+1) = 0  (forces xi in {-1,0,1})
ternary_constraints = [x_vars[i] * (x_vars[i] - 1) * (x_vars[i] + 1) for i in range(n)]

print("\nPrepared:")
print(f"  {m} public h's, total Arora–Ge polynomials: {len(pre_arora_polys)}")
print(f"  Ternary constraints: {len(ternary_constraints)}\n")

# We'll run two experiments: without ternary constraints, and with ternary constraints
for include_ternary in (False, True):
    label = "with ternary constraints" if include_ternary else "without ternary constraints"
    print("="*72)
    print(f"Running experiment: {label}")

    # Build equation set
    if include_ternary:
        all_equations = pre_arora_polys + ternary_constraints
    else:
        all_equations = pre_arora_polys

    print(f"  Total equations: {len(all_equations)}")

    # Build ideal
    I = Rpoly.ideal(all_equations)

    # Time Groebner basis computation
    t0 = perf_counter()
    try:
        G = I.groebner_basis()
        t1 = perf_counter()
        gb_time = t1 - t0
        print(f"  Gröbner basis computed in {gb_time:.6f} seconds")
        print(f"  Number of polynomials in Gröbner basis: {len(G)}")
    except Exception as e:
        t1 = perf_counter()
        gb_time = t1 - t0
        print(f"  Gröbner basis computation FAILED after {gb_time:.6f} seconds: {e}")
        G = None

    # Time variety (solution) computation
    t2 = perf_counter()
    try:
        sols = I.variety()
        t3 = perf_counter()
        var_time = t3 - t2
        print(f"  variety() computed in {var_time:.6f} seconds")
        print(f"  Number of solutions found: {len(sols)}")
    except Exception as e:
        t3 = perf_counter()
        var_time = t3 - t2
        print(f"  variety() FAILED after {var_time:.6f} seconds: {e}")
        sols = []

    # Print brief sample of Groebner basis (if available)
    if G is not None:
        max_show = 8
        print("  Sample Gröbner basis polys (up to first", max_show, "):")
        for i, gp in enumerate(G[:max_show]):
            print(f"    G[{i}] = {gp}")
        if len(G) > max_show:
            print(f"    ... (and {len(G)-max_show} more)")

    # Print up to 6 solutions (as coefficient vectors)
    if sols:
        print("  Sample solutions (up to 6):")
        for si, s in enumerate(sols[:6]):
            sol_coeffs = [int(s[x_vars[i]]) for i in range(n)]
            print(f"    sol[{si}] = {sol_coeffs}")
        if len(sols) > 6:
            print(f"    ... (and {len(sols)-6} more)")

    # Rotation check vs secret g if solutions present
    if sols:
        g_coeffs = [int(g[i]) for i in range(n)]
        rotations = [g_coeffs[i:] + g_coeffs[:i] for i in range(n)]
        matching = []
        non_matching = []
        for idx, sol in enumerate(sols):
            sol_coeffs = [int(sol[x_vars[i]]) for i in range(n)]
            matched = False
            for r_idx, rot in enumerate(rotations):
                if sol_coeffs == rot:
                    matching.append((idx, sol_coeffs, r_idx))
                    matched = True
                    break
            if not matched:
                non_matching.append((idx, sol_coeffs))

        print(f"  Matching rotations of g: {len(matching)}")
        print(f"  Non-matching solutions: {len(non_matching)}")

    # Summary timings
    total_time = gb_time + var_time
    print(f"  Timing summary (seconds): groebner = {gb_time:.6f}, variety = {var_time:.6f}, total = {total_time:.6f}\n")

print("All experiments finished.")
# End of script
