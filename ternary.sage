# multi_h_arora_ge_final.sage
# Run inside Sage

# -------------------------
# Parameters
n = 10       # ring dimension (degree)
q = 31       # modulus
B = 1        # bound on coefficients
E = [-1,0,1]    # allowed small convolution outputs
m = 4        # number of public h's to generate
# -------------------------

import random
from sage.rings.integer import Integer

# Setup Rq = Z_q[x]/(x^n - 1)
F = Integers(q)
R.<x> = PolynomialRing(F)
Rq = R.quotient(x^n - 1, 'xbar')
xbar = Rq.gen()

# Helper: random small polynomial (coeffs in 0..B)
def random_small_poly():
    coeffs = [random.randint(-B, B) for _ in range(n)]
    return sum(Integer(coeffs[i]) * xbar**i for i in range(n))

# Generate invertible g and inverse
def generate_g():
    while True:
        g = random_small_poly()
        try:
            g_inv = g.inverse_of_unit()
            return g, g_inv
        except (ZeroDivisionError, ValueError):
            continue

# Generate f
def generate_f():
    return random_small_poly()

# Build multivariate polynomial ring for x0..x_{n-1} over GF(q)
Rpoly = PolynomialRing(GF(q), [f"x{i}" for i in range(n)])
x_vars = Rpoly.gens()

# Generate secret g once
g, g_inv = generate_g()

print("Secret g (ring element):", g)
print("g coefficients:", [int(g[i]) for i in range(n)])
print("Generating", m, "public keys h_i ...\n")

# Collect all Arora-Ge polynomials from all h_i
all_arora_polys = []
all_h_coeffs = []

for idx in range(m):
    f_i = generate_f() 
    h_i = f_i * g_inv       # public key
    # extract coefficients modulo q as Python integers
    h_coeffs = [Integer(h_i[i]) for i in range(n)]
    all_h_coeffs.append(h_coeffs)
    print(f"h_{idx} coeffs:", h_coeffs)

    # build t_l for this h_i using shared Rpoly and x_vars
    for l in range(n):
        L = Rpoly(0)
        for j in range(n):
            k = (l - j) % n
            L += Integer(h_coeffs[j]) * x_vars[k]
        t_l = Rpoly(1)
        for b in E:
            t_l *= (L - Integer(b))
        all_arora_polys.append(t_l)

# Boolean constraints xi*(xi-1) = 0
bool_constraints = [x_vars[i] * (x_vars[i]-1) * (x_vars[i]+1) for i in range(n)]

g_coeffs_modq = [int(g[i]) % q for i in range(n)]

# Map modulo q coefficients back to ternary {-1,0,1}
# Using centered representative: choose the closest in {-1,0,1}
def modq_to_ternary(c, q):
    # compute minimal representative in [-floor(q/2), ceil(q/2)]
    if c > q//2:
        c -= q
    # clamp to {-1,0,1}
    if c > 1:
        return 1
    elif c < -1:
        return -1
    else:
        return c

g_coeffs = [modq_to_ternary(c, q) for c in g_coeffs_modq]

# Compute Hamming weight / sum-of-squares
h = sum([c**2 for c in g_coeffs])

print(f"\nAdding sum-of-squares constraint: sum(x_i^2) = {h}")

# Add sum-of-squares equation
sum_squares_eq = sum([x_vars[i]**2 for i in range(n)]) - h

all_equations = all_arora_polys + bool_constraints + [sum_squares_eq]


print("\nTotal equations (Arora-Ge from all h + boolean):", len(all_equations))
print(f"  ({len(all_arora_polys)} Arora-Ge polynomials from {m} h's + {len(bool_constraints)} boolean constraints)\n")

# Build ideal and compute Groebner basis (small n/m only)
I = Rpoly.ideal(all_equations)

print("Computing Groebner basis (may take time)...")
G = I.groebner_basis()

print("\n=== Gröbner basis summary ===")
print("Number of polynomials in Gröbner basis:", len(G))
for i, gp in enumerate(G):
    print(f"G[{i}] = {gp}")
    print("-" * 60)

# Attempt to get solutions (works over finite fields / small systems)
print("\nSolving variety() ...")
sols = I.variety()
print("Solutions (variable assignments):", sols)

# ---------- ROTATION CHECK AND NON-MATCHING COUNT ----------
print("\nSecret g coefficients (one-line):", [int(g[i]) for i in range(n)])

g_coeffs = [int(g[i]) for i in range(n)]

# Generate all cyclic rotations of g
rotations = [g_coeffs[i:] + g_coeffs[:i] for i in range(n)]

matching = []       # list of tuples (sol_index, sol_coeffs, rotation_index)
non_matching = []   # list of tuples (sol_index, sol_coeffs)

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

# Print summary
print("\n=== Rotation match summary ===")
print(f"Total solutions found: {len(sols)}")
print(f"Matching (rotations of g): {len(matching)}")
print(f"Non-matching: {len(non_matching)}\n")

if len(matching) > 0:
    print("Matching solutions (index in sols -> rotation index):")
    for sidx, coeffs, ridx in matching:
        print(f"  sols[{sidx}] matches rotation by {ridx} -> {coeffs}")

if len(non_matching) > 0:
    print("\nNon-matching solutions (index in sols -> coeffs):")
    for sidx, coeffs in non_matching:
        print(f"  sols[{sidx}] -> {coeffs}")

print("\nNumber of non-matching solutions:", len(non_matching))
