"""Self-test for nullspace_fast: known recurrences mod p.

1. nullspace_mod_fast: explicit rank-deficient matrix, check A v = 0 mod p.
2. Fibonacci: constant-coefficient order-2, expect (r,d)=(2,0) and the
   recovered relation proportional to Z_{n+2} - Z_{n+1} - Z_n = 0.
3. Catalan (P-finite, Legendre-type): (n+2) C_{n+1} - (4n+2) C_n = 0,
   expect (r,d)=(1,1) and coefficient match up to overall scale.
4. annihilator.py --selftest gate: green run exits 0 with ALL PASS, an
   unknown flag is rejected (exit 2), and a mutant copy with one corrupted
   series term exits 1 reporting FAIL (the gate can actually fail).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from nullspace_fast import nullspace_mod_fast, find_recurrence_fast

p = 2147483629  # < 2^31, prime


def check(name, cond):
    print(f"[{'PASS' if cond else 'FAIL'}] {name}")
    if not cond:
        sys.exit(1)


# --- 1. raw nullspace ---
A = np.array([[1, 2, 3], [2, 4, 6], [0, 1, 1]], dtype=np.int64)  # rank 2
v = nullspace_mod_fast(A, p)
check("nullspace_mod_fast returns a vector", v is not None)
res = (A @ np.array(v, dtype=object)) % p
check("A v = 0 mod p", all(int(x) == 0 for x in res))
check("v nonzero", any(x != 0 for x in v))

# --- 2. Fibonacci ---
fib = [0, 1]
for _ in range(98):
    fib.append((fib[-1] + fib[-2]) % p)
rec = find_recurrence_fast(fib, p, rmax=4, dmax=2)
check("Fibonacci recurrence found", rec is not None)
r, d, w = rec
check(f"Fibonacci (r,d)=(2,0), got ({r},{d})", (r, d) == (2, 0))
# relation c0 Z_n + c1 Z_{n+1} + c2 Z_{n+2} = 0; normalize by c2
c0, c1, c2 = w[0], w[1], w[2]
inv = pow(int(c2), p - 2, p)
check("Fibonacci coeffs = (-1,-1,1) mod p",
      (c0 * inv) % p == p - 1 and (c1 * inv) % p == p - 1)

# --- 3. Catalan (P-finite) ---
cat = [1]
for n in range(99):
    cat.append(cat[-1] * (4 * n + 2) * pow(n + 2, p - 2, p) % p)
rec = find_recurrence_fast(cat, p, rmax=3, dmax=3)
check("Catalan recurrence found", rec is not None)
r, d, w = rec
check(f"Catalan (r,d)=(1,1), got ({r},{d})", (r, d) == (1, 1))
# c0(n) = w0 + w1 n, c1(n) = w2 + w3 n; expect prop. to (-(4n+2), n+2)
inv = pow(int(w[3]), p - 2, p)  # normalize c1 linear coeff to 1
nrm = [x * inv % p for x in w]
check("Catalan coeffs = (-2,-4,2,1) mod p (i.e. -(4n+2), n+2)",
      nrm == [p - 2, p - 4, 2, 1])

# --- 4. annihilator.py --selftest gate (exit-status honesty) ---
import subprocess, tempfile

cli = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                   os.pardir, 'annihilator.py')

run = subprocess.run([sys.executable, cli, '--selftest'],
                     capture_output=True, text=True)
check("--selftest exits 0 when green", run.returncode == 0)
check("--selftest reports ALL PASS", "ALL PASS" in run.stdout)

run = subprocess.run([sys.executable, cli, '--no-such-flag'],
                     capture_output=True, text=True)
check("unknown flag rejected (exit 2)", run.returncode == 2)

# mutation control: one corrupted series term must turn the gate red
# (ladder cut to d=2 to keep the mutant run fast; a[7] += 1 shifts the
#  minimal recurrence off the expected (r,s), so the rung must FAIL)
with open(cli) as fh:
    src = fh.read()
mut = src.replace("for d in (2, 3, 4, 5):", "for d in (2,):")
mut2 = mut.replace("a = moments_sc(d, 100)",
                   "a = moments_sc(d, 100); a[7] += 1")
check("mutation anchors found in source", mut != src and mut2 != mut)
with tempfile.TemporaryDirectory() as td:
    mpath = os.path.join(td, 'annihilator_mutant.py')
    with open(mpath, 'w') as fh:
        fh.write(mut2)
    run = subprocess.run([sys.executable, mpath, '--selftest'],
                         capture_output=True, text=True)
check("corrupted-series mutant exits 1", run.returncode == 1)
check("mutant output reports FAIL", "FAIL" in run.stdout)

print("ALL PASS")
