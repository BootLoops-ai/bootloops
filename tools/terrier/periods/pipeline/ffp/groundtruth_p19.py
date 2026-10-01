# Ground truth: CDEvS 1912.06146 tab:p19 — R(T) for all phi in F_19, AESZ34.
# R(T) = 1 + a T + b p T^2 + a p^3 T^3 + p^6 T^4  (smooth fibres).
# Factored rows given as (alpha, beta) meaning (1 - alpha*p*T + p^3 T^2)(1 - beta*T + p^3 T^2).
# Singular rows given as ('sing', c) meaning (1 - p T)(1 - c T + p^3 T^2).
P = 19
ROWS = {
    1:  ('sing', 20),
    2:  ('ab',  4*P, 2),      # a = 4p, b = 2
    3:  ('ab',  -8, 242),
    4:  ('fac', -4, 60),      # (1+4pT+p^3T^2)(1-60T+p^3T^2): alpha=-4, beta=60
    5:  ('fac', -4, 60),
    6:  ('ab',   8, -318),
    7:  ('ab', -44, -238),
    8:  ('fac',  2, 80),      # (1-2pT+..)(1-80T+..)
    9:  ('fac', -4, 160),
    10: ('ab',  12, 562),
    11: ('fac', -4, 140),
    12: ('ab',  12, 82),
    13: ('ab', 178, 1082),
    14: ('ab',  12, -158),
    15: ('ab',  42, -2*P),    # 1+42T-2p^2T^2+...: b*p = -2p^2 -> b = -2p
    16: ('sing', -76),
    17: ('sing', 20),
    18: ('ab', -54, 322),
}

def ab(row):
    """Return (a, b) with R = 1 + aT + b p T^2 + a p^3 T^3 + p^6 T^4, or None if singular."""
    t = ROWS[row]
    if t[0] == 'ab':
        return t[1], t[2]
    if t[0] == 'fac':
        al, be = t[1], t[2]
        # (1 - al*p*T + p^3 T^2)(1 - be*T + p^3 T^2)
        a = -(al*P + be)
        bp = 2*P**3 + al*P*be
        assert bp % P == 0
        return a, bp // P
    return None

if __name__ == '__main__':
    for r in sorted(ROWS):
        v = ab(r)
        print(r, ROWS[r][0], v)
