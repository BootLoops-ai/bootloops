"""baller.transport — certified period/ODE transport (vendored kklt engines).

  pipe_transport  majorant-tail certified transport, parametric in operator
                  order/degree (vendor/kklt engine)
  pipe_lib        exact recurrence/tower layer (only the generic half is
                  vendored; see vendor header)
  march_lib       compact arb Taylor endpoint transport for
                  den(s)*s*v'(s) = Anum(s) v(s) — complex-step arb(acb)
                  handling, enforced input contract, outward floats. Input
                  contract: h/r as float, int, arb/acb or fmpq — NOT
                  Fraction/Decimal/str (raw TypeErrors from flint).

Engines load via the pinned vendor loader (source-hash per load, cache never
consulted); vendored sources ship untouched.
"""
from ._core import load_vendored


def configure_all(*args, **kw):
    """Configure pipe_lib AND re-snapshot pipe_transport in the required
    order (PT.configure() snapshots part of PL's state and reads
    the rest LIVE — reconfiguring PL without re-calling PT.configure() yields
    majorants matching NEITHER card). Always use this, never PL.configure
    alone, when PT is in play."""
    pl = load_vendored("pipe_lib")
    pt = load_vendored("pipe_transport")
    out = pl.configure(*args, **kw)
    pt.configure()
    return out


def __getattr__(name):
    if name in ("pipe_transport", "pipe_lib", "march_lib"):
        mod = load_vendored(name)
        globals()[name] = mod
        return mod
    raise AttributeError(f"module 'baller.transport' has no attribute {name!r}")


__all__ = ["pipe_transport", "pipe_lib", "march_lib", "configure_all"]
