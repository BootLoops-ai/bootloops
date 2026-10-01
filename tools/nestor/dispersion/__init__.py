r"""nestor.dispersion -- singularity-subtracted dispersion quadrature.

Exponentially convergent quadrature for a dispersion integral

        I = (1/pi) \int rho(w') K(w') dw'

whose spectral density rho carries KNOWN threshold non-analyticities
(power x log^m turn-ons, finite jumps at higher cuts) and whose kernel K is
analytic on each panel.  On a short sub-panel next to each threshold the
known singular model S is subtracted, so rho - S is analytic and tanh-sinh
converges exponentially; the subtracted piece \int S K is added back in
CLOSED FORM from exact power-log moments times the kernel's local Taylor
coefficients.  See disp_sub.py for the derivation and README.md (this
folder) for the measured convergence ladders.

    from nestor.dispersion import (disp_subtracted, moment_powerlog,
                                   kernel_taylor, addback_endpoint,
                                   tanhsinh_panel)

A worked kernel (the one-loop light-by-light box in 1D dilogarithmic form,
Vieta-stabilized) lives in ../examples/box1_dilog.py; the battery legs are
tests/test_dispersion_*.py and run inside `python3 -m nestor.selftest`.

mpmath only; no paths, no environment variables.
"""
from .disp_sub import (addback_endpoint, disp_subtracted, kernel_taylor,
                       moment_powerlog, tanhsinh_panel)

__all__ = ["disp_subtracted", "moment_powerlog", "kernel_taylor",
           "addback_endpoint", "tanhsinh_panel"]
