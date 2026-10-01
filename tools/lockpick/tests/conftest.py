# The ambient-dps regression imports the flat shims (mplll, pslq_gate) exactly
# as receipts cite them; resolve them from the tools/ root two levels up.
import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
