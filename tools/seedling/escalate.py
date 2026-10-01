"""Escalation state machine.

Signal -> action policy (all measured-basis):

  LEFTOVER_TARGETS / REGENERATE_LOOP  (staging-visible incompleteness)
      -> r+1 first, then d+1, then s+1 LAST. Measured: r+1/d+1 moved eqn count
         <=0.03% on the reference family; s+1 re-blew the cell (x5.78 bulk).
  GATE_FAIL_CHAIN  (two-slice FAIL_PIN whose failing rows carry ISP powers or
      super-sector chain structure — the chain-loss signature)
      -> s+1 FIRST (chain loss is an s-phenomenon), then r+1, then d+1.
  GATE_FAIL_OTHER  (two-slice FAIL_PIN, no chain signature)
      -> conservative full ladder r+1 -> d+1 -> s+1.
  RANK_DEFICIT / RANK_EXCESS  (rank check vs independent expectation)
      -> NOT a margin problem: REFER (symmetry completion / basis quotient
         needed; the missed-symmetry and retained-redundant-master exhibits
         below). Escalating margins cannot fix these.
  DE_CENSUS_UNCOVERED
      -> NOT an escalation: run the minikira certification (runner.minikira_spec).

Termination: margins are monotone, bounded by max_margins; GIVE_UP means
"fall back to unpinned Laporta on bigger iron" (sound worst case).

Every transition is a ledger record — the label stream for the learned-margin
extension (each record: margins before/after, signal, staging outcome).
"""
from dataclasses import dataclass, field

# signals
LEFTOVER_TARGETS = "LEFTOVER_TARGETS"
REGENERATE_LOOP = "REGENERATE_LOOP"
GATE_FAIL_CHAIN = "GATE_FAIL_CHAIN"
GATE_FAIL_OTHER = "GATE_FAIL_OTHER"
RANK_DEFICIT = "RANK_DEFICIT"
RANK_EXCESS = "RANK_EXCESS"
DE_CENSUS_UNCOVERED = "DE_CENSUS_UNCOVERED"

# terminal / non-margin actions
GIVE_UP = "GIVE_UP"                  # unpinned fallback (bigger iron)
REFER = "REFER"                      # rank problem: symmetry/basis cure, not margins
RUN_MINIKIRA = "RUN_MINIKIRA"        # de-census certification, no margin change

_LADDERS = {
    LEFTOVER_TARGETS: ("r", "d", "s"),
    REGENERATE_LOOP: ("r", "d", "s"),
    GATE_FAIL_CHAIN: ("s", "r", "d"),
    GATE_FAIL_OTHER: ("r", "d", "s"),
}
_IDX = {"r": 0, "s": 1, "d": 2}


@dataclass
class EscalationState:
    margins: tuple = (0, 0, 0)                # (dr, ds, dd)
    max_margins: tuple = (3, 2, 3)
    history: list = field(default_factory=list)


def next_margins(state, signal):
    """Return (action, new_margins). action in {'ESCALATE', REFER, RUN_MINIKIRA,
    GIVE_UP}. On 'ESCALATE', new_margins is the incremented tuple (recorded in
    state.history); other actions leave margins unchanged."""
    if signal in (RANK_DEFICIT, RANK_EXCESS):
        state.history.append({"signal": signal, "action": REFER,
                              "margins": state.margins})
        return REFER, state.margins
    if signal == DE_CENSUS_UNCOVERED:
        state.history.append({"signal": signal, "action": RUN_MINIKIRA,
                              "margins": state.margins})
        return RUN_MINIKIRA, state.margins
    ladder = _LADDERS.get(signal)
    if ladder is None:
        raise ValueError(f"unknown escalation signal: {signal}")
    for dim in ladder:
        i = _IDX[dim]
        if state.margins[i] < state.max_margins[i]:
            new = list(state.margins)
            new[i] += 1
            new = tuple(new)
            state.history.append({"signal": signal, "action": "ESCALATE",
                                  "dim": dim, "margins": new})
            state.margins = new
            return "ESCALATE", new
    state.history.append({"signal": signal, "action": GIVE_UP,
                          "margins": state.margins})
    return GIVE_UP, state.margins


def classify_gate_failure(fail_rows, targets_by_key=None):
    """Map two-slice FAIL_PIN failing rows to a chain/other signature.

    fail_rows: iterable of index tuples (or keys resolvable via targets_by_key).
    Chain signature (chain-loss class): any failing row carries an ISP power (a
    negative index) — the dropped equations live in the s-direction.
    """
    for row in fail_rows:
        idx = targets_by_key[row] if targets_by_key else row
        if any(n < 0 for n in idx):
            return GATE_FAIL_CHAIN
    return GATE_FAIL_OTHER
