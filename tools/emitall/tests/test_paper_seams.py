#!/usr/bin/env python3
"""paper_seams.py fixture battery — the comment-swallow seam gate must strike
the planted S1 swallow and NAME the seam (file:line), pass the cured twin, and
honour --allow on a legitimate lowercase sentence head (the head rule).

Every leg runs the tool as a subprocess from a FOREIGN cwd (a fresh temp dir),
the way a paper build calls it; two legs run it as a module
(python3 -m emitall.paper_seams) to prove both invocations. Fixtures are
committed tex trees under specs/seams/ — a master that \\inputs one section
file, so the roster walk is exercised:

  planted_s1/  a comment block dropped INSIDE a sentence ("... Every entry" on
               the comment's last line, "carries ..." live below)  -> exit 1
  cured/       the same text with the comment on its own line boundary -> exit 0
  head_rule/   a live line that legitimately begins lowercase (a tool name)
               -> exit 1 bare, exit 0 with --allow "tallyfit, the backend"
  s1c_swallowed/ a rewritten live line whose second clause survives nowhere
               (snapshot/ holds the pre-edit tree)  -> exit 1 with --snapshot, S1c named
  s1c_quoted/  the same rewrite with the superseded line quoted in a comment -> exit 0

No PDF sits beside the fixtures, so S2 is skipped by design (and says so);
these legs cover the tex rows. Run:  python3 -m pytest tools/emitall -q
"""
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)                # tools/emitall
PKG_PARENT = os.path.dirname(PKG)          # tools/  (import root for -m emitall.…)
SEAMS_PY = os.path.join(PKG, "paper_seams.py")
FIX = os.path.join(PKG, "specs", "seams")


def _run(fixture, *flags, module=False, master="main.tex"):
    """Run the gate on specs/seams/<fixture>/<master> from a foreign cwd."""
    args = [os.path.join(FIX, fixture, master), *flags] if fixture else list(flags)
    with tempfile.TemporaryDirectory() as foreign:
        if module:
            env = dict(os.environ)
            env["PYTHONPATH"] = PKG_PARENT + os.pathsep + env.get("PYTHONPATH", "")
            cmd = [sys.executable, "-m", "emitall.paper_seams", *args]
        else:
            env = None
            cmd = [sys.executable, SEAMS_PY, *args]
        return subprocess.run(cmd, cwd=foreign, env=env, capture_output=True, text=True)


class TestPaperSeamsFixtures(unittest.TestCase):

    def test_planted_s1_strikes_and_names_the_seam(self):
        p = _run("planted_s1")
        self.assertEqual(p.returncode, 1, p.stdout + p.stderr)
        self.assertIn("FAIL S1 sections/body.tex:3 ", p.stdout)   # the seam is NAMED: file + line
        self.assertIn('"carries a receipt and a checksum."', p.stdout)
        self.assertIn("paper_seams: 1 swallow(s) -- main.tex", p.stdout)
        self.assertIn("S2 skipped (no PDF or no pdftotext)", p.stdout)
        self.assertIn("roster 2 files", p.stdout)               # master + the \input section
        self.assertNotIn("S1b ", p.stdout)                      # exactly the planted row fires

    def test_cured_twin_is_clean(self):
        p = _run("cured")
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertIn("paper_seams: CLEAN -- main.tex", p.stdout)
        self.assertNotIn("FAIL", p.stdout)

    def test_head_rule_strikes_without_allow(self):
        p = _run("head_rule")
        self.assertEqual(p.returncode, 1, p.stdout + p.stderr)
        self.assertIn("FAIL S1 sections/body.tex:3 ", p.stdout)
        self.assertIn("1 swallow(s)", p.stdout)

    def test_head_rule_allow_suppresses(self):
        p = _run("head_rule", "--allow", "tallyfit, the backend")
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertIn("paper_seams: CLEAN -- main.tex", p.stdout)

    def test_allow_is_phrase_specific(self):
        # negative control: an unrelated allowed head suppresses nothing
        p = _run("head_rule", "--allow", "some other head")
        self.assertEqual(p.returncode, 1, p.stdout + p.stderr)
        self.assertIn("FAIL S1 ", p.stdout)

    def test_quiet_keeps_verdict_and_exit_code(self):
        p = _run("planted_s1", "--quiet")
        self.assertEqual(p.returncode, 1, p.stdout + p.stderr)
        self.assertNotIn("FAIL", p.stdout)
        self.assertIn("paper_seams: 1 swallow(s) -- main.tex", p.stdout)

    def test_module_invocation_planted(self):
        p = _run("planted_s1", module=True)
        self.assertEqual(p.returncode, 1, p.stdout + p.stderr)
        self.assertIn("FAIL S1 sections/body.tex:3 ", p.stdout)

    def test_module_invocation_cured(self):
        p = _run("cured", module=True)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertIn("paper_seams: CLEAN -- main.tex", p.stdout)

    def test_s1c_vanished_clause_strikes_with_snapshot(self):
        snap = os.path.join(FIX, "s1c_swallowed", "snapshot")
        r = _run("s1c_swallowed", "--snapshot", snap)
        self.assertEqual(r.returncode, 1, r.stdout + r.stderr)
        self.assertIn("S1c sections/body.tex: snapshot line 2", r.stdout)
        self.assertIn("quadrature and truncation terms", r.stdout)

    def test_s1c_quoted_clause_is_clean(self):
        snap = os.path.join(FIX, "s1c_quoted", "snapshot")
        r = _run("s1c_quoted", "--snapshot", snap)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertNotIn("S1c", r.stdout)

    def test_s1c_is_silent_without_snapshot(self):
        r = _run("s1c_swallowed")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)

    def test_usage_is_rc_2(self):
        p = _run(None)
        self.assertEqual(p.returncode, 2, p.stdout + p.stderr)
        self.assertIn("usage:", p.stdout)


if __name__ == "__main__":
    unittest.main(verbosity=2)
