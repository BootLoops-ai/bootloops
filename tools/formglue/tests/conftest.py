"""Keep engine-capability skips LOUD: pytest's quiet mode (-q) reports skips
as a bare count, but a skipped formglue engine leg must always show its named
reason (which FORM binary/version was found, what capability is missing).
This summary hook prints those reasons unconditionally for formglue tests."""


def pytest_terminal_summary(terminalreporter, exitstatus, config):
    seen = set()
    lines = []
    for rep in terminalreporter.stats.get("skipped", []):
        nodeid = getattr(rep, "nodeid", "")
        if "formglue" not in nodeid:
            continue
        lr = getattr(rep, "longrepr", None)
        reason = lr[2] if isinstance(lr, tuple) and len(lr) == 3 else str(lr)
        if reason.startswith("Skipped: "):
            reason = reason[len("Skipped: "):]
        key = (nodeid, reason)
        if key in seen:
            continue
        seen.add(key)
        lines.append(f"SKIPPED {nodeid}: {reason}")
    if lines:
        terminalreporter.section("formglue named skips")
        for ln in lines:
            terminalreporter.write_line(ln)
