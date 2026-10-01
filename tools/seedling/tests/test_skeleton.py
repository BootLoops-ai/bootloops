"""Package sanity: imports + CLI parses."""
import sys

def test_import():
    import seedling
    assert seedling.__version__

def test_cli_requires_subcommand():
    from seedling.cli import main
    try:
        main([])
        raised = False
    except SystemExit as e:
        raised = e.code != 0
    assert raised

if __name__ == "__main__":
    test_import(); test_cli_requires_subcommand(); print("skeleton tests OK")
