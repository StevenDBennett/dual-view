"""
CLI entry point for dual-view.

Usage:
    python -m dual_view           # run full demo
    python -m dual_view --quick   # skip slow demos
    python -m dual_view --test    # run test suite
    python -m dual_view --version
"""
import sys
import argparse


def main() -> None:
    parser = argparse.ArgumentParser(description="dual-view: a mathematical framework for 2-adic number systems, with diagnostics for quantized neural network weights")
    parser.add_argument("--quick", action="store_true", help="Skip slow demos")
    parser.add_argument("--test", action="store_true", help="Run test suite")
    parser.add_argument("--version", action="store_true", help="Print version")
    args = parser.parse_args()

    if args.version:
        from . import __version__
        print(f"dual-view {__version__}")
        return

    if args.test:
        import unittest
        import os
        test_dir = os.path.join(os.path.dirname(__file__), "../../tests")
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../.."))
        loader = unittest.TestLoader()
        suite = loader.discover(test_dir, pattern="test_*.py")
        runner = unittest.TextTestRunner(verbosity=2)
        result = runner.run(suite)
        sys.exit(0 if result.wasSuccessful() else 1)
        return

    from .demo import main as demo_main
    demo_main()


if __name__ == "__main__":
    main()
