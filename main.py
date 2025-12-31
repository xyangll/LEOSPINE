from __future__ import annotations

import argparse
import importlib
import sys
import traceback


def is_frozen():
    """Check if running as a frozen (PyInstaller) executable."""
    return getattr(sys, 'frozen', False) and hasattr(sys, '_MEIPASS')


def _load_main(module_name: str):
    mod = importlib.import_module(module_name)
    try:
        return getattr(mod, "main")
    except AttributeError as e:
        raise RuntimeError(f"Module {module_name!r} does not provide a main() entrypoint") from e


def run(argv: list[str] | None = None) -> int:
    """
    Launch the GUI.

    - Default is auto: try Qt (PySide6) first, fall back to Tk on failure.
    - Use --gui qt|tk|auto to force a backend.
    - In frozen (PyInstaller) mode, only Qt GUI is available.
    """
    parser = argparse.ArgumentParser(prog="LEO-SPINE", add_help=True)
    parser.add_argument(
        "--gui",
        choices=("auto", "qt", "tk"),
        default="auto",
        help="Choose GUI backend: auto (default) tries Qt then falls back to Tk; qt forces Qt; tk forces Tk",
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Print full traceback on startup failure (useful when auto falls back)",
    )
    args = parser.parse_args(argv)

    # In frozen mode, only use Qt GUI (no Tk fallback)
    if is_frozen():
        try:
            _load_main("app.gui_qt")()
            return 0
        except Exception as e:
            print(f"[LEO-SPINE] Failed to start GUI: {e}", file=sys.stderr)
            traceback.print_exc()
            input("Press Enter to exit...")  # Keep console window open
            return 1

    # Normal mode: try Qt, fallback to Tk
    if args.gui in ("auto", "qt"):
        try:
            _load_main("app.gui_qt")()
            return 0
        except Exception as e:
            if args.gui == "qt":
                raise
            # auto mode: fall back to Tk, but print the reason to stderr for debugging
            print("[LEO-SPINE] Failed to start Qt GUI; falling back to Tk GUI. Reason:", file=sys.stderr)
            if args.debug:
                traceback.print_exc()
            else:
                print(f"{type(e).__name__}: {e}", file=sys.stderr)

    _load_main("app.gui")()
    return 0


def main():
    """Main entry point for LEO-SPINE application."""
    raise SystemExit(run())


if __name__ == "__main__":
    main()