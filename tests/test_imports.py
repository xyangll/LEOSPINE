import importlib


def test_core_imports():
    importlib.import_module("core.constellation")
    importlib.import_module("core.sim_data")
    importlib.import_module("core.sat_sim")
    importlib.import_module("core.positioning")
    importlib.import_module("core.utilities")


def test_app_imports():
    importlib.import_module("app.czml")
    importlib.import_module("app.gui_qt")
    importlib.import_module("app.settings")
    importlib.import_module("app.webserver")

