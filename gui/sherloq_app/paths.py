from pathlib import Path


GUI_ROOT = Path(__file__).resolve().parents[1]
ICONS_DIR = GUI_ROOT / "icons"
MODELS_DIR = GUI_ROOT / "models"
PYEXIFTOOL_DIR = GUI_ROOT / "pyexiftool"
BUTTERAUGLI_DIR = GUI_ROOT / "butteraugli"
SSIMULACRA_DIR = GUI_ROOT / "ssimulacra"
TRUFOR_DIR = GUI_ROOT / "TruFor_main"


def icon_path(name):
    return str(ICONS_DIR / name)


def model_path(name):
    return str(MODELS_DIR / name)


def bundled_path(*parts):
    return str(GUI_ROOT.joinpath(*parts))
