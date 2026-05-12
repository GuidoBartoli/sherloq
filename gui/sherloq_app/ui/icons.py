import hashlib
import re
from pathlib import Path
from tempfile import gettempdir

from PySide6.QtWidgets import QApplication
from PySide6.QtGui import QColor, QIcon, QImage, QPalette, QPixmap

from gui.sherloq_app.paths import icon_path


_FILL_RE = re.compile(r'fill="(?!none)([^"]+)"')
_PATH_WITHOUT_FILL_RE = re.compile(r"<path\b(?![^>]*\bfill=)")


def is_dark_theme():
    app = QApplication.instance()
    if app is None:
        return False
    window_color = app.palette().color(QPalette.Window)
    return window_color.lightness() < 128


def _themed_svg_path(source):
    source_path = Path(source)
    svg = source_path.read_text(encoding="utf-8")
    fill_values = set(_FILL_RE.findall(svg))
    has_brand_colors = any(fill != "currentColor" for fill in fill_values)
    if has_brand_colors:
        return str(source_path)

    color = "#ffffff" if is_dark_theme() else "#000000"
    themed_svg = svg.replace("currentColor", color)
    themed_svg = _PATH_WITHOUT_FILL_RE.sub(f'<path fill="{color}"', themed_svg)
    if themed_svg == svg:
        return str(source_path)

    digest = hashlib.sha1(
        f"{source_path}:{color}:{themed_svg}".encode("utf-8")
    ).hexdigest()[:12]
    cache_dir = Path(gettempdir()) / "sherloq-icons"
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / f"{source_path.stem}-{digest}.svg"
    if not cache_path.exists():
        cache_path.write_text(themed_svg, encoding="utf-8")
    return str(cache_path)


def _mask_icon(source):
    image = QImage(source).convertToFormat(QImage.Format_ARGB32)
    if image.isNull():
        return QIcon(source)

    color = QColor("#ffffff" if is_dark_theme() else "#000000")
    for y in range(image.height()):
        for x in range(image.width()):
            pixel = image.pixelColor(x, y)
            if pixel.alpha() == 0:
                continue
            pixel.setRed(color.red())
            pixel.setGreen(color.green())
            pixel.setBlue(color.blue())
            image.setPixelColor(x, y, pixel)
    return QIcon(QPixmap.fromImage(image))


def themed_icon(name):
    source = icon_path(name)
    if name.lower().endswith(".svg"):
        source = _themed_svg_path(source)
    if name == "sherloq_alpha.png":
        return _mask_icon(source)
    return QIcon(source)
