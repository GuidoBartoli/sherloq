# Project Structure

Sherloq keeps the Qt application code in `gui/sherloq_app` and leaves large
runtime assets, bundled command line tools and research models in stable
top-level folders under `gui`.

```text
sherloq.py                  Compatibility launcher from the repository root
gui/
  sherloq_app/              Sherloq Qt application package
    main.py                 QApplication entry point and main window
    paths.py                Centralized paths for icons, models and bundled tools
    core/                   Shared image-processing and utility helpers
    ui/                     Reusable Qt widgets, viewers and tool tree
    tools/                  Feature widgets grouped by toolbox category
      general/
      metadata/
      inspection/
      detail/
      colors/
      noise/
      jpeg/
      tampering/
      various/
  icons/                    Qt icons and static image assets
  models/                   Lightweight bundled ML models
  noiseprint/               Bundled Noiseprint implementation and weights
  pyexiftool/               Bundled PyExifTool and ExifTool binaries
  butteraugli/              Bundled Butteraugli helper
  ssimulacra/               Bundled SSIMULACRA helper
  TruFor_main/              Optional TruFor checkout and weights
```

Prefer absolute package imports rooted at `gui.sherloq_app` for Sherloq code.
Paths to files under `gui` should go through `gui.sherloq_app.paths` instead of
depending on the current working directory.
