#  ETE3 circular tree 
import os, random, re, shutil, subprocess
from pathlib import Path
from typing import Optional
import matplotlib
matplotlib.use("Agg")  
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")  
out_pdf = Path.cwd() / "png_tree600.pdf"
def _which(name: str) -> Optional[str]:
    try:
        return shutil.which(name)
    except Exception:
        return None

def vector_to_png(src, dpi: int = 1200, out_path: Optional[str] = None) -> str:
    src_path = Path(src)
    if out_path is None:
        out_path = Path(str(src_path.with_suffix('')) + f"_{dpi}dpi.png")
    else:
        out_path = Path(out_path)
    if out_path.exists():
        return str(out_path)

    gs, inkscape = _which("gs"), _which("inkscape")

    if gs:
        cmd = [
            gs, "-dSAFER", "-dBATCH", "-dNOPAUSE",
            "-sDEVICE=pngalpha",  
            f"-r{dpi}",
            f"-sOutputFile={str(out_path)}",
            str(src_path),
        ]
        subprocess.run(cmd, check=True)
        return str(out_path)

    if inkscape:
        cmd = [
            inkscape, str(src_path),
            "--export-type=png",
            f"--export-filename={str(out_path)}",
            f"--export-dpi={dpi}",
        ]
        subprocess.run(cmd, check=True)
        return str(out_path)

    raise RuntimeError(
        "Need Ghostscript (`gs`) or Inkscape installed to convert EPS/PDF to PNG."
    )

def safe_convert(path_like, dpi: int = 1200):
    """Only convert if the file exists; never crash."""
    if path_like is None:
        return None
    p = Path(path_like)
    if not p.exists():
        print(f"[skip] {p} not found; not converting.")
        return None
    try:
        out = vector_to_png(p, dpi=dpi)
        print("[ok] converted ->", out)
        return out
    except Exception as e:
        print(f"[warn] conversion failed for {p}: {e}")
        return None

converted_pngs = [safe_convert(out_pdf, dpi=1200)]
print("Done. Converted PNGs:", [c for c in converted_pngs if c])

from PIL import Image, ImageChops

def trim_whitespace(path_in, path_out=None, bg='white'):
    im = Image.open(path_in).convert('RGB')
    bg_img = Image.new('RGB', im.size, bg)
    diff = ImageChops.difference(im, bg_img)
    bbox = diff.getbbox()
    if not bbox:  
        bbox = (0, 0, im.width, im.height)
    im.crop(bbox).save(path_out or path_in)

trim_whitespace("legend.png", "legend_cropped.png")
trim_whitespace("legend_solid_lines.png", "legend_solid_lines_cropped.png")