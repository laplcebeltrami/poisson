"""Prepare lossless ASCII GIfTI copies for MATLAB without Java.

Run from any folder with Python 3 (standard library only):
    python3 assets/readme/prepare_demo_surfaces.py
Original repository files are read only. Geometry is not changed.
"""

from array import array
import base64
from pathlib import Path
import sys
import xml.etree.ElementTree as ET
import zlib

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "build" / "readme"
OUT.mkdir(parents=True, exist_ok=True)
TYPES = {
    "NIFTI_TYPE_FLOAT32": ("f", ".9g"),
    "NIFTI_TYPE_FLOAT64": ("d", ".17g"),
    "NIFTI_TYPE_INT32": ("i", "d"),
    "NIFTI_TYPE_UINT8": ("B", "d"),
}

for name in ("pial.surfL.gii", "fs_LR.32k.L.inflated.surf.gii", "fs_LR.32k.L.flat.surf.gii"):
    tree = ET.parse(ROOT / name)
    for node in tree.getroot().findall("DataArray"):
        data = node.find("Data")
        encoding = node.attrib["Encoding"]
        if encoding == "ASCII":
            continue
        raw = base64.b64decode(data.text)
        if encoding == "GZipBase64Binary":
            raw = zlib.decompress(raw)
        elif encoding != "Base64Binary":
            raise ValueError("Unsupported GIfTI encoding: " + encoding)
        code, precision = TYPES[node.attrib["DataType"]]
        values = array(code)
        values.frombytes(raw)
        little = node.attrib["Endian"] == "LittleEndian"
        if little != (sys.byteorder == "little"):
            values.byteswap()
        data.text = "\n" + " ".join(format(v, precision) for v in values) + "\n"
        node.set("Encoding", "ASCII")
    tree.write(OUT / name, encoding="utf-8", xml_declaration=True)
    print("Prepared", OUT / name)
