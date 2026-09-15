#!/usr/bin/env python3
"""Convert public spectra with explicit provenance and strict array checks."""

import argparse
import base64
import math
import re
import struct
import subprocess
import xml.etree.ElementTree as ET
import zlib
from pathlib import Path

from provenance import atomic_json, cached_stage


def mzxml_to_mgf(source, output):
    count = 0
    with output.open("w") as out:
        for _, element in ET.iterparse(source, events=("end",)):
            if element.tag.rsplit("}", 1)[-1] != "scan":
                continue
            if element.get("msLevel") == "2":
                children = {c.tag.rsplit("}", 1)[-1]: c for c in element}
                precursor, peaks = children.get("precursorMz"), children.get("peaks")
                if precursor is None or peaks is None or not precursor.text:
                    raise ValueError("MS2 scan lacks precursor or peak array")
                precision = peaks.get("precision")
                if precision not in ("32", "64"):
                    raise ValueError("Unsupported mzXML precision")
                order = peaks.get("byteOrder", "network")
                if order not in ("network", "big", "little"):
                    raise ValueError("Unknown mzXML byte order")
                if peaks.get("pairOrder", "m/z-int") != "m/z-int":
                    raise ValueError("Unsupported mzXML pair order")
                raw = base64.b64decode("".join((peaks.text or "").split()), validate=True)
                compression = peaks.get("compressionType", "none")
                if compression == "zlib":
                    raw = zlib.decompress(raw)
                elif compression != "none":
                    raise ValueError("Unsupported mzXML compression")
                code = ("<" if order == "little" else ">") + ("ff" if precision == "32" else "dd")
                if len(raw) != int(element.get("peaksCount")) * struct.calcsize(code):
                    raise ValueError("mzXML peak count disagrees with binary array")
                rt = element.get("retentionTime", "")
                match = re.fullmatch(r"PT([0-9.]+)(S|M)", rt)
                if not match:
                    raise ValueError(f"Unsupported retention time: {rt}")
                seconds = float(match[1]) * (60 if match[2] == "M" else 1)
                mz = float(precursor.text)
                if not math.isfinite(mz) or mz <= 0:
                    raise ValueError("Invalid precursor m/z")
                out.write(f"BEGIN IONS\nTITLE=scan={element.get('num')}\nSCANS={element.get('num')}\nPEPMASS={mz:.12g}\nRTINSECONDS={seconds:.12g}\n")
                charge = int(precursor.get("precursorCharge", "0"))
                if charge > 0:
                    out.write(f"CHARGE={charge}+\n")
                for mz, intensity in struct.iter_unpack(code, raw):
                    if not math.isfinite(mz) or not math.isfinite(intensity) or mz <= 0 or intensity < 0:
                        raise ValueError("Invalid peak")
                    out.write(f"{mz:.12g} {intensity:.12g}\n")
                out.write("END IONS\n\n")
                count += 1
            element.clear()
    if not count:
        raise ValueError("No MS2 spectra converted")
    return count


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("source", type=Path)
    p.add_argument("--root", type=Path, required=True)
    args = p.parse_args()
    source = args.source.resolve(strict=True)
    destination = args.root / "converted" / source.parent.name
    destination.mkdir(parents=True, exist_ok=True)
    is_xml = source.suffix.lower() == ".mzxml"
    output = destination / (source.stem + (".mgf" if is_xml else ".mzML"))
    inputs = [source, Path(__file__).resolve()]
    command = None
    if not is_xml:
        if source.suffix.lower() != ".raw":
            raise ValueError("Only mzXML and RAW input are accepted")
        converter = args.root / "tools/thermo/ThermoRawFileParser"
        inputs.extend(p for p in converter.parent.iterdir() if p.is_file())
        command = [str(converter), f"-i={source}", f"-b={output}", "-f=1", "-m=0"]
    with cached_stage(output.with_suffix(".conversion.json"), "scientific-conversion-v1",
                      inputs, [output], {"command": command, "mzxml_ms2_only": is_xml}) as hit:
        if not hit:
            if is_xml:
                count = mzxml_to_mgf(source, output)
                atomic_json(output.with_suffix(".spectra.json"), {"ms2_spectra": count, "ms1_retained": False})
            else:
                with output.with_suffix(".conversion.log").open("w") as log:
                    subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=1800)
    print(output, flush=True)


if __name__ == "__main__":
    main()
