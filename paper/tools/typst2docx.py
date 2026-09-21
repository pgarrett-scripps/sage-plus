#!/usr/bin/env python3
"""Convert the Typst-exported HTML of paper.typ into a Word .docx.

Invoked by `just docx`; not usually run by hand. Two problems sit between Typst's
HTML export and a usable Word file, and this script fixes both.

1. Inline math. Typst's HTML export silently DROPS every equation, gutting
   sentences into things like "the threshold's __ must be recalibrated ...
   ranges __ to __". paper.typ works around this in docx mode by
   wrapping each equation in `html.frame()`, which preserves it as an SVG. But
   Typst emits that frame as a BLOCK-level sibling, so a sentence containing inline
   math is shattered into separate <p> fragments:

       <p>The threshold</p> <svg/> <p>must be recalibrated ...</p>

   We stitch those fragments back into one paragraph, then rasterize each math SVG
   to an inline PNG. PNG rather than SVG because pandoc cannot size a raw SVG
   without `rsvg-convert` (absent here) and older Word versions will not render it.

2. Everything else already survives the round trip: headings become real Word
   heading styles, tables stay tables, and figures embed as images.

Usage: typst2docx.py <input.html> <output.docx>
"""

from __future__ import annotations

import base64
import io
import re
import sys
from pathlib import Path
from xml.dom import minidom
from zipfile import ZipFile

import cairosvg
import pypandoc
from PIL import Image

# Rasterization factor for math. 4x keeps symbols crisp on a high-DPI screen and
# when printed; the <img> is then sized back down to the glyph's natural size.
SCALE = 4

# What pandoc's default reference document sets as the Word body size
# (docDefaults, measured from the produced file). Typst frames carry their own
# em -- an 11pt manuscript emits 11pt-per-em frames -- so without rescaling,
# every equation lands ~9% smaller than the Word text around it. Each frame's
# em is read from its own style attribute, so a project that changes its text
# size stays correctly scaled.
WORD_BODY_PT = 12.0


def _em_pt(svg: str) -> float:
    """pt per em for this frame, from its width in pt vs its width in em."""
    pt = _pt(svg, "width")
    em = re.search(r'width:\s*([0-9.]+)em', svg)
    if pt and em and float(em.group(1)):
        return pt / float(em.group(1))
    return WORD_BODY_PT                      # no ratio available: no rescale

# Only stitch a frame back into the surrounding paragraph if it is glyph-sized.
# Inline math runs ~7-15pt tall; a display equation is taller and deserves to stay
# its own block. Raise this if a tall inline construct (a fraction, a stacked
# subscript) is being wrongly promoted to its own paragraph.
MAX_INLINE_PT = 30.0

FRAME_RE = re.compile(r'<svg class="typst-frame".*?</svg>', re.S)
# </p>  <frame>  <p>   ->   one paragraph, frame inline
MERGE_RE = re.compile(r'</p>\s*(<svg class="typst-frame".*?</svg>)\s*<p>', re.S)


def _pt(svg: str, attr: str) -> float:
    m = re.search(rf'\b{attr}="([0-9.]+)pt"', svg)
    return float(m.group(1)) if m else 0.0


# Typst sizes a math frame to its LINE BOX and then relies on `overflow: visible`
# to paint the parts that stick out (subscripts, descenders, italic overhang).
# Browsers honour that; cairosvg does not -- it rasterizes to a canvas of exactly
# width x height and silently shears off anything outside, which cropped every
# subscripted equation to a uniform 44px band. So grow the canvas before
# rendering, then crop back to the ink. PAD is a multiple of the frame height and
# only has to be big enough to contain the overhang; the crop removes the slack.
PAD = 1.0


def _pad_frame(svg: str) -> str:
    """Enlarge an SVG's canvas symmetrically so nothing is clipped when raster-
    ized. Keeps the user-unit scale identical, so the glyph renders unchanged."""
    w, h = _pt(svg, "width"), _pt(svg, "height")
    vb = re.search(r'viewBox="([-\d.]+) ([-\d.]+) ([-\d.]+) ([-\d.]+)"', svg)
    if not (w and h and vb):
        return svg
    x0, y0, vw, vh = (float(g) for g in vb.groups())
    pad = PAD * h
    head_end = svg.index(">") + 1
    head, body = svg[:head_end], svg[head_end:]
    head = re.sub(r'\bwidth="[0-9.]+pt"', f'width="{w + 2 * pad}pt"', head)
    head = re.sub(r'\bheight="[0-9.]+pt"', f'height="{h + 2 * pad}pt"', head)
    head = re.sub(
        r'viewBox="[^"]+"',
        f'viewBox="{x0 - pad} {y0 - pad} {vw + 2 * pad} {vh + 2 * pad}"',
        head,
    )
    return head + body


def _render(svg: str) -> tuple[bytes, float]:
    """Rasterize one math frame and return (png bytes, ink width in pt)."""
    png = cairosvg.svg2png(bytestring=_pad_frame(svg).encode(), scale=SCALE)
    im = Image.open(io.BytesIO(png))
    box = im.getbbox()  # bounding box of the non-transparent pixels
    if box:
        im = im.crop(box)
    buf = io.BytesIO()
    im.save(buf, format="PNG")
    return buf.getvalue(), im.width / SCALE


def merge_inline_frames(html: str) -> tuple[str, int]:
    """Rejoin paragraphs that Typst split around an inline math frame."""

    def repl(m: re.Match) -> str:
        svg = m.group(1)
        if _pt(svg, "height") > MAX_INLINE_PT:
            return m.group(0)  # tall enough to be display math; leave it alone
        return f" {svg} "

    return MERGE_RE.subn(repl, html)


def rasterize_frames(html: str) -> tuple[str, int]:
    """Replace each math SVG with an inline PNG <img> at its natural size."""

    def repl(m: re.Match) -> str:
        svg = m.group(0)
        png, ink_pt = _render(svg)
        b64 = base64.b64encode(png).decode()
        # pt -> CSS px so Word lays the glyph out at the size Typst intended,
        # rescaled from the manuscript's em to Word's, so an 11pt paper's math
        # is not 9% small beside 12pt Word text. Measured off the cropped image
        # rather than the declared width, since the declared box understates
        # glyphs that overhang it.
        w = ink_pt * (WORD_BODY_PT / _em_pt(svg)) * 96 / 72
        size = f' width="{w:.1f}"' if w else ""
        return f'<img src="data:image/png;base64,{b64}"{size} alt="equation" />'

    return FRAME_RE.subn(repl, html)


def keep_table_captions(path: Path) -> None:
    """Keep compact tables and figures with complete captions in Word."""
    with ZipFile(path) as source:
        entries = [(item, source.read(item.filename)) for item in source.infolist()]
    document = next(data for item, data in entries if item.filename == "word/document.xml")
    dom = minidom.parseString(document)
    namespace = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"

    # Pandoc wraps a table figure in a one-cell table. Remove that wrapper so
    # the data table uses the page width and its caption remains a sibling.
    for table in list(dom.getElementsByTagName("w:tbl")):
        styles = table.getElementsByTagName("w:tblStyle")
        nested = table.getElementsByTagName("w:tbl")
        if styles and styles[0].getAttribute("w:val") == "FigureTable" and len(nested) == 1:
            inner = nested[0]
            inner.parentNode.removeChild(inner)
            table.parentNode.replaceChild(inner, table)

    def property_on(paragraph, name):
        properties = next((node for node in paragraph.childNodes if getattr(node, "tagName", "") == "w:pPr"), None)
        if properties is None:
            properties = dom.createElementNS(namespace, "w:pPr")
            paragraph.insertBefore(properties, paragraph.firstChild)
        if not properties.getElementsByTagName(name):
            properties.appendChild(dom.createElementNS(namespace, name))

    bibliography = False
    for paragraph in dom.getElementsByTagName("w:p"):
        text = "".join(node.firstChild.data for node in paragraph.getElementsByTagName("w:t") if node.firstChild)
        if text == "References":
            bibliography = True
        if text == "Supporting Information":
            bibliography = False
            property_on(paragraph, "w:pageBreakBefore")
        if bibliography:
            for numbering in list(paragraph.getElementsByTagName("w:numPr")):
                numbering.parentNode.removeChild(numbering)
        if re.match(r"^(Figure|Table)\s", text):
            property_on(paragraph, "w:keepLines")
        if paragraph.getElementsByTagName("w:drawing"):
            property_on(paragraph, "w:keepNext")
    for table in dom.getElementsByTagName("w:tbl"):
        cells = table.getElementsByTagName("w:tc")
        headers = ["".join(node.firstChild.data for node in cell.getElementsByTagName("w:t") if node.firstChild) for cell in cells[:6]]
        development_table = headers[:3] == ["Area", "Sage Plus change", "Evaluation in this chapter"]
        scaling_table = headers == ["Engine", "Workers", "Seconds", "Time range", "Peak MiB", "PSM range"]
        mass_offset_table = headers == ["Configuration", "Search", "Wall", "Peak RSS", "DB peptides", "Fragments"]
        if mass_offset_table:
            widths = (1250, 850, 850, 1200, 1400, 1550, 1100, 1160)
            for column, width in zip(table.getElementsByTagName("w:gridCol"), widths):
                column.setAttributeNS(namespace, "w:w", str(width))
        if development_table:
            # The narrative inventory needs wider change columns and left alignment.
            widths = (1800, 3300, 3900)
            properties = table.getElementsByTagName("w:tblPr")[0]
            for name, attributes in (
                ("w:tblInd", {"w:w": "0", "w:type": "dxa"}),
                ("w:jc", {"w:val": "left"}),
                ("w:tblW", {"w:w": "9000", "w:type": "dxa"}),
            ):
                elements = properties.getElementsByTagName(name)
                element = elements[0] if elements else dom.createElementNS(namespace, name)
                for attribute, value in attributes.items():
                    element.setAttributeNS(namespace, attribute, value)
                if not elements:
                    properties.appendChild(element)
            grid = table.getElementsByTagName("w:gridCol")
            for column, width in zip(grid, widths):
                column.setAttributeNS(namespace, "w:w", str(width))
            for index, cell in enumerate(cells):
                for width in cell.getElementsByTagName("w:tcW"):
                    width.setAttributeNS(namespace, "w:w", str(widths[index % len(widths)]))
                    width.setAttributeNS(namespace, "w:type", "dxa")
                for paragraph in cell.getElementsByTagName("w:p"):
                    property_on(paragraph, "w:jc")
                    paragraph.getElementsByTagName("w:jc")[0].setAttributeNS(namespace, "w:val", "left")
                    property_on(paragraph, "w:spacing")
                    spacing = paragraph.getElementsByTagName("w:spacing")[0]
                    spacing.setAttributeNS(namespace, "w:before", "0")
                    spacing.setAttributeNS(namespace, "w:after", "40")
        # Compact table text prevents header and numeric-range fragmentation.
        for run in table.getElementsByTagName("w:r"):
            properties = next((node for node in run.childNodes if getattr(node, "tagName", "") == "w:rPr"), None)
            if properties is None:
                properties = dom.createElementNS(namespace, "w:rPr")
                run.insertBefore(properties, run.firstChild)
            for name in ("w:sz", "w:szCs"):
                elements = properties.getElementsByTagName(name)
                element = elements[0] if elements else dom.createElementNS(namespace, name)
                element.setAttributeNS(namespace, "w:val", "18" if development_table or scaling_table or mass_offset_table else "20")
                if not elements:
                    properties.appendChild(element)
        following = table.nextSibling
        while following is not None and following.nodeType != following.ELEMENT_NODE:
            following = following.nextSibling
        if following is None or following.tagName != "w:p":
            continue
        caption = "".join(node.firstChild.data for node in following.getElementsByTagName("w:t") if node.firstChild)
        if not re.match(r"^Table\s", caption):
            continue
        # Word and LibreOffice reliably keep a preceding caption with a table.
        # A keep-next flag inside the final cell does not cross that boundary.
        table.parentNode.insertBefore(following, table)
        property_on(following, "w:keepNext")
        rows = [node for node in table.childNodes if getattr(node, "tagName", "") == "w:tr"]
        if not rows:
            continue
        for row in rows:
            properties = next((node for node in row.childNodes if getattr(node, "tagName", "") == "w:trPr"), None)
            if properties is None:
                properties = dom.createElementNS(namespace, "w:trPr")
                row.insertBefore(properties, row.firstChild)
            if not properties.getElementsByTagName("w:cantSplit"):
                properties.appendChild(dom.createElementNS(namespace, "w:cantSplit"))
        # The manuscript's short tables fit within one page. Longer future
        # tables may paginate without splitting individual rows.
        held_rows = rows[:-1] if len(rows) <= 21 else []
        for row in held_rows:
            for paragraph in row.getElementsByTagName("w:p"):
                property_on(paragraph, "w:keepNext")
    output = dom.toxml(encoding="UTF-8")
    with ZipFile(path, "w") as target:
        for item, data in entries:
            target.writestr(item, output if item.filename == "word/document.xml" else data)


def manuscript_style(path: Path) -> None:
    """Keep Word headings black and make short comparison tables legible."""
    with ZipFile(path) as archive:
        entries = [(item, archive.read(item.filename)) for item in archive.infolist()]
    data = dict((item.filename, content) for item, content in entries)
    namespace = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"

    def child(dom, parent, name):
        existing = next((node for node in parent.childNodes
                         if getattr(node, "tagName", "") == name), None)
        if existing is None:
            existing = dom.createElementNS(namespace, name)
            if name in ("w:tcPr", "w:trPr", "w:tblPr", "w:pPr", "w:rPr") and parent.tagName != "w:style":
                parent.insertBefore(existing, parent.firstChild)
            else:
                parent.appendChild(existing)
        return existing

    styles = minidom.parseString(data["word/styles.xml"])
    for style in styles.getElementsByTagName("w:style"):
        name = style.getAttribute("w:styleId")
        if name in ("Title", "Subtitle") or name.startswith("Heading"):
            properties = child(styles, style, "w:rPr")
            color = child(styles, properties, "w:color")
            for attribute in ("w:themeColor", "w:themeTint", "w:themeShade"):
                if color.hasAttribute(attribute):
                    color.removeAttribute(attribute)
            color.setAttribute("w:val", "000000")
            child(styles, properties, "w:u").setAttribute("w:val", "none")
            for border in list(style.getElementsByTagName("w:pBdr")):
                border.parentNode.removeChild(border)
    document = minidom.parseString(data["word/document.xml"])
    first = document.getElementsByTagName("w:body")[0].getElementsByTagName("w:p")[0]
    props = child(document, first, "w:pPr")
    child(document, props, "w:pStyle").setAttribute("w:val", "Title")
    section = document.getElementsByTagName("w:sectPr")[-1]
    page = child(document, section, "w:pgSz")
    page.setAttribute("w:w", "12240")
    page.setAttribute("w:h", "15840")
    margins = child(document, section, "w:pgMar")
    for name, value in (("top", 1440), ("bottom", 1440), ("left", 1440),
                        ("right", 1440), ("header", 720), ("footer", 720), ("gutter", 0)):
        margins.setAttribute(f"w:{name}", str(value))
    for table in document.getElementsByTagName("w:tbl"):
        props = child(document, table, "w:tblPr")
        child(document, props, "w:tblLayout").setAttribute("w:type", "fixed")
        width = child(document, props, "w:tblW")
        width.setAttribute("w:type", "dxa")
        width.setAttribute("w:w", "9360")
        columns = table.getElementsByTagName("w:gridCol")
        original = [int(column.getAttribute("w:w") or "1") for column in columns]
        widths = [round(9360 * value / sum(original)) for value in original]
        widths[-1] += 9360 - sum(widths)
        for column, value in zip(columns, widths):
            column.setAttribute("w:w", str(value))
        for row in table.getElementsByTagName("w:tr"):
            for cell, value in zip(row.getElementsByTagName("w:tc"), widths):
                cell_width = child(document, child(document, cell, "w:tcPr"), "w:tcW")
                cell_width.setAttribute("w:type", "dxa")
                cell_width.setAttribute("w:w", str(value))
        borders = child(document, props, "w:tblBorders")
        for edge in ("top", "bottom", "left", "right", "insideH", "insideV"):
            border = child(document, borders, f"w:{edge}")
            for key, value in (("val", "single"), ("sz", "4"), ("color", "D9D9D9")):
                border.setAttribute(f"w:{key}", value)
        rows = table.getElementsByTagName("w:tr")
        if rows:
            child(document, child(document, rows[0], "w:trPr"), "w:tblHeader")
            for cell in rows[0].getElementsByTagName("w:tc"):
                shade = child(document, child(document, cell, "w:tcPr"), "w:shd")
                shade.setAttribute("w:fill", "F2F2F2")
    data["word/styles.xml"] = styles.toxml(encoding="UTF-8")
    data["word/document.xml"] = document.toxml(encoding="UTF-8")
    with ZipFile(path, "w") as archive:
        for item, content in entries:
            archive.writestr(item, data[item.filename])


def main() -> int:
    if len(sys.argv) != 3:
        print(__doc__.strip().splitlines()[-1], file=sys.stderr)
        return 2
    src, out = Path(sys.argv[1]), Path(sys.argv[2])

    html = src.read_text(encoding="utf-8")
    html, merged = merge_inline_frames(html)
    html, drawn = rasterize_frames(html)
    print(f"math: {drawn} equations rasterized, {merged} paragraphs rejoined")

    staged = src.with_suffix(".staged.html")
    staged.write_text(html, encoding="utf-8")
    try:
        pypandoc.convert_file(
            str(staged), "docx", outputfile=str(out),
            extra_args=["--standalone",
                        "--shift-heading-level-by=-1"],
        )
    finally:
        staged.unlink(missing_ok=True)

    keep_table_captions(out)
    manuscript_style(out)
    print(f"wrote {out.name} ({out.stat().st_size / 1e6:.1f} MB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
