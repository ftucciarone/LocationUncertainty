#!/usr/bin/env python3
"""
Draw noise-pipeline-railroad.svg. Run from this folder: python make_railroad.py

Railroad libraries make every track rejoin. This one has tracks that stop, so it is drawn by hand.
"""
R, ROW = 10, 44
EDGE, SPLIT, JOIN, TURN, W = 20, 60, 1000, 1040, 1060
BLUE, ORANGE, GREY = "#e8f1f8", "#fdf1dc", "#f2f2f2"
tracks, bars, boxes = [], [], []


def line(*points):
    r"""
    line: add a track through the given points, with rounded corners.

        Parameters:
            *points (tuple): (x, y) points, in order

        Returns:
            None; the path is added to tracks

        Formula:
            \text{SVG path through the points, each corner replaced by a quadratic curve of radius } R

        See also:
            branch
    """
    d = f"M{points[0][0]} {points[0][1]}"
    for (x0, y0), (x1, y1), (x2, y2) in zip(points, points[1:], points[2:]):
        ax, ay = x1 - R*((x1 > x0) - (x1 < x0)), y1 - R*((y1 > y0) - (y1 < y0))
        bx, by = x1 + R*((x2 > x1) - (x2 < x1)), y1 + R*((y2 > y1) - (y2 < y1))
        d += f"L{ax} {ay}Q{x1} {y1} {bx} {by}"
    d += f"L{points[-1][0]} {points[-1][1]}"
    tracks.append(d)


def text(x, y, s, anchor="middle", weight="normal"):
    r"""
    text: add a label.

        Parameters:
            x      (float): Position of the anchor
            y      (float): Vertical centre
            s        (str): Text
            anchor   (str): SVG text-anchor
            weight   (str): SVG font-weight

        Returns:
            None; the label is added to boxes

        Formula:
            \text{SVG text at } (x, y + 5)

        See also:
            box
    """
    boxes.append(f'<text x="{x}" y="{y + 5}" font-family="monospace" font-size="13" '
                 f'text-anchor="{anchor}" fill="#111111" font-weight="{weight}">{s}</text>')


def width(s):
    r"""
    width: the width of a box holding a text.

        Parameters:
            s (str): Text

        Returns:
            w (int): Width, in pixels

        Formula:
            w = 8\,|s| + 24

        See also:
            box
    """
    return 8*len(s) + 24


def box(x, y, s, fill=BLUE, rx=10, dash=None):
    r"""
    box: add a box with a label; boxes hide the track behind them.

        Parameters:
            x    (float): Left edge
            y    (float): Vertical centre
            s      (str): Label
            fill   (str): Fill colour
            rx   (float): Corner radius
            dash   (str): SVG dash pattern of the outline, or None

        Returns:
            right (float): Right edge of the box

        Formula:
            \text{rectangle } [x, x + w] \times [y - 13, y + 13],\ w = \text{width}(s)

        See also:
            width, text
    """
    d = f' stroke-dasharray="{dash}"' if dash else ""
    boxes.append(f'<rect x="{x}" y="{y - 13}" width="{width(s)}" height="26" rx="{rx}" fill="{fill}" '
                 f'stroke="#333333" stroke-width="1.5"{d}/>')
    text(x + width(s)/2, y, s)
    return x + width(s)


def bar(x, y):
    r"""
    bar: add a short vertical bar, marking where a track starts or stops.

        Parameters:
            x (float): Position
            y (float): Vertical centre

        Returns:
            None

        Formula:
            \text{segment } (x, y - 10) \to (x, y + 10)

        See also:
            line
    """
    bars.append(f'<path d="M{x} {y - 10}V{y + 10}" stroke="#333333" stroke-width="3"/>')


def branch(y0, y, x0=SPLIT, x1=JOIN):
    r"""
    branch: add a track that leaves the main line for another row and comes back.

        Parameters:
            y0 (float): Main line
            y  (float): Row of the branch
            x0 (float): Where it leaves
            x1 (float): Where it rejoins

        Returns:
            None

        Formula:
            (x_0, y_0) \to (x_0, y) \to (x_1, y) \to (x_1, y_0)

        See also:
            line, section
    """
    line((x0 - R, y0), (x0, y0), (x0, y), (x0 + 20, y), (x1 - 20, y), (x1, y), (x1, y0), (x1 + R, y0))


def section(y0, title, previous=None):
    r"""
    section: start the track of a section, coming from the end of the previous one.

        Parameters:
            y0       (float): Line of this section
            title      (str): Section title
            previous (float): Line of the previous section, or None for the first

        Returns:
            None

        Formula:
            \text{main line at } y_0, \text{ joined from the previous section's end}

        See also:
            branch
    """
    text(SPLIT, y0 - 26, title, anchor="start", weight="bold")
    if previous is None:
        bar(EDGE, y0)
        line((EDGE, y0), (JOIN, y0))
    else:
        mid = y0 - 50
        line((JOIN, previous), (TURN, previous), (TURN, mid), (EDGE, mid), (EDGE, y0), (JOIN, y0))


# Stage 0. Raw u, v goes on to the noise, or is preprocessed into a new source, where the track stops.
y0 = 50
section(y0, "Stage 0: source")
a1, a4, a2, surface, a3, w, pre = [y0 + i*ROW for i in range(1, 8)]
tap, merge, block = 300, 745, 760

box(SPLIT + 20, y0, "raw u, v: V1 V2 D1 D2")
line((tap - R, y0), (tap, y0), (tap, a1), (block, a1))
end = box(tap + 20, a1, "PP1 barotropic, baroclinic")
line((end + 20, a1), (end + 40, a1), (end + 40, a4), (end + 60, a4), (merge - 20, a4), (merge, a4), (merge, a1), (block, a1))
box(end + 60, a4, "PP4 spatial filter")
line((tap, a1 - R), (tap, a2), (merge - 20, a2), (merge, a2), (merge, a1), (block, a1))
box(tap + 20, a2, "PP2 w from divergence")
end = box(block, a1, "new source, preprocessed", fill=ORANGE, rx=2)
line((end, a1), (end + 14, a1))
bar(end + 14, a1)

for y, label in ((surface, "raw surface, SSH: V3"), (w, "raw vertical velocity: W1"), (pre, "preprocessed field: V4 V5 W1 S1")):
    branch(y0, y)
    box(SPLIT + 20, y, label)

line((tap - R, surface), (tap, surface), (tap, a3), (block, a3))
box(tap + 20, a3, "PP3 geostrophic u, v")
end = box(block, a3, "new source, preprocessed", fill=ORANGE, rx=2)
line((end, a3), (end + 14, a3))
bar(end + 14, a3)

# Conditioning
y1 = pre + 92
section(y1, "Conditioning: optional", y0)
branch(y1, y1 + ROW)
box(SPLIT + 20, y1 + ROW, "smooth, coastal taper: V2 D1 D2 W1")

# B1
y2 = y1 + ROW + 92
section(y2, "B1: noise quantity", y1)
end = box(SPLIT + 20, y2, "statistics: V1-V5")
branch(y2, y2 + ROW, end + 40, end + 240)
box(end + 60, y2 + ROW, "low-pass: V3 V4")
branch(y2, y2 + 2*ROW)
box(SPLIT + 20, y2 + 2*ROW, "length scale: V5")
branch(y2, y2 + 3*ROW)
end = box(SPLIT + 20, y2 + 3*ROW, "running variance: D1 D2")
box(end + 22, y2 + 3*ROW, "DMD or kernel DMD")
branch(y2, y2 + 4*ROW)
box(SPLIT + 20, y2 + 4*ROW, "DMD: W1 S1")

# B2 and B3
y3 = y2 + 4*ROW + 92
mid = y3 - 50
text(SPLIT, y3 - 26, "B2: sampling", anchor="start", weight="bold")
text(420, y3 - 26, "B3", anchor="start", weight="bold")
line((JOIN, y2), (TURN, y2), (TURN, mid), (EDGE, mid), (EDGE, y3), (420, y3))
box(SPLIT + 20, y3, "eORCA1")
branch(y3, y3 + ROW, SPLIT, 350)
box(SPLIT + 20, y3 + ROW, "Med (planned)", fill=GREY, dash="5 3")
end = box(420, y3, "file for NEMO")
line((end, y3), (end + 14, y3))
bar(end + 14, y3)

height = y3 + ROW + 30
with open("noise-pipeline-railroad.svg", "w") as f:
    f.write(f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{height}" viewBox="0 0 {W} {height}">'
            f'<rect width="100%" height="100%" fill="#ffffff"/>'
            # One path for all tracks, so that shared stretches are not drawn darker
            f'<path d="{"".join(tracks)}" fill="none" stroke="#333333" stroke-width="2"/>' + "\n".join(bars + boxes) + "</svg>\n")
