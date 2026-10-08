"""Approved transfer-block glyphs in prototype coordinates.

Placement applies one uniform affine transform; geometry and layering are fixed.
"""
from __future__ import annotations
import re
import numpy as np
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch, Polygon, Rectangle
from matplotlib.transforms import Affine2D

HX_PATH_D = r"""M 44.985885 192.86006
C 45.484725 192.86007 45.888672 193.26402 45.888672 193.76285
L 45.888672 198.62974
L 45.890739 198.62974
C 45.896611 203.44106 53.116138 203.43929 53.116138 198.62405
L 53.116138 185.07552
C 53.117138 183.87227 54.923273 183.87227 54.922229 185.07552
L 54.922229 198.62405
C 54.922229 203.44125 62.148145 203.44125 62.148145 198.62405
L 62.148145 185.07552
C 62.149745 183.87286 63.953113 183.87286 63.954753 185.07552
L 63.954753 198.62405
C 63.954753 203.44125 71.180668 203.44125 71.180668 198.62405
L 71.180668 185.07552
C 71.181668 183.87227 72.987803 183.87227 72.986759 185.07552
L 72.986759 198.62405
C 72.986759 203.44125 80.212675 203.44125 80.212675 198.62405
L 80.212675 193.76027
C 80.212678 193.26144 80.617144 192.85697 81.115979 192.85696
L 87.266508 192.85696
L 87.266508 190.14757
L 81.115979 190.14757
C 79.120629 190.14757 77.50328 191.76492 77.50328 193.76027
L 77.50328 198.62405
C 77.50228 199.8273 75.695628 199.8273 75.696672 198.62405
L 75.696672 185.07552
C 75.693272 180.26172 68.467358 180.26172 68.470756 185.07552
L 68.470756 198.62405
C 68.47011 199.82771 66.663502 199.82771 66.664148 198.62405
L 66.664148 185.07552
C 66.663598 182.66749 64.856871 181.46386 63.050415 181.46386
C 61.243959 181.46386 59.438183 182.66749 59.438749 185.07552
L 59.438749 198.62405
C 59.438103 199.82771 57.631495 199.82771 57.632141 198.62405
L 57.632141 185.07552
C 57.628741 180.26172 50.402828 180.26172 50.406226 185.07552
L 50.406226 198.62405
C 50.405226 199.8273 48.598574 199.8273 48.599618 198.62405
L 48.598584 198.62405
L 48.598584 193.76285
C 48.598584 191.7675 46.981235 190.15015 44.985885 190.15015
Z"""


def _tokenize_svg_path(d: str):
    return re.findall(r"[MLCZmlcz]|[-+]?(?:\d*\.\d+|\d+)(?:[eE][-+]?\d+)?", d)


def _svg_path_to_mpl(d: str) -> MplPath:
    tokens = _tokenize_svg_path(d)
    vertices, codes = [], []
    i = 0
    cmd = None
    current = (0.0, 0.0)
    start = (0.0, 0.0)

    def point(k):
        return float(tokens[k]), float(tokens[k + 1])

    while i < len(tokens):
        if re.fullmatch(r"[MLCZmlcz]", tokens[i]):
            cmd = tokens[i]
            i += 1

        if cmd in ("M", "m"):
            x, y = point(i); i += 2
            if cmd == "m":
                x += current[0]; y += current[1]
            current = (x, y)
            start = current
            vertices.append(current); codes.append(MplPath.MOVETO)
            cmd = "L" if cmd == "M" else "l"

        elif cmd in ("L", "l"):
            x, y = point(i); i += 2
            if cmd == "l":
                x += current[0]; y += current[1]
            current = (x, y)
            vertices.append(current); codes.append(MplPath.LINETO)

        elif cmd in ("C", "c"):
            pts = []
            for _ in range(3):
                x, y = point(i); i += 2
                if cmd == "c":
                    x += current[0]; y += current[1]
                pts.append((x, y))
            vertices.extend(pts)
            codes.extend([MplPath.CURVE4] * 3)
            current = pts[-1]

        elif cmd in ("Z", "z"):
            vertices.append(start)
            codes.append(MplPath.CLOSEPOLY)
            current = start
            cmd = None

        else:
            raise ValueError(f"Unsupported path command: {cmd!r}")

    return MplPath(np.asarray(vertices, float), codes)


HX_PATH = _svg_path_to_mpl(HX_PATH_D)
HX_BBOX = HX_PATH.get_extents()
SRC_PIPE_H = 2.709395
SRC_TRI_W = 36.683549 - 20.965645
SRC_TRI_H = 201.99904 - 181.14036
SRC_BAR_W = 39.392944 - 36.683549
PIPE_H = 0.12
SCALE = PIPE_H / SRC_PIPE_H
TRI_W = SRC_TRI_W * SCALE
TRI_H = SRC_TRI_H * SCALE
BAR_W = SRC_BAR_W * SCALE
# Keep the source bar thickness. Its triangle-facing edge reaches the
# intersection of the sloped triangle edge and the conduit edge (similar
# triangles); this overlap scales uniformly with the whole transfer block.
VALVE_BAR_WIDTH = BAR_W
VALVE_BAR_OVERLAP = TRI_W * PIPE_H / TRI_H
HX_VISIBLE_X0_SRC = HX_BBOX.x0
HX_VISIBLE_X1_SRC = 81.115979
HX_VISIBLE_W_SRC = HX_VISIBLE_X1_SRC - HX_VISIBLE_X0_SRC
HX_W = HX_VISIBLE_W_SRC * SCALE
HX_H = HX_BBOX.height * SCALE
HX_PORT_CENTER_Y_SRC = 0.25 * (190.15015 + 192.86006 + 190.14757 + 192.85696)
SEG_MID = 0.40
SEG_OUTER = 0.40
HX_Y_OFFSET = 0.004
BLOCK_WIDTH = 2 * SEG_OUTER + SEG_MID + HX_W + TRI_W + BAR_W

def mix(c0, c1, t):
    a = np.asarray(c0, float)
    b = np.asarray(c1, float)
    return tuple((1 - t) * a + t * b)


def gradient_rect(ax, x0, x1, y, h, c0, c1, zorder=2):
    """Horizontal gradient conduit with perfectly square ends."""
    if x1 <= x0:
        return
    n = 512
    arr = np.zeros((2, n, 3), float)
    for i in range(n):
        arr[:, i, :] = mix(c0, c1, i / (n - 1))
    ax.imshow(
        arr,
        extent=(x0, x1, y - h / 2, y + h / 2),
        origin="lower",
        aspect="auto",
        interpolation="bicubic",
        zorder=zorder,
    )


def gradient_polygon(ax, vertices, c0, c1, zorder=5):
    xs = [p[0] for p in vertices]
    ys = [p[1] for p in vertices]
    xmin, xmax = min(xs), max(xs)
    ymin, ymax = min(ys), max(ys)

    n = 512
    arr = np.zeros((2, n, 3), float)
    for i in range(n):
        arr[:, i, :] = mix(c0, c1, i / (n - 1))

    image = ax.imshow(
        arr,
        extent=(xmin, xmax, ymin, ymax),
        origin="lower",
        aspect="auto",
        interpolation="bicubic",
        zorder=zorder,
    )
    clip = Polygon(vertices, closed=True, facecolor="none", edgecolor="none")
    ax.add_patch(clip)
    image.set_clip_path(clip)


def add_exchanger(ax, x_left, y_center, color, mirror=False):
    """Clip only the long cylinder-side stub; align on the port centreline."""
    transform = Affine2D().translate(-HX_VISIBLE_X0_SRC, -HX_PORT_CENTER_Y_SRC)
    if mirror:
        transform = transform.scale(-1.0, 1.0).translate(HX_VISIBLE_W_SRC, 0.0)
    transform = transform.scale(SCALE, SCALE)
    transform = transform.translate(x_left, y_center + HX_Y_OFFSET)
    patch = PathPatch(
        HX_PATH, transform=transform + ax.transData,
        facecolor=color, edgecolor="none", zorder=1,
    )
    clip_rect = Rectangle(
        (x_left, y_center - 2.0 * HX_H), HX_W, 4.0 * HX_H,
        facecolor="none", edgecolor="none",
    )
    ax.add_patch(clip_rect)
    patch.set_clip_path(clip_rect)
    ax.add_patch(patch)


def valve_geometry(x_left, y, direction):
    """
    x_left is the left edge of the whole valve component.

    The transverse bar is on the TRIANGLE TIP side, exactly as in the
    original SVG diode. Closed = triangle + bar. Open = triangle only.
    """
    if direction == "right":
        base_x = x_left
        tip_x = base_x + TRI_W
        bar_x0 = tip_x
        vertices = [
            (base_x, y - TRI_H / 2),
            (base_x, y + TRI_H / 2),
            (tip_x, y),
        ]
        component_right = bar_x0 + BAR_W

    elif direction == "left":
        # Mirrored original: bar first, then triangle point, then triangle base.
        bar_x0 = x_left
        tip_x = bar_x0 + BAR_W
        base_x = tip_x + TRI_W
        vertices = [
            (base_x, y - TRI_H / 2),
            (base_x, y + TRI_H / 2),
            (tip_x, y),
        ]
        component_right = base_x

    else:
        raise ValueError(direction)

    return vertices, base_x, tip_x, bar_x0, component_right


def draw_valve(ax, x_left, y, direction, is_open, c_left, c_right):
    vertices, base_x, tip_x, bar_x0, x_right = valve_geometry(x_left, y, direction)

    bar_left = (tip_x + VALVE_BAR_OVERLAP - VALVE_BAR_WIDTH
                if direction == "left" else tip_x - VALVE_BAR_OVERLAP)

    if 0 < is_open < 1:
        # Only continuous model openings receive a visual opacity transition.
        # Keep the approved shapes, layering and conduit gradients unchanged.
        gradient_polygon(ax, vertices, c_left, c_right, zorder=6)
        ax.add_patch(Polygon(vertices, closed=True, facecolor="black",
                             edgecolor="none", zorder=6, alpha=1-float(is_open)))
        ax.add_patch(Rectangle((bar_left, y - TRI_H / 2), VALVE_BAR_WIDTH, TRI_H,
                               facecolor="black", edgecolor="none", zorder=7,
                               alpha=1-float(is_open)))
    elif is_open:
        gradient_polygon(ax, vertices, c_left, c_right, zorder=6)
    else:
        ax.add_patch(
            Polygon(vertices, closed=True, facecolor="black", edgecolor="none", zorder=6)
        )
        # Geometric bar thickness scales with the entire transfer block.
        ax.add_patch(
            Rectangle(
                (bar_left, y - TRI_H / 2),
                VALVE_BAR_WIDTH,
                TRI_H,
                facecolor="black",
                edgecolor="none",
                zorder=7,
            )
        )

    return base_x, tip_x, x_right


def tip_overlap():
    """Approved similar-triangle penetration with the visual safety factor."""
    return 1.60 * TRI_W * PIPE_H / TRI_H


def draw_row(
    ax,
    *,
    y,
    name,
    order,
    valve_direction,
    valve_open,
    left_color,
    right_color,
    hx_color,
    mirror_hx=False,
):
    """
    Exactly:
      conduit 1 + element 1 + conduit 2 + element 2 + conduit 3

    Outer and middle conduits have equal prototype lengths.
    No extra cylinder-side 'stub' is added.
    """
    x = 1.20
    x_start = x

    # component widths
    valve_w = TRI_W + BAR_W

    # Segment 1
    s1 = (x, x + SEG_OUTER)
    x = s1[1]

    # Element 1
    e1_w = HX_W if order[0] == "hx" else valve_w
    e1 = (x, x + e1_w)
    x = e1[1]

    # Segment 2
    s2 = (x, x + SEG_MID)
    x = s2[1]

    # Element 2
    e2_w = HX_W if order[1] == "hx" else valve_w
    e2 = (x, x + e2_w)
    x = e2[1]

    # Segment 3
    s3 = (x, x + SEG_OUTER)
    x_end = s3[1]

    hx = e1 if order[0] == "hx" else e2
    valve = e1 if order[0] == "valve" else e2

    # ---- Thermal gradients in conduits ---------------------------------
    # Interpolate from each cylinder up to the exchanger.
    if order == ("hx", "valve"):
        # Left: only s1 from S to HX
        gradient_rect(ax, *s1, y, PIPE_H, left_color, hx_color, zorder=2)

        # Right: s2 + valve + s3 from HX to L
        long_x0 = s2[0]
        long_x1 = s3[1]

        def right_c(xp):
            t = (xp - long_x0) / (long_x1 - long_x0)
            return mix(hx_color, right_color, np.clip(t, 0, 1))

        verts, base_x, tip_x, _, _ = valve_geometry(valve[0], y, valve_direction)
        overlap = tip_overlap()

        if valve_direction == "left":
            # tip faces HX: s2 penetrates under tip; s3 touches triangle base.
            gradient_rect(
                ax, s2[0], tip_x + overlap, y, PIPE_H,
                right_c(s2[0]), right_c(tip_x + overlap), zorder=3
            )
            gradient_rect(
                ax, base_x, s3[1], y, PIPE_H,
                right_c(base_x), right_c(s3[1]), zorder=3
            )
        else:
            gradient_rect(
                ax, s2[0], base_x, y, PIPE_H,
                right_c(s2[0]), right_c(base_x), zorder=3
            )
            gradient_rect(
                ax, tip_x - overlap, s3[1], y, PIPE_H,
                right_c(tip_x - overlap), right_c(s3[1]), zorder=3
            )

        cva0 = right_c(valve[0])
        cva1 = right_c(valve[1])

    elif order == ("valve", "hx"):
        # Left: s1 + valve + s2 from S to HX
        long_x0 = s1[0]
        long_x1 = s2[1]

        def left_c(xp):
            t = (xp - long_x0) / (long_x1 - long_x0)
            return mix(left_color, hx_color, np.clip(t, 0, 1))

        verts, base_x, tip_x, _, _ = valve_geometry(valve[0], y, valve_direction)
        overlap = tip_overlap()

        if valve_direction == "right":
            gradient_rect(
                ax, s1[0], base_x, y, PIPE_H,
                left_c(s1[0]), left_c(base_x), zorder=3
            )
            gradient_rect(
                ax, tip_x - overlap, s2[1], y, PIPE_H,
                left_c(tip_x - overlap), left_c(s2[1]), zorder=3
            )
        else:
            gradient_rect(
                ax, s1[0], tip_x + overlap, y, PIPE_H,
                left_c(s1[0]), left_c(tip_x + overlap), zorder=3
            )
            gradient_rect(
                ax, base_x, s2[1], y, PIPE_H,
                left_c(base_x), left_c(s2[1]), zorder=3
            )

        # Right: only s3 from HX to L
        gradient_rect(ax, *s3, y, PIPE_H, hx_color, right_color, zorder=2)

        cva0 = left_c(valve[0])
        cva1 = left_c(valve[1])

    else:
        raise ValueError(order)

    # Exchanger below conduits; valve above conduits.
    add_exchanger(ax, hx[0], y, hx_color, mirror=mirror_hx)
    draw_valve(
        ax,
        valve[0],
        y,
        valve_direction,
        valve_open,
        cva0,
        cva1,
    )

    return x_end


LAYOUTS = {
    "UU": {
        "Ho": {"order": ("hx", "valve"), "mirror_hx": True},
        "Hi": {"order": ("valve", "hx"), "mirror_hx": False},
    },
    "DD": {
        "Ho": {"order": ("valve", "hx"), "mirror_hx": False},
        "Hi": {"order": ("hx", "valve"), "mirror_hx": True},
    },
}
# Each placement is independently supported by the production hydraulic model.
LAYOUTS["UD"] = {"Ho": LAYOUTS["UU"]["Ho"], "Hi": LAYOUTS["DD"]["Hi"]}
LAYOUTS["DU"] = {"Ho": LAYOUTS["DD"]["Ho"], "Hi": LAYOUTS["UU"]["Hi"]}


class PlacedAxes:
    """Apply the same uniform transform to images, glyphs and clipping patches."""
    def __init__(self, axes, left_head, right_head, row_y, *, zorder_offset=0):
        self.axes = axes
        self.zorder_offset = zorder_offset
        self.factor = (right_head - left_head) / BLOCK_WIDTH
        if self.factor <= 0:
            raise ValueError("Cylinder inner faces must be ordered left to right.")
        self.placement = Affine2D().scale(self.factor).translate(
            left_head - 1.20 * self.factor, row_y)
        self.transData = self.placement + axes.transData

    def imshow(self, *args, **kwargs):
        kwargs["zorder"] = kwargs.get("zorder", 0) + self.zorder_offset
        return self.axes.imshow(*args, transform=self.transData, **kwargs)

    def add_patch(self, patch):
        patch.set_zorder(patch.get_zorder() + self.zorder_offset)
        if not patch.is_transform_set():
            patch.set_transform(self.transData)
        return self.axes.add_patch(patch)
