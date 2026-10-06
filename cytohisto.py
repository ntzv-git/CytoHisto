#!/usr/bin/env python3
"""CytoHisto: publication-ready flow cytometry histograms (version: see __version__).

Graphical tool to overlay one or more .fcs files, style each curve and save the
figure (PNG, TIFF, PDF, SVG). Works with any FCS 2.0/3.0/3.1 list-mode file.
"""

import csv
import json
import logging
import math
import os
import re
import sys
import tkinter as tk
from tkinter import colorchooser, filedialog, messagebox, ttk

import matplotlib
matplotlib.use("TkAgg")
import matplotlib.colors as mcolors
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.lines import Line2D
from matplotlib.ticker import FixedLocator, FuncFormatter

try:  # drag and drop (optional: the app still works without it)
    from tkinterdnd2 import COPY, DND_FILES, TkinterDnD
    BaseTk = TkinterDnD.Tk
except Exception:  # pragma: no cover
    COPY, DND_FILES, BaseTk = None, None, tk.Tk

# Sans-serif font: Arial on Windows, closest equivalent elsewhere (no warning if missing)
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Liberation Sans", "DejaVu Sans"]
logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)

APP = "CytoHisto"
__version__ = "1.1.1"
PALETTE = ["#000000", "#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
CHANNEL_CHOICES = ["256", "512", "1024", "2048", "4096", "8192", "16384", "32768", "65536",
                   "262144", "1048576"]
LINETYPES = {"solid": "-", "dashed": (0, (4, 2)), "dotted": (0, (1, 1.5)),
             "dash-dot": (0, (4, 1.5, 1, 1.5)), "long dash": (0, (8, 2))}
MM = 72 / 25.4        # points per mm
# independent text sizes (pt, for the reference width; they scale with the image width)
FONTS = [("labels", "Labels", 9), ("xtitle", "X title", 9), ("ytitle", "Y title", 9),
         ("xvalues", "X values", 7), ("yvalues", "Y values", 7)]
REF_WIDTH = 170.0     # figure width (mm) the default sizes are designed for
MAX_TICKS, MAX_LABELS = 400, 60   # beyond this, user steps are too small: automatic steps used


# ---------------------------------------------------------------------------
# FCS reader (handles parameters with different bit widths)
# ---------------------------------------------------------------------------
def read_fcs(path):
    """Return (columns, ranges): {name: values}, {name: $PnR}."""
    with open(path, "rb") as fh:
        raw = fh.read()
    header = raw[:58].decode("ascii", "replace")
    pos = [header[10:18], header[18:26], header[26:34], header[34:42]]
    pos = [int(p) if p.strip().isdigit() else 0 for p in pos]

    text = raw[pos[0]:pos[1] + 1].replace(b"\x00", b"").decode("latin-1")
    delim = text[0]
    text = text[1:].replace(delim * 2, "\x01")          # escaped delimiter
    fields = [f.replace("\x01", delim) for f in text.split(delim)]
    if fields and fields[-1] == "":
        fields = fields[:-1]
    kw = {fields[i].strip().upper(): fields[i + 1].strip() for i in range(0, len(fields) - 1, 2)}

    start = pos[2] or int(kw["$BEGINDATA"])
    npar, ntot = int(kw["$PAR"]), int(kw["$TOT"])
    names = [kw.get(f"$P{i}N", f"P{i}") for i in range(1, npar + 1)]
    bits = [int(kw[f"$P{i}B"]) for i in range(1, npar + 1)]
    ranges = [float(kw[f"$P{i}R"]) for i in range(1, npar + 1)]
    order = "<" if kw.get("$BYTEORD", "1,2,3,4") in ("1,2,3,4", "1,2") else ">"
    dtype = kw.get("$DATATYPE", "I")
    data = raw[start:]

    if dtype in ("F", "D"):
        dt = np.dtype(order + ("f4" if dtype == "F" else "f8"))
        m = np.frombuffer(data, dtype=dt, count=npar * ntot).reshape(ntot, npar)
        cols = {names[i]: m[:, i].astype(float) for i in range(npar)}
    elif dtype == "I":
        nbytes = [b // 8 for b in bits]
        row = sum(nbytes)
        m = np.frombuffer(data, dtype=np.uint8, count=row * ntot).reshape(ntot, row)
        cols, c = {}, 0
        for i in range(npar):
            part = m[:, c:c + nbytes[i]].astype(np.uint64)
            if order == ">":
                part = part[:, ::-1]
            v = np.zeros(ntot, dtype=np.uint64)
            for b in range(nbytes[i]):
                v |= part[:, b] << np.uint64(8 * b)
            used = math.ceil(math.log2(ranges[i])) if ranges[i] > 1 else bits[i]
            if used < bits[i]:                            # bits actually used ($PnR)
                v &= np.uint64((1 << used) - 1)
            cols[names[i]] = v.astype(float)
            c += nbytes[i]
    else:
        raise ValueError(f"$DATATYPE '{dtype}' is not supported")
    return cols, dict(zip(names, ranges))


def default_channel(names):
    """DNA if present, otherwise the first parameter that is not time."""
    for n in names:
        if n.upper() == "DNA":
            return n
    for n in names:
        if "TIME" not in n.upper():
            return n
    return names[0]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
CSS_EXTRA = {"rebeccapurple": "#663399"}


def css_color(x):
    """HTML/CSS colour -> hex (#RGB, #RRGGBB, rgb(r,g,b), CSS names)."""
    y = x.strip().lower().replace(" ", "")
    m = re.fullmatch(r"rgb\((\d+),(\d+),(\d+)\)", y)
    if m:
        return "#%02x%02x%02x" % tuple(min(255, int(v)) for v in m.groups())
    if re.fullmatch(r"#[0-9a-f]{3}", y):
        y = "#" + "".join(ch * 2 for ch in y[1:])
    y = CSS_EXTRA.get(y, y)
    try:
        return mcolors.to_hex(y)
    except ValueError:
        raise ValueError(f"Unknown colour: {x}")


def nice_step(span, n_max):
    """Round step (1, 2, 2.5 or 5 x 10^n) giving at most n_max intervals."""
    k = 10 ** math.floor(math.log10(span / n_max))
    for c in (1, 2, 2.5, 5, 10):
        if span / (k * c) <= n_max:
            return k * c
    return 10 * k


def number(txt, default=None):
    txt = str(txt).strip().replace(",", ".")
    return default if txt == "" else float(txt)


def numbers(txt, n_max):
    """'200,100,20', '200;100;20' or '200 100 20' -> [200, 100, 20] (1 to n_max positive numbers)."""
    txt = str(txt).strip()
    if not txt:
        return None
    v = [float(p) for p in re.split(r"[,;/\s]+", txt) if p]
    if not 1 <= len(v) <= n_max or min(v) <= 0:
        raise ValueError(f"'{txt}': expected 1 to {n_max} positive numbers")
    return v


# ---------------------------------------------------------------------------
# Computation and drawing
# ---------------------------------------------------------------------------
def prepare(files, g, cache):
    """Compute curves and statistics. Returns (curves, stats, context)."""
    out = float(g["out_channels"])
    xlim = sorted(g["xlim"]) if g["xlim"] else [0, out]

    def inside_edges(x):          # ignore debris (low end) and the saturation channel (high end)
        return (x > 0.05 * out) & (x < 0.95 * out)

    curves, stats = [], []
    for f in files:
        cols, ranges = cache[f["path"]]
        channel = f["channel"] if f["channel"] in cols else default_channel(list(cols))
        in_range = f.get("in_range") or ranges.get(channel) or 65536   # this file's own range
        v = cols[channel] * out / in_range                # input range -> output channels
        smooth = max(1.0, float(f["smooth"]))
        w = smooth
        nb = math.ceil(out / w)
        cnt = np.bincount(np.clip((v // w).astype(int), 0, nb - 1), minlength=nb)
        x = np.arange(nb) * w + (w - 1) / 2
        y = cnt / w                                       # events per channel
        vis = (x >= xlim[0]) & (x <= xlim[1])
        curves.append((x[vis], y[vis]))
        ok = inside_edges(x)
        s = {"file": os.path.basename(f["path"]), "label": f["name"], "events": len(v),
             "peak": round(float(x[ok][np.argmax(y[ok])]), 1) if ok.any() else "",
             "gate": "", "% in gate": "", "mean": "", "CV %": ""}
        if f.get("gate"):
            a, b = f["gate"]
            inside = v[(v >= a) & (v <= b)]
            s["gate"] = f"{a:g}-{b:g}"
            s["% in gate"] = round(100 * len(inside) / len(v), 1) if len(v) else ""
            if len(inside) > 1:
                s["mean"] = round(float(inside.mean()), 1)
                s["CV %"] = round(float(100 * inside.std(ddof=1) / inside.mean()), 2)
        stats.append(s)
    return curves, stats, {"xlim": xlim, "inside_edges": inside_edges}


def draw(files, g, cache, dpi=100):
    """Build the Matplotlib figure. Returns (figure, stats, info)."""
    k = g["width"] / REF_WIDTH
    fs = {key: g["fonts"][key] * k for key, _, _ in FONTS}
    axis_lw = 9 * k / 22 * 2.845     # axis line width (same rule as ggplot2)
    curves, stats, ctx = prepare(files, g, cache)
    xlim = ctx["xlim"]

    fig = Figure(figsize=(g["width"] / 25.4, g["height"] / 25.4), dpi=dpi)
    fig.patch.set_facecolor("white")
    ax = fig.add_subplot(111)

    # --- Y axis: steps estimated from the highest visible peak
    ymax = max([float(y.max()) for _, y in curves if len(y)] + [1.0])
    warnings = []
    if g["y_steps"]:
        y_major = g["y_steps"][0]
        y_minor = g["y_steps"][1] if len(g["y_steps"]) > 1 else y_major / 5
        if ymax / y_minor > MAX_TICKS or ymax / y_major > MAX_LABELS:
            warnings.append(f"Y steps {g['y_steps']} too small for a peak of {ymax:.0f}: automatic steps used")
            y_major = None
    if not g["y_steps"] or y_major is None:
        y_major = nice_step(ymax, 8)
        y_minor = y_major / (4 if round(y_major / 10 ** math.floor(math.log10(y_major)), 6) == 2 else 5)
    top = math.ceil(ymax * (1.12 if g["labels"] else 1.05) / y_minor) * y_minor
    y_maj = np.arange(0, top + 1e-9, y_major)
    y_min = np.arange(0, top + 1e-9, y_minor)

    # --- X axis: labels / grid + long ticks / short ticks, from the visible range
    xs, span = None, xlim[1] - xlim[0]
    if g["x_steps"]:               # labels, grid + long ticks, short ticks (missing ones derived)
        xs = list(g["x_steps"])
        if len(xs) == 1:
            xs.append(xs[0] / 2)
        if len(xs) == 2:
            xs.append(xs[1] / 5)
        if span / xs[2] > MAX_TICKS or span / xs[1] > MAX_TICKS / 2 or span / xs[0] > MAX_LABELS:
            warnings.append(f"X steps {g['x_steps']} too small for {span:g} channels: automatic steps used")
            xs = None
    if xs is None:
        p = nice_step(xlim[1] - xlim[0], 6)
        xs = [p, p / 2, p / 10]
    on = lambda step: np.arange(math.ceil(xlim[0] / step - 1e-9) * step, xlim[1] + 1e-9, step)
    x_grid, x_ticks = on(xs[1]), on(xs[2])

    # --- grid, gates, curves
    grey = "#d9d9d9"
    for yv in y_maj[1:]:
        ax.axhline(yv, color=grey, lw=0.25 * 2.845 * k, zorder=0)
    for xv in x_grid:
        if xlim[0] < xv <= xlim[1]:
            ax.axvline(xv, color=grey, lw=0.25 * 2.845 * k, zorder=0)

    for f in files:
        if f.get("gate"):
            a, b = f["gate"]
            if b > xlim[0] and a < xlim[1]:
                ax.axvspan(max(a, xlim[0]), min(b, xlim[1]), color=f["color"], alpha=0.15,
                           lw=0, zorder=1)

    for f, (x, y) in zip(files, curves):
        ax.plot(x, y, color=f["color"], lw=f["linewidth"] * 2.845 * k,
                ls=LINETYPES.get(f["linetype"], "-"), zorder=3, solid_joinstyle="round")

    # --- names: labels above the peaks, or a legend
    if g["labels"]:
        for f, (x, y) in zip(files, curves):
            zone = ctx["inside_edges"](x)
            if not zone.any():
                continue
            j = np.flatnonzero(zone)[np.argmax(y[zone])]
            ax.text(x[j], y[j] + 0.02 * top, f["name"], color=f["color"], ha="center",
                    va="bottom", fontsize=fs["labels"], zorder=4)
    elif len(files) > 1:
        handles = [Line2D([], [], color=f["color"], lw=f["linewidth"] * 2.845 * k,
                          ls=LINETYPES.get(f["linetype"], "-")) for f in files]
        ax.legend(handles, [f["name"] for f in files], loc="upper right", frameon=False,
                  prop={"size": fs["labels"]}, handlelength=3)

    # --- axes and ticks
    ax.set_ylim(0, top)
    ax.yaxis.set_major_locator(FixedLocator(y_maj))
    ax.yaxis.set_minor_locator(FixedLocator(y_min))
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
    ax.set_xlim(*xlim)
    ax.xaxis.set_major_locator(FixedLocator(x_grid))
    ax.xaxis.set_minor_locator(FixedLocator(x_ticks))
    ax.xaxis.set_major_formatter(FuncFormatter(
        lambda v, _: f"{v:g}" if abs(v / xs[0] - round(v / xs[0])) < 1e-9 else ""))

    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(axis_lw)
    ax.tick_params(axis="both", which="both", direction="out", width=axis_lw, colors="#333333",
                   labelcolor="black", pad=2 * k)
    ax.tick_params(axis="x", which="major", labelsize=fs["xvalues"])
    ax.tick_params(axis="y", which="major", labelsize=fs["yvalues"])
    ax.tick_params(axis="y", which="major", length=2.5 * MM * k)
    ax.tick_params(axis="y", which="minor", length=1.2 * MM * k)
    ax.tick_params(axis="x", which="major", length=1.5 * MM * k)
    ax.tick_params(axis="x", which="minor", length=0.8 * MM * k)
    ax.set_xlabel(g["xlab"], fontsize=fs["xtitle"], labelpad=3 * k)
    ax.set_ylabel(g["ylab"], fontsize=fs["ytitle"], labelpad=3 * k)
    fig.tight_layout(pad=0.6 * k)

    info = {"warnings": warnings, "ymax": ymax, "y_steps": (y_major, y_minor), "x_steps": xs,
            "ax": ax, "curves": [(f["name"], f["color"], x, y) for f, (x, y) in zip(files, curves)]}
    return fig, stats, info


# ---------------------------------------------------------------------------
# Graphical interface
# ---------------------------------------------------------------------------
class App(BaseTk):
    def __init__(self):
        super().__init__()
        self.title(f"{APP} {__version__} - flow cytometry histograms")
        self.geometry(f"{min(1350, self.winfo_screenwidth() - 40)}x"
                      f"{min(820, self.winfo_screenheight() - 80)}+10+10")
        self.minsize(900, 500)
        if sys.platform.startswith("win"):
            try:
                self.state("zoomed")
            except tk.TclError:
                pass
        self.cache = {}       # path -> (columns, ranges)
        self.files = []       # per-file settings
        self._loading = False
        self._timer = None
        self.canvas = None
        self.last_stats = []

        style = ttk.Style(self)
        for theme in ("vista", "clam"):
            if theme in style.theme_names():
                style.theme_use(theme)
                break

        panes = ttk.PanedWindow(self, orient="horizontal")
        panes.pack(fill="both", expand=True)
        left = self._scrollable_left(panes)
        right = ttk.Frame(panes, padding=6)
        panes.add(right, weight=1)

        # shared settings first, then the files and their own settings
        self._build_plot(left)
        self._build_image(left)
        self._build_files(left)
        self._build_curve(left)
        self._build_toolbar(right)
        self._build_preview(right)

        hint = ("Drag and drop .fcs files anywhere in this window, paste them (Ctrl+V) or click 'Add...'."
                if DND_FILES else "Click 'Add...' to load .fcs files (drag and drop unavailable: "
                "install the 'tkinterdnd2' package).")
        self.status = tk.StringVar(value=hint)
        ttk.Label(self, textvariable=self.status, anchor="w", relief="sunken", padding=(6, 2)).pack(
            fill="x", side="bottom")
        self._enable_curve_form(False)
        self._warnings = []
        self.bind_all("<Control-v>", self.paste_files, add="+")
        if DND_FILES:
            # register once the window is really on screen: registering earlier makes some
            # Linux desktops ignore the first drop
            self.after(400, self._register_drop)

    def _register_drop(self):
        if DND_FILES:
            self.update_idletasks()
            self.drop_target_register(DND_FILES)
            # accept the drop from the moment the pointer enters: without these handlers
            # some systems refuse the first drop and only the second one gets through
            self.dnd_bind("<<DropEnter>>", self._on_drop_enter)
            self.dnd_bind("<<DropPosition>>", lambda e: COPY)
            self.dnd_bind("<<Drop>>", self._on_drop)

    # ----- layout -----
    def _scrollable_left(self, panes):
        holder = ttk.Frame(panes)
        bg = ttk.Style(self).lookup("TFrame", "background") or "#f0f0f0"
        canvas = tk.Canvas(holder, highlightthickness=0, width=480, background=bg)
        bar = ttk.Scrollbar(holder, orient="vertical", command=canvas.yview)
        canvas.configure(yscrollcommand=bar.set)
        bar.pack(side="right", fill="y")
        canvas.pack(side="left", fill="both", expand=True)
        inner = ttk.Frame(canvas, padding=6)
        win = canvas.create_window((0, 0), window=inner, anchor="nw")
        inner.bind("<Configure>", lambda e: (canvas.configure(scrollregion=canvas.bbox("all")),
                                              canvas.configure(width=inner.winfo_reqwidth())))
        canvas.bind("<Configure>", lambda e: canvas.itemconfigure(win, width=e.width))

        # mouse wheel scrolls the panel only while the pointer is over it
        def wheel(e):
            canvas.yview_scroll(-1 if (getattr(e, "delta", 0) > 0 or e.num == 4) else 1, "units")

        def enter(_):
            for ev in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
                self.bind_all(ev, wheel)

        def leave(_):
            for ev in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
                self.unbind_all(ev)

        holder.bind("<Enter>", enter)
        holder.bind("<Leave>", leave)
        panes.add(holder, weight=0)
        return inner

    def _row(self, frame, r, text, var=None, width=10, widget=None, hint=None):
        ttk.Label(frame, text=text).grid(row=r, column=0, sticky="w", pady=1)
        w = widget or ttk.Entry(frame, textvariable=var, width=width)
        w.grid(row=r, column=1, sticky="we", pady=1)
        if hint:
            ttk.Label(frame, text=hint, foreground="#666").grid(row=r, column=2, sticky="w", padx=4)
        return w

    def _pair(self, parent, v1, v2, sep=" to ", w=7):
        fr = ttk.Frame(parent)
        ttk.Entry(fr, textvariable=v1, width=w).pack(side="left")
        ttk.Label(fr, text=sep).pack(side="left")
        ttk.Entry(fr, textvariable=v2, width=w).pack(side="left")
        return fr

    def _build_plot(self, parent):
        fr = ttk.LabelFrame(parent, text="Plot settings (all curves)", padding=6)
        fr.pack(fill="x")
        fr.columnconfigure(1, weight=1)
        self.v_xlab = tk.StringVar(value="Fluorescence intensity (channels)")
        self.v_ylab = tk.StringVar(value="Count")
        self.v_out = tk.StringVar(value="1024")
        self.v_xmin, self.v_xmax = tk.StringVar(), tk.StringVar()
        self.v_labels = tk.BooleanVar(value=True)
        self.v_xsteps, self.v_ysteps = tk.StringVar(), tk.StringVar()

        self._row(fr, 0, "X axis title", self.v_xlab, 30)
        self._row(fr, 1, "Y axis title", self.v_ylab, 30)
        cb = ttk.Combobox(fr, textvariable=self.v_out, width=10, values=CHANNEL_CHOICES)
        self._row(fr, 2, "Output channels", widget=cb, hint="X axis resolution (all curves)")
        self._row(fr, 4, "X zoom", widget=self._pair(fr, self.v_xmin, self.v_xmax), hint="empty = all")
        self._row(fr, 5, "X steps", self.v_xsteps, 12, hint="labels, grid, ticks (empty = auto)")
        self._row(fr, 6, "Y steps", self.v_ysteps, 12, hint="labels, ticks (empty = auto)")
        cf = ttk.Frame(fr)
        ttk.Checkbutton(cf, text="Names above the peaks (otherwise legend)",
                        variable=self.v_labels).pack(anchor="w")
        cf.grid(row=7, column=0, columnspan=3, sticky="w", pady=(2, 0))
        for v in (self.v_xlab, self.v_ylab, self.v_out, self.v_xmin, self.v_xmax,
                  self.v_labels, self.v_xsteps, self.v_ysteps):
            v.trace_add("write", lambda *a: self.schedule())

    def _build_image(self, parent):
        fr = ttk.LabelFrame(parent, text="Image", padding=6)
        fr.pack(fill="x", pady=(6, 0))
        self.v_pxw, self.v_pxh = tk.StringVar(value="2000"), tk.StringVar(value="1000")
        r1 = self._pair(fr, self.v_pxw, self.v_pxh, " x ", 6)
        ttk.Label(r1, text=" px").pack(side="left")
        self._row(fr, 0, "Size", widget=r1)
        # one size per text element, laid out on two lines
        self.v_fonts = {key: tk.StringVar(value=f"{size:g}") for key, _, size in FONTS}
        grid = ttk.Frame(fr)
        for n, (key, label, _) in enumerate(FONTS):
            r, c = divmod(n, 3)
            ttk.Label(grid, text=label).grid(row=r, column=2 * c, sticky="w", padx=(0 if c == 0 else 8, 2))
            ttk.Entry(grid, textvariable=self.v_fonts[key], width=4).grid(row=r, column=2 * c + 1, pady=1)
        self._row(fr, 1, "Font (pt)", widget=grid)
        ttk.Label(fr, text="Text sizes, lines and ticks scale with the image width.",
                  foreground="#666").grid(row=2, column=0, columnspan=3, sticky="w")
        for v in (self.v_pxw, self.v_pxh, *self.v_fonts.values()):
            v.trace_add("write", lambda *a: self.schedule())

    def _build_files(self, parent):
        fr = ttk.LabelFrame(parent, text="FCS files (overlaid)", padding=6)
        fr.pack(fill="x", pady=(10, 0))
        cols = ("file", "label", "input")
        self.tree = ttk.Treeview(fr, columns=cols, show="tree headings", height=4, selectmode="browse")
        self.tree.heading("#0", text="")
        self.tree.column("#0", width=52, minwidth=52, stretch=False, anchor="w")
        for c, t, w in zip(cols, ("File", "Label", "Input ch."), (160, 140, 70)):
            self.tree.heading(c, text=t)
            self.tree.column(c, width=w, anchor="w")
        self.tree.pack(fill="x")
        self.tree.bind("<<TreeviewSelect>>", lambda e: self._fill_curve_form())
        self._swatches = {}      # colour -> small square image shown in the list
        b = ttk.Frame(fr)
        b.pack(fill="x", pady=(4, 0))
        ttk.Button(b, text="Add...", command=self.add_dialog).pack(side="left")
        ttk.Button(b, text="Remove", command=self.remove).pack(side="left", padx=4)
        ttk.Button(b, text="Up", width=4, command=lambda: self.move(-1)).pack(side="left")
        ttk.Button(b, text="Down", width=5, command=lambda: self.move(1)).pack(side="left", padx=4)
        ttk.Label(b, text="or drop files / Ctrl+V" if DND_FILES else "or Ctrl+V",
                  foreground="#666").pack(side="left", padx=6)

    def _build_curve(self, parent):
        fr = ttk.LabelFrame(parent, text="Selected curve", padding=6)
        fr.pack(fill="x", pady=(6, 0))
        fr.columnconfigure(1, weight=1)
        self.v_name, self.v_color = tk.StringVar(), tk.StringVar()
        self.v_lw, self.v_lt = tk.StringVar(), tk.StringVar()
        self.v_smooth, self.v_channel = tk.StringVar(), tk.StringVar()
        self.v_gmin, self.v_gmax = tk.StringVar(), tk.StringVar()
        self.curve_widgets = []

        self.curve_widgets.append(self._row(fr, 0, "Label", self.v_name, 22, hint="name on the plot"))
        cf = ttk.Frame(fr)
        e = ttk.Entry(cf, textvariable=self.v_color, width=12)
        e.pack(side="left", fill="x", expand=True)
        self.swatch = tk.Label(cf, width=3, relief="solid", bd=1)
        self.swatch.pack(side="left", padx=4)
        bc = ttk.Button(cf, text="Pick...", command=self.pick_color)
        bc.pack(side="left")
        self._row(fr, 1, "Colour", widget=cf)
        self.curve_widgets += [e, bc]
        sp = ttk.Spinbox(fr, textvariable=self.v_lw, from_=0.1, to=5, increment=0.1, width=8)
        self.curve_widgets.append(self._row(fr, 2, "Line width", widget=sp))
        cb = ttk.Combobox(fr, textvariable=self.v_lt, values=list(LINETYPES), state="readonly", width=12)
        self.curve_widgets.append(self._row(fr, 3, "Line type", widget=cb))
        ss = ttk.Spinbox(fr, textvariable=self.v_smooth, from_=1, to=20, increment=1, width=8)
        self.curve_widgets.append(self._row(fr, 4, "Smoothing", widget=ss, hint="channels per bin"))
        gf = ttk.Frame(fr)
        e1 = ttk.Entry(gf, textvariable=self.v_gmin, width=7)
        e2 = ttk.Entry(gf, textvariable=self.v_gmax, width=7)
        e1.pack(side="left"); ttk.Label(gf, text=" to ").pack(side="left"); e2.pack(side="left")
        self._row(fr, 5, "Gate", widget=gf, hint="empty = none")
        self.curve_widgets += [e1, e2]
        self.cb_channel = ttk.Combobox(fr, textvariable=self.v_channel, state="readonly", width=12)
        self.curve_widgets.append(self._row(fr, 6, "Parameter", widget=self.cb_channel))
        self.v_in = tk.StringVar()
        cb_in = ttk.Combobox(fr, textvariable=self.v_in, width=12, values=CHANNEL_CHOICES)
        self._editable_combos = {cb_in}          # typed values allowed (not read-only)
        self.curve_widgets.append(self._row(fr, 7, "Input channels", widget=cb_in,
                                            hint="read from this file ($PnR)"))
        for v in (self.v_name, self.v_color, self.v_lw, self.v_lt, self.v_smooth, self.v_channel,
                  self.v_gmin, self.v_gmax, self.v_in):
            v.trace_add("write", lambda *a: self._curve_form_changed())

    def _build_toolbar(self, parent):
        bar = ttk.Frame(parent)
        bar.pack(fill="x", pady=(0, 6))
        ttk.Button(bar, text="Save figure...", command=self.save_figure).pack(side="left")
        ttk.Button(bar, text="Statistics (CSV)...", command=self.save_stats).pack(side="left", padx=4)
        ttk.Button(bar, text="Save settings...", command=self.save_settings).pack(side="right")
        ttk.Button(bar, text="Open settings...", command=self.open_settings).pack(side="right", padx=4)
        self.readout = tk.StringVar(value="")
        ttk.Label(parent, textvariable=self.readout, foreground="#333", anchor="w").pack(fill="x")

    def _build_preview(self, parent):
        self.area = ttk.Frame(parent, relief="sunken", borderwidth=1)
        self.area.pack(fill="both", expand=True)
        self.area.bind("<Configure>", lambda e: self.schedule())
        fr = ttk.LabelFrame(parent, text="Statistics", padding=4)
        fr.pack(fill="x", pady=(6, 0))
        cols = ("file", "label", "events", "peak", "gate", "% in gate", "mean", "CV %")
        self.table = ttk.Treeview(fr, columns=cols, show="headings", height=4)
        for c in cols:
            self.table.heading(c, text=c[0].upper() + c[1:])
            self.table.column(c, width=130 if c in ("file", "label") else 70, minwidth=50,
                              anchor="center", stretch=False)
        hs = ttk.Scrollbar(fr, orient="horizontal", command=self.table.xview)
        self.table.configure(xscrollcommand=hs.set)
        self.table.pack(fill="x")
        hs.pack(fill="x")

    # ----- files -----
    def add_dialog(self):
        paths = filedialog.askopenfilenames(title="FCS files",
                                            filetypes=[("FCS files", "*.fcs *.FCS"), ("All files", "*.*")])
        self.add_files(paths)

    def _on_drop_enter(self, event):
        print("[CytoHisto] drag entered the window", flush=True)
        return COPY

    def _on_drop(self, event):
        print(f"[CytoHisto] drop received, data = {event.data!r}", flush=True)
        paths = [p for p in self.tk.splitlist(event.data) if p.lower().endswith(".fcs")]
        if not paths:
            self.status.set("Drop ignored: no .fcs file in it (try again, or copy the files "
                            "and press Ctrl+V).")
            return COPY
        self.after(20, lambda: self._after_drop(paths))  # finish the drop before reading files
        return COPY

    def _after_drop(self, paths):
        self.add_files(paths)
        # some desktops (e.g. Wayland) do not repaint a window that does not have the focus:
        # bring it to the front and redraw now, so the first drop is visible immediately
        self.lift()
        self.focus_force()
        self.update_idletasks()
        self.preview()

    def paste_files(self, _=None):
        """Ctrl+V: add the .fcs files copied in the file manager (or paths copied as text)."""
        if isinstance(self.focus_get(), (tk.Entry, ttk.Entry, ttk.Combobox, ttk.Spinbox)):
            return None                      # normal paste inside a text field
        try:
            text = self.clipboard_get()
        except tk.TclError:
            text = ""
        paths = []
        for line in text.replace("\r", "\n").split("\n"):
            line = line.strip()
            if line.startswith("file://"):
                from urllib.parse import unquote, urlparse
                line = unquote(urlparse(line).path)
                if re.match(r"^/[A-Za-z]:/", line):   # Windows file URI
                    line = line[1:]
            if line.lower().endswith(".fcs") and os.path.isfile(line):
                paths.append(line)
        if paths:
            self.add_files(paths)
        else:
            self.status.set("Nothing to paste: copy .fcs files (or their paths) first.")
        return "break"

    def add_files(self, paths):
        added, skipped = 0, []
        for p in paths:
            if any(f["path"] == p for f in self.files):   # same file dropped twice
                skipped.append(os.path.basename(p))
                continue
            try:
                if p not in self.cache:
                    self.cache[p] = read_fcs(p)
            except Exception as e:
                messagebox.showerror(APP, f"Cannot read {os.path.basename(p)}:\n{e}")
                continue
            cols, ranges = self.cache[p]
            channel = default_channel(list(cols))
            self.files.append({
                "path": p, "name": os.path.splitext(os.path.basename(p))[0],
                "color": PALETTE[len(self.files) % len(PALETTE)], "linewidth": 0.5,
                "linetype": "solid", "smooth": 1, "gate": None, "channel": channel,
                "in_range": ranges.get(channel) or 65536})   # each file keeps its own range
            added += 1
        if added:
            self._refresh_tree(len(self.files) - 1)
            self.schedule()
        if skipped:
            self.status.set("Already loaded: " + ", ".join(skipped))

    def remove(self):
        i = self._index()
        if i is None:
            return
        del self.files[i]
        self._refresh_tree(min(i, len(self.files) - 1))
        self.schedule()

    def move(self, d):
        i = self._index()
        if i is None or not 0 <= i + d < len(self.files):
            return
        self.files[i], self.files[i + d] = self.files[i + d], self.files[i]
        self._refresh_tree(i + d)
        self.schedule()

    def _index(self):
        sel = self.tree.selection()
        return int(sel[0]) if sel else None

    def _swatch(self, color):
        if color not in self._swatches:
            img = tk.PhotoImage(width=22, height=12)
            img.put("#444444", to=(0, 0, 22, 12))
            img.put(color, to=(1, 1, 21, 11))
            self._swatches[color] = img
        return self._swatches[color]

    def _row_values(self, f):
        return dict(image=self._swatch(f["color"]),
                    values=(os.path.basename(f["path"]), f["name"], f"{self._in_range(f):g}"))

    def _refresh_tree(self, select=None):
        self.tree.delete(*self.tree.get_children())
        for i, f in enumerate(self.files):
            self.tree.insert("", "end", iid=str(i), text="", **self._row_values(f))
        if select is not None and 0 <= select < len(self.files):
            self.tree.selection_set(str(select))
        else:
            self._loading = True             # empty the form without touching any file
            for v in (self.v_name, self.v_color, self.v_lw, self.v_lt, self.v_smooth,
                      self.v_channel, self.v_gmin, self.v_gmax, self.v_in):
                v.set("")
            self._loading = False
            self.swatch.configure(bg=self.cget("bg"))
            self._enable_curve_form(False)

    def _enable_curve_form(self, on):
        for w in self.curve_widgets:
            try:
                editable = not isinstance(w, ttk.Combobox) or w in self._editable_combos
                w.configure(state=("normal" if editable else "readonly") if on else "disabled")
            except tk.TclError:
                pass

    def _fill_curve_form(self):
        i = self._index()
        if i is None:
            self._enable_curve_form(False)
            return
        f = self.files[i]
        self._loading = True
        self._enable_curve_form(True)
        self.cb_channel.configure(values=list(self.cache[f["path"]][0]))
        self.v_name.set(f["name"]); self.v_color.set(f["color"])
        self.v_lw.set(f"{f['linewidth']:g}"); self.v_lt.set(f["linetype"])
        self.v_smooth.set(f"{f['smooth']:g}"); self.v_channel.set(f["channel"])
        self.v_gmin.set(f"{f['gate'][0]:g}" if f["gate"] else "")
        self.v_gmax.set(f"{f['gate'][1]:g}" if f["gate"] else "")
        self.v_in.set(f"{self._in_range(f):g}")
        self.swatch.configure(bg=f["color"])
        self._loading = False

    def _in_range(self, f):
        """Input range of a file: set by the user, otherwise read from the file ($PnR)."""
        return f.get("in_range") or self.cache[f["path"]][1].get(f["channel"]) or 65536

    def _curve_form_changed(self):
        if self._loading:
            return
        i = self._index()
        if i is None:
            return
        f = self.files[i]
        f["name"] = self.v_name.get()
        try:
            f["color"] = css_color(self.v_color.get())
            self.swatch.configure(bg=f["color"])
        except ValueError:
            pass
        for key, var in (("linewidth", self.v_lw), ("smooth", self.v_smooth)):
            try:
                val = number(var.get())
                if val and val > 0:
                    f[key] = val
            except ValueError:
                pass
        f["linetype"] = self.v_lt.get() or "solid"
        new_channel = self.v_channel.get() or f["channel"]
        if new_channel != f["channel"]:          # other parameter: take its own range from the file
            f["channel"] = new_channel
            f["in_range"] = self.cache[f["path"]][1].get(new_channel) or 65536
            self._loading = True
            self.v_in.set(f"{f['in_range']:g}")
            self._loading = False
        else:
            try:
                val = number(self.v_in.get())
                if val and val > 1:
                    f["in_range"] = val
            except ValueError:
                pass
        try:
            a, b = number(self.v_gmin.get()), number(self.v_gmax.get())
            f["gate"] = sorted([a, b]) if a is not None and b is not None else None
        except ValueError:
            pass
        self.tree.item(str(i), **self._row_values(f))
        self.schedule()

    def pick_color(self):
        i = self._index()
        if i is not None:
            c = colorchooser.askcolor(color=self.files[i]["color"], title="Curve colour")[1]
            if c:
                self.v_color.set(c)

    # ----- shared settings -----
    def settings(self):
        g = {"xlab": self.v_xlab.get(), "ylab": self.v_ylab.get(), "labels": self.v_labels.get()}
        g["out_channels"] = number(self.v_out.get(), 1024)
        if g["out_channels"] <= 1:
            raise ValueError("Channel numbers must be positive.")
        xmin, xmax = number(self.v_xmin.get()), number(self.v_xmax.get())
        g["xlim"] = [xmin, xmax] if xmin is not None and xmax is not None and xmin != xmax else None
        # an invalid step is ignored (automatic steps) instead of blocking the preview
        self._warnings = []
        for key, var, n_max, label in (("x_steps", self.v_xsteps, 3, "X steps"),
                                       ("y_steps", self.v_ysteps, 2, "Y steps")):
            try:
                g[key] = numbers(var.get(), n_max)
            except ValueError:
                g[key] = None
                self._warnings.append(f"{label} ignored ('{var.get()}': use e.g. 200,100,20)")
        # the figure is laid out at REF_WIDTH mm wide, then rendered at the requested pixels
        pxw, pxh = number(self.v_pxw.get(), 2000), number(self.v_pxh.get(), 1000)
        if not 100 <= pxw <= 20000 or not 100 <= pxh <= 20000:
            raise ValueError("Image size must be between 100 and 20000 pixels.")
        g["width"], g["height"] = REF_WIDTH, REF_WIDTH * pxh / pxw
        g["px"], g["px_height"] = pxw, pxh
        g["fonts"] = {}
        for key, label, default in FONTS:
            try:
                size = number(self.v_fonts[key].get(), default)
            except ValueError:
                size = default
            g["fonts"][key] = size if size and 1 <= size <= 72 else default
        if g["width"] <= 0 or g["height"] <= 0:
            raise ValueError("Width and height must be positive.")
        return g

    # ----- preview -----
    def schedule(self):
        if self._timer:
            self.after_cancel(self._timer)
        self._timer = self.after(600, self.preview)

    def _clear_preview(self):
        if self.canvas:
            self.canvas.get_tk_widget().destroy()
            self.canvas = None
        self.last_stats = []
        self.table.delete(*self.table.get_children())
        self.readout.set("")

    def preview(self):
        self._timer = None
        if not self.files:
            self._clear_preview()
            return
        try:
            g = self.settings()
            aw, ah = max(self.area.winfo_width() - 10, 100), max(self.area.winfo_height() - 10, 100)
            # real figure size, scaled to the area: the preview has the exact proportions
            dpi = min(aw / (g["width"] / 25.4), ah / (g["height"] / 25.4))
            fig, stats, info = draw(self.files, g, self.cache, dpi=dpi)
        except Exception as e:
            self.status.set(f"Error: {e}")
            return
        if self.canvas:
            self.canvas.get_tk_widget().destroy()
        # matplotlib gives the keyboard focus to every new canvas, which would interrupt
        # typing in a field: disable that while the canvas is created
        grab = tk.Canvas.focus_set
        tk.Canvas.focus_set = lambda *a: None
        try:
            self.canvas = FigureCanvasTkAgg(fig, master=self.area)
        finally:
            tk.Canvas.focus_set = grab
        self.canvas.draw()
        self.canvas.get_tk_widget().place(relx=0.5, rely=0.5, anchor="center")
        self._setup_tracking(info)
        self.last_stats = stats
        self.table.delete(*self.table.get_children())
        for s in stats:
            self.table.insert("", "end", values=[s[c] for c in self.table["columns"]])
        xs = info["x_steps"]
        warnings = self._warnings + info["warnings"]
        if warnings:
            self.status.set("Warning: " + "; ".join(warnings))
            return
        self.status.set(f"Y axis: highest peak {info['ymax']:.0f} -> steps {info['y_steps'][0]:g} / "
                        f"{info['y_steps'][1]:g}" + (f"   |   X axis: labels {xs[0]:g}, grid {xs[1]:g}, "
                                                     f"ticks {xs[2]:g}" if xs else ""))

    # ----- coordinates under the mouse -----
    def _setup_tracking(self, info):
        """Vertical cursor line + readout of X, Y and each curve's value at X."""
        ax, canvas = info["ax"], self.canvas
        self._track = {"ax": ax, "curves": info["curves"], "bg": None,
                       "line": ax.axvline(np.nan, color="#888888", lw=0.8, ls=(0, (3, 2)),
                                          zorder=5, animated=True)}
        self.readout.set("Move the mouse over the plot to read the coordinates.")

        def grab_background(_=None):
            self._track["bg"] = canvas.copy_from_bbox(ax.bbox)

        def move(e):
            t = self._track
            if e.inaxes is not ax or e.xdata is None:
                leave()
                return
            self.readout.set(f"x = {e.xdata:.1f}     y = {e.ydata:.0f}")
            if t["bg"] is None:
                return
            canvas.restore_region(t["bg"])
            t["line"].set_xdata([e.xdata, e.xdata])
            ax.draw_artist(t["line"])
            canvas.blit(ax.bbox)

        def leave(_=None):
            t = self._track
            self.readout.set("Move the mouse over the plot to read the coordinates.")
            if t["bg"] is not None:
                canvas.restore_region(t["bg"])
                canvas.blit(ax.bbox)

        canvas.mpl_connect("draw_event", grab_background)
        canvas.mpl_connect("motion_notify_event", move)
        canvas.mpl_connect("axes_leave_event", leave)
        canvas.mpl_connect("figure_leave_event", leave)
        grab_background()

    # ----- saving -----
    def save_figure(self):
        if not self.files:
            messagebox.showinfo(APP, "Add at least one .fcs file first.")
            return
        formats = {"PNG": ".png", "TIFF": ".tiff", "PDF": ".pdf", "SVG": ".svg"}
        chosen = tk.StringVar(self, value="PNG")   # format picked in the dialog's type list
        path = filedialog.asksaveasfilename(
            title="Save figure", initialfile="histogram", typevariable=chosen,
            filetypes=[(name, "*" + ext) for name, ext in formats.items()])
        if not path:
            return
        # Tk's Linux dialog ignores the selected type when adding an extension: use the
        # extension typed by the user if it is a known one, otherwise the selected format
        if os.path.splitext(path)[1].lower() not in (".png", ".tif", ".tiff", ".pdf", ".svg"):
            path += formats.get(chosen.get(), ".png")
        try:
            g = self.settings()
            dpi = g["px"] / (g["width"] / 25.4)
            fig, _, _ = draw(self.files, g, self.cache, dpi=dpi)
            ext = os.path.splitext(path)[1].lower()
            opts = {"pil_kwargs": {"compression": "tiff_lzw"}} if ext in (".tif", ".tiff") else {}
            fig.savefig(path, dpi=dpi, facecolor="white", **opts)
        except Exception as e:
            messagebox.showerror(APP, f"Could not save the figure:\n{e}")
            return
        size = ""
        if ext in (".png", ".tif", ".tiff"):
            size = f" ({round(g['width'] / 25.4 * dpi)} x {round(g['height'] / 25.4 * dpi)} px)"
        self.status.set(f"Figure saved: {path}{size}")

    def save_stats(self):
        if not self.last_stats:
            return
        path = filedialog.asksaveasfilename(title="Statistics", defaultextension=".csv",
                                            filetypes=[("CSV", "*.csv")], initialfile="statistics")
        if not path:
            return
        with open(path, "w", newline="", encoding="utf-8-sig") as fh:
            w = csv.DictWriter(fh, fieldnames=list(self.last_stats[0]))
            w.writeheader()
            w.writerows(self.last_stats)
        self.status.set(f"Statistics saved: {path}")

    def _shared_vars(self):
        return {"xlab": self.v_xlab, "ylab": self.v_ylab, 
                "out_channels": self.v_out, "xmin": self.v_xmin, "xmax": self.v_xmax,
                "labels": self.v_labels, "x_steps": self.v_xsteps,
                "y_steps": self.v_ysteps, "px_width": self.v_pxw, "px_height": self.v_pxh,
                **{f"font_{key}": var for key, var in self.v_fonts.items()}}

    def save_settings(self):
        path = filedialog.asksaveasfilename(title="Save settings", defaultextension=".json",
                                            filetypes=[("Settings", "*.json")], initialfile="settings")
        if not path:
            return
        data = {"shared": {k: v.get() for k, v in self._shared_vars().items()}, "files": self.files}
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(data, fh, ensure_ascii=False, indent=2)
        self.status.set(f"Settings saved: {path}")

    def open_settings(self):
        path = filedialog.askopenfilename(title="Open settings", filetypes=[("Settings", "*.json")])
        if not path:
            return
        try:
            with open(path, encoding="utf-8") as fh:
                data = json.load(fh)
            for k, v in data.get("shared", {}).items():
                if k in self._shared_vars():
                    self._shared_vars()[k].set(v)
            files = []
            for f in data.get("files", []):
                if f["path"] not in self.cache:
                    self.cache[f["path"]] = read_fcs(f["path"])
                files.append(f)
            self.files = files
        except Exception as e:
            messagebox.showerror(APP, f"Unreadable settings file:\n{e}")
            return
        self._refresh_tree(0 if self.files else None)
        self.schedule()


if __name__ == "__main__":
    if "--version" in sys.argv[1:]:
        print(f"{APP} {__version__}")
        sys.exit(0)
    App().mainloop()
