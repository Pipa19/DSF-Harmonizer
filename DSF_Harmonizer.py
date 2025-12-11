# PARA CORRER ESTE SCRIPT IR A SU DIRECTORIO EN LA TERMINAL Y CORRER:
    # (El .gdsf debe de estar en el mismo directorio que el script (o si no especificas el path)!)

    # > > > python3 dsf_step_fixer.py ARCHIVO_DSF.gdsf < < <

import sys
import os
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from tkinter import font as tkfont
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib import gridspec

# try Savitzky–Golay; fallback to simple moving average if not available
try:
    from scipy.signal import savgol_filter as _savgol
except Exception:
    _savgol = None

APP_TITLE = "DSF Harmonizer"


# ------------------------ helpers (robust stats & smoothing) ------------------------
def robust_mad_sigma(diffs):
    if len(diffs) == 0:
        return 0.0
    diffs = np.asarray(diffs, dtype=float)
    med = np.nanmedian(diffs)
    mad = np.nanmedian(np.abs(diffs - med))
    return 1.4826 * mad


def _odd(n):
    n = int(max(3, n))
    return n if n % 2 == 1 else n + 1


def smooth_signal(y, strength=0):
    """Smooth y with Savitzky–Golay if available; else symmetric moving average.
    strength: 0..100 (0 = off). Higher => wider window.
    """
    y = np.asarray(y, dtype=float)
    n = len(y)
    if n < 3 or strength <= 0:
        return y.copy()

    # Map 0..100 -> window fraction ~0.03..0.25 of n (clamped & odd)
    frac = max(0.0, min(1.0, float(strength) / 100.0))
    w_target = int(round(0.03 * n + 0.22 * frac * n))
    w = _odd(max(5, min(w_target, 101, n - (1 - n % 2))))  # clamp
    if w >= n:
        w = _odd(max(5, n - 1))
    if w < 5:
        return y.copy()

    if _savgol is not None and w >= 5:
        poly = min(3, max(2, w - 2))
        try:
            return _savgol(y, window_length=w, polyorder=poly, mode="interp")
        except Exception:
            pass

    # Fallback: symmetric moving average with edge reflection
    k = max(5, min(w, n - 1))
    kernel = np.ones(k) / k
    ypad = np.r_[y[k-1:0:-1], y, y[-2:-k-1:-1]]
    ys = np.convolve(ypad, kernel, mode="valid")
    start = (len(ys) - n) // 2
    return ys[start:start+n]


def smooth_derivative(d, strength=20):
    """Light smoothing for derivative trace; strength is 0..100 like smooth_signal."""
    return smooth_signal(d, strength=max(0, int(strength)))


# ------------------------ main app ------------------------
class DSF_Harmonizer(tk.Tk):
    def __init__(self, path=None):
        super().__init__()
        self.title(APP_TITLE)
        self.geometry("1380x980")

        # data / state
        self.df_orig = None
        self.df_work = None
        self.wells = []
        self.per_well_orig = {}
        self.per_well_work = {}
        self.current_well = None
        self.selected_idx = None
        self.max_idx = 0

        # masks / trimming
        self.trim_ranges = {}

        # per-well Tm cache & Tm-outlier list
        self.tm_values = {}
        self.tm_outlier_wells = []
        self.tm_sorted_wells = []

        # Umbral y Tm de referencia para definir outliers
        self.tm_thr_var = tk.StringVar(value="20")
        self.tm_ref_var = tk.StringVar(value="")
        self.tm_mean_var = tk.StringVar(value="Current mean Tm: n/a")

        # settings
        self.op_var = tk.StringVar(value="auto")
        self.status_var = tk.StringVar(value="")
        self.abs_thr_var = tk.StringVar(value="0")
        self.kdisp_var = tk.StringVar(value="6")
        self.disp_method = tk.StringVar(value="MAD")
        self.iterative_var = tk.BooleanVar(value=False)
        self.show_deriv_var = tk.BooleanVar(value=True)
        self.animate_var = tk.BooleanVar(value=True)
        self.multi_jump_var = tk.BooleanVar(value=True)

        # smoothing controls
        self.smooth_on_var = tk.BooleanVar(value=False)
        self.smooth_strength_var = tk.IntVar(value=35)

        # T-range for analysis
        self.tmin_var = tk.DoubleVar(value=0.0)
        self.tmax_var = tk.DoubleVar(value=0.0)

        # Expected Tm range (only used for auto-trim)
        self.tm_win_min_var = tk.StringVar(value="")
        self.tm_win_max_var = tk.StringVar(value="")

        # suspected / corrected / DELETED tracking
        self.suspected_wells = []
        self.corrected_wells = []
        self.deleted_wells = set()          # pocillos eliminados lógicamente
        self.auto_trimmed_wells = set()     # pozos que han sufrido auto-trim
        self.auto_suspect_index = {}

        # undo/redo per well (fluorescencia + trimming + flag auto-trim)
        self.history = {}
        self.redo_history = {}
        self.trim_history = {}
        self.trim_redo = {}
        self.auto_trim_history = {}
        self.auto_trim_redo = {}

        # Cache para cálculos de Tm
        self._tm_cache = {}
        self._last_tm_params = None
        self._cached_smooth_value = 25

        # UI
        self._build_ui()
        self._bind_shortcuts()

        # load file
        if path is not None and os.path.exists(path):
            self._load_gdsf(path)
        else:
            self._ask_and_load()

    # ------------------------ UI ------------------------
    def _build_ui(self):
        # ===== Toolbar row 1 =====
        bar1 = ttk.Frame(self)
        bar1.pack(side=tk.TOP, fill=tk.X, padx=8, pady=(8,4))
        ttk.Button(bar1, text="Open .gdsf", command=self._ask_and_load).pack(side=tk.LEFT)
        ttk.Button(bar1, text="Export corrected", command=self._export_corrected).pack(side=tk.LEFT, padx=(6,0))
        ttk.Button(bar1, text="Export corrected + smoothed", command=self._export_corrected_smoothed).pack(side=tk.LEFT, padx=(6,0))
        ttk.Button(bar1, text="Export Tm table", command=self._export_tm_table).pack(side=tk.LEFT, padx=(6,0))
        ttk.Separator(bar1, orient="vertical").pack(side=tk.LEFT, fill=tk.Y, padx=8)
        self.undo_btn = ttk.Button(bar1, text="Undo (Ctrl/Cmd+Z)", command=self._undo_current_well, state="disabled")
        self.undo_btn.pack(side=tk.LEFT, padx=(0,6))
        self.redo_btn = ttk.Button(bar1, text="Redo (Ctrl+Y / Ctrl+Shift+Z)", command=self._redo_current_well, state="disabled")
        self.redo_btn.pack(side=tk.LEFT, padx=(0,6))
        ttk.Checkbutton(bar1, text="Animate correction", variable=self.animate_var).pack(side=tk.LEFT, padx=(8,0))

        # ===== Toolbar row 2 =====
        bar2 = ttk.Frame(self)
        bar2.pack(side=tk.TOP, fill=tk.X, padx=8, pady=(0,4))
        ttk.Label(bar2, text="Step correction:").pack(side=tk.LEFT, padx=(0,12))
        ttk.Label(bar2, text="Abs threshold:").pack(side=tk.LEFT)
        self.abs_thr_entry = self._create_validated_entry(
            bar2, self.abs_thr_var,
            lambda p, s: self._validate_float(p, s, 0, 1000),
            width=7
        )
        self.abs_thr_entry.pack(side=tk.LEFT, padx=(4,12))
        ttk.Label(bar2, text="k·disp:").pack(side=tk.LEFT)
        self.kdisp_entry = self._create_validated_entry(
            bar2, self.kdisp_var,
            lambda p, s: self._validate_float(p, s, 0, 100),
            width=5
        )
        self.kdisp_entry.pack(side=tk.LEFT, padx=(4,12))
        ttk.Label(bar2, text="Dispersion:").pack(side=tk.LEFT)
        self.disp_combo = ttk.Combobox(bar2, width=6, state="readonly", values=["MAD","STD"], textvariable=self.disp_method)
        self.disp_combo.pack(side=tk.LEFT, padx=(4,12))
        self.iter_chk = ttk.Checkbutton(bar2, text="Iterative", variable=self.iterative_var)
        self.iter_chk.pack(side=tk.LEFT, padx=(0,12))
        ttk.Checkbutton(bar2, text="Show derivative & Tm", variable=self.show_deriv_var, command=self._draw_current).pack(side=tk.LEFT)
        ttk.Separator(bar2, orient="vertical").pack(side=tk.LEFT, fill=tk.Y, padx=8)
        self.scan_btn = ttk.Button(bar2, text="Scan suspects", command=self._scan_suspects, state="disabled")
        self.scan_btn.pack(side=tk.LEFT)
        self.correct_all_btn = ttk.Button(bar2, text="Correct all suspects", command=self._correct_all_suspects, state="disabled")
        self.correct_all_btn.pack(side=tk.LEFT, padx=(6,0))
        ttk.Checkbutton(bar2, text="Multi-jump", variable=self.multi_jump_var).pack(side=tk.LEFT, padx=(12,0))

        # ===== Toolbar row 3 (Review mode) =====
        bar3 = ttk.Frame(self)
        bar3.pack(side=tk.TOP, fill=tk.X, padx=8, pady=(0,4))
        ttk.Label(bar3, text="Review mode:").pack(side=tk.LEFT)
        ttk.Button(bar3, text="Prev Suspected  (Shift+S)", command=lambda: self._jump_review(list_name="suspected", step=-1)).pack(side=tk.LEFT, padx=(6,4))
        ttk.Button(bar3, text="Next Suspected  (s)", command=lambda: self._jump_review(list_name="suspected", step=+1)).pack(side=tk.LEFT, padx=(0,10))
        ttk.Button(bar3, text="Prev Corrected  (Shift+C)", command=lambda: self._jump_review(list_name="corrected", step=-1)).pack(side=tk.LEFT, padx=(0,4))
        ttk.Button(bar3, text="Next Corrected  (c)", command=lambda: self._jump_review(list_name="corrected", step=+1)).pack(side=tk.LEFT)
        ttk.Button(
            bar3,
            text="Prev Tm outlier",
            command=lambda: self._jump_tm_outlier(step=-1)
        ).pack(side=tk.LEFT, padx=(12,4))

        ttk.Button(
            bar3,
            text="Next Tm outlier",
            command=lambda: self._jump_tm_outlier(step=+1)
        ).pack(side=tk.LEFT, padx=(0,4))

        # ===== Toolbar row 4 (Smoothing) =====
        bar4 = ttk.Frame(self)
        bar4.pack(side=tk.TOP, fill=tk.X, padx=8, pady=(0,6))
        ttk.Label(bar4, text="Smoothing:").pack(side=tk.LEFT)
        ttk.Checkbutton(bar4, text="Smooth on/off", variable=self.smooth_on_var, command=self._on_smooth_change).pack(side=tk.LEFT, padx=(6,8))
        ttk.Label(bar4, text="Strength:").pack(side=tk.LEFT)
        self.smooth_slider = ttk.Scale(bar4, from_=0, to=100, value=self.smooth_strength_var.get(), orient=tk.HORIZONTAL,
                                       command=lambda _=None: self._on_smooth_change())
        self.smooth_slider.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(6,8))
        self.smooth_entry = ttk.Entry(bar4, width=4)
        self.smooth_entry.insert(0, str(self.smooth_strength_var.get()))
        self.smooth_entry.bind("<Return>", lambda e: self._set_smooth_from_entry())
        self.smooth_entry.pack(side=tk.LEFT)

        # ===== Main area =====
        main = ttk.Frame(self)
        main.pack(side=tk.TOP, fill=tk.BOTH, expand=True, padx=8, pady=(0,8))

        # Left column: All wells + Tm outliers + shortcuts
        left = ttk.Frame(main, width=320)
        left.pack(side=tk.LEFT, fill=tk.Y)

        all_frame = ttk.LabelFrame(left, text="All wells (with data)")
        all_frame.pack(fill=tk.BOTH, expand=True)
        self.well_list = tk.Listbox(all_frame, exportselection=False)
        self.well_list.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)
        self.well_list.bind("<<ListboxSelect>>", self._on_select_well)
        self.well_list.bind("<Double-1>", self._goto_selected_well)

        # Frame para Tm outliers con contador dinámico
        self.tm_outlier_frame = ttk.LabelFrame(left, text="Tm outliers (0)")
        self.tm_outlier_frame.pack(fill=tk.BOTH, expand=True, pady=(6,0))
        self.tm_outlier_list = tk.Listbox(self.tm_outlier_frame, exportselection=False)
        self.tm_outlier_list.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)
        self.tm_outlier_list.bind("<<ListboxSelect>>", self._on_select_tm_outlier)
        self.tm_outlier_list.bind("<Double-1>", self._goto_tm_outlier)

        # Umbral para considerar una Tm como outlier (en ºC)
        tm_thr_row = ttk.Frame(left)
        tm_thr_row.pack(fill=tk.X, padx=4, pady=(2,0))
        ttk.Label(tm_thr_row, text="Tm outlier threshold (°C):").pack(side=tk.LEFT)
        self.tm_thr_entry = self._create_validated_entry(
            tm_thr_row, self.tm_thr_var,
            lambda p, s: self._validate_float(p, s, 0.1, 100),
            width=6
        )
        self.tm_thr_entry.pack(side=tk.LEFT, padx=(4,0))
        self.tm_thr_entry.bind("<Return>", lambda e: self._on_tm_thr_change())

        # Tm de referencia conocida por el usuario
        tm_ref_row = ttk.Frame(left)
        tm_ref_row.pack(fill=tk.X, padx=4, pady=(2,0))
        ttk.Label(tm_ref_row, text="I know my Tm (°C):").pack(side=tk.LEFT)
        self.tm_ref_entry = self._create_validated_entry(
            tm_ref_row, self.tm_ref_var,
            lambda p, s: self._validate_float(p, s, -50, 150),
            width=6
        )
        self.tm_ref_entry.pack(side=tk.LEFT, padx=(4,0))
        self.tm_ref_entry.bind("<Return>", lambda e: self._on_tm_thr_change())

        # Label con la media actual de las Tm calculadas
        self.tm_mean_label = ttk.Label(left, textvariable=self.tm_mean_var, justify="left")
        self.tm_mean_label.pack(anchor="w", padx=8, pady=(2,0))

        # --- Small shortcuts panel at the bottom of the left column ---
        shortcuts_frame = ttk.LabelFrame(left, text="Shortcuts")
        shortcuts_frame.pack(side=tk.BOTTOM, fill=tk.X, padx=0, pady=(6,0))
        small_font = tkfont.Font(size=8)
        shortcuts_text = (
            "Undo: Ctrl/Cmd+Z   •   Redo: Ctrl+Y / Ctrl+Shift+Z\n"
            "Next suspected: s   •   Prev suspected: Shift+S\n"
            "Next corrected: c   •   Prev corrected: Shift+C\n"
            "Move breakpoint: ← →   •   Go: Enter"
        )
        lbl_short = ttk.Label(shortcuts_frame, text=shortcuts_text, justify="left")
        lbl_short.pack(anchor="w", padx=6, pady=4)
        try:
            lbl_short.configure(font=small_font)
        except tk.TclError:
            pass

        # Middle column: corrected + suspected + delete/recover buttons
        mid = ttk.Frame(main, width=260)
        mid.pack(side=tk.LEFT, fill=tk.Y, padx=(8,8))
        corr_frame = ttk.LabelFrame(mid, text="Corrected")
        corr_frame.pack(fill=tk.BOTH, expand=True)
        self.corrected_list = tk.Listbox(corr_frame, exportselection=False)
        self.corrected_list.pack(fill=tk.BOTH, expand=True, padx=6, pady=6)
        self.corrected_list.bind("<<ListboxSelect>>", self._on_select_corrected)

        susp_frame = ttk.LabelFrame(mid, text="Suspected")
        susp_frame.pack(fill=tk.BOTH, expand=True, pady=(8,0))
        self.suspected_list = tk.Listbox(susp_frame, exportselection=False)
        self.suspected_list.pack(fill=tk.BOTH, expand=True, padx=6, pady=6)
        self.suspected_list.bind("<<ListboxSelect>>", self._on_select_suspected)

        # Botones Delete/Recover debajo de Suspected
        btn_frame = ttk.Frame(mid)
        btn_frame.pack(side=tk.TOP, pady=(6,0))
        self.delete_well_btn = ttk.Button(btn_frame, text="Delete Well", command=self._delete_selected_well)
        self.delete_well_btn.pack(side=tk.LEFT)
        self.recover_well_btn = ttk.Button(btn_frame, text="Recover Well", command=self._recover_selected_well)
        self.recover_well_btn.pack(side=tk.LEFT, padx=(6,0))

        # right: plots + per-well controls
        right = ttk.Frame(main)
        right.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

        perwell = ttk.LabelFrame(right, text="Per-well controls")
        perwell.pack(fill=tk.X, padx=4, pady=(0,6))
        ttk.Label(perwell, text="Breakpoint index:").pack(side=tk.LEFT, padx=(8,4))
        self.idx_slider = ttk.Scale(perwell, from_=0, to=1, value=0, orient=tk.HORIZONTAL, command=self._on_slider)
        self.idx_slider.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0,8), pady=8)
        self.idx_left_btn = ttk.Button(perwell, text="◀", width=3, command=lambda: self._nudge_index(-1))
        self.idx_left_btn.pack(side=tk.LEFT, padx=(0,2))
        self.idx_entry = ttk.Entry(perwell, width=8)
        self.idx_entry.pack(side=tk.LEFT, padx=(0,2))
        self.idx_entry.bind("<Return>", lambda e: self._set_index_from_entry())
        self.idx_right_btn = ttk.Button(perwell, text="▶", width=3, command=lambda: self._nudge_index(+1))
        self.idx_right_btn.pack(side=tk.LEFT, padx=(2,4))
        ttk.Button(perwell, text="Go", command=self._set_index_from_entry).pack(side=tk.LEFT)
        self.idx_label = ttk.Label(perwell, text="Index: –/–   |   Temp: –   |   Tm: –")
        self.idx_label.pack(side=tk.LEFT, padx=8)

        # Operation & Apply correction
        perwell2 = ttk.Frame(right)
        perwell2.pack(fill=tk.X, padx=4, pady=(0,8))
        ttk.Label(perwell2, text="Operation:").pack(side=tk.LEFT)
        ttk.Radiobutton(perwell2, text="Auto", variable=self.op_var, value="auto").pack(side=tk.LEFT, padx=(6,0))
        ttk.Radiobutton(perwell2, text="Add", variable=self.op_var, value="add").pack(side=tk.LEFT, padx=(6,0))
        ttk.Radiobutton(perwell2, text="Subtract", variable=self.op_var, value="sub").pack(side=tk.LEFT, padx=(6,12))
        self.apply_btn = ttk.Button(perwell2, text="Apply correction at breakpoint", command=self._apply_correction, state="disabled")
        self.apply_btn.pack(side=tk.LEFT)

        # plots (main + derivative)
        self.fig = Figure(figsize=(6,4), dpi=100)
        gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.25)
        self.ax = self.fig.add_subplot(gs[0])
        self.axd = self.fig.add_subplot(gs[1], sharex=self.ax)

        self.canvas = FigureCanvasTkAgg(self.fig, master=right)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.canvas.mpl_connect("button_press_event", self._on_plot_click)

        # ===== T-range & Expected T controls (below derivative plot) =====
        range_frame = ttk.LabelFrame(right, text="Temperature ranges")
        range_frame.pack(side=tk.TOP, fill=tk.X, padx=4, pady=(4,0))

        # Row 1: analysis T range + Remove button
        r1 = ttk.Frame(range_frame)
        r1.pack(side=tk.TOP, fill=tk.X, pady=(2,0))
        ttk.Label(r1, text="Analysis T range:").pack(side=tk.LEFT)
        self.tmin_scale = ttk.Scale(
            r1, from_=0, to=100, orient=tk.HORIZONTAL,
            command=lambda _=None: self._on_tmin_slider()
        )
        self.tmin_scale.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(6,2))
        self.tmin_entry = ttk.Entry(r1, width=6)
        self.tmin_entry.pack(side=tk.LEFT, padx=(0,4))
        self.tmin_entry.bind("<Return>", lambda e: self._set_tmin_from_entry())

        self.tmax_scale = ttk.Scale(
            r1, from_=0, to=100, orient=tk.HORIZONTAL,
            command=lambda _=None: self._on_tmax_slider()
        )
        self.tmax_scale.pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(4,2))
        self.tmax_entry = ttk.Entry(r1, width=6)
        self.tmax_entry.pack(side=tk.LEFT, padx=(0,8))
        self.tmax_entry.bind("<Return>", lambda e: self._set_tmax_from_entry())

        self.remove_data_btn = ttk.Button(r1, text="Remove data", command=self._remove_data_outside_range)
        self.remove_data_btn.pack(side=tk.LEFT, padx=(4,0))

        # Row 2: expected Tm range (used only for auto-trim)
        r2 = ttk.Frame(range_frame)
        r2.pack(side=tk.TOP, fill=tk.X, pady=(4,2))
        ttk.Label(r2, text="Expected Tm range (°C):").pack(side=tk.LEFT)
        self.tm_min_entry = ttk.Entry(r2, width=6, textvariable=self.tm_win_min_var)
        self.tm_min_entry.pack(side=tk.LEFT, padx=(6,2))
        ttk.Label(r2, text="to").pack(side=tk.LEFT)
        self.tm_max_entry = ttk.Entry(r2, width=6, textvariable=self.tm_win_max_var)
        self.tm_max_entry.pack(side=tk.LEFT, padx=(2,4))
        # No recalcula Tm, sólo redibuja por si quieres feedback visual
        self.tm_min_entry.bind("<KeyRelease>", lambda e: self._on_tm_range_change())
        self.tm_max_entry.bind("<KeyRelease>", lambda e: self._on_tm_range_change())
        ttk.Button(
            r2,
            text="Auto-trim to expected range",
            command=self._auto_trim_to_expected_range
        ).pack(side=tk.LEFT, padx=(8,0))

        # status bar
        status_bar = ttk.Frame(self)
        status_bar.pack(side=tk.BOTTOM, fill=tk.X)
        ttk.Label(status_bar, textvariable=self.status_var, anchor="w").pack(side=tk.LEFT, fill=tk.X, expand=True, padx=8, pady=4)

        # light theme defaults
        self._apply_light_theme()

    # ------------------------ VALIDACIÓN DE ENTRADAS ------------------------
    def _create_validated_entry(self, parent, variable, command, width=6):
        """Crea entry con validación numérica."""
        vcmd = (self.register(command), '%P', '%s')
        entry = ttk.Entry(parent, width=width, textvariable=variable,
                          validate="key", validatecommand=vcmd)
        return entry

    def _validate_float(self, new_value, old_value, min_val=None, max_val=None):
        """Validador genérico para números flotantes."""
        if new_value == "" or new_value == "-":
            return True
        try:
            val = float(new_value)
            if min_val is not None and val < min_val:
                return False
            if max_val is not None and val > max_val:
                return False
            return True
        except ValueError:
            return False

    # ------------------------ light theme (default) ------------------------
    def _apply_light_theme(self):
        style = ttk.Style()
        try:
            style.theme_use("clam")
        except Exception:
            pass
        style.configure(".", background="#f7f7f7", foreground="#111111", fieldbackground="#ffffff")
        style.configure("TFrame", background="#f7f7f7")
        style.configure("TLabelframe", background="#f7f7f7", foreground="#111111")
        for lb in [getattr(self, "well_list", None),
                   getattr(self, "corrected_list", None),
                   getattr(self, "suspected_list", None),
                   getattr(self, "tm_outlier_list", None)]:
            if lb:
                lb.config(bg="#ffffff", fg="#111111", highlightbackground="#f7f7f7",
                          selectbackground="#c7d2fe", selectforeground="#000000")
        # matplotlib light
        face, grid, txt = "#ffffff", "#cccccc", "#111111"
        for ax in [getattr(self, "ax", None), getattr(self, "axd", None)]:
            if ax:
                ax.set_facecolor(face)
                ax.tick_params(colors=txt)
                ax.xaxis.label.set_color(txt)
                ax.yaxis.label.set_color(txt)
                ax.grid(True, color=grid, alpha=0.35)
        if hasattr(self, "fig"):
            self.fig.patch.set_facecolor(face)

    # ------------------------ bindings ------------------------
    def _bind_shortcuts(self):
        # Undo/Redo
        self.bind_all("<Control-z>", lambda e: self._undo_current_well())
        self.bind_all("<Command-z>", lambda e: self._undo_current_well())
        self.bind_all("<Control-Shift-Z>", lambda e: self._redo_current_well())
        self.bind_all("<Control-y>", lambda e: self._redo_current_well())
        self.bind_all("<Command-Shift-Z>", lambda e: self._redo_current_well())

        # arrows for index
        self.bind_all("<Left>", lambda e: self._nudge_index(-1))
        self.bind_all("<Right>", lambda e: self._nudge_index(+1))

        # 🔹 arrows for wells navigation: solo cuando la lista tiene el foco
        if hasattr(self, "well_list"):
            def _on_up(e):
                self._move_well_selection(-1)
                return "break"  # evita que el Listbox mueva otra vez

            def _on_down(e):
                self._move_well_selection(+1)
                return "break"

            self.well_list.bind("<Up>", _on_up)
            self.well_list.bind("<Down>", _on_down)

        # review mode
        self.bind_all("s", lambda e: self._jump_review("suspected", +1))
        self.bind_all("S", lambda e: self._jump_review("suspected", -1))
        self.bind_all("c", lambda e: self._jump_review("corrected", +1))
        self.bind_all("C", lambda e: self._jump_review("corrected", -1))

    # ------------------------ helpers ------------------------
    @staticmethod
    def _well_sortkey(w):
        w = (w or "").strip().upper()
        if not w:
            return ("Z", 999)
        row = w[0]
        try:
            col = int("".join(ch for ch in w[1:] if ch.isdigit()))
        except Exception:
            col = 999
        return (row, col)

    def _clamp_index(self, i):
        if self.current_well is None:
            return 0
        i = max(0, i)
        i = min(i, self.max_idx)
        return i

    def _get_thresholds(self):
        try:
            abs_thr = float(self.abs_thr_var.get())
        except ValueError:
            abs_thr = 0.0
        try:
            k = float(self.kdisp_var.get())
        except ValueError:
            k = 0.0
        return max(0.0, abs_thr), max(0.0, k), self.disp_method.get()

    def _update_undo_redo_state(self):
        w = self.current_well
        hist = self.history.get(w, []) if w else []
        redo = self.redo_history.get(w, []) if w else []
        self.undo_btn.config(state=("normal" if hist else "disabled"))
        self.redo_btn.config(state=("normal" if redo else "disabled"))

    # ------------------------ smoothing helpers ------------------------
    def _get_smooth_strength(self):
        try:
            return int(self.smooth_strength_var.get())
        except Exception:
            return 0

    def _get_smoothing_for_derivative(self):
        """Smoothing para derivada: siempre base=25 + extra si smooth_on activado."""
        SMOOTH_BASE = 25
        if not self.smooth_on_var.get():
            return SMOOTH_BASE
        try:
            user = int(self.smooth_strength_var.get())
            return max(SMOOTH_BASE, user)
        except:
            return SMOOTH_BASE

    def _maybe_smooth(self, y):
        if self.smooth_on_var.get():
            return smooth_signal(y, self._get_smooth_strength())
        return np.asarray(y, dtype=float)

    def _on_smooth_change(self):
        try:
            val = float(self.smooth_slider.get())
        except Exception:
            val = 0
        val = max(0, min(100, int(val)))
        self.smooth_strength_var.set(val)
        self.smooth_entry.delete(0, tk.END)
        self.smooth_entry.insert(0, str(val))

        # Invalidar cache de Tm solo si smoothing afecta derivada
        old_smooth = self._cached_smooth_value
        new_smooth = self._get_smoothing_for_derivative()
        if old_smooth != new_smooth:
            self._tm_cache.clear()
            self._cached_smooth_value = new_smooth
            self._recompute_tm_all_wells()
            self._refresh_all_wells_with_tm() # refrescar también las Tm mostradas en "All wells"

        self._draw_current()

    def _set_smooth_from_entry(self):
        try:
            val = int(self.smooth_entry.get())
        except Exception:
            val = 0
        val = max(0, min(100, val))
        self.smooth_strength_var.set(val)
        self.smooth_slider.set(val)
        self._on_smooth_change()

    # ------------------------ Expected Tm range behaviour ------------------------
    def _on_tm_range_change(self):
        """
        Se llama cuando el usuario cambia el rango esperado de Tm.
        IMPORTANTE: esto ya NO afecta al cálculo de Tm.
        El rango sólo se usa en el botón 'Auto-trim to expected range'.
        """
        # Sólo redibujamos por si quieres ver algo actualizado
        self._draw_current()

    # ------------------------ analysis T-range helpers ------------------------
    def _init_t_range_for_well(self):
        """Initialize analysis T range sliders for current well."""
        if self.current_well is None:
            return
        g = self.per_well_work[self.current_well]
        if g.empty:
            return
        tmin_data = float(g["Temperature"].min())
        tmax_data = float(g["Temperature"].max())
        if not np.isfinite(tmin_data) or not np.isfinite(tmax_data) or tmin_data >= tmax_data:
            return

        if self.current_well in self.trim_ranges:
            tmin, tmax = self.trim_ranges[self.current_well]
        else:
            tmin, tmax = tmin_data, tmax_data

        tmin = max(tmin_data, min(tmin, tmax_data))
        tmax = max(tmin_data, min(tmax, tmax_data))
        if tmin > tmax:
            tmin, tmax = tmax, tmin

        self.tmin_var.set(tmin)
        self.tmax_var.set(tmax)
        self.tmin_scale.configure(from_=tmin_data, to=tmax_data)
        self.tmax_scale.configure(from_=tmin_data, to=tmax_data)
        self.tmin_scale.set(tmin)
        self.tmax_scale.set(tmax)
        self.tmin_entry.delete(0, tk.END)
        self.tmin_entry.insert(0, f"{tmin:.2f}")
        self.tmax_entry.delete(0, tk.END)
        self.tmax_entry.insert(0, f"{tmax:.2f}")

    def _get_current_t_range(self):
        """Return (tmin, tmax) from sliders/entries, clamped to data range."""
        if self.current_well is None:
            return (None, None)
        g = self.per_well_work[self.current_well]
        if g.empty:
            return (None, None)
        tmin_data = float(g["Temperature"].min())
        tmax_data = float(g["Temperature"].max())
        try:
            tmin = float(self.tmin_var.get())
        except Exception:
            tmin = tmin_data
        try:
            tmax = float(self.tmax_var.get())
        except Exception:
            tmax = tmax_data
        tmin = max(tmin_data, min(tmin, tmax_data))
        tmax = max(tmin_data, min(tmax, tmax_data))
        if tmin > tmax:
            tmin, tmax = tmax, tmin
        return (tmin, tmax)

    def _on_tmin_slider(self):
        if self.current_well is None:
            return
        val = float(self.tmin_scale.get())
        self.tmin_var.set(val)
        self.tmin_entry.delete(0, tk.END)
        self.tmin_entry.insert(0, f"{val:.2f}")
        self._draw_current()

    def _on_tmax_slider(self):
        if self.current_well is None:
            return
        val = float(self.tmax_scale.get())
        self.tmax_var.set(val)
        self.tmax_entry.delete(0, tk.END)
        self.tmax_entry.insert(0, f"{val:.2f}")
        self._draw_current()

    def _set_tmin_from_entry(self):
        if self.current_well is None:
            return
        g = self.per_well_work[self.current_well]
        if g.empty:
            return
        try:
            val = float(self.tmin_entry.get())
        except Exception:
            return
        tmin_data = float(g["Temperature"].min())
        tmax_data = float(g["Temperature"].max())
        val = max(tmin_data, min(val, tmax_data))
        self.tmin_var.set(val)
        self.tmin_scale.configure(from_=tmin_data, to=tmax_data)
        self.tmin_scale.set(val)
        self.tmin_entry.delete(0, tk.END)
        self.tmin_entry.insert(0, f"{val:.2f}")
        self._draw_current()

    def _set_tmax_from_entry(self):
        if self.current_well is None:
            return
        g = self.per_well_work[self.current_well]
        if g.empty:
            return
        try:
            val = float(self.tmax_entry.get())
        except Exception:
            return
        tmin_data = float(g["Temperature"].min())
        tmax_data = float(g["Temperature"].max())
        val = max(tmin_data, min(val, tmax_data))
        self.tmax_var.set(val)
        self.tmax_scale.configure(from_=tmin_data, to=tmax_data)
        self.tmax_scale.set(val)
        self.tmax_entry.delete(0, tk.END)
        self.tmax_entry.insert(0, f"{val:.2f}")
        self._draw_current()

    def _remove_data_outside_range(self):
        """Store analysis T range [tmin, tmax] for current well (reversible via undo)."""
        if self.current_well is None:
            return
        g = self.per_well_work[self.current_well].copy()
        if g.empty:
            return
        tmin, tmax = self._get_current_t_range()
        if tmin is None or tmax is None or tmin >= tmax:
            messagebox.showinfo("Invalid range", "Set a valid analysis T range first.")
            return
        mask = (g["Temperature"] >= tmin) & (g["Temperature"] <= tmax)
        if mask.sum() < 3:
            messagebox.showinfo("Too few points", "Range would leave fewer than 3 points.")
            return

        # snapshot estado antes de cambiar el trimming (para undo/redo)
        self._push_history(self.current_well)
        self.trim_ranges[self.current_well] = (tmin, tmax)

        # Invalidar cache de Tm para este pozo
        self._tm_cache.pop(self.current_well, None)
        self._recompute_tm_all_wells()
        self._refresh_all_wells_with_tm()
        self._draw_current()
        self.status_var.set(f"{self.current_well}: analysis restricted to [{tmin:.2f}, {tmax:.2f}] °C (reversible).")

    # ------------------------ trimming helpers ------------------------
    def _get_visible_df(self, well):
        """Return dataframe for well after applying any trim range."""
        if well in self.deleted_wells:
            # Si está eliminado, devolver DataFrame vacío para que no se exporte ni se calcule
            return pd.DataFrame(columns=["Well","Temperature","Fluorescence"])

        g = self.per_well_work.get(well)
        if g is None:
            return None
        if well not in self.trim_ranges:
            return g
        tmin, tmax = self.trim_ranges[well]
        mask = (g["Temperature"] >= tmin) & (g["Temperature"] <= tmax)
        g2 = g.loc[mask].reset_index(drop=True)
        if len(g2) < 3:
            # Fall back to full data if overly aggressive trim
            return g
        return g2

    # ------------------------ file I/O ------------------------
    def _ask_and_load(self):
        fpath = filedialog.askopenfilename(
            title="Select a .gdsf file",
            filetypes=[("GDSF text", "*.gdsf"), ("TSV", "*.tsv"), ("Text", "*.txt"), ("All files", "*.*")]
        )
        if fpath:
            self._load_gdsf(fpath)

    def _load_gdsf(self, path):
        try:
            df = pd.read_csv(path, sep="\t", header=None, names=["Well","Temperature","Fluorescence"])
        except Exception as e:
            messagebox.showerror("Read error", "Could not read file:\n" + str(e))
            return
        for col in ["Well","Temperature","Fluorescence"]:
            if col not in df.columns:
                messagebox.showerror("Invalid format", "The .gdsf must have 3 columns: Well, Temperature, Fluorescence")
                return
        df["Well"] = df["Well"].astype(str).str.upper().str.strip()
        df["Temperature"] = pd.to_numeric(df["Temperature"], errors="coerce")
        df["Fluorescence"] = pd.to_numeric(df["Fluorescence"], errors="coerce")
        df = df.dropna(subset=["Temperature","Fluorescence"])

        self.df_orig = df.copy()
        self.df_work = df.copy()

        # Elegimos wells que tengan alguna fluorescencia distinta de 0
        wells_sorted = sorted(
            [w for w, g in df.groupby("Well") if (g["Fluorescence"] != 0).any()],
            key=self._well_sortkey
        )

        self.per_well_orig = {
            w: self.df_orig[self.df_orig["Well"]==w].copy().sort_values("Temperature").reset_index(drop=True)
            for w in wells_sorted
        }
        self.per_well_work = {
            w: self.df_work[self.df_work["Well"]==w].copy().sort_values("Temperature").reset_index(drop=True)
            for w in wells_sorted
        }

        self.wells = wells_sorted
        self.suspected_wells = []
        self.corrected_wells = []
        self.deleted_wells = set()
        self.auto_suspect_index = {}
        self.history = {w: [] for w in self.wells}
        self.redo_history = {w: [] for w in self.wells}
        self.trim_history = {w: [] for w in self.wells}
        self.trim_redo = {w: [] for w in self.wells}
        self.auto_trim_history = {w: [] for w in self.wells}
        self.auto_trim_redo = {w: [] for w in self.wells}
        self.trim_ranges = {}
        self.tm_values = {}
        self.tm_outlier_wells = []
        self._tm_cache = {}
        self.auto_trimmed_wells = set()
        self._last_tm_params = None
        self._cached_smooth_value = 25

        self._populate_lists()
        self.status_var.set(f"Loaded: {os.path.basename(path)} | Wells with data: {len(self.wells)}")

        if self.wells:
            self.well_list.selection_clear(0, tk.END)
            self.well_list.selection_set(0)
            self._on_select_well()

        self.scan_btn.config(state="normal")
        self.correct_all_btn.config(state="normal")
        self._update_undo_redo_state()

    def _build_export_df_corrected(self):
        """Build corrected dataframe, respecting trims and deleted wells."""
        frames = []
        for w in self.wells:
            g = self._get_visible_df(w)
            if len(g) > 0:  # Solo exportar si no está vacío (no eliminado)
                frames.append(g[["Well","Temperature","Fluorescence"]])
        out_df = pd.concat(frames, ignore_index=True)[["Well","Temperature","Fluorescence"]]
        return out_df

    def _export_corrected(self):
        if self.df_work is None:
            messagebox.showinfo("Nothing to save", "Load a .gdsf first.")
            return
        out_df = self._build_export_df_corrected()
        fpath = filedialog.asksaveasfilename(
            title="Export corrected .gdsf",
            defaultextension=".gdsf",
            filetypes=[("GDSF text","*.gdsf"), ("TSV","*.tsv"), ("Text","*.txt"), ("All files","*.*")]
        )
        if not fpath:
            return
        try:
            out_df.to_csv(fpath, sep="\t", header=False, index=False, float_format="%.10g")
            self.status_var.set(f"Saved corrected: {os.path.basename(fpath)}")
        except Exception as e:
            messagebox.showerror("Save error", "Could not save:\n" + str(e))

    def _export_corrected_smoothed(self):
        if self.df_work is None:
            messagebox.showinfo("Nothing to save", "Load a .gdsf first.")
            return
        strength = self._get_smooth_strength()
        frames = []
        for w in self.wells:
            g = self._get_visible_df(w)
            if len(g) > 0:  # Solo exportar si no está vacío
                g = g.copy().sort_values("Temperature").reset_index(drop=True)
                y = g["Fluorescence"].values.astype(float)
                ys = smooth_signal(y, strength)
                g["Fluorescence"] = ys
                frames.append(g[["Well","Temperature","Fluorescence"]])
        out_df = pd.concat(frames, ignore_index=True)
        fpath = filedialog.asksaveasfilename(
            title="Export corrected+smoothed .gdsf",
            defaultextension=".gdsf",
            filetypes=[("GDSF text","*.gdsf"), ("TSV","*.tsv"), ("Text","*.txt"), ("All files","*.*")]
        )
        if not fpath:
            return
        try:
            out_df.to_csv(fpath, sep="\t", header=False, index=False, float_format="%.10g")
            self.status_var.set(f"Saved corrected+smoothed: {os.path.basename(fpath)} (strength={strength})")
        except Exception as e:
            messagebox.showerror("Save error", "Could not save:\n" + str(e))

    # ------------------------ Tm helpers ------------------------
    def _get_tm_window(self):
        """Return (lo, hi) for expected Tm range, or None if not set/valid."""
        s_lo = (self.tm_win_min_var.get() or "").strip()
        s_hi = (self.tm_win_max_var.get() or "").strip()
        if not s_lo or not s_hi:
            return None
        try:
            lo = float(s_lo)
            hi = float(s_hi)
        except ValueError:
            return None
        if not np.isfinite(lo) or not np.isfinite(hi) or lo >= hi:
            return None
        return (lo, hi)

    def _get_tm_reference(self):
        """Obtiene Tm de referencia (usuario o media)."""
        ref_str = (self.tm_ref_var.get() or "").strip()
        if ref_str:
            try:
                return float(ref_str)
            except:
                pass
        valid_tms = [tm for tm in self.tm_values.values() if tm is not None and np.isfinite(tm)]
        return np.mean(valid_tms) if valid_tms else 0.0

    def _get_tm_threshold(self):
        """Obtiene umbral de outlier."""
        try:
            thr = float(self.tm_thr_var.get())
            return max(0.1, thr)
        except:
            return 20.0

    def _tm_from_xy(self, x, y, smooth=False, strength=None):
        """
        Versión simple: calcula Tm como máximo global de la derivada orientada hacia arriba.
        NO usa el Expected Tm range.
        """
        base_strength = self._get_smooth_strength() if strength is None else strength
        s = max(25, int(base_strength)) if smooth else max(25, 0)
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        order = np.argsort(x)
        x = x[order]; y = y[order]
        if np.any(np.diff(x) == 0):
            eps = np.linspace(0, 1e-9, num=len(x))
            x = x + eps
        ys = smooth_signal(y, s)
        d = np.gradient(ys, x)
        if d.size == 0 or not np.isfinite(d).any():
            return None

        max_abs = np.nanmax(np.abs(d))
        if not np.isfinite(max_abs) or max_abs <= 0:
            return None

        # Orient derivative so the dominant peak is upward
        i_dom = int(np.nanargmax(np.abs(d)))
        orient = 1.0 if d[i_dom] > 0 else -1.0
        d_oriented = d * orient

        i_tm = int(np.nanargmax(d_oriented))
        return float(x[i_tm]) if np.isfinite(d_oriented[i_tm]) else None

    def _compute_tm_for_xy(self, x, y):
        """Internal helper: compute Tm and normalized derivative for arbitrary x,y.
        IMPORTANTE: Tm siempre se calcula globalmente (no se restringe por Expected Tm range).
        """
        # Verificar si está en cache
        cache_key = f"{len(x)}_{hash(x.tobytes())}_{hash(y.tobytes())}"
        current_params = (self._get_smoothing_for_derivative(), None)  # ventana ya no se usa

        if cache_key in self._tm_cache and self._last_tm_params == current_params:
            return self._tm_cache[cache_key]

        # Parámetros de smoothing
        s = self._get_smoothing_for_derivative()

        # convertir datos
        x = np.asarray(x, dtype=float)
        y = np.asarray(y, dtype=float)
        if x.size < 3:
            result = (None, x, None)
            self._tm_cache[cache_key] = result
            return result

        # ordenar por temperatura
        order = np.argsort(x)
        x = x[order]
        y = y[order]

        # evitar duplicados exactos
        if np.any(np.diff(x) == 0):
            eps = np.linspace(0, 1e-9, num=len(x))
            x = x + eps

        # aplicar smoothing para derivada
        ys = smooth_signal(y, s)
        d = np.gradient(ys, x)
        if d.size == 0 or not np.isfinite(d).any():
            result = (None, x, None)
            self._tm_cache[cache_key] = result
            return result

        max_abs = np.nanmax(np.abs(d))
        if not np.isfinite(max_abs) or max_abs <= 0:
            result = (None, x, None)
            self._tm_cache[cache_key] = result
            return result

        # orientar derivada hacia arriba
        i_dom = int(np.nanargmax(np.abs(d)))
        orient = 1.0 if d[i_dom] > 0 else -1.0
        d_up = d * orient

        # Tm = máximo global de la derivada orientada
        i_tm = int(np.nanargmax(d_up))
        tm = float(x[i_tm]) if np.isfinite(d_up[i_tm]) else None

        # derivada normalizada para el plot
        maxpos = np.nanmax(d_up) if np.isfinite(d_up).any() else 0.0
        dplot = d_up / maxpos if maxpos > 0 else d_up

        result = (tm, x, dplot)
        self._tm_cache[cache_key] = result
        self._last_tm_params = current_params
        return result

    def _compute_tm(self, well):
        g = self._get_visible_df(well)
        if g is None or len(g) < 3:
            return (None, None, None)
        return self._compute_tm_for_xy(
            g["Temperature"].values,
            g["Fluorescence"].values
        )

    def _export_tm_table(self):
        if self.df_work is None:
            messagebox.showinfo("Nothing to export", "Load a .gdsf first.")
            return
        strength = self._get_smooth_strength()
        rows = []
        for w in self.wells:
            g = self._get_visible_df(w)
            if len(g) > 0:  # Solo incluir si no está eliminado
                g = g.copy().sort_values("Temperature").reset_index(drop=True)
                x = g["Temperature"].values
                y = g["Fluorescence"].values.astype(float)
                
                tm_raw, _, _ = self._compute_tm_for_xy(x, y)  # como ahora

                # si quisieras un Tm con smoothing más fuerte, algo así:
                tm_smooth, _, _ = self._compute_tm_for_xy(
                    x, smooth_signal(y, max(25, int(strength)))
                )

                rows.append({
                    "Well": w,
                    "Tm_corrected": tm_raw,
                    "Tm_smoothed": tm_smooth,
                    "Smooth_strength": max(25, int(strength))
                })
        out = pd.DataFrame(rows, columns=["Well","Tm_corrected","Tm_smoothed","Smooth_strength"])
        fpath = filedialog.asksaveasfilename(
            title="Export Tm table (.tsv)",
            defaultextension=".tsv",
            filetypes=[("TSV","*.tsv"), ("CSV","*.csv"), ("All files","*.*")]
        )
        if not fpath:
            return
        try:
            sep = "\t" if fpath.lower().endswith(".tsv") else ","
            out.to_csv(fpath, sep=sep, index=False, float_format="%.6g")
            self.status_var.set(f"Saved Tm table: {os.path.basename(fpath)}")
        except Exception as e:
            messagebox.showerror("Save error", "Could not save:\n" + str(e))

    def _on_tm_thr_change(self):
        """Se llama cuando el usuario cambia el umbral o la Tm de referencia."""
        self._recompute_tm_all_wells()
        self._refresh_all_wells_with_tm()

    # ------------------------ Tm cache & outliers ------------------------
    def _recompute_tm_all_wells(self):
        """Calcula Tm para cada pozo, IGNORANDO los eliminados."""
        self.tm_values = {}
        valid = []  # (well, tm) con Tm finita

        for w in self.wells:
            if w in self.deleted_wells:
                continue  # Ignorar pocillos eliminados
            g = self._get_visible_df(w)
            if g is None or len(g) < 3:
                tm = None
            else:
                tm = self._compute_tm(w)[0]
            self.tm_values[w] = tm
            if tm is not None and np.isfinite(tm):
                valid.append((w, tm))

        # reset listas
        self.tm_outlier_wells = []
        self.tm_sorted_wells = []

        if not valid:
            self.tm_mean_var.set("Current mean Tm: n/a")
            self._refresh_tm_outlier_list()
            return

        # media de Tm de todos los pozos válidos
        arr = np.asarray([tm for _, tm in valid], dtype=float)
        mean = float(np.mean(arr))
        self.tm_mean_var.set(f"Current mean Tm: {mean:.2f} °C (n={len(arr)})")

        # umbral y Tm de referencia
        thr = self._get_tm_threshold()
        ref_tm = self._get_tm_reference()

        # Calcular desviaciones y ordenar
        outliers_with_dev = []
        for w, tm in valid:
            dev = abs(tm - ref_tm)
            if dev >= thr:
                outliers_with_dev.append((w, tm, dev))

        # Ordenar por nombre de pocillo (A1, A2, A3... B1, B2, etc.)
        outliers_with_dev.sort(key=lambda x: self._well_sortkey(x[0]))
        self.tm_outlier_wells = [w for w, _, _ in outliers_with_dev]
        self.tm_sorted_wells = self.tm_outlier_wells.copy()

        # Actualizar lista con desviaciones
        self._refresh_tm_outlier_list(outliers_with_dev)

    def _format_well_label(self, w):
        tm = self.tm_values.get(w)
    
        if w in self.deleted_wells:
            base = f"{w} [DELETED]"
        else:
            if tm is None or not np.isfinite(tm):
                base = w
            else:
                base = f"{w} — Tm={tm:.2f} °C"
    
            # ✂ tijerita si el pozo está recortado por cualquier método
            if self._well_is_trimmed(w):
                base += " ✂"
    
        return base


    def _well_is_trimmed(self, w):
        """True si el pozo tiene trimming manual o auto-trim."""
        return (w in self.trim_ranges) or (w in self.auto_trimmed_wells)


    def _refresh_all_wells_with_tm(self):
        """Refresh All wells list text (with Tm) while preserving selection & scroll."""
        if not self.wells:
            self.well_list.delete(0, tk.END)
            return

        sel = self.well_list.curselection()
        sel_index = sel[0] if sel else None
        try:
            first_visible = self.well_list.nearest(0)
        except Exception:
            first_visible = 0

        self.well_list.delete(0, tk.END)
        for w in self.wells:
            self.well_list.insert(tk.END, self._format_well_label(w))

        if sel_index is not None and 0 <= sel_index < len(self.wells):
            self.well_list.selection_set(sel_index)
            self.well_list.see(sel_index)
        else:
            self.well_list.see(first_visible)

        self._paint_all_wells_list()

    def _refresh_tm_outlier_list(self, outliers_with_dev=None):
        """Rellena la lista 'Tm outliers' SOLO con los pozos marcados, mostrando desviación."""
        self.tm_outlier_list.delete(0, tk.END)
        count = len(self.tm_outlier_wells)
        self.tm_outlier_frame.config(text=f"Tm outliers ({count})" if count > 0 else "Tm outliers")

        if outliers_with_dev is None:
            outliers_with_dev = []
            ref_tm = self._get_tm_reference()
            for w in self.tm_sorted_wells:
                tm = self.tm_values.get(w)
                if tm is not None:
                    dev = abs(tm - ref_tm)
                    outliers_with_dev.append((w, tm, dev))

        for w, tm, dev in outliers_with_dev:
            if w not in self.deleted_wells:  # No mostrar outliers eliminados
                trim_mark = " ✂" if self._well_is_trimmed(w) else ""
                label = f"{w} — Tm={tm:.2f}°C (Δ={dev:.1f}){trim_mark}"
                self.tm_outlier_list.insert(tk.END, label)

    # ------------------------ lists & color mapping ------------------------
    def _populate_lists(self):
        self._recompute_tm_all_wells()
        self.well_list.delete(0, tk.END)
        for w in self.wells:
            self.well_list.insert(tk.END, self._format_well_label(w))
        self._refresh_corrected_list()
        self._refresh_suspected_list()
        self._paint_all_wells_list()

    def _refresh_corrected_list(self):
        self.corrected_list.delete(0, tk.END)
        uniq, seen = [], set()
        for w in self.corrected_wells:
            if w not in self.deleted_wells and w not in seen:  # No mostrar eliminados
                uniq.append(w)
                seen.add(w)

        # 🔹 Ordenar siempre por orden de pozo (A1, A2, ... B1...)
        self.corrected_wells = sorted(uniq, key=self._well_sortkey)

        for w in self.corrected_wells:
            label = f"{w} ✂" if self._well_is_trimmed(w) else w
            self.corrected_list.insert(tk.END, label)

        self._paint_all_wells_list()

    def _refresh_suspected_list(self):
        self.suspected_list.delete(0, tk.END)
        filtered = [
            w for w in self.suspected_wells
            if w not in set(self.corrected_wells) and w not in self.deleted_wells
        ]

        # 🔹 Ordenar también por orden de pozo
        self.suspected_wells = sorted(filtered, key=self._well_sortkey)

        for w in self.suspected_wells:
            self.suspected_list.insert(tk.END, w)
        self._paint_all_wells_list()

    def _paint_all_wells_list(self):
        """Color mapping en main well list (light theme)."""
        if not self.wells:
            return

        # Colores
        col_norm_bg = "#ffffff"
        col_norm_fg = "#111111"
        col_susp_bg = "#fff4cc"
        col_susp_fg = "#000000"
        col_corr_bg = "#e6f4ea"
        col_corr_fg = "#0b3d20"
        col_out_fg = "#b00000"
        col_deleted_bg = "#e0e0e0"  # Gris claro para eliminados
        col_deleted_fg = "#888888"

        susp = set(self.suspected_wells)
        corr = set(self.corrected_wells)
        out = set(self.tm_outlier_wells)
        deleted = self.deleted_wells

        self.well_list.update_idletasks()
        for i, w in enumerate(self.wells):
            bg = col_norm_bg
            fg = col_norm_fg
            if w in deleted:
                bg, fg = col_deleted_bg, col_deleted_fg
            elif w in corr:
                bg, fg = col_corr_bg, col_corr_fg
            elif w in susp:
                bg, fg = col_susp_bg, col_susp_fg
            if w in out and w not in deleted:
                fg = col_out_fg
            self.well_list.itemconfig(i, bg=bg, fg=fg)

    # ------------------------ suspects / auto ------------------------
    def _scan_suspects(self):
        abs_thr, k, method = self._get_thresholds()
        suspects, auto_idx = [], {}
        corrected_set = set(self.corrected_wells)
        for w in self.wells:
            if w in corrected_set or w in self.deleted_wells:  # No escanear eliminados
                continue
            g = self._get_visible_df(w)
            if g is None or len(g) < 2:
                continue
            y = g["Fluorescence"].values.astype(float)
            diffs = np.diff(y)
            if len(diffs) == 0:
                continue
            disp = robust_mad_sigma(diffs) if method == "MAD" else (np.std(diffs, ddof=1) if len(diffs) > 1 else 0.0)
            i_star = int(np.argmax(np.abs(diffs)))
            maxjump = float(np.abs(diffs[i_star]))
            cond_abs = maxjump > abs_thr if abs_thr > 0 else False
            cond_k = (disp > 0) and (maxjump > k * disp) if k > 0 else False
            if cond_abs or cond_k:
                suspects.append(w)
                auto_idx[w] = i_star
        self.suspected_wells = sorted(list(set(suspects) - corrected_set), key=self._well_sortkey)
        self.auto_suspect_index = auto_idx
        self._refresh_suspected_list()
        self.status_var.set(f"Suspected: {len(self.suspected_wells)} wells (method={method}, abs>{abs_thr}, k={k})")

    # ------------------------ multi-jump correction engine ------------------------
    def _find_step_indices(self, y, abs_thr, k, method):
        diffs = np.diff(y)
        if diffs.size == 0:
            return []
        disp = robust_mad_sigma(diffs) if method == "MAD" else (np.std(diffs, ddof=1) if len(diffs) > 1 else 0.0)
        thr_rel = (k * disp) if k > 0 and disp > 0 else -np.inf
        thr_abs = abs_thr if abs_thr > 0 else -np.inf
        thr = max(thr_rel, thr_abs)
        if not np.isfinite(thr) or thr <= 0:
            return []
        idx = np.where(np.abs(diffs) > thr)[0]
        return idx.tolist()

    def _apply_multi_jump(self, well, iterative=False):
        abs_thr, k, method = self._get_thresholds()
        g = self.per_well_work[well].copy()
        y0 = g["Fluorescence"].values.astype(float)
        changed = False
        max_loops = 20 if iterative else 1
        y = y0.copy()
        for _ in range(max_loops):
            idxs = self._find_step_indices(y, abs_thr, k, method)
            if not idxs:
                break
            changed = True
            diffs = np.diff(y)
            deltas = diffs[idxs]
            n = len(y)
            adjust = np.zeros(n, dtype=float)
            for i, dlt in zip(idxs, deltas):
                adjust[i+1] -= dlt
            adjust = np.cumsum(adjust)
            y = y + adjust
            if not iterative:
                break
        if changed:
            g["Fluorescence"] = y
            self.per_well_work[well] = g
        return changed

    # ------------------------ batch correct ------------------------
    def _correct_all_suspects(self):
        if not self.suspected_wells:
            self._scan_suspects()
            if not self.suspected_wells:
                messagebox.showinfo("Nothing to correct", "No suspected wells under current thresholds.")
                return

        iterative = self.iterative_var.get()
        op = self.op_var.get()
        use_multi = self.multi_jump_var.get()

        corrected_now = []
        for w in list(self.suspected_wells):
            if w in self.deleted_wells:  # Saltar eliminados
                continue
            if use_multi:
                self._push_history(w)
                changed = self._apply_multi_jump(w, iterative=iterative)
            else:
                changed = False
                while True:
                    g_vis = self._get_visible_df(w)
                    y = g_vis["Fluorescence"].values.astype(float)
                    diffs = np.diff(y)
                    if len(diffs) == 0:
                        break
                    abs_thr, k, method = self._get_thresholds()
                    disp = robust_mad_sigma(diffs) if method == "MAD" else (np.std(diffs, ddof=1) if len(diffs) > 1 else 0.0)
                    i_star = int(np.argmax(np.abs(diffs)))
                    maxjump = float(np.abs(diffs[i_star]))
                    cond_abs = maxjump > abs_thr if abs_thr > 0 else False
                    cond_k = (disp > 0) and (maxjump > k * disp) if k > 0 else False
                    if not (cond_abs or cond_k):
                        break

                    # Map index in visible data to full dataframe index
                    g_full = self.per_well_work[w].copy()
                    # visible mask
                    if w in self.trim_ranges:
                        tmin, tmax = self.trim_ranges[w]
                        mask = (g_full["Temperature"] >= tmin) & (g_full["Temperature"] <= tmax)
                        idxs = np.where(mask)[0]
                        if len(idxs) <= i_star:
                            break
                        full_idx = idxs[i_star]
                    else:
                        full_idx = i_star

                    self._push_history(w)
                    y_full = g_full["Fluorescence"].values.astype(float)
                    delta = float(y_full[full_idx+1] - y_full[full_idx])
                    adj = -delta if op in ("auto","sub") else +delta
                    y_full[full_idx+1:] = y_full[full_idx+1:] + adj
                    g_full["Fluorescence"] = y_full
                    self.per_well_work[w] = g_full
                    changed = True
                    if not iterative:
                        break
            if changed:
                corrected_now.append(w)

        for w in corrected_now:
            if w not in self.corrected_wells:
                self.corrected_wells.append(w)

        self._recompute_tm_all_wells()
        self._refresh_all_wells_with_tm()
        self._refresh_corrected_list()
        self._scan_suspects()
        self._update_undo_redo_state()
        self.status_var.set(f"Corrected {len(corrected_now)} wells. Remaining suspects: {len(self.suspected_wells)}")
        self._draw_current()

    # ------------------------ AUTO-TRIM TO EXPECTED RANGE ------------------------
    def _restore_main_focus(self):
        """Intenta devolver el foco y los eventos de teclado a la ventana principal."""
        try:
            # Forzar foco a la ventana principal
            self.focus_force()
        except Exception:
            try:
                self.focus_set()
            except Exception:
                pass
    
        # Dar foco a algo razonable para navegación con teclado (lista principal)
        try:
            self.well_list.focus_set()
        except Exception:
            pass
    
        # Asegura que Tk actualiza el estado de focus/eventos
        try:
            self.update_idletasks()
        except Exception:
            pass

    def _auto_trim_single_well(self, w, tm_lo, tm_hi):
        """
        Devuelve una propuesta de trimming para el pozo w:
        - Sólo actúa si la Tm actual está fuera del rango [tm_lo, tm_hi].
        - Quita puntos de un lado u otro, de uno en uno, hasta meter la Tm en el rango.
        - Minimiza el número de puntos eliminados y respeta mínimo 3 puntos.
        """
        if w in self.deleted_wells:
            return None

        g = self._get_visible_df(w)
        if g is None or len(g) < 3:
            return None

        x = g["Temperature"].values.astype(float)
        y = g["Fluorescence"].values.astype(float)

        tm0, _, _ = self._compute_tm_for_xy(x, y)
        if tm0 is None or not np.isfinite(tm0):
            return None

        # Si ya está dentro del rango, no se toca
        if tm_lo <= tm0 <= tm_hi:
            return None

        n = len(x)
        left = 0
        right = n - 1
        removed_low = 0
        removed_high = 0
        tm_current = tm0

        max_remove = n - 3  # mínimo 3 puntos

        last_good = None
        orig_min = float(x[0])
        orig_max = float(x[-1])

        while (removed_low + removed_high) < max_remove:
            if tm_current is None or not np.isfinite(tm_current):
                break
            if tm_lo <= tm_current <= tm_hi:
                last_good = (left, right, tm_current, removed_low, removed_high)
                break

            if tm_current < tm_lo:
                # Tm demasiado baja -> quitar puntos por abajo
                if left >= right:
                    break
                left += 1
                removed_low += 1
            elif tm_current > tm_hi:
                # Tm demasiado alta -> quitar puntos por arriba
                if right <= left:
                    break
                right -= 1
                removed_high += 1
            else:
                # Ya estaría dentro, pero lo habríamos capturado arriba
                break

            if (right - left + 1) < 3:
                break

            x_seg = x[left:right+1]
            y_seg = y[left:right+1]
            tm_new, _, _ = self._compute_tm_for_xy(x_seg, y_seg)
            tm_current = tm_new

            if tm_current is not None and np.isfinite(tm_current) and tm_lo <= tm_current <= tm_hi:
                last_good = (left, right, tm_current, removed_low, removed_high)
                break

        if last_good is None:
            return None

        left, right, tm_final, removed_low, removed_high = last_good
        total_removed = removed_low + removed_high
        if total_removed == 0:
            return None

        new_tmin = float(x[left])
        new_tmax = float(x[right])

        # grados eliminados por abajo y por arriba
        deg_low = max(0.0, new_tmin - orig_min)
        deg_high = max(0.0, orig_max - new_tmax)

        return {
            "tm_before": float(tm0),
            "tm_after": float(tm_final),
            "removed_low": int(removed_low),
            "removed_high": int(removed_high),
            "removed_total": int(total_removed),
            "new_tmin": new_tmin,
            "new_tmax": new_tmax,
            "deg_low": float(deg_low),
            "deg_high": float(deg_high),
        }

    def _auto_trim_to_expected_range(self):
        """
        Botón 'Auto-trim to expected range' con diálogo tipo tabla y checkboxes.

        IMPORTANTE:
        - El Expected Tm range NO restringe el cálculo de Tm; sólo se usa para decidir el recorte.
        - Aquí sí metemos el auto-trim en el flujo de undo/redo:
          * Antes de aplicar el recorte a un pozo, guardamos su estado en history (self._push_history).
          * El recorte se hace sobre self.per_well_work[w], de modo que Undo/Redo lo pueda revertir.
        """
        tm_win = self._get_tm_window()
        if tm_win is None:
            messagebox.showinfo(
                "Expected Tm range",
                "Please set a valid Expected Tm range (min and max in °C) first."
            )
            return

        lo, hi = tm_win

        try:
            changes = {}
            for w in self.wells:
                if w in self.deleted_wells:
                    continue
                res = self._auto_trim_single_well(w, lo, hi)
                if res is not None:
                    changes[w] = res

            if not changes:
                messagebox.showinfo(
                    "Auto-trim",
                    "No wells required trimming or no valid auto-trim solution was found."
                )
                return

            # ------- Crear diálogo con tabla + checkboxes -------
            dlg = tk.Toplevel(self)
            dlg.title("Auto-trim to expected Tm range")
            dlg.transient(self)
            dlg.grab_set()  # modal

            info_txt = (
                f"Review auto-trim proposals to keep Tm in [{lo:.2f}, {hi:.2f}] °C.\n"
                "Tick the wells you want to apply and click 'Apply selected'."
            )
            ttk.Label(dlg, text=info_txt, justify="left").pack(fill=tk.X, padx=8, pady=(8, 4))

            container = ttk.Frame(dlg)
            container.pack(fill=tk.BOTH, expand=True, padx=8, pady=4)

            canvas = tk.Canvas(container, height=320)
            canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
            vsb = ttk.Scrollbar(container, orient="vertical", command=canvas.yview)
            vsb.pack(side=tk.RIGHT, fill=tk.Y)
            canvas.configure(yscrollcommand=vsb.set)

            inner = ttk.Frame(canvas)
            canvas.create_window((0, 0), window=inner, anchor="nw")

            def _on_configure(event):
                canvas.configure(scrollregion=canvas.bbox("all"))

            inner.bind("<Configure>", _on_configure)

            # ---------- Cabecera de la tabla ----------
            headers = [
                "Apply",
                "Well",
                "Δlow (°C)",
                "Δhigh (°C)",
                "Tm before (°C)",
                "Tm after (°C)",
                "ΔTm (°C)",
            ]
            for col, h in enumerate(headers):
                lbl = ttk.Label(inner, text=h, font=("TkDefaultFont", 9, "bold"))
                lbl.grid(row=0, column=col, padx=4, pady=2, sticky="w")

            # ---------- Filas por pozo ----------
            check_vars = {}
            sorted_wells = sorted(changes.keys(), key=self._well_sortkey)
            for r, w in enumerate(sorted_wells, start=1):
                res = changes[w]
                var = tk.BooleanVar(value=True)
                check_vars[w] = var

                deg_low = float(res.get("deg_low", 0.0) or 0.0)
                deg_high = float(res.get("deg_high", 0.0) or 0.0)
                tm_before = res.get("tm_before", None)
                tm_after = res.get("tm_after", None)
                if tm_before is not None and tm_after is not None:
                    dtm = tm_after - tm_before
                else:
                    dtm = None

                cb = ttk.Checkbutton(inner, variable=var)
                cb.grid(row=r, column=0, padx=4, pady=2, sticky="w")

                ttk.Label(inner, text=w).grid(row=r, column=1, padx=4, pady=2, sticky="w")
                ttk.Label(inner, text=f"{deg_low:.2f}").grid(row=r, column=2, padx=4, pady=2, sticky="w")
                ttk.Label(inner, text=f"{deg_high:.2f}").grid(row=r, column=3, padx=4, pady=2, sticky="w")

                txt_before = f"{tm_before:.2f}" if tm_before is not None else "n/a"
                ttk.Label(inner, text=txt_before).grid(row=r, column=4, padx=4, pady=2, sticky="w")

                txt_after = f"{tm_after:.2f}" if tm_after is not None else "n/a"
                ttk.Label(inner, text=txt_after).grid(row=r, column=5, padx=4, pady=2, sticky="w")

                txt_dtm = f"{dtm:+.2f}" if dtm is not None else "n/a"
                ttk.Label(inner, text=txt_dtm).grid(row=r, column=6, padx=4, pady=2, sticky="w")

            # ---------- Botones inferior ----------
            btn_frame = ttk.Frame(dlg)
            btn_frame.pack(fill=tk.X, padx=8, pady=(4, 8))

            def on_select_none():
                for v in check_vars.values():
                    v.set(False)

            def on_select_all():
                for v in check_vars.values():
                    v.set(True)

            def on_apply():
                applied = 0
                for w in sorted_wells:
                    if not check_vars[w].get():
                        continue

                    res = changes[w]
                    tmin_new = float(res.get("new_tmin"))
                    tmax_new = float(res.get("new_tmax"))
                    
                    g_full = self.per_well_work[w].copy().sort_values("Temperature").reset_index(drop=True)
                    x_full = g_full["Temperature"].values.astype(float)
                    
                    mask = (x_full >= tmin_new) & (x_full <= tmax_new)
                    if mask.sum() < 3:
                        continue
                    
                    # Guardar estado anterior (DF + trim + flag auto-trim) para Undo/Redo
                    self._push_history(w)
                    
                    # 🔹 NO tocamos per_well_work: mantenemos todos los puntos
                    # Solo guardamos el rango de análisis como trimming lógico
                    self.trim_ranges[w] = (tmin_new, tmax_new)
                    
                    # Marcamos que este pozo tiene auto-trim
                    self.auto_trimmed_wells.add(w)
                    
                    if w not in self.corrected_wells:
                        self.corrected_wells.append(w)
                    
                    # Invalidar cache de Tm para este pozo
                    self._tm_cache.pop(w, None)
                    applied += 1

                if applied > 0:
                    self._recompute_tm_all_wells()
                    self._refresh_all_wells_with_tm()
                    self._refresh_corrected_list()
                    self._refresh_suspected_list()
                
                    # 🔹 Sincronizar sliders de Analysis T range del pozo actual
                    if self.current_well is not None:
                        self._init_t_range_for_well()
                
                    self._draw_current()
                    self.status_var.set(f"Auto-trim applied to {applied} wells.")
                else:
                    self.status_var.set("Auto-trim: no changes applied.")

                # 🔹 SOLUCIÓN: liberar el grab y devolver foco a la ventana principal
                try:
                    dlg.grab_release()
                except Exception:
                    pass
                dlg.destroy()
                self._restore_main_focus()

            def on_cancel():
                self.status_var.set("Auto-trim cancelled.")
                # 🔹 Igual que en apply: liberar grab, destruir y restaurar foco
                try:
                    dlg.grab_release()
                except Exception:
                    pass
                dlg.destroy()
                self._restore_main_focus()

            ttk.Button(btn_frame, text="Select none", command=on_select_none).pack(side=tk.LEFT, padx=(0, 4))
            ttk.Button(btn_frame, text="Select all", command=on_select_all).pack(side=tk.LEFT, padx=(0, 4))
            ttk.Button(btn_frame, text="Apply selected", command=on_apply).pack(side=tk.RIGHT, padx=(4, 0))
            ttk.Button(btn_frame, text="Cancel", command=on_cancel).pack(side=tk.RIGHT, padx=(4, 0))

            dlg.bind("<Return>", lambda e: on_apply())
            dlg.bind("<Escape>", lambda e: on_cancel())

        finally:
            # Aseguramos que los entries se quedan normales
            self.abs_thr_entry.config(state="normal")
            self.kdisp_entry.config(state="normal")

    # ------------------------ events ------------------------
    def _goto_selected_well(self, event=None):
        self._on_select_well()

    def _parse_well_label_to_name(self, label):
        """Extract well name from label that may include Tm text, [DELETED] or ✂."""
        if "—" in label:
            return label.split("—", 1)[0].strip()
        if "[DELETED]" in label:
            return label.split("[DELETED]", 1)[0].strip()
        if "✂" in label:
            return label.split("✂", 1)[0].strip()
        if "-" in label:
            return label.split("-", 1)[0].strip()
        return label.strip()

    def _on_select_well(self, event=None):
        sel = self.well_list.curselection()
        if not sel:
            self.current_well = None
            self._clear_plot()
            self.apply_btn.config(state="disabled")
            self.undo_btn.config(state="disabled")
            self.redo_btn.config(state="disabled")
            return

        label = self.well_list.get(sel[0])
        w = self._parse_well_label_to_name(label)
        self.current_well = w

        g_work = self.per_well_work[w].sort_values("Temperature").reset_index(drop=True)
        self.per_well_work[w] = g_work

        n = len(g_work)
        if n < 2:
            self.idx_slider.configure(from_=0, to=0)
            self.idx_slider.set(0)
            self.apply_btn.config(state="disabled")
            self.max_idx = 0
        else:
            self.idx_slider.configure(from_=0, to=max(0, n - 2))
            self.idx_slider.set((n - 2) / 2)
            self.apply_btn.config(state="normal" if w not in self.deleted_wells else "disabled")
            self.max_idx = n - 2

        if w in self.auto_suspect_index:
            i_star = self._clamp_index(self.auto_suspect_index[w])
            self.idx_slider.set(i_star)

        self.idx_entry.delete(0, tk.END)
        self.idx_entry.insert(0, str(int(round(self.idx_slider.get()))))

        self._init_t_range_for_well()
        self._update_selected_idx_label()
        self._draw_current()
        self._update_undo_redo_state()

    def _on_select_corrected(self, event=None):
        sel = self.corrected_list.curselection()
        if not sel:
            return
        w = self.corrected_list.get(sel[0])
        # Puede venir con ✂, limpiamos
        w = self._parse_well_label_to_name(w)
        try:
            idx = self.wells.index(w)
            self.well_list.selection_clear(0, tk.END)
            self.well_list.selection_set(idx)
            self.well_list.see(idx)
            self._on_select_well()
        except ValueError:
            pass

    def _on_select_suspected(self, event=None):
        sel = self.suspected_list.curselection()
        if not sel:
            return
        w = self.suspected_list.get(sel[0])
        try:
            idx = self.wells.index(w)
            self.well_list.selection_clear(0, tk.END)
            self.well_list.selection_set(idx)
            self.well_list.see(idx)
            self._on_select_well()
        except ValueError:
            pass

    def _on_select_tm_outlier(self, event=None):
        sel = self.tm_outlier_list.curselection()
        if not sel:
            return
        label = self.tm_outlier_list.get(sel[0])
        w = self._parse_well_label_to_name(label)
        try:
            idx = self.wells.index(w)
            self.well_list.selection_clear(0, tk.END)
            self.well_list.selection_set(idx)
            self.well_list.see(idx)
            self._on_select_well()
        except ValueError:
            pass

    def _goto_tm_outlier(self, event=None):
        self._on_select_tm_outlier()

    def _on_slider(self, *args):
        try:
            val = float(self.idx_slider.get())
        except Exception:
            val = 0
        i = self._clamp_index(int(round(val)))
        self.selected_idx = i
        self.idx_entry.delete(0, tk.END)
        self.idx_entry.insert(0, str(i))
        self._update_selected_idx_label()
        self._draw_current()

    def _nudge_index(self, step):
        if self.current_well is None:
            return
        try:
            i = int(self.idx_entry.get())
        except ValueError:
            i = int(round(self.idx_slider.get()))
        i = self._clamp_index(i + step)
        self.idx_slider.set(i)
        self.selected_idx = i
        self.idx_entry.delete(0, tk.END)
        self.idx_entry.insert(0, str(i))
        self._update_selected_idx_label()
        self._draw_current()

    def _move_well_selection(self, step):
        """Mueve la selección en la lista de All wells con las flechas ↑ / ↓."""
        if not self.wells:
            return

        sel = self.well_list.curselection()
        if sel:
            idx = sel[0]
        else:
            # Si no hay nada seleccionado, empezamos desde el principio o el final
            idx = 0 if step > 0 else len(self.wells) - 1

        idx = max(0, min(len(self.wells) - 1, idx + step))

        self.well_list.selection_clear(0, tk.END)
        self.well_list.selection_set(idx)
        self.well_list.see(idx)
        self._on_select_well()

    def _set_index_from_entry(self):
        try:
            i = int(self.idx_entry.get())
        except ValueError:
            messagebox.showinfo("Invalid index", "Enter an integer (0..n-2).")
            return
        i = self._clamp_index(i)
        self.idx_slider.set(i)
        self.selected_idx = i
        self._update_selected_idx_label()
        self._draw_current()

    def _on_plot_click(self, event):
        if self.current_well is None or event.xdata is None:
            return
        g = self.per_well_work[self.current_well].sort_values("Temperature").reset_index(drop=True)
        x = g["Temperature"].values
        i = int(np.argmin(np.abs(x - event.xdata)))
        i = self._clamp_index(i)
        self.idx_slider.set(i)
        self.selected_idx = i
        self._update_selected_idx_label()
        self._draw_current()

    # ------------------------ correction / undo / redo / animation ------------------------
    def _push_history(self, well):
        """
        Guarda en la pila de undo el estado actual del pozo:
        - DataFrame de trabajo (per_well_work[well])
        - trim_ranges[well] (o None si no tiene)
        - flag de auto_trim (si el pozo está en auto_trimmed_wells)
        Y limpia las pilas de redo correspondientes.
        """
        if not well or well not in self.per_well_work:
            return

        # Aseguramos estructuras
        if well not in self.history:
            self.history[well] = []
        if well not in self.redo_history:
            self.redo_history[well] = []
        if well not in self.trim_history:
            self.trim_history[well] = []
        if well not in self.trim_redo:
            self.trim_redo[well] = []
        if well not in self.auto_trim_history:
            self.auto_trim_history[well] = []
        if well not in self.auto_trim_redo:
            self.auto_trim_redo[well] = []

        # Guardar snapshot del DF de trabajo
        self.history[well].append(self.per_well_work[well].copy())
        # Guardar snapshot del trimming actual (o None)
        self.trim_history[well].append(self.trim_ranges.get(well, None))
        # Guardar si estaba auto-trimmeado
        self.auto_trim_history[well].append(well in self.auto_trimmed_wells)

        # Borrar las pilas de redo (al hacer un cambio nuevo, el redo previo ya no tiene sentido)
        self.redo_history[well].clear()
        self.trim_redo[well].clear()
        self.auto_trim_redo[well].clear()

    def _animate_transition(self, x, y_from, y_to, steps=12, duration_ms=240):
        if not self.animate_var.get():
            return
        steps = max(3, int(steps))
        delay = max(10, int(duration_ms // steps))
        for k in range(1, steps + 1):
            alpha = k / steps
            yk = y_from * (1 - alpha) + y_to * alpha
            self._draw_current_override(x, yk)
            self.update_idletasks()
            self.after(delay)

    def _draw_current_override(self, x_override, y_override):
        self.ax.clear()
        self.axd.clear()
        w = self.current_well
        if w is None:
            self.canvas.draw_idle()
            return
        g0 = self.per_well_orig[w]
        x0, y0 = g0["Temperature"].values, g0["Fluorescence"].values
        self.ax.plot(x0, y0, alpha=0.4, linewidth=1, label="Original")

        # visible / trimmed data
        g_vis = self._get_visible_df(w)
        x_vis = g_vis["Temperature"].values
        y_plot = self._maybe_smooth(y_override)
        self.ax.plot(x_vis, y_plot, linewidth=1.6, label="Corrected")

        tm, xd, dplot = self._compute_tm(w)

        if tm is not None:
            self.ax.axvline(tm, linestyle=":", alpha=0.4)
        if self.show_deriv_var.get() and (xd is not None) and (dplot is not None):
            self.axd.plot(xd, dplot, linewidth=1)
            if tm is not None:
                self.axd.axvline(tm, linestyle=":", alpha=0.8)

        tmin, tmax = self._get_current_t_range()
        if tmin is not None and tmax is not None and tmin < tmax:
            self.ax.axvline(tmin, linestyle="--", alpha=0.3)
            self.ax.axvline(tmax, linestyle="--", alpha=0.3)
            self.axd.axvline(tmin, linestyle="--", alpha=0.3)
            self.axd.axvline(tmax, linestyle="--", alpha=0.3)

        if tm is not None and np.isfinite(tm):
            is_out = (w in self.tm_outlier_wells)
            txt_color = "#b00000" if is_out else "#111111"
            self.ax.text(
                0.02,
                0.95,
                f"Tm = {tm:.2f} °C",
                transform=self.ax.transAxes,
                ha="left",
                va="top",
                fontsize=10,
                color=txt_color,
                bbox=dict(boxstyle="round,pad=0.25", fc="#ffffff", ec="none", alpha=0.75),
            )

        self.ax.set_ylabel("Fluorescence")
        self.axd.set_xlabel("Temperature")
        self.axd.set_ylabel("-dF/dT (oriented up)")
        self.ax.legend(loc="upper right")
        self.canvas.draw_idle()

    def _apply_correction(self):
        if self.current_well is None or self.current_well in self.deleted_wells:
            return
        g_full = self.per_well_work[self.current_well].copy()
        n = len(g_full)
        if n < 2:
            return
        i = self._clamp_index(int(round(self.idx_slider.get())))
        y_full = g_full["Fluorescence"].values.astype(float)

        if self.multi_jump_var.get():
            self._push_history(self.current_well)
            y_from = y_full.copy()
            changed = self._apply_multi_jump(self.current_well, iterative=self.iterative_var.get())
            y_to = self.per_well_work[self.current_well]["Fluorescence"].values.astype(float)
            if changed and self.animate_var.get():
                self._animate_transition(g_full["Temperature"].values, y_from, y_to)
        else:
            diffs = np.diff(y_full)
            if len(diffs) == 0:
                return
            delta = float(y_full[i + 1] - y_full[i])
            op = self.op_var.get()
            adj = -delta if op in ("auto", "sub") else +delta
            self._push_history(self.current_well)
            y_from = y_full.copy()
            y_full[i + 1 :] = y_full[i + 1 :] + adj
            g_full["Fluorescence"] = y_full
            self.per_well_work[self.current_well] = g_full
            if self.animate_var.get():
                self._animate_transition(g_full["Temperature"].values, y_from, y_full)

        if self.current_well not in self.corrected_wells:
            self.corrected_wells.append(self.current_well)
        if self.current_well in self.suspected_wells:
            self.suspected_wells.remove(self.current_well)

        # Invalidar cache de Tm para este pozo
        self._tm_cache.pop(self.current_well, None)
        self._recompute_tm_all_wells()
        self._refresh_all_wells_with_tm()
        self._refresh_corrected_list()
        self._refresh_suspected_list()
        self._update_undo_redo_state()
        self._draw_current()
        self.status_var.set(f"{self.current_well}: correction applied.")

    def _undo_current_well(self):
        w = self.current_well
        if not w or w not in self.history or not self.history[w]:
            return

        # Asegurar estructuras de redo para este pozo
        if w not in self.redo_history:
            self.redo_history[w] = []
        if w not in self.trim_redo:
            self.trim_redo[w] = []
        if w not in self.auto_trim_redo:
            self.auto_trim_redo[w] = []

        # Guardar estado actual en redo
        self.redo_history[w].append(self.per_well_work[w].copy())
        self.trim_redo[w].append(self.trim_ranges.get(w, None))
        self.auto_trim_redo[w].append(w in self.auto_trimmed_wells)

        # Recuperar último estado de undo
        prev_df = self.history[w].pop()
        prev_trim = None
        if w in self.trim_history and self.trim_history[w]:
            prev_trim = self.trim_history[w].pop()
        prev_auto_flag = None
        if w in self.auto_trim_history and self.auto_trim_history[w]:
            prev_auto_flag = self.auto_trim_history[w].pop()

        # Animación SEGURA (solo si longitudes coinciden)
        if self.animate_var.get():
            x_cur = self.per_well_work[w]["Temperature"].values
            y_from = self.per_well_work[w]["Fluorescence"].values.astype(float)
            y_to = prev_df["Fluorescence"].values.astype(float)
            if y_from.shape == y_to.shape:
                self._animate_transition(x_cur, y_from, y_to)

        # Restaurar DF
        self.per_well_work[w] = prev_df

        # Restaurar trimming
        if prev_trim is None:
            self.trim_ranges.pop(w, None)
        else:
            self.trim_ranges[w] = prev_trim

        # Restaurar flag de auto-trim
        if prev_auto_flag is not None:
            if prev_auto_flag:
                self.auto_trimmed_wells.add(w)
            else:
                self.auto_trimmed_wells.discard(w)

        # Ver si el DF es igual al original -> sacar de corrected_wells
        try:
            is_original = prev_df.equals(self.per_well_orig[w])
        except Exception:
            is_original = False
        if is_original and w in self.corrected_wells:
            self.corrected_wells.remove(w)

        # Invalidar cache de Tm para este pozo y refrescar todo
        self._tm_cache.pop(w, None)
        self._recompute_tm_all_wells()
        self._refresh_all_wells_with_tm()
        self._refresh_corrected_list()
        self._scan_suspects()
        self._update_undo_redo_state()
        self.status_var.set(f"{w}: undo applied.")
        self._on_select_well()

    def _redo_current_well(self):
        w = self.current_well
        if not w or w not in self.redo_history or not self.redo_history[w]:
            return

        # Asegurar estructuras de undo para este pozo
        if w not in self.history:
            self.history[w] = []
        if w not in self.trim_history:
            self.trim_history[w] = []
        if w not in self.auto_trim_history:
            self.auto_trim_history[w] = []

        # Guardar estado actual en undo
        self.history[w].append(self.per_well_work[w].copy())
        self.trim_history[w].append(self.trim_ranges.get(w, None))
        self.auto_trim_history[w].append(w in self.auto_trimmed_wells)

        # Recuperar estado desde redo
        next_df = self.redo_history[w].pop()
        next_trim = None
        if w in self.trim_redo and self.trim_redo[w]:
            next_trim = self.trim_redo[w].pop()
        next_auto_flag = None
        if w in self.auto_trim_redo and self.auto_trim_redo[w]:
            next_auto_flag = self.auto_trim_redo[w].pop()

        # Animación SEGURA (solo si longitudes coinciden)
        if self.animate_var.get():
            x_cur = self.per_well_work[w]["Temperature"].values
            y_from = self.per_well_work[w]["Fluorescence"].values.astype(float)
            y_to = next_df["Fluorescence"].values.astype(float)
            if y_from.shape == y_to.shape:
                self._animate_transition(x_cur, y_from, y_to)

        # Restaurar DF
        self.per_well_work[w] = next_df

        # Restaurar trimming
        if next_trim is None:
            self.trim_ranges.pop(w, None)
        else:
            self.trim_ranges[w] = next_trim

        # Restaurar flag de auto-trim
        if next_auto_flag is not None:
            if next_auto_flag:
                self.auto_trimmed_wells.add(w)
            else:
                self.auto_trimmed_wells.discard(w)

        # Si hay redo es que ha habido alguna corrección -> aseguramos que está en corrected
        try:
            is_original = next_df.equals(self.per_well_orig[w])
        except Exception:
            is_original = False
        if is_original:
            if w in self.corrected_wells:
                self.corrected_wells.remove(w)
        else:
            if w not in self.corrected_wells:
                self.corrected_wells.append(w)

        # Invalidar cache de Tm para este pozo y refrescar todo
        self._tm_cache.pop(w, None)
        self._recompute_tm_all_wells()
        self._refresh_all_wells_with_tm()
        self._refresh_corrected_list()
        self._scan_suspects()
        self._update_undo_redo_state()
        self.status_var.set(f"{w}: redo applied.")
        self._on_select_well()


    # ------------------------ plotting & Tm ------------------------
    def _clear_plot(self):
        self.ax.clear()
        self.axd.clear()
        self.canvas.draw_idle()

    def _draw_current(self):
        self.ax.clear()
        self.axd.clear()

        if self.current_well is None:
            self.canvas.draw_idle()
            return

        g0 = self.per_well_orig[self.current_well]
        x0, y0 = g0["Temperature"].values, g0["Fluorescence"].values
        self.ax.plot(x0, y0, alpha=0.5, linewidth=1, label="Original")

        # Si NO está eliminado, dibujar el trazo corregido
        if self.current_well not in self.deleted_wells:
            g1 = self._get_visible_df(self.current_well)
            x1, y1 = g1["Temperature"].values, g1["Fluorescence"].values
            y1_plot = self._maybe_smooth(y1)
            self.ax.plot(x1, y1_plot, linewidth=1.6, label="Corrected")

            if len(x1) >= 2:
                i = self._clamp_index(int(round(self.idx_slider.get())))
                i = max(0, min(i, len(x1) - 1))
                self.ax.axvline(x1[i], linestyle="--", alpha=0.6)
                y_marker = y1_plot[i] if i < len(y1_plot) else y1[i]
                self.ax.scatter([x1[i]], [y_marker], s=40, zorder=5, label=f"i={i}")

        tm, xd, dplot = (None, None, None)
        if self.show_deriv_var.get() and self.current_well not in self.deleted_wells:
            tm, xd, dplot = self._compute_tm(self.current_well)
            if xd is not None and dplot is not None:
                self.axd.plot(xd, dplot, linewidth=1)
            if tm is not None:
                self.axd.axvline(tm, linestyle=":", alpha=0.8)
                self.ax.axvline(tm, linestyle=":", alpha=0.4)

        tmin, tmax = self._get_current_t_range()
        if tmin is not None and tmax is not None and tmin < tmax:
            self.ax.axvline(tmin, linestyle="--", alpha=0.3)
            self.ax.axvline(tmax, linestyle="--", alpha=0.3)
            if self.show_deriv_var.get() and self.current_well not in self.deleted_wells:
                self.axd.axvline(tmin, linestyle="--", alpha=0.3)
                self.axd.axvline(tmax, linestyle="--", alpha=0.3)

        if tm is not None and np.isfinite(tm):
            # Tm en rojo si el pozo actual es outlier y no está eliminado
            is_out = (self.current_well in self.tm_outlier_wells)
            txt_color = "#b00000" if is_out else "#111111"
            self.ax.text(
                0.02,
                0.95,
                f"Tm = {tm:.2f} °C",
                transform=self.ax.transAxes,
                ha="left",
                va="top",
                fontsize=10,
                color=txt_color,
                bbox=dict(boxstyle="round,pad=0.25", fc="#ffffff", ec="none", alpha=0.75),
            )

        self.ax.set_ylabel("Fluorescence")
        self.axd.set_xlabel("Temperature")
        self.axd.set_ylabel("-dF/dT (oriented up)")
        self.ax.legend(loc="upper right")
        self.canvas.draw_idle()
        self._update_undo_redo_state()
        self._update_selected_idx_label()

    def _update_selected_idx_label(self):
        if self.current_well is None:
            self.idx_label.config(text="Index: –/–   |   Temp: –   |   Tm: –")
            return
        g = self._get_visible_df(self.current_well)
        n = len(g)
        if n < 2:
            self.idx_label.config(text="(Curve too short)")
            return
        try:
            i = int(round(float(self.idx_slider.get())))
        except Exception:
            i = 0
        i = max(0, min(i, n - 1))
        t = g.loc[i, "Temperature"]
        tm_val = self._compute_tm(self.current_well)[0] if self.current_well not in self.deleted_wells else None
        tm_txt = f"{tm_val:.2f} °C" if tm_val is not None else "–"
        self.idx_label.config(text=f"Index: {i}/{n-2}   |   Temp: {t:.2f} °C   |   Tm: {tm_txt}")

    # ------------------------ Delete / Recover well ------------------------
    def _get_selected_well_any_list(self):
        # Priority: All wells, Tm outliers, Corrected, Suspected
        sel = self.well_list.curselection()
        if sel:
            label = self.well_list.get(sel[0])
            return self._parse_well_label_to_name(label)
        sel = self.tm_outlier_list.curselection()
        if sel:
            label = self.tm_outlier_list.get(sel[0])
            return self._parse_well_label_to_name(label)
        sel = self.corrected_list.curselection()
        if sel:
            w = self.corrected_list.get(sel[0])
            return self._parse_well_label_to_name(w)
        sel = self.suspected_list.curselection()
        if sel:
            return self.suspected_list.get(sel[0])
        return None

    def _delete_selected_well(self):
        w = self._get_selected_well_any_list()
        if not w:
            messagebox.showinfo("No well selected", "Select a well in any list first.")
            return
        if w not in self.per_well_work:
            messagebox.showinfo("Unknown well", f"Well {w} not found.")
            return

        self.deleted_wells.add(w)  # Marcar como eliminado
        self.status_var.set(f"{w}: well deleted (will not be exported).")

        # Invalidar caché y recalcular
        self._tm_cache.clear()
        self._recompute_tm_all_wells()
        self._refresh_all_wells_with_tm()
        self._refresh_corrected_list()
        self._refresh_suspected_list()
        self._draw_current()

    def _recover_selected_well(self):
        w = self._get_selected_well_any_list()
        if not w:
            messagebox.showinfo("No well selected", "Select a well in any list first.")
            return
        if w not in self.deleted_wells:
            messagebox.showinfo("Not deleted", f"Well {w} is not deleted.")
            return

        self.deleted_wells.remove(w)  # Recuperar
        self.status_var.set(f"{w}: well recovered.")

        # Recalcular todo
        self._tm_cache.clear()
        self._recompute_tm_all_wells()
        self._refresh_all_wells_with_tm()
        self._refresh_corrected_list()
        self._refresh_suspected_list()
        self._draw_current()

    # ------------------------ review mode ------------------------
    def _jump_review(self, list_name="suspected", step=+1):
        arr = self.suspected_wells if list_name == "suspected" else self.corrected_wells
        if not arr:
            return
        w_cur = self.current_well
        if w_cur in arr:
            idx = arr.index(w_cur)
        else:
            idx = -1 if step > 0 else 0
        idx = (idx + step) % len(arr)
        w = arr[idx]
        try:
            pos = self.wells.index(w)
            self.well_list.selection_clear(0, tk.END)
            self.well_list.selection_set(pos)
            self.well_list.see(pos)
            self._on_select_well()
        except ValueError:
            pass

    def _jump_tm_outlier(self, step=+1):
        """Moverse entre pozos marcados como Tm outliers con feedback."""
        # Filtrar outliers eliminados
        valid_outliers = [w for w in self.tm_sorted_wells if w not in self.deleted_wells]
        if not valid_outliers:
            self.status_var.set("No Tm outliers to review")
            return

        w_cur = self.current_well
        if w_cur in valid_outliers:
            idx = valid_outliers.index(w_cur)
        else:
            idx = -1 if step > 0 else 0

        idx = (idx + step) % len(valid_outliers)
        w = valid_outliers[idx]

        # Actualizar UI con posición
        self.status_var.set(f"Tm outlier {idx+1}/{len(valid_outliers)}: {w}")

        try:
            pos = self.wells.index(w)
            self.well_list.selection_clear(0, tk.END)
            self.well_list.selection_set(pos)
            self.well_list.see(pos)
            self._on_select_well()
        except ValueError:
            pass


# ------------------------ run ------------------------
def main():
    path = sys.argv[1] if len(sys.argv) > 1 else None
    app = DSF_Harmonizer(path)
    app.mainloop()


if __name__ == "__main__":
    main()

