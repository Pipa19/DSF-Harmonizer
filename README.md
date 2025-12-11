# DSF Harmonizer — Complete User Guide

**DSF Harmonizer** is an interactive tool for cleaning, correcting, and analyzing DSF curves (`.gdsf`).
Although the window displays the title *DSF Harmonizer*, it is the same tool traditionally known as your *DSF Step Fixer*.

## How to Run

From the terminal, navigate to the script directory and run:

```bash
# The .gdsf file must be in the same directory (or provide a full path)
python3 dsf_step_fixer.py YOUR_FILE.gdsf
```

Expected `.gdsf` format (tab-separated, no header):

Well (tab) Temp (tab) Fluor.
A1	20.0001	16588.30664
A1	20.2411	16380.11035
A1	20.4821	16248.78516
A1	20.7232	15978.29785
A1	20.9642	16185.03418
A1	21.2052	15990.87012
…	…	…
A2	20.0001	25971.25586
A2	20.2411	25679.30859
…	…	…
…	…	…
P24	99.7791	5478.471191


# DSF Harmonizer — Complete User Guide

**DSF Harmonizer** is an interactive tool for cleaning, correcting, and analyzing DSF curves (`.gdsf`).
Although the window displays the title *DSF Harmonizer*, it is the same tool traditionally known as your *DSF Step Fixer*.

---

## How to Run

From the terminal, navigate to the script directory and run:

```bash
# The .gdsf file must be in the same directory (or provide a full path)
python3 dsf_step_fixer.py YOUR_FILE.gdsf
```

Expected `.gdsf` format (tab-separated, no header):

```
Well    Temperature    Fluorescence
```

---

# Main Features

* Manual or automatic correction of DSF curve jumps (“steps”).
* Automatic detection of suspicious wells using MAD or STD.
* Multi-jump engine for wells with multiple discontinuities.
* Real-time smoothing (Savitzky–Golay if available, otherwise adaptive moving average).
* Tm calculation as the global maximum of the upward-oriented derivative (−dF/dT).
* Interactive trimming of temperature data using sliders (reversible; per-well).
* Auto-trimming based on an expected Tm range defined by the user.
* Identification of Tm outliers with configurable thresholds and optional reference Tm.
* Quick navigation through:

  * Suspected wells
  * Corrected wells
  * Tm outliers
* Ability to delete wells (logical deletion) and later recover them.
* Export options for corrected curves, smoothed curves, and Tm tables.
* Graphical indicators for wells that are corrected, suspected, trimmed, outliers, or deleted.

---

# Manual Workflow (Recommended)

1. **Open file**
   Use the *Open .gdsf* button.

2. **Select a well**
   From the left panel (*All wells (with data)*).

3. **Navigate the breakpoint**
   Use:

   * The slider
   * Arrow buttons
   * Index input + *Go*
   * Clicking directly on the plot (nearest temperature)

4. **Read point information**
   Displayed as:
   `Index i/n | Temp: XX.XX °C | Tm: YY.YY °C`

5. **Choose operation**

   * **Auto** → program decides Add/Subtract
   * **Add** → shifts curve upward from index+1
   * **Subtract** → shifts curve downward from index+1

6. **Apply correction**
   Applies the shift from index+1 onward (animated if enabled).

7. **Repeat if needed**
   You may continue fixing steps in the same well or use Review Mode to jump between wells.

8. **Export results**

   * *Export corrected*
   * *Export corrected + smoothed*
   * *Export Tm table*

---

# Undo / Redo

* **Undo:** Ctrl+Z / Cmd+Z
* **Redo:** Ctrl+Y / Ctrl+Shift+Z / Cmd+Shift+Z

Undo/Redo is *per well* and supports:

* Step corrections
* Remove data trims
* Auto-trim changes

**Note:** Delete Well / Recover Well do NOT participate in Undo/Redo.

---

# Automatic Step Detection & Correction

In the "Step correction" toolbar:

* Abs threshold
* k·disp threshold
* Dispersion method: MAD or STD
* Iterative mode
* Multi-jump
* "Show derivative & Tm"
* Buttons:

  * **Scan suspects**
  * **Correct all suspects**

**Typical pipeline:**

1. Adjust Abs threshold, k, and dispersion method.
2. Run **Scan suspects**.
3. Run **Correct all suspects**:

   * Uses Multi-jump if enabled
   * Applies multiple rounds if *Iterative* is on
4. Review corrected wells and use Undo if necessary.

---

# Smoothing

Controls:

* Toggle: **Smooth ON/OFF**
* Slider: Strength (0–100)

Behavior:

### Smooth OFF

* Corrected curve displayed raw.
* Derivative always uses a minimal fixed smoothing.
* Slider has no effect.

### Smooth ON

* Corrected curve smoothed using slider value.
* Derivative smoothing = max(base smoothing, slider).
* Affects:

  * Displayed derivative
  * Tm value
  * `Tm_smoothed` in export tables

Internally:

* Uses Savitzky–Golay if available.
* Falls back to adaptive moving average.

---

# Temperature Trimming

Below the plots there is a **Temperature ranges** panel.

## 1. Analysis T range

Adjust Tmin and Tmax using sliders or text boxes.
Press **Remove data** to set the active analysis range for that well.

Effects:

* Limits Tm calculation to that range.
* Affects export.
* Fully reversible with Undo/Redo.

Useful for removing:

* Noisy low-temperature regions (<30 °C)
* Saturated high-temperature regions (>90 °C)

---

## 2. Auto-trim to Expected Tm Range

Fields:

* Expected Tm range (min, max)
  Button: **Auto-trim to expected range**

Important:

* Expected range **does not** restrict Tm directly.
* It is only used to *propose* trims that would move Tm into the desired range.

Workflow:

1. Enter expected Tm range (e.g., 50–65 °C).
2. Press **Auto-trim**.
3. A proposal table appears with suggested crops.
4. Select wells to apply; press **Apply selected**.
5. Changes are Undo-compatible and wells get the ✂ mark.

Deleted wells are excluded.

---

# Derivative & Tm Calculation

For each non-deleted well:

1. Use corrected + trimmed curve.
2. Sort by temperature.
3. Apply smoothing (base + optional).
4. Compute derivative and orient it upward.
5. Normalize derivative for plotting.
6. **Tm = temperature where derivative reaches its global maximum.**
7. Store Tm for:

   * All wells list
   * Tm outliers list
   * Export table

On screen:

* Derivative shown in lower panel
* Vertical dashed line at Tm
* "Tm = XX.XX °C" text (black = normal, red = outlier)

---

# Tm Outliers

## All wells (with data)

Each entry shows:

* `A01 — Tm=57.32 °C`
* Deleted wells → `A01 [DELETED]` (gray)
* Trimmed wells → ✂
* Outliers → red text

## Tm outliers list

Includes wells where:

```
|Tm − Tm_ref| ≥ threshold
```

Tm_ref:

* User-specified ("I know my Tm"), or
* Mean Tm of valid wells (if field is empty)

Controls:

* Threshold (°C)
* Reference Tm (optional)
* Current mean Tm

---

# Review Mode

Navigation for:

* Prev/Next Suspected
* Prev/Next Corrected
* Prev/Next Tm outlier

Status bar indicates progress (e.g., `Tm outlier 3/7: B05`).

---

# Color System

* Corrected → light green
* Suspected → amber
* Normal → white
* Tm outlier → red text
* Deleted → light gray + `[DELETED]`
* Trimmed (Remove data or Auto-trim) → ✂

On plot:

* Black Tm label = normal
* Red Tm label = outlier

---

# Delete Well / Recover Well

**Delete Well:**

* Marks well as deleted
* Removes it from:

  * Tm calculation
  * Tm outliers
  * Exports
* Original curve remains visible
* Corrected curve is hidden

**Recover Well:**

* Restores previous state
* Re-enables corrected curve
* Well becomes eligible for Tm, outliers, export

Delete/Recover do not use Undo/Redo.

---

# Keyboard Shortcuts

* Undo: Ctrl+Z / Cmd+Z
* Redo: Ctrl+Y / Ctrl+Shift+Z / Cmd+Shift+Z
* Next suspected: `s`
* Prev suspected: `Shift+S`
* Next corrected: `c`
* Prev corrected: `Shift+C`
* Move breakpoint: ← →
* Move well selection: ↑ ↓
* Apply typed index: Enter

---

# Export Options

## 1. Export corrected

* Corrected curves (no export-time smoothing)
* Respects trims and step corrections
* Excludes deleted wells
* Formats: `.gdsf`, `.tsv`, `.txt`

## 2. Export corrected + smoothed

* Same as above but applies smoothing at export time
* Uses current slider value even if smoothing is OFF
* Ideal for machine learning / fitting workflows

## 3. Export Tm table

Includes:

* Well
* `Tm_corrected`
* `Tm_smoothed`
* `Smooth_strength`
  Excludes deleted wells
  Formats: `.tsv`, `.csv`

---

# MAD, k, and Thresholds (Concepts)

For each well:

```
diffs[i] = y[i+1] - y[i]
```

Dispersion:

* **MAD (recommended)**
  Robust; σᵣ = 1.4826 × MAD(diffs)
* **STD**

Jumps are flagged if:

* diff > Abs threshold
* diff > k × dispersion

Typical values:

* **k ≈ 6 with MAD** works very well.

---

# Practical Tips

* Fix steps first; refine Tm later.
* If too many suspects → increase k or Abs threshold.
* If jumps are missed → lower k.
* For noisy derivatives → trim temperature range or enable smoothing.
* For multi-phase melting → use Auto-trim to isolate relevant region.
* Use Tm outliers list to review problematic wells quickly.
* If expected Tm is known → fill "I know my Tm" and adjust threshold.

---

# Troubleshooting

**Derivative changes while smoothing OFF?**
Derivative always uses a fixed minimal smoothing; slider does nothing unless smoothing is ON.

**Tm peak appears downward?**
Derivative is auto-oriented upward.

**Tm jumps to unrealistic values (e.g., 20 or 95 °C)?**
Check step corrections, trimming, or outlier status.

**Undo not working?**
Undo is per well and only works after an action.

**Exported curves have fewer points?**
Due to applied Remove data or Auto-trim trims.
