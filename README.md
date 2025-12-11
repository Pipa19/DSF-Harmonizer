# DSF-Harmonizer
A small script I created to pre-process and polish Differential Scanning Fluorimetry curves. It might be useful when dealing with unexpected steps in the graphs, tricky regions of the curves that affect the first derivative and Tm estimation, or noisy data that needs softening.

To run this script, go to its directory in the terminal and execute:

python3 dsf_step_fixer.py YOUR_DSF_FILE.gdsf

# The .gdsf file must be in the same directory (or you must specify the path). This .gdsf file is meant to be compatible with the HTSDSF Explorer. It has to have this data structure (consider making an script to change your data to this format):

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

🧪 What does this program do?
Interactive tool to clean and analyze DSF curves (.gdsf)..

Main features:
	•	Fixes curve jumps (steps) manually or automatically.
	•	Automatically detects suspicious wells using MAD or STD.
	•	Corrects multiple jumps per well with a multi-jump engine.
	•	Real-time curve smoothing (Savitzky–Golay if available, otherwise moving average).
	•	Computes Tm as the global maximum of the upward-oriented −dF/dT.
	•	Interactive temperature-range trimming with sliders (Remove data, reversible).
	•	Auto-trimming so each well’s Tm falls within a user-defined Expected Tm range.
	•	Marks wells with unusual Tm as Tm outliers, with a configurable threshold and optional reference Tm.
	•	Fast navigation through: ◦ Suspected ◦ Corrected ◦ Tm outliers
	•	Allows deleting entire wells (Delete Well) and restoring them (Recover Well), excluding them from: ◦ Tm calculation ◦ Outlier lists ◦ Export
	•	Marks wells with trimmed data (manual or auto-trim) with a ✂ icon in lists.
	•	Exports: ◦ Corrected curves ◦ Corrected + smoothed curves ◦ Tm table (corrected + smoothed)
Expected .gdsf file format (tab-separated, no header):
Well    Temperature    Fluorescence

🚀 Recommended workflow (manual)
	1	Open file Click Open .gdsf → select the file.
	2	Select a well In the left panel (All wells (with data)) choose a well.
	3	Navigate the breakpoint In Breakpoint index (row “Per-well controls”):
	◦	Move the slider.
	◦	Use ◀ ▶ buttons.
	◦	Type an index + press Go.
	◦	Click on the plot → picks the nearest temperature.
	4	Read the info of the point Below the slider: Index: i/n | Temp: XX.XX °C | Tm: YY.YY °C
	5	Choose operation (step correction) In Operation:
	◦	Auto → program decides add or subtract.
	◦	Add → shifts upward from i+1.
	◦	Subtract → shifts downward from i+1.
	6	Apply correction Click Apply correction at breakpoint:
	◦	Adjusts all points from i+1 onward.
	◦	If animation is enabled, you’ll see the smooth transition.
	7	Repeat if there are more jumps
	◦	Stay in the same well, or
	◦	Use Review mode to hop through suspected / corrected / Tm outliers.
	8	Export when the plate looks good:
	◦	Export corrected
	◦	Export corrected + smoothed
	◦	Export Tm table

⚡ Undo / Redo
	•	Undo: Ctrl+Z / Cmd+Z
	•	Redo: Ctrl+Y / Ctrl+Shift+Z / Cmd+Shift+Z
Characteristics:
	•	Per-well, each well has its own history.
	•	Affects: ◦ Step corrections ◦ Remove data trims ◦ Auto-trim changes
	•	If animation is ON, transitions between states are shown (if curves have same length).
Note: Delete Well / Recover Well are not part of Undo/Redo; they’re controlled by their own buttons.

🤖 Automatic pipeline (step detection & correction) Second top bar (Step correction):
	•	Step correction (group label)
	•	Abs threshold
	•	k·disp
	•	Dispersion: MAD or STD
	•	Iterative
	•	Show derivative & Tm
	•	Multi-jump
	•	Buttons: ◦ Scan suspects ◦ Correct all suspects
Typical flow:
	1	Adjust Abs threshold, k, and Dispersion (MAD/STD).
	2	Click Scan suspects → fills the Suspected list.
	3	Click Correct all suspects → batch correction:
	◦	With Multi-jump ON: uses the multi-step engine.
	◦	With Iterative ON: can run several rounds until no more jumps appear.
	4	Review Corrected wells:
	◦	Use Review mode or click on the list.
	◦	If something looks wrong → Undo for that well.

📈 Smoothing Fourth top bar:
	•	Smooth on/off
	•	Strength (slider 0–100 + numeric box)
Behavior:
	•	Smooth OFF: ◦ Corrected curve displayed raw. ◦ Derivative always uses a minimum smoothing (base e.g. 25). ◦ Slider does nothing while OFF.
	•	Smooth ON: ◦ Corrected curve is smoothed using the slider value. ◦ Derivative smoothing uses max(base, slider). ◦ Active smoothing affects: ▪ Displayed derivative ▪ Tm calculation (for All wells, Tm outliers, and per-well plot) ▪ Tm_smoothed column of the Tm table
Internally:
	•	Tries Savitzky–Golay (scipy.signal.savgol_filter).
	•	If unavailable, uses adaptive symmetric moving average.

🔍 Interactive trimming (Analysis T range + Remove data)
Below the plots there is a Temperature ranges block.
1️⃣ Analysis T range
Top row:
	•	Sliders: Tmin and Tmax
	•	Boxes: editable Tmin/Tmax
	•	Button: Remove data
Purpose:
	•	Restrict the analyzed/exported range per well.
	•	Useful for removing: ◦ Cold noisy regions (<30 °C) ◦ Hot saturated regions (>90 °C)
What Remove data does:
	•	Saves the range [Tmin, Tmax] as the active analysis range for that well.
	•	Doesn’t delete original data; just uses a cropped working view.
	•	Affects: ◦ Tm calculation ◦ Export (corrected / smoothed)
	•	Fully reversible with Undo/Redo.

✂ Auto-trim to an expected Tm range
2️⃣ Expected Tm range (for auto-trim)
Bottom row of Temperature ranges:
	•	Boxes: Expected Tm range min / max
	•	Button: Auto-trim to expected range
⚠️ Important:
	•	Expected Tm range NO LONGER directly restricts Tm calculation.
	•	Tm is always computed as the global max of the derivative (within trimmed data).
	•	Expected Tm range is only used to propose auto-trim suggestions.
How Auto-trim works:
	1	Enter a reasonable Tm range (e.g., 50–65 °C).
	2	Press Auto-trim to expected range.
	3	Program:
	◦	Computes current Tm per well (ignoring deleted wells).
	◦	If Tm outside [min, max], it proposes trims that: ▪ Remove points from cold side, hot side, or both ▪ Keep at least 3 points ▪ Force Tm into range
	4	A dialog opens with the proposals: Columns:
	◦	Apply (checkbox)
	◦	Well
	◦	Δlow (°C)
	◦	Δhigh (°C)
	◦	Tm before
	◦	Tm after
	◦	ΔTm
	5	Controls:
	◦	Select none
	◦	Select all
	◦	Apply selected
	◦	Cancel
	6	When applying selected proposals:
	◦	Saves previous state (Undo compatible).
	◦	Trims curve to [Tmin_new, Tmax_new].
	◦	Marks well as auto-trimmed (adds ✂).
	◦	Ensures the well enters Corrected.
	◦	Recalculates Tm and list statuses.
	7	Undo per well fully reverts the auto-trim.
Deleted wells are excluded from proposals.

🌡️ Derivative & Tm — real mechanics For each non-deleted well:
	1	Use its corrected and trimmed curve (Remove data and/or Auto-trim).
	2	Sort points by temperature.
	3	Apply minimum smoothing to derivative + optional extra if Smooth ON.
	4	Compute derivative and force it upright:
	◦	Find the strongest peak
	◦	Flip sign if needed so that it becomes positive
	5	Normalize derivative (max ≈ 1).
	6	Tm = temperature at which the oriented derivative is maximal.
	7	Store Tm for:
	◦	All wells list
	◦	Tm outliers list
	◦	Tm export table
On screen:
	•	Derivative plot in the lower panel
	•	Dotted vertical line at Tm
	•	In the main plot: “Tm = XX.XX °C”
	◦	Black if normal
	◦	Red if marked as Tm outlier

🔥 Tm outlier handling
1️⃣ All wells (with data)
Each line looks like:
	•	A01 — Tm=57.32 °C
	•	If deleted:
	◦	A01 [DELETED] (light gray)
	•	If trimmed (Remove data or Auto-trim):
	◦	✂ appears at the end
	•	If Tm outlier:
	◦	Text in red
2️⃣ Tm outliers
Lists wells satisfying: |Tm − Tm_ref| ≥ threshold
	•	Tm_ref: ◦ If "I know my Tm" is filled → that value ◦ If empty → mean Tm of valid wells
	•	Deleted wells are excluded
	•	Trimmed wells show ✂
	•	Header shows count: Tm outliers (n)
Interactions:
	•	Click → jump to well
	•	Double-click → same
	•	Use Review mode (Prev/Next Tm outlier)
3️⃣ Tm controls
	•	Tm outlier threshold (°C)
	•	I know my Tm (°C)
	•	Current mean Tm: XX.XX °C (n=…)

🎨 Review Mode Third top bar:
	•	Prev/Next Suspected
	•	Prev/Next Corrected
	•	Prev/Next Tm outlier
Lets you cycle only through the relevant wells. Status bar shows e.g.: Tm outlier 3/7: B05

🟩 Color system and marks
In All wells:
	•	Corrected → light green
	•	Suspected → amber
	•	Normal → white
	•	Tm outlier → red text
	•	Deleted → light gray + [DELETED]
Extra marks:
	•	✂ → trimmed (Remove data or Auto-trim)
In the plot:
	•	Black Tm label → normal
	•	Red Tm label → Tm outlier

🗑️ Delete Well / Recover Well
Under the Suspected list:
	•	Delete Well
	•	Recover Well
Behavior:
	•	Delete Well: ◦ Logically removes the well ◦ Displays as A01 [DELETED] ◦ Effects: ▪ No Tm calculation ▪ No outlier status ▪ Excluded from exports ▪ Original curve still visible ▪ Corrected curve disappears
	•	Recover Well: ◦ Restores the well ◦ Included in Tm / outliers / exports again ◦ Restores its previous corrected curve
Delete/Recover are not Undo/Redo actions.

⌨️ Useful shortcuts
	•	Undo: Ctrl+Z / Cmd+Z
	•	Redo: Ctrl+Y / Ctrl+Shift+Z / Cmd+Shift+Z
	•	Next suspected: s
	•	Prev suspected: Shift+S
	•	Next corrected: c
	•	Prev corrected: Shift+C
	•	Move breakpoint: ← / →
	•	Move well selection: ↑ / ↓
	•	Apply typed index: Enter in index box

💾 Export options
	1	Export corrected
	◦	Corrected curves without export-time smoothing
	◦	Respects: ▪ Step corrections ▪ Remove data trims ▪ Auto-trim trims ▪ Excludes deleted wells
	◦	Formats: .gdsf, .tsv, .txt
	2	Export corrected + smoothed
	◦	Same as above but with export-smoothing
	◦	Applies smoothing using current slider value even if Smooth OFF
	◦	Ideal for ML / downstream fitting
	3	Export Tm table Includes: ▪ Well ▪ Tm_corrected ▪ Tm_smoothed ▪ Smooth_strength Excludes deleted wells Formats: .tsv or .csv

🧠 Key concepts (MAD, k, thresholds) For each well:
	•	Compute diffs: diffs[i] = y[i+1] - y[i]
	•	Dispersion: ◦ MAD (recommended): ▪ Robust to outliers ▪ σᵣ = 1.4826 × MAD(diffs) ◦ STD: standard deviation
	•	Thresholds: ◦ Abs threshold → jump > fixed value ◦ k·disp → jump > k × σ
Typically k ≈ 6 with MAD works well.

🛠️ Practical tips
	•	Fix steps first; refine Tm later.
	•	MAD + k=6 = great starting point.
	•	Too many suspects → increase k or Abs threshold.
	•	Missing clear jumps → lower k.
	•	Noisy derivative: ◦ Adjust Analysis T range ◦ Enable Smooth and increase strength
	•	Two melting phases: ◦ Use Auto-trim to focus on the relevant phase
	•	Use Tm outliers to find problematic wells fast.
	•	If you know the “correct” Tm (e.g., 57 °C): ◦ Put 57 in I know my Tm (°C) ◦ Adjust outlier threshold

🐛 Troubleshooting
	•	❌ “Derivative changes even when smoothing is OFF” → Not anymore: with Smooth OFF derivative uses fixed minimal smoothing.
	•	❌ “Tm peak appears downward” → Program auto-orients derivative upward.
	•	❌ “Tm jumps to 20 °C or 95 °C randomly” Check: ◦ Missing step corrections ◦ Bad Analysis T range ◦ Over-aggressive auto-trim ◦ Whether it's a Tm outlier (likely)
	•	❌ “Undo does nothing” → Undo only works after some action has been done in that well.
	•	❌ “Export corrected loses points” → Because Remove data or Auto-trim trimmed them. Export respects the working state.
