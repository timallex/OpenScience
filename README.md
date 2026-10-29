# T96 Nonlinear Wave-Particle Interaction Suite

Given a satellite probe, time range, selected time, electron energy, wave-normal angle, 
this module downloads THEMIS and OMNI data, traces a magnetic field line using the T96 model at the location 
of the satellite at the selected time on the selected date, and calculates numerous global 
quantities and quantities along the field line for each $\Omega$ ratio

## 1. Requirements

- Python 3.11 or newer
- An internet connection (data are downloaded on the first run and cached locally, so later runs are faster)

## 2. Install

All required libraries can be downloaded by opening a terminal and running:

```bash
pip install pyspedas geopack spacepy numpy pandas matplotlib tkcalendar
```

## 3. Run it



```bash
python t96_pipeline.py
```

An input window opens, already filled in with an example case (THEMIS-D, 2019-09-29 at 14:58 UTC).
Change any values you like, optionally click **Browse...** to choose where to save the results as a CSV (leave it blank to skip saving), 
then click **Run Pipeline**. The first run may take longer due to the necessary downloading of data.

If the required data aren't available for the chosen date and time, the program stops with an error such as
"Required data are missing for the selected time". Try a different date, time or probe.

## 4. Optional: run from your own script
 
To skip the window, for example when running many cases, call the pipeline directly:

```python
import datetime as dt
from t96_pipeline import run_pipeline

full_df, global_values = run_pipeline(
  probe="d",  # THEMIS probe: a, b, c, d or e
  time_range=["2019-09-29", "2019-09-30"],  # days of data to download
  selected_time=dt.datetime(2019, 9, 29, 14, 58, tzinfo=dt.UTC),
  energy_kev=200,  # electron energy (keV)
  theta=0.0,  # wave-normal angle (radians)
  output_csv=None,  # optional: save results
)
```



## 5. Read the results
 
- `global_values` holds values that are constant along the field line: L-shell, equatorial field strength, equatorial cyclotron frequency, equatorial density and the Lorentz factor.
- `full_df` (and `results.csv`, if you asked for it) has one row per point along the field line for each wave frequency. Key columns:
  - `s`: distance from the magnetic equator (m)
  - `MLAT`: magnetic latitude (deg)
  - `omega`: wave frequency (rad/s)
  - `alpha_eq`: resonant equatorial pitch angle (deg)
  - `S`: Omura's S factor

- Since 'full_df' contains the data for each wave frequency, the file contains 4 trace lines

- By default the code analyzes wave frequencies at 1/2, 1/5, 1/10 and 1/20 of the equatorial electron cyclotron frequency.
- To change them, edit `OMEGA_RATIOS` in `T96pipeline3.py`.


## Files
 
| File | What it does                                        |
|---|-----------------------------------------------------|
| `T96pipeline3.py` | Main script; start here                             |
| `T96inputHelper.py` | Default inputs and the input pop-up window          |
| `T96DataHelper.py` | Downloads THEMIS and OMNI data, prepares T96 inputs |
| `T96WaveHelper.py` | Turns measured wave spectra into wave amplitudes    |
| `T96TraceLineHelper.py` | Traces the field line and builds the results table  |
| `T96EquationHelper.py` | The physics equations                               |
| `T96plotsHelper.py` | Plot helpers                                        |


## Main functions

The pipeline runs these in order. Some steps read data that an earlier step stored inside PySPEDAS; this is noted under **Needs**.

### `get_inputs(auto)` — `T96inputHelper.py`
Gets the inputs from a pop-up window (`auto=False`) or uses preset example values (`auto=True`).
- **Returns:** `probe`, `time_range`, `selected_time`, `energy_kev`, `theta`, `max_step`, `output_csv`, ready to pass to `run_pipeline`

### `run_pipeline(probe, time_range, selected_time, energy_kev, theta, max_step=0.05, plot_enabled=False, no_update=False, output_csv=None)` — `t96_pipeline.py`
Runs the whole program.
- **Inputs:**
  - `probe`: THEMIS probe, `"a"`, `"b"`, `"c"`, `"d"` or `"e"`
  - `time_range`: start and end dates of data to download, `["YYYY-MM-DD", "YYYY-MM-DD"]`
  - `selected_time`: time of interest as a UTC `datetime`; the nearest spacecraft position is used
  - `energy_kev`: electron energy (keV)
  - `theta`: wave-normal angle (radians)
  - `max_step`: optional, field-line tracing step in Earth radii (default 0.05)
  - `plot_enabled`: optional, show plots (default `False`)
  - `no_update`: optional, use already-downloaded data only (default `False`)
  - `output_csv`: optional, path to save the results as a CSV (default: don't save)
- **Returns:** `full_df` (results table), `global_values` (dict of field-line constants)

### `load_required_data(probe, time_range, no_update)` — `T96DataHelper.py`
Downloads THEMIS position, THEMIS SCM wave spectra and OMNI data.
- **Inputs:** `probe`, `time_range` (`["YYYY-MM-DD", "YYYY-MM-DD"]`), `no_update` (bool)
- **Returns:** `scm2_times`, `scm2_psd`, `frequencies_hz`, `scm3_times`, `scm3_psd`, `scm3_frequencies` (times in Unix seconds, PSD in nT²/Hz, frequencies in Hz)

### `prepare_trace_seed(probe, selected_time_epoch, plot_enabled=False)` — `T96DataHelper.py`
Builds the T96 model inputs and picks the spacecraft position closest to the chosen time.
- **Inputs:** `probe`, `selected_time_epoch` (Unix seconds), `plot_enabled`
- **Returns:** `seed_time_epoch` (Unix seconds), `seed_parameters` (T96 model inputs at that time: solar wind pressure, Dst, IMF By and Bz)
- **Needs:** data loaded by `load_required_data`

### `get_bw_values(probe, seed_time_epoch, time_range, no_update, plot_enabled, scm2_times, scm2_psd, frequencies_hz, scm3_times, scm3_psd, scm3_frequencies)` — `T96WaveHelper.py`
Turns the two SCM spectra at the chosen time into wave amplitudes.
- **Inputs:** `seed_time_epoch` from `prepare_trace_seed`; the SCM data from `load_required_data`
- **Returns:** `frequencies_hz` (Hz), `measured_bw_values` (nT), `measured_omega_values` (rad/s)

### `trace_field_line(seed_parameters, max_step, plot_enabled=False)` — `T96TraceLineHelper.py`
Traces the T96 field line from the spacecraft to the equator and both ionospheres.
- **Inputs:** `seed_parameters` from `prepare_trace_seed`, `max_step` (R_E)
- **Returns:** `s_full` (km from equator), `r_full_gsm` (km), `b_full` (nT), `bmag_full` (nT)
- **Needs:** the seed position `seed_pos_gsm` stored by `prepare_trace_seed`

### `process_trace_results(s_full, r_full_gsm, b_full, bmag_full, plot_enabled=False)` — `T96TraceLineHelper.py`
Computes magnetic latitude and L-shell and builds the field-line table.
- **Inputs:** the four outputs of `trace_field_line`
- **Returns:** `df` (field-line table: `s` in m, `MLAT` in deg, positions in R_E, field in nT), `equator_time` (UTC datetime)
- **Needs:** the traces stored by `trace_field_line`

### `calculate_global_values(l_shell, equatorial_bmag, energy_kev)` — `T96EquationHelper.py`
Computes values that are constant along the field line.
- **Inputs:** `l_shell`, `equatorial_bmag` (nT), `energy_kev`
- **Returns:** `f_ce_eq` (Hz), `omega_ce_eq` (rad/s), `n_pe_eq` (cm⁻³), `gamma`

### `select_wave_values(measured_omega_values, measured_bw_values, omega_ce_eq, omega_ratios)` — `T96WaveHelper.py`
Picks wave frequencies as fractions of the equatorial cyclotron frequency and interpolates their amplitudes.
- **Inputs:** `measured_omega_values`, `measured_bw_values` from `get_bw_values`; `omega_ce_eq` from `calculate_global_values`; `omega_ratios` (e.g. `(1/2, 1/5, 1/10, 1/20)`)
- **Returns:** `omega_values` (rad/s), `bw_values` (nT)

### `add_trace_equation_values(df, n_pe_eq)` — `T96EquationHelper.py`
Calculates values that depend on field line position but not omega or Bw.
- **Inputs:** `df` from `process_trace_results`, `n_pe_eq` from `calculate_global_values`
- **Returns:** `df` with added columns `f_ce`, `omega_ce`, `dOce_ds`, `f_pe`, `omega_pe`

### `calculate_wave_results(df, omega_values, bw_values, theta, gamma, omega_ce_eq)` — `T96EquationHelper.py`
Computes values that depend on Bw and/or omega, including S, for each wave frequency.
- **Inputs:** `df` from `add_trace_equation_values`; `omega_values`, `bw_values` from `select_wave_values`; `gamma`, `omega_ce_eq` from `calculate_global_values`; `theta` (rad)
- **Returns:** `full_df`, one copy of the field line per frequency with added columns `k`, `N`, `alpha`, `alpha_eq`, `v_perp`, `v_R`, `S`

## Time Statistics

Estimated Runtime

Average time to download new data: 6.27 seconds
Average time to load previously downloaded data: 5.03 seconds
Time to complete program after downloading or loading the data: 4.31 seconds

Therefore, the program should take approximately 9.34–10.58 seconds to run fully, excluding the time taken to enter inputs
