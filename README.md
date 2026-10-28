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
python T96pipeline3.py
```

An input window opens, already filled in with an example case (THEMIS-D, 2019-09-29 at 14:58 UTC). 
Change any values you like, optionally click **Browse...** to choose where to save the results as a CSV (leave it blank to skip saving), 
then click **Run Pipeline**. The first run may take longer due to the necessary downloading of data.


## 4. Optional: run from your own script
 
To skip the window, for example when running many cases, call the pipeline directly:
 
```python
import datetime as dt
from T96pipeline3 import run_pipeline
 
full_df, global_values = run_pipeline(
    probe="d",                                   # THEMIS probe: a, b, c, d or e
    time_range=["2019-09-29", "2019-09-30"],     # days of data to download
    selected_time=dt.datetime(2019, 9, 29, 14, 58, tzinfo=dt.UTC),
    energy_kev=200,                              # electron energy (keV)
    theta=0.0,                                   # wave-normal angle (radians)
    output_csv="results.csv",                    # optional: save results
)
```



## 5. Read the results
 
- `global_values` holds values that are constant along the field line: L-shell, equatorial field strength, equatorial cyclotron frequency, equatorial density and the Lorentz factor.
- `full_df` (and `results.csv`, if you asked for it) has one row per point along the field line for each wave frequency. Key columns:
  - `s`: distance from the magnetic equator (m)
  - `MLAT`: magnetic latitude (deg)
  - `omega`: wave frequency (rad/s)
  - `alpha_eq`: resonant equatorial pitch angle (deg)
  - `S`: inhomogeneity factor. |S| < 1 means nonlinear interaction is possible; |S| ≫ 1 means quasi-linear.

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


## Time Statistics

Estimated Runtime

Average time to download new data: 6.27 seconds
Average time to load previously downloaded data: 5.03 seconds
Time to complete program after downloading or loading the data: 4.31 seconds

Therefore, the program should take approximately 9.34–10.58 seconds to run fully, excluding the time taken to enter inputs
