
import datetime as dt
import tkinter as tk
from tkinter import ttk, filedialog
from tkcalendar import DateEntry


def get_inputs(auto):
    """Retrieve pipeline inputs using default values or user input.

    Parameters
    ----------
    auto : bool
        Whether to use predefined input values (True) or prompt
        the user to enter values through the GUI (False).

    Returns
    -------
    probe : str
        THEMIS probe identifier: "a", "b", "c", "d", or "e".
    trange : list of str
        Start and end dates for data retrieval in YYYY-MM-DD format.
    time : datetime.datetime
        Selected UTC date and time for the field line calculation.
    E : float
        Electron kinetic energy in keV.
    theta : float
        Wave-normal angle in radians.
    maxstep : float
        Maximum field-line tracing step size in Earth radii.
    output_csv : str or None
        File location of final CSV file
    """

    if auto:
        time_range = ['2019-09-29', '2019-09-30']
        max_step = .05
        probe = 'd'
        selected_time = dt.datetime(2019,9,29,14,58,0,tzinfo=dt.UTC)
        theta = 0  # rad
        energy_kev = 200  # keV

        output_csv = None

        return probe,time_range,selected_time, energy_kev, theta, max_step, output_csv


    else:
        inputs = promptUser()

        if not inputs:
            raise SystemExit("Input window closed; pipeline not run.")

        start_date = dt.datetime.strptime(inputs["date"], "%Y-%m-%d").date()
        end_date = start_date + dt.timedelta(days=1)

        time_range = [start_date.strftime("%Y-%m-%d"), end_date.strftime("%Y-%m-%d")]

        date = dt.datetime(start_date.year, start_date.month, start_date.day, inputs["hour"],inputs["minute"],inputs["second"],tzinfo=dt.UTC)

        csv_location = inputs["CSV"]

        return inputs["probe"], time_range, date, inputs["E"],inputs["theta"], inputs["maxstep"], csv_location






def promptUser():
    """Display a GUI to collect user inputs for the T96 pipeline.

    Returns
    -------
    inputs : dict
        Dictionary containing the selected pipeline parameters:
        - probe : str
            THEMIS probe identifier: "a", "b", "c", "d", or "e".
        - date : str
            Selected date in YYYY-MM-DD format.
        - hour : int
            Selected UTC hour (0-23).
        - minute : int
            Selected UTC minute (0-59).
        - second : int
            Selected UTC second (0-59).
        - maxstep : float
            Maximum field-line tracing step size in Earth radii.
        - theta : float
            Wave-normal angle in radians.
        - E : float
            Electron kinetic energy in keV.
    """

    inputs = {}

    def submit():
        inputs["probe"] = probe_var.get()
        inputs["date"] = date_var.get()
        inputs["hour"] = int(hour_var.get())
        inputs["minute"] = int(minute_var.get())
        inputs["second"] = int(second_var.get())
        inputs["maxstep"] = float(maxstep_var.get())
        inputs["theta"] = float(theta_var.get())
        inputs["E"] = float(E_var.get())
        inputs["CSV"] = csv_var.get().strip() or None
        root.destroy()

    def browse_csv():
        path = filedialog.asksaveasfilename(
            parent=root,
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv")],
            initialfile="results.csv",
        )
        if path:
            csv_var.set(path)

    root = tk.Tk()
    root.title("T96 Pipeline Inputs")

    probe_var = tk.StringVar(value="d")
    date_var = tk.StringVar(value="2019-09-29")
    hour_var = tk.StringVar(value="14")
    minute_var = tk.StringVar(value="58")
    second_var = tk.StringVar(value="00")
    maxstep_var = tk.StringVar(value="0.05")
    theta_var = tk.StringVar(value="0")
    E_var = tk.StringVar(value="50")
    csv_var = tk.StringVar(value="")


    fields = [
        ("Probe", probe_var),
        ("Date", date_var),
        ("Chosen time UTC", None),
        ("Trace maxstep", maxstep_var),
        ("theta [rad]", theta_var),
        ("E [keV]", E_var),
        ("CSV file location", csv_var)
    ]

    for row, (label, var) in enumerate(fields):
        ttk.Label(root, text=label).grid(row=row, column=0, padx=8, pady=4, sticky="e")
        if label == "Probe":
            ttk.Combobox(
                root,
                textvariable=var,
                values=["a", "b", "c", "d", "e"],
                state="readonly",
                width=22
            ).grid(row=row, column=1, padx=8, pady=4)
        elif label == "Date":
            date_entry = DateEntry(
                root,
                textvariable=date_var,
                date_pattern="yyyy-mm-dd",
                width=22
            )
            date_entry.set_date(dt.date(2019, 9, 29))
            date_entry.grid(row=row, column=1, padx=8, pady=4)
        elif label == "Chosen time UTC":
            time_frame = ttk.Frame(root)
            time_frame.grid(row=row, column=1, padx=8, pady=4)

            ttk.Spinbox(
                time_frame,
                from_=0,
                to=23,
                textvariable=hour_var,
                width=4,
                format="%02.0f"
            ).grid(row=0, column=0)

            ttk.Label(time_frame, text=":").grid(row=0, column=1)

            ttk.Spinbox(
                time_frame,
                from_=0,
                to=59,
                textvariable=minute_var,
                width=4,
                format="%02.0f"
            ).grid(row=0, column=2)

            ttk.Label(time_frame, text=":").grid(row=0, column=3)

            ttk.Spinbox(
                time_frame,
                from_=0,
                to=59,
                textvariable=second_var,
                width=4,
                format="%02.0f"
            ).grid(row=0, column=4)

        elif label == "CSV file location":
            csv_frame = ttk.Frame(root)
            csv_frame.grid(row=row, column=1, padx=8, pady=4)
            ttk.Entry(csv_frame, textvariable=csv_var, width=18).grid(row=0, column=0)
            ttk.Button(csv_frame, text="Browse...", command=browse_csv).grid(row=0, column=1)
        else:
            ttk.Entry(root, textvariable=var, width=25).grid(row=row, column=1, padx=8, pady=4)

    ttk.Button(root, text="Run Pipeline", command=submit).grid(
        row=len(fields), column=0, columnspan=2, pady=10
    )

    root.mainloop()
    return inputs