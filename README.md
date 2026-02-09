
# Computing_Accum — Finding accumulation rates from OIB Snow Radar

This README describes how to: 

   (1) download NASA Operation IceBridge (OIB) Snow Radar echograms, \
   (2) preprocess them in MATLAB, \
   (3) pick layers, and \
   (4) compute relative accumulation rates from the picked layers.

Currently, the full code is available upon request.

---

## Overview of workflow

1. **Download OIB Snow Radar data** (NetCDF echograms) from NSIDC for a given flight date.
2. **Convert and preprocess** the echograms in MATLAB:
   - NetCDF → MAT
   - Surface normalization
   - Along-track smoothing
   - Merge into ~25 km segments
3. **Pick layers** (manual/assisted UI).
4. **Compute relative accumulation** and inspect outputs.

---

## Requirements

- MATLAB (tested with standard MATLAB toolboxes used by the scripts)
- An **Earthdata Login** (required by NSIDC to download data)
- This codebase on your machine, with the included `SAMPLE_DATA/` directory structure as a template

---

## Step 1 — Download OIB Snow Radar data (before opening MATLAB)

Download OIB Snow Radar echogram files from NSIDC:

`https://n5eil01u.ecs.nsidc.org/ICEBRIDGE/IRSNO1B.002/{DATE}`

- Replace `{DATE}` with the flight date folder used by NSIDC (example: `2009.10.18`)
- You will need an Earthdata login.


---

## Step 2 — Preprocess data in MATLAB

### 2.1 Define ROOTPATH and date
**ROOTPATH** should point to the *parent directory containing date folders*, not the date folder itself.

Example:
- If your files are in:
  `/Users/yourname/Desktop/Computing_Accum/SAMPLE_DATA/2009.10.18/`
- Then set:
  `ROOTPATH = '/Users/yourname/Desktop/Computing_Accum/SAMPLE_DATA/'`

### 2.2 Run preprocessing commands
In the MATLAB command window (or in a `.m` script), run:

```matlab
ROOTPATH = '/Users/yourname/Desktop/Computing_Accum/SAMPLE_DATA/';
% Location of OIB NetCDF echogram files on your computer

date = '2009.10.18';
% Flight date folder you are processing

rad_nc2mat(date, ROOTPATH);
% Converts NetCDF echogram files to MATLAB .mat files

radar_processing(date, 100, 9, ROOTPATH);
% Normalizes surface to 0 and applies along-track smoothing
% Arguments here reflect:
%   100 = along-track spacing (m) used downstream
%   9   = smoothing parameter 

comb_echo(date, ROOTPATH);
% Merges processed echograms into ~25 km segments
````

### 2.3 What to check after Step 2

After these steps, you should see processed outputs and 25 km “patched” echograms created under the date directory (exact subfolders depend on your functions, but typically a `comb/` directory is created later).

**If you do not see new files:**

* Confirm `ROOTPATH` ends with `/SAMPLE_DATA/` (or your equivalent parent folder).
* Confirm the date string exactly matches the folder name.
* Confirm the `.nc` files are actually present under `{ROOTPATH}/{date}/`.

---

## Step 3 — Pick layers (manual UI)

Layer picking is the most time-consuming step.

### 3.1 Launch the picker

In MATLAB command window:

```matlab
tmpl_picker_man_auto_brooke
```

Follow the prompts. The automated option only works well when the layer is extremely clear; manual picking is often more reliable.

**Note:** If your project requires automated multi-layer picking, your workflow will likely diverge here (this pipeline assumes interactive picking of layers).

---

## Step 4 — Calculate relative accumulation

Once layer picks exist, compute relative accumulation:

```matlab
ROOTPATH = '/Users/yourname/Desktop/Computing_Accum/SAMPLE_DATA/';
date = '2009.10.18';

cal_rel_accum(date, ROOTPATH);
```

---

## Outputs: where to find results and what they contain

You should now have processed `.mat` files for each ~25 km segment located in:

```
{ROOTPATH}/{DATE}/comb/
```

For each segment file, you should see variables such as:

* `lay`
  A 1-D array of **two-way travel times (TWTT)** for your layer picks
  Typically one pick per **100 m** along-track sample.
* `accR`
  Accumulation rates associated with those picks (relative accumulation output from `cal_rel_accum`).


---

## Common pitfalls

* **Wrong ROOTPATH**: ROOTPATH should be the parent folder containing the date folder.
* **Date mismatch**: The `date` string must exactly match the directory name.
* **Missing credentials for NSIDC downloads**: verify Earthdata login and cookie/`wget` configuration.
* **Picker issues**: if the IRH is discontinuous, the UI may be frustrating; consider choosing a clearer horizon or excluding problematic segments.

---

## Citation

Dattler, M. E., Lenaerts, J. T., & Medley, B. (2019). Significant spatial variability in radar‐derived west Antarctic accumulation linked to surface winds and topography. Geophysical Research Letters, 46(22), 13126-13134.

