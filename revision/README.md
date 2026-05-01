# Revisions Readme

This directory contains the data, workflows, and scripts used to evaluate synthesis route ranking using the SyngPS library across a curated set of target molecules.

---

## Prerequisites

- [Conda](https://docs.conda.io/en/latest/miniconda.html) (Miniconda or Anaconda)
- Python 3.12 (managed by the conda environment)

---

## Environment Setup

A conda environment file is provided in `workflows/environment.yml`. It installs all required dependencies, including JupyterLab, RDKit, NetworkX, pandas, psutil, and a local editable install of the `syngps` library.

```bash
# From the revision/ directory
conda env create -f workflows/environment.yml

# Activate the environment
conda activate syngps
```

To update an existing environment after changes to `environment.yml`:

```bash
conda env update -f workflows/environment.yml --prune
```

---

## Downloads

Some of the following data has been included in the repository, but some of it is too large to include and must be downloaded separately. Please reach out to the team if you need access to any of the data.

- `simpretro_smiles.txt` — [https://github.com/catalystforyou/SimpRetro/blob/master/SMILES.txt](https://github.com/catalystforyou/SimpRetro/blob/master/SMILES.txt)
- `retrostar_subset_tms.tsv` — Generated TBD
- `retrostar_tms.txt` — TBD

To download the eMolecules sample inventory file (~1 GB), run the provided script from the `revision/` directory:

```bash
bash download_inventory_data.sh
```

This downloads `emolecules_inv.txt` into `revision/data/`.

---

## Data

Synthesis graphs for each target molecule are stored under `data/graphs/` in two sub-folders:

| Folder | Description |
|---|---|
| `leaves_as_sm/` | Graphs where all leaf nodes are assumed to be commercially available starting materials. All leaf intermediates are treated as in-inventory, so routes are rarely disqualified. |
| `with_inv_file/` | Graphs where availability is provided as a separate inventory file. Leaf intermediates are only considered available if they appear in the inventory with `available: true`. |

Each JSON file corresponds to a single target molecule and contains:

- `search_params` — the query parameters used to generate the graph
- `time_info` — server-side timing metadata
- `summary` — graph topology summary (node/edge counts, roles)
- `availability` — per-substance inventory lookup results (InChIKey keyed)
- `nodes` / `edges` — the flat synthesis graph in node-link format

---

## Understanding Viable Routes

A synthesis route is considered **viable** only if every leaf node in the route (substances with `srole='im'` and `in_degree=0`) is available. Availability is determined by checking:

1. `inventory.available: true` — found in the provided inventory file
2. `commercial_availability.available: true` — commercially purchasable

If **all** leaf intermediates in a route candidate are unavailable, the route is pruned down to just the target molecule and marked as **"Route Candidate Only Contains Target Molecule"** — it is disqualified and excluded from top-N ranking.

### Common disqualification causes

- The graph was generated with `leaves_as_sm: false` (i.e., `with_inv_file` graphs), but none of the required leaf intermediates appear in the inventory. This is the expected behaviour when the inventory does not cover the molecules needed by any reaction path.
- A combination graph is empty due to incompatible reactions being removed upstream.

### `leaves_as_sm` vs `with_inv_file`

| | `leaves_as_sm` | `with_inv_file` |
|---|---|---|
| Availability source | Leaf role assumed in-stock | External inventory file |
| Expected viable routes | Higher — leaf nodes always pass | Lower — depends on inventory coverage |
| Use case | Theoretical route enumeration | Realistic purchasability-constrained ranking |

---

## Running the Notebooks

Start JupyterLab from the `revision/` directory with the `syngps` environment active:

```bash
conda activate syngps
jupyter lab
```

Then open any notebook from the `workflows/` folder:

| Notebook | Description |
|---|---|
| `process_graphs_with_leaves_as_sm.ipynb` | Loads graphs from `data/graphs/leaves_as_sm/`, identifies routes, ranks by aggregated yield, and sends to Cytoscape. |
| `process_graphs_with_inv.ipynb` | Runs the full `find_top_n_routes` pipeline across all graphs in `data/graphs/with_inv_file/`, benchmarks each file, and produces a pandas DataFrame of topology, route, and performance metrics. |

### Benchmark output (`process_graphs_with_inv.ipynb`)

Cell 3 iterates all JSON files and collects per-file metrics into a DataFrame (cell 4) with the following columns:

| Column | Description |
|---|---|
| `nodes` / `edges` | Graph size |
| `leaf_ims_total` / `available` / `unavailable` | Leaf intermediate availability breakdown |
| `combination_graphs` | Number of combination graphs generated |
| `route_candidates` | Total route candidates evaluated |
| `viable_routes` | Routes passing availability check |
| `top_n_returned` | Routes returned after yield ranking (≤ `TOP_N`) |
| `build_time_s` | Time to load JSON and build `SynthGraph` |
| `identify_time_s` | Time for `identify_individual_synthesis_routes` |
| `find_top_n_time_s` | Time for `find_top_n_routes` (includes yield projection) |
| `total_time_s` | End-to-end wall-clock time per file |
| `cpu_user_s` / `cpu_system_s` | CPU time consumed |
| `peak_alloc_mb` | Peak heap allocation (tracemalloc) |
| `rss_delta_mb` | Change in process RSS during processing |

Cell 6 provides summary analysis: timing phase breakdown, viable vs. disqualified counts, size-vs-time correlation, and the 10 slowest graphs.
