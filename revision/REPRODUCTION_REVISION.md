# This is the Source Code and Data Repository of the SynGPS and BOYA Algorithms

1.  [Cite Us](#cite-us)
2.  [Prerequisite](#prerequisite)
3.  [Installation](#installation)
4.  [Reproducing the Experiments of Revision](#reproducing-the-experiments-of-revision)
   
## Cite Us

Zahoránszky-Kőhalmi G, Walker B, Cole W, Miller N, Yang B, Vorontcov I, et al. Synthesis Route Identification and Prioritization in Reaction Knowledge Graphs. ChemRxiv. 2025; doi:10.26434/chemrxiv-2025-0s3jp This content is a preprint and has not been peer-reviewed.

[Link to preprint](https://chemrxiv.org/engage/chemrxiv/article-details/683f330a1a8f9bdab576d9ea)

## Prerequisite

Conda Python environment manager:

- For installation on your specific platform please refer to [Conda official documentation](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).


- [Git](https://git-scm.com/)

- [Git LFS](https://git-lfs.com/) — required to download the MongoDB dump and Memgraph snapshot (large binary files stored via LFS)

## Installation



### Clone this Git repository

```bash
git clone https://github.com/ncats/syngps
```

Enter into the subdirectory containing the source code.

```bash
cd syngps
```

```bash
git checkout publication
```

### Download Large Database Artifacts (Git LFS)

The MongoDB dump and Memgraph snapshot are stored in Git LFS. After cloning, pull them with:

```bash
git lfs install
git lfs pull
```

This will populate:

- `revision/data/input/mongo_dump/syngps-mongo.gz`
- `revision/data/input/memgraph_snapshot/`



### Download Required Input Files



```bash
cd revision/src
```

```bash
bash 01_download_and_process_input.sh.sh
```

The script will put the input files into `revision/data/input` (relative to the repository root).



### Create Environments

The SynGPS code can be installed as a Python package using the following command:


```bash
conda env create -f environment_1.yml
```

```bash
conda env create -f environment_2.yml
```

```bash
conda env create -f environment_3.yml
```



### Run the App

In order to reproduce the experiments, we need to run the evidence-based route search engine powered by SynGPS and BOYA algorithms, i.e.: "AICP light", AICP: ASPIRE Integrated Computational Platform. The following command will start the app and services in Docker containers.

```bash
# Copy env file
cp app/.env.sample app/.env

# Build custom app
docker compose -f docker-compose.yml build

# Run app and services in Docker containers
docker compose -f docker-compose.yml up
```



## Reproducing the Experiments of Revision

At this point we assume the input files have been downloaded, all conda environments have been created and the app is running as described above.
<BR></BR>
The experiments related to revision can be reproduced via a series of Jupyter Notebooks located at `revision/src/` (relative to the repository root). The notebooks are numbered incrementally, starting from 02. They need to be executed in the order of the prefix, i.e. 02, 03, 04 ...

<BR></BR>
Due to conflicting dependencies, we needed to create three environemnts (as described above) so that all experiments can be performed. NOTE, that at the top of the import section of each Jupyetr NoteBook we indicate in a comment line which environment needs to be used in connection with that specific JupyterNotebook.

Results of the experiments will be placed into `revision/data/output` (relative to the repository root). We provide a brief description of each output file in `revision/data/MANIFEST.md` (relative to the repository root).
