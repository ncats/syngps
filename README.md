# This is the Source Code and Data Repository of the SynGPS and BOYA Algorithms

### Preprint

Citation...

1.  [Prerequisite](#prerequisite)
2.  [Installation](#installation)
3.  [Case Study Reproduction](#case-study-reproduction)
4.  [Post Publication Notes](#post-publication-notes)
5.  [For Contributors](#for-contributors)
6.  [License Related Section](#license-related-section)

## Prerequisite

Conda Python environment manager:

- For installation on your specific platform please refer to [Conda official documentation](https://conda.io/projects/conda/en/latest/user-guide/install/index.html).

Optional (not required, but nice to have):

- [Git](https://git-scm.com/)

## Installation

### Clone this Git repository

```bash
git clone https://github.com/ncats/syngps
```

Enter into the subdirectory containing the source code.

```bash
cd syngps
```

If you want to utilize the publication branch of this repository, you can check out the publication branch using the following commands:

```bash
git fetch origin
git checkout -b publication origin/publication
```

To ensure you have the latest updates for this repository, you can use the following command:

```bash
git pull
```

### Create the required Python environment via Conda

```bash
cd src
```

```bash
conda env create -f environment.yml
```

## Case Study Reproduction

Code and data required to reproduce the case studies is stored under subdirectories `syngps/src`, and `syngps/data`, respectively.

Run the workflow in Jupyter Notebook in JupyterLab using the above created `syngps` [Conda environment](#create-the-required-python-environment-via-conda):

```bash
src/syngps_cs1.ipynb
src/case_study2.ipynb
src/use_case_3_syngps_from_file.ipynb
```

## Post Publication Notes

To view the repository code at the time of publication view the `publication` branch.

## For Contributors

We welcome contributions to the SynGPS and BOYA algorithms! Here's how you can contribute:

1. **Create a Feature Branch**:

   - Start by creating a personal feature branch (e.g., `feat_new_1`) from the `feature` branch.

2. **Develop and Test**:

   - Implement your feature or fix.
   - Ensure your code is thoroughly tested and works as expected.

3. **Submit a Pull Request**:

   - Once your code is complete and tested, create a pull request (PR) from your feature branch to the `feature` branch of the SynGPS repository.

4. **Approval and Merge**:
   - Approved features will be merged into the `main` branch.

Thank you for contributing to the SynGPS repository!

## License Related Section

This repository contains source code, Jupyter notebook, data and results files which are organized into various subdirectories.

- Source code subdirectories

    `src/`

- Data and results subdirectory

    `data/`

- Worklows

    The workflow to reproduce Case Study 3 in the form of Jupyter Notebook can be found at `src/use_case_3_syngps_from_file.ipynb`.

- Source Code License of SynGPS Repository

    The applicable license to source code can be found under filename: [src/LICENSE](src/LICENSE) (license type: [MIT License](https://opensource.org/licenses/MIT)). This license is applicable to all files recursively in the source code subdirectories as defined above. The files [src/NOTES](src/NOTES) list source code modules that were utilized and their respective licenses. These modules have their own licenses which might be different from the Source Code License of this repository, and they need to be respected accordingly.

- Data License of SynGPS Repository

    The applicable license to data and results can be found under: [data/LICENSE](data/LICENSE) that is a [Creative Commons Attribution 4.0 International Public License](https://creativecommons.org/licenses/by/4.0/legalcode.txt). This license is applicable to all files recursively in the data and results subdirectory as defined above. The files listed in [data/NOTES](data/NOTES) and [plots/NOTES](plots/NOTES) lists input files and resources utilized to perform the experiments and are considered as derivative work of those resources. These input files and resources have their own licenses which might be different from the Data License of this repository, and they need to be respected accordingly. In the same file we also list which results files can be considered as derivative works, and we also list the the respective ascendent data source(s).

- Jupyter Notebook License of SynGPS Repository

    Jupyter Notebboks are special in the sense that they are comprised of source code, but they can also contain embedded data and plot (graph) sections. This duality is resolved via dual-licensing as follows. The code sections of Jupter Notebooks fall under the same license as source codes, i.e. the MIT License [src/LICENSE](src/LICENSE)whereas data and plot sections embedded into the Jupyter Notebooks fall under the same license as data and result files, i.e. Creative Commons Attribution 4.0 International Public License CC-BY 4.0 International License [data/LICENSE](data/LICENSE). Remarks enclosed in the [src/NOTES](src/NOTES) file are also valid for code section of the Jupyter Notebooks. Remarks enclosed in the [data/NOTES](data/NOTES) files are also valid for the embedded data and plots of the Jupyter Notebook files.

- Links to licenses

  - MIT License: [https://opensource.org/licenses/MIT](https://opensource.org/licenses/MIT)

  - Creative Commons Attribution 4.0 International Public License (CC-BY 4.0 International): [https://creativecommons.org/licenses/by/4.0/legalcode.txt](https://creativecommons.org/licenses/by/4.0/legalcode.txt)

  - Creative Commons Attribution-ShareAlike 3.0 Unported License: [https://creativecommons.org/licenses/by-sa/3.0/](https://creativecommons.org/licenses/by-sa/3.0/)

  - GPL-2 License: [https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

  - GPL-3 License: [https://www.gnu.org/licenses/gpl-3.0.en.html](https://www.gnu.org/licenses/gpl-3.0.en.html)

  - Apache License 2.0: [https://www.apache.org/licenses/LICENSE-2.0](https://www.apache.org/licenses/LICENSE-2.0)

  - 3-Clause BSD license: [https://opensource.org/licenses/BSD-3-Clause](https://opensource.org/licenses/BSD-3-Clause), [https://github.com/scikit-learn/scikit-learn/blob/master/COPYING](https://github.com/scikit-learn/scikit-learn/blob/master/COPYING)
