# Conda environment for running the parsers

This repository follows advice outlined in [“Reproducible and upgradable Conda environments”.](https://pythonspeed.com/articles/conda-dependency-management/). In short, there are two environment files:
* [`environment.yml`](environment.yml) is maintained manually and contains only direct dependencies of the modules in this repository, with their versions pinned.
* [`environment-lock.yml`](environment-lock.yml) contains the complete locked environment, including all transitive dependencies, with their versions pinned. It is generated from the first file and must never be edited directly. This file can be used to reproduce the exact environment which was present at a particular point in repository history.

All commands below should be run from within `envs/`, the directory which contains this README file.

See also additional discussion in https://github.com/opentargets/platform/issues/1628.



## To reproduce the existing environment
This recreates the exact environment which was present during the time of a particular commit in the repository. Make sure to substitute the appropriate lock file for your platform (Linux/OSX).
```bash
conda env remove -n evidence_datasource_parsers_lock
conda env create --file environment-lock.yml
conda activate evidence_datasource_parsers_lock
```



## To update the environment
To ensure reproducibility, this must be done on a Linux machine only, because `conda env export` files are not cross-platform compatible. Use a local Linux installation or a Google Cloud instance for this.

Add or update any direct dependencies to [`environment.yml`](environment.yml) (Conda).
* Do not add any transitive dependencies.
* Do not forget to remove any dependencies which are no longer required.
* When a package is available through both Conda and PIP, install via Conda.

Recreate and activate the environment:
```bash
if [[ "${CONDA_DEFAULT_ENV}" == "evidence_datasource_parsers" ]]; then conda deactivate; fi
conda env remove -n evidence_datasource_parsers
conda env create -f environment.yml
conda activate evidence_datasource_parsers
```
* Resolve the version conflicts if any arise.
* Make sure that the new environment works not only for the functionality which you are adding, but for the unchanged modules as well.
* Re-run the test suite.
* Ideally, run some or all of the modules using the old and the new environment and compare the results.
* Never install any packages directly to this environment without also adding them to `environment.yml`.

Once the testing is complete, create a locked environment file:
```bash
conda env export --no-builds --file environment-lock.yml
sed -i '/^prefix/d' environment-lock.yml
sed -i 's/name: evidence_datasource_parsers/name: evidence_datasource_parsers_lock/' environment-lock.yml
```

Commit the changes to both `environment.yml` and `environment-lock.yml` alongside with your pull request.
