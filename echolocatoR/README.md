
### Import environment
#### For the same platform (e.g. Linux --> Linux)
`conda env create -f echolocatoR/echo_env.yml`
#### Across platforms (e.g. Linux --> MacOS)
NOTE: You need to latest version of conda to use the `--from-history` flag.
`conda env export --from-history`

### Export environment
`conda env export --no-builds > echolocatoR/echo_env.yml`


# Installation on Linux

```
# create new conda environment with all non-R dependencies
conda env create -f echolocator_no_R_packages.yml

# needed to get R working
conda install -c conda-forge ncurses=6.1 

# install R package

Rscript install_echolocatoR_R_packages.R
```
