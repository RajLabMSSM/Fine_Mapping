
### Import environment
#### For the same platform (e.g. Linux --> Linux)
`conda env create -f echolocatoR/echo_env.yml`
#### Across platforms (e.g. Linux --> MacOS)
NOTE: You need to latest version of conda to use the `--from-history` flag.
`conda env export --from-history`

### Export environment
`conda env export --no-builds > echolocatoR/echo_env.yml`
