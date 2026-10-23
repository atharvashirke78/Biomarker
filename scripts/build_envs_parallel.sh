#!/usr/bin/env bash
# Same as build_envs.sh but every env gets its OWN package cache, so builds do not serialize on the shared pkgs lock.
B=/workspace/Biomarker_pipeline
for n in "$@"; do
  ( export MAMBA_ROOT_PREFIX=$B/env/mamba CONDA_PKGS_DIRS=$B/env/mamba/pkgs_$n TMPDIR=$B/tmp
    start=$(date +%s)
    if $B/env/bin/micromamba create -y -p $B/env/prefix/$n -f $B/env/$n.yml > $B/logs/env_build/$n.log 2>&1; then
      $B/env/bin/micromamba list -p $B/env/prefix/$n > $B/docs/env_lock/$n.txt; echo "$n built in $(( $(date +%s) - start )) s"
    else echo "$n FAILED (see logs/env_build/$n.log)"; fi ) &
done
wait
