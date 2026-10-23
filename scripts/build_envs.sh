#!/usr/bin/env bash
# Build every pinned env/*.yml into env/prefix/<name> in parallel (creating them lazily one by one inside Nextflow is very slow on this FUSE disk).
# Resolved package lists are saved to docs/env_lock/<name>.txt.
B=/workspace/Biomarker_pipeline
export MAMBA_ROOT_PREFIX=$B/env/mamba
mkdir -p $B/env/prefix $B/logs/env_build $B/docs/env_lock
build() {
  y=$1; n=$(basename $y .yml); B=/workspace/Biomarker_pipeline
  if [ -d $B/env/prefix/$n/conda-meta ]; then echo "$n exists"; return; fi
  start=$(date +%s)
  $B/env/bin/micromamba create -y -p $B/env/prefix/$n -f $y > $B/logs/env_build/$n.log 2>&1 && \
    $B/env/bin/micromamba list -p $B/env/prefix/$n > $B/docs/env_lock/$n.txt && echo "$n built in $(( $(date +%s) - start )) s" || echo "$n FAILED (see logs/env_build/$n.log)"
}
export -f build
ls $B/env/*.yml | xargs -P 8 -I{} bash -c 'build {}'
