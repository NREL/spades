#!/usr/bin/env bash

set -e
cmd() {
  echo "+ $@"
  eval "$@"
}

echo "Run starting at $(date '+%Y-%m-%d %H:%M:%S')"

EXECDIR="../Build"
EXECNAME="spades"
EXECPATH="${EXECDIR}/${EXECNAME}"

FNAME="phold.inp"
FPATH="${EXECDIR}/${FNAME}"

LOGDIR="logs"
cmd "rm -rf ${LOGDIR}"
cmd "mkdir -p ${LOGDIR}"

ARGS="max_step=10000 amr.n_cell=128 128 spades.write_particles=f amr.max_grid_size=128 amr.plot_int=-1 amr.chk_int=-1"

SEED=42
RANDOM=$SEED

NITERS=10000
for ((i = 0 ; i < NITERS ; i++ )); do
    SPADES_SEED=$RANDOM
    LOGPATH="${LOGDIR}/log_${SPADES_SEED}.out"
    cmd "${EXECPATH} ${FPATH} ${ARGS} spades.seed=${SPADES_SEED} > ${LOGPATH}"
done

echo "Run ended successfully at $(date '+%Y-%m-%d %H:%M:%S')"
