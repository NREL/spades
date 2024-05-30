#!/usr/bin/env bash

set -e
cmd() {
  echo "+ $@"
  eval "$@"
}

unset -v EXECPATH
unset -v FPATH

while getopts e:f: flag
do
    case "${flag}" in
        e)
            EXECPATH=${OPTARG}
            ;;
        f)
            FPATH=${OPTARG}
            ;;
        '?')
            echo "INVALID OPTION -- ${OPTARG}" >&2
            exit 1
            ;;
        ':')
            echo "MISSING ARGUMENT for option -- ${OPTARG}" >&2
            exit 1
            ;;
        *)
            echo "UNIMPLEMENTED OPTION -- ${flag}" >&2
            exit 1
            ;;
    esac
done

if [ -z "${EXECPATH}" ] ; then
    echo "Missing -e argument" >&2
    exit 1
fi

if [ ! -f "${EXECPATH}" ]; then
    echo "Executable ${EXECPATH} not found!"
    exit 1
fi

if [ -z "${FPATH}" ] ; then
    echo "Missing -f argument" >&2
    exit 1
fi

if [ ! -f "${FPATH}" ]; then
    echo "Input file ${FPATH} not found!"
    exit 1
fi

echo "Run starting at $(date '+%Y-%m-%d %H:%M:%S')"

LOGDIR="logs"
cmd "rm -rf ${LOGDIR}"
cmd "mkdir -p ${LOGDIR}"

SEED=42
RANDOM=$SEED

NITERS=2
for ((i = 0 ; i < NITERS ; i++ )); do
    SPADES_SEED=$RANDOM
    LOGPATH="${LOGDIR}/log_${SPADES_SEED}.out"
    cmd "${EXECPATH} ${FPATH} spades.seed=${SPADES_SEED} > ${LOGPATH}"
done

echo "Run ended successfully at $(date '+%Y-%m-%d %H:%M:%S')"
