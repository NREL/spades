#!/usr/bin/env bash

unset -v EXECPATH
unset -v FPATH
NITERS=1000

set -e
cmd() {
  echo "+ $@"
  eval "$@"
}

usage () {
    cat <<-END
Usage:
------
   -e
     executable path (required)
   -f
     input file path (required)
   -n
     number of seeds to test (optional, default $NITERS)
   -h
     Print the usage information

   test_seeds.sh -e <executable path> -f <input file> -n <number of seeds>
END
}

while getopts e:f:n:h flag
do
    case "${flag}" in
        e)
            EXECPATH=${OPTARG}
            ;;
        f)
            FPATH=${OPTARG}
            ;;
        n)
            NITERS=${OPTARG}
            ;;
        h)
            usage
            exit 0
            ;;
        '?')
            echo "INVALID OPTION -- ${OPTARG}" >&2
            usage
            exit 1
            ;;
        ':')
            echo "MISSING ARGUMENT for option -- ${OPTARG}" >&2
            usage
            exit 1
            ;;
        *)
            echo "UNIMPLEMENTED OPTION -- ${flag}" >&2
            usage
            exit 1
            ;;
    esac
done

if [ -z "${EXECPATH}" ] ; then
    echo "Missing -e argument" >&2
    usage
    exit 1
fi

if [ ! -f "${EXECPATH}" ]; then
    echo "Executable ${EXECPATH} not found!"
    usage
    exit 1
fi

if [ -z "${FPATH}" ] ; then
    echo "Missing -f argument" >&2
    usage
    exit 1
fi

if [ ! -f "${FPATH}" ]; then
    echo "Input file ${FPATH} not found!"
    usage
    exit 1
fi

echo "Run starting at $(date '+%Y-%m-%d %H:%M:%S')"

LOGDIR="logs"
cmd "rm -rf ${LOGDIR}"
cmd "mkdir -p ${LOGDIR}"

SEED=42
RANDOM=$SEED

for ((i = 0 ; i < NITERS ; i++ )); do
    SPADES_SEED=$RANDOM
    LOGPATH="${LOGDIR}/out_${SPADES_SEED}.log"
    cmd "${EXECPATH} ${FPATH} spades.seed=${SPADES_SEED} > ${LOGPATH}"
done

echo "Run ended successfully at $(date '+%Y-%m-%d %H:%M:%S')"
