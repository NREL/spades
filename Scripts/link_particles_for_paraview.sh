#!/usr/bin/env bash

usage() {
    cat <<-END
                Usage:
                ------
                   -n <name_of_particle_directory>
                     Link <name_of_particle_directory> to "particles"
                   -c
                     Remove links
                   -h
                     Link particle directories to "particles" for easy viewing in ParaView

                   link_particles_for_paraview.sh -n <name_of_particle_directory>
        END
}

unset -v PARTICLE_DIR_NAME
CWD="$(pwd)"
PARAVIEW_DIR_NAME="particles"
CLEAN="false"

while getopts n::ch flag; do
    case "${flag}" in
    n)
        PARTICLE_DIR_NAME="${OPTARG}"
        ;;
    c)
        CLEAN="true"
        ;;
    h)
        usage
        exit
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

if ${CLEAN}; then
    echo "Cleaning links"
    for DIRNAME in "${CWD}/plt"*; do
        echo "Unlinking ${DIRNAME}/${PARAVIEW_DIR_NAME}"
        rm "${DIRNAME}/${PARAVIEW_DIR_NAME}"
    done
    exit
fi

if [ -z "${PARTICLE_DIR_NAME}" ]; then
    echo "Missing -n"
    exit 1
fi

for DIRNAME in "${CWD}/plt"*; do
    if [ -d "${DIRNAME}/${PARAVIEW_DIR_NAME}" ]; then
        rm "${DIRNAME}/${PARAVIEW_DIR_NAME}"
    fi
    echo "Linking ${DIRNAME}/${PARTICLE_DIR_NAME} to ${DIRNAME}/${PARAVIEW_DIR_NAME}"
    ln -s "${DIRNAME}/${PARTICLE_DIR_NAME}" "${DIRNAME}/${PARAVIEW_DIR_NAME}"
done
