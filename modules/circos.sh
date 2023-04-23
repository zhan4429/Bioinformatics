#!/bin/bash

VERSION=0.69-9
PACKAGE=circos
TOOL=circos
DIRECTORY=$(dirname "$0")

STORAGES=(/ocean /hive /bil /local)

OPTIONS=""
for STORAGE in "${STORAGES[@]}"
do
  	if [ -d "$STORAGE" ]; then
                OPTIONS=$OPTIONS" -B $STORAGE"
        fi
done

singularity exec $OPTIONS $DIRECTORY/singularity-$PACKAGE-$VERSION.sif $TOOL "$@"