#!/bin/bash
set -e -x

BASEDIR=/raid/users/imai/Work/dev/tagged_moisture_model/
EXPNAME=test
USEDATA=era5
DATADIR="$BASEDIR"/Pre/"$USEDATA"/"$EXPNAME"

OUTDATA="$DATADIR"/merged.nc

cdo merge "$DATADIR"/tcwv.nc "$DATADIR"/prec.nc "$DATADIR"/evap.nc "$DATADIR"/viwve.nc "$DATADIR"/viwvn.nc $OUTDATA

echo "Pre-processing complete: $OUTPUT is ready."
