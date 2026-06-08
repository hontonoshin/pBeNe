#!/usr/bin/env bash
set -euo pipefail

cd ~/Downloads/pBN

mkdir -p macros results_scan logs_scan

NPROTONS=1000000

ENERGIES=(
  3.1
  4.0
  5.0
  6.0
  7.0
  8.0
  9.0
  10.0
  12.0
  14.0
  16.0
  18.0
  20.0
)

for E in "${ENERGIES[@]}"
do
    ESAFE=${E//./p}

    MACRO="macros/run_${ESAFE}MeV.mac"
    ROOTOUT="../results_scan/pBN_${ESAFE}MeV_${NPROTONS}p.root"
    LOGOUT="../logs_scan/pBN_${ESAFE}MeV_${NPROTONS}p.log"

    sed \
      -e "s/ENERGY_MEV/${E}/g" \
      -e "s/NPROTONS/${NPROTONS}/g" \
      macros/run_template.mac > "${MACRO}"

    echo "============================================================"
    echo "Running E = ${E} MeV, N = ${NPROTONS}"
    echo "ROOT: results_scan/pBN_${ESAFE}MeV_${NPROTONS}p.root"
    echo "LOG : logs_scan/pBN_${ESAFE}MeV_${NPROTONS}p.log"
    echo "============================================================"

    cd build

    PBN_OUTPUT_FILE="${ROOTOUT}" ./pBN "../${MACRO}" | tee "${LOGOUT}"

    cd ..
done

echo "Scan finished."
