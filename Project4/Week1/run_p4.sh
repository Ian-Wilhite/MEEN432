#!/usr/bin/env bash
# run_p4.sh — run P4 simulation then generate all figures
#
# Usage:
#   ./run_p4.sh           # full pipeline: sim + figures
#   ./run_p4.sh --figs    # figures only (loads existing simdata.mat, skips sim)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

MATLAB="matlab"
NOISE_FILTER='Warning|serialize|Permission denied|Make sure|write-protected|Error opening|for mode|In Stateflow|In construct_module|In generate_code|In sfc|In targetman|In autobuild_|In slsf|_jitprj'

run_matlab() {
    local script="$1"
    local label="$2"
    echo ""
    echo "══════════════════════════════════════════"
    echo "  $label"
    echo "══════════════════════════════════════════"
    $MATLAB -batch "$script" 2>&1 | grep -Ev "$NOISE_FILTER" | grep -v '^$' | grep -v '^\['
    echo "  ✓ done"
}

FIGS_ONLY=false
for arg in "$@"; do
    [[ "$arg" == "--figs" ]] && FIGS_ONLY=true
done

START=$(date +%s)

if [ "$FIGS_ONLY" = false ]; then
    run_matlab "p4_runsim" "Running simulation (p4_runsim.m)"
fi

run_matlab "gen_all_figs" "Generating figures (gen_all_figs.m)"

END=$(date +%s)
ELAPSED=$(( END - START ))
echo ""
echo "══════════════════════════════════════════"
printf "  Total wall time: %dm %ds\n" $(( ELAPSED/60 )) $(( ELAPSED%60 ))
echo "  Figures in: $SCRIPT_DIR/figures/"
echo "══════════════════════════════════════════"
