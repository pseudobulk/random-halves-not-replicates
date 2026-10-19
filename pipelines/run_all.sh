#!/usr/bin/env bash
# The whole analysis, in the order it has to run. Paths come from config.env at the repository root.
#   bash pipelines/run_all.sh                 all stages in order
#   bash pipelines/run_all.sh figures         one or more named stages
# Stages (README, "Order of runs"):
#   extract_cd4  extract_tahoe  cd4  cd4_pseudobulk  celleval  tahoe  replogle  sim  collect  figures
# Every stage reads only the outputs of earlier stages; "figures" needs nothing but the small files in FIGDATA (or
# data/ in this repository, see README "Figures only").
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CFG="${RHNR_CONFIG:-$REPO/config.env}"
[ -f "$CFG" ] || { echo "missing $CFG: copy config.example.env to config.env and edit it"; exit 1; }
_env="$(export -p)"; set -a; . "$CFG"; set +a; eval "$_env"   # variables already exported take precedence
export RHNR_CONFIG="$CFG"
cd "$REPO/code"
TAHOE_ROOT="${TAHOE_ROOT:-$DATA_ROOT/tahoe}"; export TAHOE_ROOT
PY="${PYTHON:-python3}"; D="$DATA_ROOT"; W="$WORK"; N="${WORKERS:-16}"; RES="$W/results"
mkdir -p "$W/cd4" "$RES"
log() { echo "$(date '+%F %T') $*"; }

CD4_FILES="D1_Rest D1_Stim48hr D2_Rest D2_Stim48hr"

extract_cd4() {   # per-cell aggregates of the CD4 single cells: random halves, guide halves, 8-way cells, non-targeting guides
  for n in $CD4_FILES; do
    h5="$D/cd4_sc/$n.assigned_guide.h5ad"
    log "extract $n"
    $PY cd4_sc_extract.py       "$h5" "$W/cd4/${n}_random" "$N"
    $PY cd4_sc_extract_guide.py "$h5" "$W/cd4/${n}_guide"  "$N"
    $PY cd4_sc_extract_cells.py "$h5" "$W/cd4/${n}_cells"  "$N"
    $PY cd4_ntc_extract.py      "$h5" "$W/cd4/${n}_ntc"    "$N"
  done
}

extract_tahoe() { # per line x sample x random half aggregates over all Tahoe-100M shards (map, then reduce)
  $PY tahoe_extract.py map "$TAHOE_ROOT/Tahoe-100M/data" "$W/tahoe_agg" "$N"
  $PY tahoe_extract.py reduce "$W/tahoe_agg"
}

cd4() {           # single-cell analyses, the same for donors 1 and 2
  for d in D1 D2; do
    local R="$W/cd4/${d}_Rest_random" S="$W/cd4/${d}_Stim48hr_random" RG="$W/cd4/${d}_Rest_guide" SG="$W/cd4/${d}_Stim48hr_guide"
    local CR="$W/cd4/${d}_Rest_cells" CS="$W/cd4/${d}_Stim48hr_cells" O="$RES/$d" sfx=""
    [ "$d" = D2 ] && sfx="_D2"
    log "cd4 $d"
    $PY cd4_inert2.py "$R" "$S" "$O/cell_split"
    $PY cd4_inert2.py "$RG" "$SG" "$O/guide_split"
    $PY robust_cd4_tahoe.py cd4 "$R" "$S" "$RG" "$SG" "$O/robust"
    $PY cd4_kd_matched.py "$R" "$S" "$RG" "$SG" "$O/kd_matched"
    $PY readout_inert_line.py cd4 "$R" "$S" "$O/inert_line_cell"
    $PY readout_inert_line.py cd4 "$RG" "$SG" "$O/inert_line_guide"
    $PY ceiling_compare.py cd4 "$R" "$RG" "rest$sfx" "$O/ceiling_rest"
    $PY ceiling_compare.py cd4 "$S" "$SG" "stim48$sfx" "$O/ceiling_stim48"
    $PY clone_tests_cd4.py "$CR" rest "$CS" stim48 "$O/clone_tests"
    $PY ci_bootstrap.py cd4 "$CR" "$CS" "$O/ci"
    $PY export_fig_data.py "$CR" "$CS" "$O/cd4_$d"
    [ "$d" = D1 ] && $PY verify_headline.py cd4 "$CR" "$CS" "$O/verify"
  done
  mkdir -p "$RES/ntc"
  $PY cd4_ntc_null.py "$RES/ntc/metrics_ntc_null.json" \
    "D1_rest:$W/cd4/D1_Rest_ntc:$W/cd4/D1_Rest_cells"   "D1_stim:$W/cd4/D1_Stim48hr_ntc:$W/cd4/D1_Stim48hr_cells" \
    "D2_rest:$W/cd4/D2_Rest_ntc:$W/cd4/D2_Rest_cells"   "D2_stim:$W/cd4/D2_Stim48hr_ntc:$W/cd4/D2_Stim48hr_cells"
  $PY cd4_ntc_ratio.py "$RES/ntc/cd4_ntc_ratio.json" \
    "D1_rest:$W/cd4/D1_Rest_ntc:$W/cd4/D1_Rest_cells" "D2_rest:$W/cd4/D2_Rest_ntc:$W/cd4/D2_Rest_cells"
  mkdir -p "$RES/trbc"
  $PY cd4_trbc_state_keff.py "$RES/trbc/trbc_state_keff.json" \
    "D1:$W/cd4/D1_Rest_cells:$W/cd4/D1_Stim48hr_cells" "D2:$W/cd4/D2_Rest_cells:$W/cd4/D2_Stim48hr_cells"
}

cd4_pseudobulk() { # released four-donor pseudobulk: donor halves, same vs different guide across donors
  local PB="$D/cd4/GWCD4i.pseudobulk_merged.h5ad"
  $PY cd4_donor_split.py "$PB" "$RES/donor_split"
  $PY cd4_guide_donor.py "$PB" "$RES/guide_donor"
}

celleval() {      # cell-eval 0.8.2 itself on a donor-1 resting subset (separate environment)
  mkdir -p "$RES/celleval"
  $PY cd4_celleval_subset.py "$D/cd4_sc/D1_Rest.assigned_guide.h5ad" "$W/cd4/D1_Rest_cells" "$RES/celleval/cd4_d1_rest_subset.npz" "$N"
  "${CELLEVAL_PYTHON:-python3}" run_celleval_ceiling.py "$RES/celleval/cd4_d1_rest_subset.npz" "$RES/celleval/full" 48
}

tahoe() {         # replicate plates 6 / 14, DMSO wells, ceilings, predictor
  local T="$RES/tahoe"; mkdir -p "$T"
  $PY tahoe_plate_overlap.py > "$T/plate_overlap.txt"
  $PY tahoe_plate_recompute.py "$W/tahoe_agg" "$D/hpa/rna_celline.tsv" "$T/plate_recompute"
  $PY tahoe_dmso_null.py "$W/tahoe_agg" "$T/dmso"
  $PY ceiling_compare.py tahoe "$T/plate_recompute/D_arrays.npz" "$T/ceiling"
  $PY robust_cd4_tahoe.py tahoe "$T/plate_recompute/D_arrays.npz" "$T/robust"
  $PY ci_bootstrap.py tahoe "$T/plate_recompute/D_arrays.npz" "$T/ci"
  $PY verify_headline.py tahoe "$W/tahoe_agg" "$T/verify"
}

replogle() {      # K562 / RPE1 (Replogle et al. 2022) and Jurkat / HepG2 (Nadig et al. 2025)
  local K6="$D/replogle/K562_essential_raw_singlecell_01.h5ad" K8="$D/replogle/K562_gwps_raw_singlecell_01.h5ad"
  local RP="$D/replogle/rpe1_raw_singlecell_01.h5ad" O="$RES/replogle"
  $PY replogle_indep_readout.py "$K6" "$K8" "$RP" "$O/indep"
  $PY k562_inert_screens.py "$O/indep/responses.npz" "$K6" "$K8" "$O/k562_inert_screens"
  $PY k562_inert.py "$K8" "$O/k562_template"
  $PY readout_inert_line.py h5ad "$K6" "$RP" k562 rpe1 "$O/inert_line_k562_rpe1"
  mkdir -p "$O/cl_ntc"
  WORKERS="$N" $PY cl_ntc_null.py "$O/cl_ntc/cl_ntc_replogle.json" "k562_ess:$K6" "rpe1:$RP" "k562_gw1:$K8"
  WORKERS="$N" $PY cl_ntc_null.py "$O/cl_ntc/cl_ntc_nadig.json" "jurkat:$D/nadig/jurkat.h5ad" "hepg2:$D/nadig/hepg2.h5ad"
  WORKERS="$N" $PY cl_ntc_null.py "$O/cl_ntc/cl_ntc_h1.json" "h1:$D/vcc2025/adata_Training.h5ad"
  $PY readout_inert_line.py h5ad "$D/nadig/jurkat.h5ad" "$D/nadig/hepg2.h5ad" jurkat hepg2 "$O/inert_line_jh"
}

sim() {           # three-term unit model (simulation, no data)
  mkdir -p "$RES/sim"
  $PY clone_model_sim.py "$RES/sim/sim_v6.json"
}

collect() {       # copy the small result files the figure and table scripts read into FIGDATA
  local F="${FIGDATA:-$W/figdata}"; mkdir -p "$F"
  while read -r dst src; do
    [ -z "$dst" ] && continue
    cp "$RES/$src" "$F/$dst"
  done <<EOF
cd4_D1.npz              D1/cd4_D1.npz
cd4_D2.npz              D2/cd4_D2.npz
ceiling_D1_rest.json    D1/ceiling_rest/metrics_ceiling_compare.json
ceiling_D2_rest.json    D2/ceiling_rest/metrics_ceiling_compare.json
ceiling_tahoe.json      tahoe/ceiling/metrics_ceiling_compare.json
celleval.json           celleval/full/metrics_celleval_ceiling.json
ci_cd4.json             D1/ci/metrics_ci_bootstrap.json
ci_cd4_D2.json          D2/ci/metrics_ci_bootstrap.json
ci_tahoe.json           tahoe/ci/metrics_ci_bootstrap.json
clone_D1.json           D1/clone_tests/metrics_clone_tests_cd4.json
clone_D2.json           D2/clone_tests/metrics_clone_tests_cd4.json
dmso.json               tahoe/dmso/metrics_tahoe_dmso_null.json
donor_split.json        donor_split/metrics_cd4_donor_split.json
guide_donor.json        guide_donor/metrics_cd4_guide_donor.json
k562_inert_screens.json replogle/k562_inert_screens/metrics_k562_inert_screens.json
k562_template.json      replogle/k562_template/metrics_k562_inert.json
kd_D1.json              D1/kd_matched/metrics_cd4_kd_matched.json
kd_D2.json              D2/kd_matched/metrics_cd4_kd_matched.json
ntc_null.json           ntc/metrics_ntc_null.json
cd4_ntc_ratio.json      ntc/cd4_ntc_ratio.json
trbc_state_keff.json    trbc/trbc_state_keff.json
cl_ntc_replogle.json    replogle/cl_ntc/cl_ntc_replogle.json
cl_ntc_nadig.json       replogle/cl_ntc/cl_ntc_nadig.json
cl_ntc_h1.json          replogle/cl_ntc/cl_ntc_h1.json
pass_cd4_D1_cell.json   D1/inert_line_cell/metrics_readout_inert_line.json
pass_cd4_D1_guide.json  D1/inert_line_guide/metrics_readout_inert_line.json
pass_cd4_D2_cell.json   D2/inert_line_cell/metrics_readout_inert_line.json
pass_cd4_D2_guide.json  D2/inert_line_guide/metrics_readout_inert_line.json
pass_k562_rpe1.json     replogle/inert_line_k562_rpe1/metrics_readout_inert_line.json
pass_jh.json            replogle/inert_line_jh/metrics_readout_inert_line.json
replogle_indep.json     replogle/indep/metrics_replogle_indep_readout.json
robust_D1.json          D1/robust/metrics_robust.json
robust_D2.json          D2/robust/metrics_robust.json
robust_tahoe.json       tahoe/robust/metrics_robust.json
sim_v5.json             sim/sim_v6.json
sim_v6.json             sim/sim_v6.json
tahoe_recompute.json    tahoe/plate_recompute/metrics_tahoe_plate_recompute.json
EOF
  log "collected into $F"
}

figures() {       # Figures 1-5, Supplementary Figures S1-S3 and the Supplementary Information text (Tables S1-S6)
  local F="${FIGDATA:-$W/figdata}" O="${FIGURES:-$W/figures}"
  mkdir -p "$O/supp"
  $PY plot_schematics.py "$O"
  $PY plot_main_figs.py "$F" "$O"
  $PY make_supplement.py "$F" "$REPO/paper/supplementary_v1.template.md" "$O/supplementary_v1.md" "$O/supp"
}

STAGES="${*:-extract_cd4 extract_tahoe cd4 cd4_pseudobulk celleval tahoe replogle sim collect figures}"
for s in $STAGES; do
  declare -F "$s" >/dev/null || { echo "unknown stage $s"; exit 1; }
  log "stage $s"; "$s"
done
log "done: $STAGES"
