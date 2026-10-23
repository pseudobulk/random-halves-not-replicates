#!/usr/bin/env bash
# The analyses added in the second version of the preprint, in the order they have to run. Paths come from config.env
# at the repository root, as for run_all.sh, which must have run first for the stages extract_cd4, cd4 and celleval
# (this pipeline reads their per-cell extractions and the cell-eval subset).
#   bash pipelines/run_next.sh                 all stages in order
#   bash pipelines/run_next.sh figures_next    one or more named stages
# Stages (README, "Second version"):
#   extract_next  cd4_next  unified  generality  organoid  jin2020  wholebrain  stats
#   model_inputs  model_gpu  model_pilot  model_post  collect_next  figures_next
# The original runs were spread over two servers (extraction on a CPU server, model training on a GPU server, files
# copied between them); here every command runs on one machine with the same arguments. GPU stages need GPU, STATE_BIN,
# STATE_PYTHON and the model environments of code/model_eval/mkenv_models.sh.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CFG="${RHNR_CONFIG:-$REPO/config.env}"
[ -f "$CFG" ] || { echo "missing $CFG: copy config.example.env to config.env and edit it"; exit 1; }
_env="$(export -p)"; set -a; . "$CFG"; set +a; eval "$_env"   # variables already exported take precedence
export RHNR_CONFIG="$CFG"
cd "$REPO/code"
PY="${PYTHON:-python3}"; D="$DATA_ROOT"; W="$WORK"; N="${WORKERS:-16}"; RES="$W/results"; NX="$RES/next"; G="${GPU:-0}"
ME="$W/model_eval"; MC="$REPO/code/model_eval"
export WORK STATE_BIN STATE_PYTHON
mkdir -p "$W/cd4" "$NX" "$ME"
log() { echo "$(date '+%F %T') $*"; }
h5() { echo "$D/cd4_sc/$1.assigned_guide.h5ad"; }

extract_next() {  # per-cell aggregates of donors 3 and 4; 16-code aggregates (depth grid) and TRBC classes, donors 1 and 2
  for n in D3_Rest D3_Stim48hr D4_Rest D4_Stim48hr; do
    log "extract cells $n"; $PY cd4_sc_extract_cells.py "$(h5 $n)" "$W/cd4/${n}_cells" "$N"
  done
  for n in D1_Rest D1_Stim48hr D2_Rest D2_Stim48hr; do
    log "extract grid and TRBC classes $n"
    $PY cd4_sc_extract_grid.py "$(h5 $n)" "$W/cd4/${n}_cells" "$W/cd4/${n}_grid16" 8
    $PY cd4_trbc_split.py "$(h5 $n)" "$W/cd4/${n}_cells/cells.npz" "$W/cd4/${n}_trbc" 12
  done
}

cd4_next() {      # Fig. 2d depth grid; Fig. 3d constant-region (TRBC) halves; TRBC clone tests of donors 3 and 4
  mkdir -p "$NX/depth_grid" "$NX/trbc_split"
  for d in D1 D2; do
    $PY depth_grid_cd4.py "$W/cd4/${d}_Rest_cells" "$W/cd4/${d}_Rest_grid16" rest \
      "$W/cd4/${d}_Stim48hr_cells" "$W/cd4/${d}_Stim48hr_grid16" stim48 "$NX/depth_grid/depth_grid_$d.json"
    $PY trbc_split_test.py "$W/cd4/${d}_Rest_trbc" ${d}_Rest "$NX/trbc_split/${d}_Rest.json"
    $PY trbc_split_test.py "$W/cd4/${d}_Stim48hr_trbc" ${d}_Stim48 "$NX/trbc_split/${d}_Stim48.json"
  done
  for d in D3 D4; do
    $PY clone_tests_cd4.py "$W/cd4/${d}_Rest_cells" rest "$W/cd4/${d}_Stim48hr_cells" stim48 "$RES/$d/clone_tests"
  done
}

unified() {       # one clone-strength estimator on the screens with several null sgRNAs (Fig. 4c, Supplementary Table S3)
  local U="$NX/unified" R="$D/replogle"; mkdir -p "$U"
  export OMP_NUM_THREADS=4
  $PY model_eval/generality/unified_null.py "$R/K562_essential_raw_singlecell_01.h5ad" sgID_AB gene gem_group "$U/k562ess.json"
  $PY model_eval/generality/unified_null.py "$R/rpe1_raw_singlecell_01.h5ad" sgID_AB gene gem_group "$U/rpe1.json"
  $PY model_eval/generality/unified_null.py "$R/K562_gwps_raw_singlecell_01.h5ad" sgID_AB gene gem_group "$U/k562gw.json"
  $PY model_eval/generality/unified_null.py "$D/nadig/hepg2.h5ad" sgID_AB gene gem_group "$U/hepg2.json"
  $PY model_eval/generality/unified_null.py "$D/nadig/jurkat.h5ad" sgID_AB gene gem_group "$U/jurkat.json"
  $PY model_eval/generality/unified_null.py "$D/vcc2025/adata_Training.h5ad" guide_id target_gene batch "$U/vcc2025_h1.json"
  for t in Control IFNγ Co-culture; do
    $PY model_eval/generality/unified_null.py "$D/frangieh/FrangiehIzar2021_RNA.h5ad" sgRNA perturbation - "$U/frangieh_$t.json" "keep=perturbation_2:$t"
  done
  for d in D1 D2 D3 D4; do for s in Rest Stim48hr; do
    $PY model_eval/generality/unified_null.py "$(h5 ${d}_$s)" guide_id @guide_type lane_id "$U/cd4_${d}_$s.json" qc=low_quality "drop=guide_group:multi sgRNA"
  done; done
  unset OMP_NUM_THREADS
}

generality() {    # the same estimator on in vivo and primary-cell screens with their own loaders (Fig. 4c, Table S3)
  local O="$NX/generality"; mkdir -p "$O"
  $PY model_eval/generality/run_generality.py zhou "$D/zhou" "$O/zhou.json"
  $PY model_eval/generality/run_generality.py knudsen "$D/knudsen" "$O/knudsen.json"
  $PY model_eval/generality/run_generality.py knudsen_lineage "$D/knudsen" "$O/knudsen_lineage.json"
  $PY model_eval/generality/run_generality.py fagerberg "$D/fitdb" "$O/fagerberg.json"
  $PY model_eval/generality/run_generality.py schmidt_resting "$D/fitdb" "$O/schmidt_resting.json"
  $PY model_eval/generality/run_generality.py schmidt_restim "$D/fitdb" "$O/schmidt_restim.json"
  $PY model_eval/generality/run_generality.py cao_invitro "$D/cao" "$O/cao_invitro.json"
  $PY model_eval/generality/run_generality.py cao_invivo "$D/cao" "$O/cao_invivo.json"
}

organoid() {      # CHOOSE organoids: clone barcodes, clone halves and the formal tests (Fig. 3, Fig. 4a-b, Table N1)
  local C="$W/choose"; mkdir -p "$C/cbc"
  "${RSCRIPT:-Rscript}" model_eval/choose/export_object.R "$D/choose/CHOOSE_ASD_Modules.rds" "$C/export"
  # clone-barcode reads streamed from ENA (86 runs; columns: run, R1 url, R1 bytes, R2 url, R2 bytes)
  while read -r run u1 b1 u2 b2; do
    [ -f "$C/cbc/$run.ok" ] && continue
    $PY model_eval/choose/cbc_parse.py "$run" "$u1" "$b1" "$u2" "$b2" "$C/cbc/$run" && touch "$C/cbc/$run.ok"
  done < "$REPO/data/next/inputs/choose_cbc_runs.txt"
  $PY model_eval/choose/choose_clone.py "$C/export" "$C/cbc" "$REPO/data/next/inputs/choose_cbc_libs.tsv" "$NX/choose_clone.json"
  $PY model_eval/choose/choose_stats.py "$C/export" "$C/cbc" "$REPO/data/next/inputs/choose_cbc_libs.tsv" "$NX/choose_stats.json" 200 1000
}

jin2020() {       # embryonic mouse cortex, library halves (Fig. 4c, Table S3); <dir> holds the GEO files and barcodes.csv
  local J="$W/jin2020"; mkdir -p "$J"
  for f in "$D"/jin2020/*; do ln -sfn "$f" "$J/$(basename "$f")"; done
  cp "$REPO/data/next/jin2020_barcodes.csv" "$J/barcodes.csv"
  $PY model_eval/jin2020/jin_clone.py "$J" "$NX/jin2020_clone.json"
}

wholebrain() {    # mouse whole brain, AAV in non-dividing neurons (Fig. 4c, Table S3)
  local B="$W/wholebrain"; mkdir -p "$B"
  $PY model_eval/wholebrain/wb_meta.py "$D/wholebrain/all_obs.parquet" "$B/wb_meta.npz"
  $PY model_eval/wholebrain/wb_aggregate.py "$B/wb_meta.npz" "$D/wholebrain/h5ads" "$B/wb_agg.npz" 7
  OMP_NUM_THREADS=8 $PY model_eval/wholebrain/wb_analyze.py "$B/wb_agg.npz" "$NX/wholebrain_wb_results.json" 0
  $PY model_eval/wholebrain/wb_mouse_halves.py "$B/wb_meta.npz" "$B/wb_meta_mouse.npz"
  $PY model_eval/wholebrain/wb_aggregate.py "$B/wb_meta_mouse.npz" "$D/wholebrain/h5ads" "$B/wb_agg_mouse.npz" 7
  OMP_NUM_THREADS=8 $PY model_eval/wholebrain/wb_split.py "$B/wb_agg.npz" "$B/wb_agg_mouse.npz" "$NX/wholebrain_wb_split.json" 0
}

stats() {         # Supplementary Note 1: symbolic checks (sympy) and simulation recovery
  mkdir -p "$NX/stats"
  $PY stats/clone_theory_check.py "$NX/stats/clone_theory_check.json"
}

model_inputs() {  # cell subsets for the transfer models (CPU): CD4 donors 1-4, K562 -> RPE1, Frangieh folds
  for d in D1 D2 D3 D4; do
    local L="$ME/sc_$d"; mkdir -p "$L"
    for s in Rest Stim48hr; do
      ln -sfn "$W/cd4/${d}_${s}_cells" "$L/cells_${d}_${s}_v2"; ln -sfn "$(h5 ${d}_$s)" "$L/${d}_$s.assigned_guide.h5ad"
    done
    local SEL="$ME/select_$d.npz"
    if [ $d = D1 ]; then SEL="$RES/celleval/cd4_d1_rest_subset.npz"   # donor 1: the targets of the cell-eval comparison
    else $PY model_eval/me_select_targets.py "$W/cd4/${d}_Rest_cells" "$SEL"; fi
    PILOT_DONOR=$d $PY cd4_pilot_extract.py "$L" "$SEL" "$ME/pilot_$d" "$N"
  done
  PILOT_DONOR=D1 PILOT_ALL_GENES=1 $PY cd4_pilot_extract.py "$ME/sc_D1" "$RES/celleval/cd4_d1_rest_subset.npz" "$ME/pilot_D1" "$N"
  RATIO=0 OMP_NUM_THREADS=8 $PY model_eval/me_replogle_extract.py "$D/replogle/K562_essential_raw_singlecell_01.h5ad" \
    "$D/replogle/rpe1_raw_singlecell_01.h5ad" "$ME/pilot_K562_RPE1"
  $PY model_eval/me_frangieh_extract.py "$D/frangieh/FrangiehIzar2021_RNA.h5ad" "$ME/frangieh"
}

model_gpu() {     # baselines, STATE, scoring, prediction export and fit check per run (code/model_eval/run_transfer.sh)
  for d in D1 D2 D3 D4; do PILOT_DONOR=$d bash "$MC/run_transfer.sh" "$ME/pilot_$d" $d "$G"; done
  PILOT_STATES=K562,RPE1 bash "$MC/run_transfer.sh" "$ME/pilot_K562_RPE1" K562_RPE1 "$G"
  PILOT_STATES=Control,IFNg bash "$MC/run_transfer.sh" "$ME/frangieh/IFNg/fold0" fr_IFNg_f0 "$G" local
  PILOT_STATES=Control,IFNg bash "$MC/run_transfer.sh" "$ME/frangieh/IFNg/fold1" fr_IFNg_f1 "$G" local
  PILOT_STATES=Control,Coculture bash "$MC/run_transfer.sh" "$ME/frangieh/Coculture/fold0" fr_Cocu_f0 "$G" baselines
  PILOT_STATES=Control,Coculture bash "$MC/run_transfer.sh" "$ME/frangieh/Coculture/fold1" fr_Cocu_f1 "$G" baselines
}

model_pilot() {   # donor 1 only: the other models of Fig. 7a and Supplementary Table S6 (score files in $ME/pilot_out)
  local DD="$ME/pilot_D1" O="$ME/pilot_out" SP="${STATE_PYTHON:-python3}"; mkdir -p "$O"
  cp "$ME/out_D1/baselines.json" "$ME/out_D1/state_st_scores.json" "$O/"
  for m in perturb_mean pseudobulk; do
    bash "$MC/run_state_model.sh" $m "$G"
    (cd "$ME/state_D1" && CUDA_VISIBLE_DEVICES=$G "${STATE_BIN:-state}" tx predict --output-dir runs/cd4_d1_$m --checkpoint final.ckpt --predict-only)
    PYTHONPATH="$MC" $SP "$MC/pilot_state_score.py" "$DD" "$ME/state_D1/runs/cd4_d1_$m/eval_final.ckpt/adata_pred.h5ad" "$O/arc_${m}_scores.json"
  done
  # State Embedding (SE-600M weights in SE_DIR), then STATE on the embeddings
  mkdir -p "$ME/state_se"
  CUDA_VISIBLE_DEVICES=$G PYTHONPATH="$MC" $SP "$MC/pilot_se_embed.py" "$DD" "$SE_DIR" "$ME/state_se"
  sed "s#$ME/state_D1/data/#$ME/state_se/data/#" "$ME/state_D1/cd4_fewshot.toml" > "$ME/state_se/cd4_fewshot_se.toml"
  bash "$MC/run_state_se.sh" "$G"
  (cd "$ME/state_se" && CUDA_VISIBLE_DEVICES=$G "${STATE_BIN:-state}" tx predict --output-dir runs/cd4_d1_se --checkpoint final.ckpt --predict-only)
  PYTHONPATH="$MC" $SP "$MC/pilot_state_score.py" "$DD" "$ME/state_se/runs/cd4_d1_se/eval_final.ckpt/adata_pred.h5ad" "$O/state_se_scores.json"
  # PIE: inputs, training (run_pie.sh), then `pie eval` with save_predictions=true on checkpoint best_auprc, and scoring
  PYTHONPATH="$MC" $PY "$MC/pilot_pie_prep.py" "$DD" "$ME/state_D1/prep.json" "${PIE_DIR:-$ME/pie_data}/cd4"
  bash "$MC/run_pie.sh" "$G"
  echo "PIE: run 'pie eval' (save_predictions=true, checkpoint best_auprc) and then"
  echo "  pilot_pie_score.py $DD <predictions.parquet> $O/pie_scores.json"
  # Stack (zero-shot), biolord and CPA (tutorial / published settings and the adapted variants)
  CUDA_VISIBLE_DEVICES=$G PYTHONPATH="$MC" "$W/envs/stackenv/bin/python" "$MC/pilot_stack.py" "$DD" "$W/envs/stack_weights" "$O/stack"
  CUDA_VISIBLE_DEVICES=$G PYTHONPATH="$MC" "$W/envs/biolordenv/bin/python" "$MC/pilot_biolord.py" "$DD" "$O/biolord"
  CUDA_VISIBLE_DEVICES=$G BIOLORD_CAT_DIM=32 BIOLORD_NO_ES=1 PYTHONPATH="$MC" "$W/envs/biolordenv/bin/python" "$MC/pilot_biolord.py" "$DD" "$O/biolord_adapted" 200
  CUDA_VISIBLE_DEVICES=$G PYTHONPATH="$MC" "$W/envs/cpaenv/bin/python" "$MC/pilot_cpa.py" "$DD" "$O/cpa"
  echo "CPA adapted ($O/cpa_adapted): pilot_cpa.py with CPA_NO_ES=1; its epoch setting is not recorded (README)"
  PYTHONPATH="$MC" $PY "$MC/pilot_summary.py" "$O"
}

model_post() {    # clonality, cross-donor control, expression clonality, matched cells, summary table
  local DS="" CD=""
  for d in D1 D2 D3 D4; do
    DS="$DS $d:$W/cd4/${d}_Rest_cells:$W/cd4/${d}_Stim48hr_cells:$ME/pilot_$d/pilot_design.json"
    CD="$CD $d:$ME/out_$d/preds_$d.npz:$ME/pilot_$d/pilot_design.json:$W/cd4/${d}_Stim48hr_cells"
  done
  $PY model_eval/me_clonality.py "$NX/clonality.json" $DS
  $PY model_eval/me_cross_donor.py "$NX/cross_donor.json" $CD
  mkdir -p "$ME/step3"
  for d in D1 D2 D3 D4; do
    PILOT_DONOR=$d $PY model_eval/me_vc_expr.py "$ME/pilot_$d" inert "$NX/vc_expr_$d.json"
    PILOT_DONOR=$d $PY model_eval/me_downsample.py "$ME/pilot_$d" "$ME/step3/ds_$d" 48 75
    PILOT_DONOR=$d PYTHONPATH="$MC" $PY model_eval/pilot_baselines.py "$ME/step3/ds_$d" "$NX/baselines_ds48_75_$d.json"
    rm -f "$ME/step3/ds_$d"/pilot_*.npz
  done
  PILOT_STATES=Control,IFNg $PY model_eval/me_vc_expr.py "$ME/frangieh/IFNg/fold0" null "$NX/vc_expr_fr_IFNg.json"
  PILOT_STATES=Control,Coculture $PY model_eval/me_vc_expr.py "$ME/frangieh/Coculture/fold0" null "$NX/vc_expr_fr_Coculture.json"
}

collect_next() {  # gather the small result files into FIGDATA/next under the names in data/next
  local F="${FIGDATA:-$W/figdata}/next"; mkdir -p "$F/pilot"
  cp -r "$NX"/. "$F/"
  for r in D1 D2 D3 D4 K562_RPE1 fr_IFNg_f0 fr_IFNg_f1 fr_Cocu_f0 fr_Cocu_f1; do
    [ -f "$ME/out_$r/baselines.json" ] && cp "$ME/out_$r/baselines.json" "$F/baselines_$r.json"
    [ -f "$ME/out_$r/detail.json" ] && cp "$ME/out_$r/detail.json" "$F/detail_$r.json"
    [ -f "$ME/out_$r/fitcheck.json" ] && cp "$ME/out_$r/fitcheck.json" "$F/fitcheck_$r.json"
  done
  for d in D1 D2 D3 D4; do
    cp "$ME/pilot_$d/pilot_design.json" "$F/pilot_design_$d.json"
    cp "$RES/$d/clone_tests/metrics_clone_tests_cd4.json" "$F/clone_tests_$d.json"
  done
  cp "$ME/pilot_K562_RPE1/pilot_design.json" "$F/pilot_design_K562_RPE1.json"
  for t in IFNg Coculture; do for k in 0 1; do
    cp "$ME/frangieh/$t/fold$k/pilot_design.json" "$F/pilot_design_fr_${t}_f$k.json"
  done; done
  cp "$REPO/data/next/jin2020_barcodes.csv" "$F/"
  cp -r "$ME/pilot_out"/. "$F/pilot/"; cp "$ME/pilot_out/pilot_summary.md" "$F/model_pilot_summary_table.md"
  $PY model_eval/me_summary.py "$F" "$W/figures_next/summary"
}

figures_next() {  # Figures 1-7 of the second version from data/ (first version) and FIGDATA/next
  local F="${FIGDATA:-$W/figdata}"
  $PY plot_next_figs.py "$F" "$F/next" "$W/figures_next"
}

STAGES=("$@"); [ ${#STAGES[@]} -eq 0 ] && STAGES=(extract_next cd4_next unified generality organoid jin2020 wholebrain stats model_inputs model_gpu model_pilot model_post collect_next figures_next)
for s in "${STAGES[@]}"; do log "== $s"; "$s"; done
log "done"
