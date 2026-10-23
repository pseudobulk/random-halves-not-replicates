# Result files read by the figure and table scripts

These are the small outputs of the analysis scripts (2.3 MB in total), copied here under short names so that the
figures and supplementary tables can be redrawn without the raw data (README, "Quick start"). `bash
pipelines/run_all.sh collect` recreates the same set from a full run; the third column is the output it copies,
relative to `${WORK}/results`.

Provenance: on 2026-10-19 the whole analysis was rerun from the raw data with the code in this repository
(`REPRODUCTION.md`, `repro/comparison.tsv`). Every file here equals its rerun counterpart either bit for bit or up to
float32 accumulation noise that changes no printed value. Before that, each file had also been compared with the
original output on the compute server: all are byte-identical to it, with three exceptions. `celleval.json` is the
rerun with the released pipeline (the earlier run, used up to preprint v11, had been built from a cell subset made
before the guide-parity fix; v12 uses this file). `sim_v6.json` was produced on a workstation (line endings
normalized to LF here); rerunning `clone_model_sim.py` reproduces all of its values exactly. `sim_v5.json` is the
earlier simulation output; its values equal the corresponding entries of `sim_v6.json` (largest difference
2 × 10⁻¹⁶). The files added for preprints v13 and v14 (`cl_ntc_replogle.json`, `cl_ntc_nadig.json`, `cl_ntc_h1.json`,
`cd4_ntc_ratio.json`, `trbc_state_keff.json`) hold only the input file name in their `file` fields.

Each JSON file stores the producing script's docstring under `rules`, so the definitions travel with the numbers.
`cd4_D*.npz` hold per-target arrays (targets, unexpressed / expressed flags, cell counts, cross-half signals and
Pearson correlations under random and guide halves, TRBC1-fraction deviations); see `code/export_fig_data.py`.

| File | Script | Output copied by `collect` | Read by | SHA-256 |
|---|---|---|---|---|
| `cd4_D1.npz` | `export_fig_data.py` | `D1/cd4_D1.npz` | plot_main_figs.py | `e0ee7c7c4387d4d2df4e12ae2961d9fac48d9da915ef1eb176ab46ce33cf9f4c` |
| `cd4_D2.npz` | `export_fig_data.py` | `D2/cd4_D2.npz` | plot_main_figs.py | `7857d47b586da4035d7847cba7383dcfa9e40504a87e79ad2b4b3feb35bf4b00` |
| `cd4_ntc_ratio.json` | `cd4_ntc_ratio.py` | `ntc/cd4_ntc_ratio.json` | plot_main_figs.py | `f249935c05ea855ccaeb740f1c45b258045a2a0a916226d1bbfc898a5a994f70` |
| `ceiling_D1_rest.json` | `ceiling_compare.py cd4` | `D1/ceiling_rest/metrics_ceiling_compare.json` | text only | `cd62cfbdc3b2d54e85e679d516ff14279f66c01ca347863a96bfe8bb4315ba95` |
| `ceiling_D2_rest.json` | `ceiling_compare.py cd4` | `D2/ceiling_rest/metrics_ceiling_compare.json` | text only | `b39b6a6070436a35af85fe44093c8f726cca333feab24a43faf8e724d6140633` |
| `ceiling_tahoe.json` | `ceiling_compare.py tahoe` | `tahoe/ceiling/metrics_ceiling_compare.json` | plot_main_figs.py | `25583a6db8a934f35e8fe480745d87a03ce3664dc3aae2d135346b96590346df` |
| `celleval.json` | `run_celleval_ceiling.py` | `celleval/full/metrics_celleval_ceiling.json` | make_supplement.py, plot_main_figs.py | `629127984b01375fe95ca274cb5b1aa39ab6129cc48e839724b8ea25fb56f1ee` |
| `ci_cd4.json` | `ci_bootstrap.py cd4` | `D1/ci/metrics_ci_bootstrap.json` | plot_main_figs.py | `625bba012a8589d13218ecbf1dcf07b001e89147f50b12592ae144e2bce96249` |
| `ci_cd4_D2.json` | `ci_bootstrap.py cd4` | `D2/ci/metrics_ci_bootstrap.json` | plot_main_figs.py | `1a8261e22e96814cd90c18c906b79fb2ab4c686e25549f590fab07093fe80c3f` |
| `ci_tahoe.json` | `ci_bootstrap.py tahoe` | `tahoe/ci/metrics_ci_bootstrap.json` | plot_main_figs.py | `7149cb0ab09ff117df8a30b81a008057a7a8797894d4cc28c674bb5f8376aa72` |
| `cl_ntc_h1.json` | `cl_ntc_null.py` | `replogle/cl_ntc/cl_ntc_h1.json` | make_supplement.py, plot_main_figs.py | `f39685440563d07349f48351f9b561ecc2d306b476e72b619e27e51356d92112` |
| `cl_ntc_nadig.json` | `cl_ntc_null.py` | `replogle/cl_ntc/cl_ntc_nadig.json` | make_supplement.py, plot_main_figs.py | `d7f2475034fd2993f22fc94451efe31898e9c183299ee34edbe58bfff42d6ae8` |
| `cl_ntc_replogle.json` | `cl_ntc_null.py` | `replogle/cl_ntc/cl_ntc_replogle.json` | make_supplement.py, plot_main_figs.py | `a66f39e9234c1fb86a9246c7a90b5125f685c0069eca0d0a3e389c9d15a53c3d` |
| `clone_D1.json` | `clone_tests_cd4.py` | `D1/clone_tests/metrics_clone_tests_cd4.json` | plot_main_figs.py | `e9266834c46290a34ea14abc189fc3f13aa98f265b80d875c6480bb2ce3a81f3` |
| `clone_D2.json` | `clone_tests_cd4.py` | `D2/clone_tests/metrics_clone_tests_cd4.json` | plot_main_figs.py | `ed8f309d5cb76c08e7d96d3e37a8715eb1bb8bc0715fe44ff6d69f2b96c898e6` |
| `dmso.json` | `tahoe_dmso_null.py` | `tahoe/dmso/metrics_tahoe_dmso_null.json` | make_supplement.py, plot_main_figs.py | `549876a7353cd2cf3eb83e22a04100aa1649a0259b3bfe902842d1cc35b624f7` |
| `donor_split.json` | `cd4_donor_split.py` | `donor_split/metrics_cd4_donor_split.json` | plot_main_figs.py | `dbdc436a43a68b1e513782bfad5f5e489511951d3bea0b401b02a88f9d507f81` |
| `guide_donor.json` | `cd4_guide_donor.py` | `guide_donor/metrics_cd4_guide_donor.json` | make_supplement.py, plot_main_figs.py | `5ba09303dbb738bd553893739e81ee2e8aaf33f874c0eb1baa96a64a093b5185` |
| `k562_inert_screens.json` | `k562_inert_screens.py` | `replogle/k562_inert_screens/metrics_k562_inert_screens.json` | make_supplement.py, plot_main_figs.py | `1e9732911db328532def6637a8bf3d14fca3180e01b84b8d5baa33161bb13730` |
| `k562_template.json` | `k562_inert.py` | `replogle/k562_template/metrics_k562_inert.json` | make_supplement.py | `8d0ec94202ae54e575eaf9e251a7f4e3227a8080ef6586e4ebea7d4a96ccfbed` |
| `kd_D1.json` | `cd4_kd_matched.py` | `D1/kd_matched/metrics_cd4_kd_matched.json` | plot_main_figs.py | `1045b61bd93bd33c3045abfd5499799769f0f6809e3f37688c5e19d2d6074c65` |
| `kd_D2.json` | `cd4_kd_matched.py` | `D2/kd_matched/metrics_cd4_kd_matched.json` | plot_main_figs.py | `58dffd889119d0178b252965be6a61638995cda419c32e9cabda3ad822717b2a` |
| `ntc_null.json` | `cd4_ntc_null.py` | `ntc/metrics_ntc_null.json` | make_supplement.py | `caa0a952722c86ee1faff755d6d0a06a362e3ab413b08edc0e3763a674becdf6` |
| `pass_cd4_D1_cell.json` | `readout_inert_line.py cd4` | `D1/inert_line_cell/metrics_readout_inert_line.json` | make_supplement.py | `d892f7a574ceafb53d28e0d6ce35d24f4c6075d8db9c9f1727f6e9f0f51f5737` |
| `pass_cd4_D1_guide.json` | `readout_inert_line.py cd4` | `D1/inert_line_guide/metrics_readout_inert_line.json` | make_supplement.py | `e590c04469397cfe66b38ffb8a777a08d8f302c2681bc3120716db116ef4b2a1` |
| `pass_cd4_D2_cell.json` | `readout_inert_line.py cd4` | `D2/inert_line_cell/metrics_readout_inert_line.json` | make_supplement.py | `0fcb7a844dbb79c2eefdcb83525e1d0987bbc9968994d498251950b884c6b436` |
| `pass_cd4_D2_guide.json` | `readout_inert_line.py cd4` | `D2/inert_line_guide/metrics_readout_inert_line.json` | make_supplement.py | `0e8b3c7204b8c038d470d72f20d4eb42827a120f1f7121d310a116423cf22874` |
| `pass_jh.json` | `readout_inert_line.py h5ad` | `replogle/inert_line_jh/metrics_readout_inert_line.json` | make_supplement.py | `de557b7221504501af46350099e11ca5e2e684d3b2f3e8aebc4eefdb1764f7f3` |
| `pass_k562_rpe1.json` | `readout_inert_line.py h5ad` | `replogle/inert_line_k562_rpe1/metrics_readout_inert_line.json` | make_supplement.py | `6c18c45b4b48ee99a76f062e6119c1270ce43e8fe1d742ee6e1589d9092eccc0` |
| `replogle_indep.json` | `replogle_indep_readout.py` | `replogle/indep/metrics_replogle_indep_readout.json` | plot_main_figs.py | `320485793d023df4e747a54fa5bca32b1dc26ffb14bc49026764439fdb9ccb9b` |
| `robust_D1.json` | `robust_cd4_tahoe.py cd4` | `D1/robust/metrics_robust.json` | make_supplement.py | `9394fa756d6e7a027b8038f5b7abb69224236ff0d8c2040dc2f371c4575a225b` |
| `robust_D2.json` | `robust_cd4_tahoe.py cd4` | `D2/robust/metrics_robust.json` | make_supplement.py | `397d4ae670efd7a888250b8516f6fb51f84224a998b574179257777289dbf346` |
| `robust_tahoe.json` | `robust_cd4_tahoe.py tahoe` | `tahoe/robust/metrics_robust.json` | make_supplement.py | `f04916ce24bc33f2ad9f3b0716e87dcf70a4c6a116609849f86ffbb1cd5ac32e` |
| `sim_v5.json` | `clone_model_sim.py` | `sim/sim_v6.json` | make_supplement.py | `87746650c842ce802b2d8a42403b67310ccc73aec0e11e50f47675d840fd5a55` |
| `sim_v6.json` | `clone_model_sim.py` | `sim/sim_v6.json` | plot_main_figs.py | `dcf2758372d30512c2c4946d311657d649fea1122645081616a1fd4981ce9499` |
| `tahoe_recompute.json` | `tahoe_plate_recompute.py` | `tahoe/plate_recompute/metrics_tahoe_plate_recompute.json` | text only | `192904646cc5f975910303e7df0d974e8f43e5fab6e20bc5973ae62098466103` |
| `trbc_state_keff.json` | `cd4_trbc_state_keff.py` | `trbc/trbc_state_keff.json` | text only | `fe8594774e7847a07fbbb3123c5aa9d81b569d33e1321d6ed8f0fb665f2fb2fd` |

## Second version: data/next/

The result files of the second version of the preprint (6.5 MB) are in `data/next/`, under the names that
`plot_next_figs.py` and `me_summary.py` read; `bash pipelines/run_next.sh collect_next` recreates the set from a full run
(`FIGDATA/next`). The README section "Second version of the preprint" lists the script and manuscript item behind each
file. `data/next/pilot/` holds the score files of the donor-1 model pilot (Fig. 7a, Supplementary Table S6), from which
`model_eval/pilot_summary.py` writes `model_pilot_summary_table.md`; `data/next/inputs/` holds the run list and library
table of the organoid clone-barcode reads and the file list of the whole-brain data. The `file` fields of
`data/next/unified/*.json` hold only the input file name.
