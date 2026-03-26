# Cofolding Evaluation Code

## Code to postprocess and run ost on cofolded models

Install/Activate the environment for prep scripts
```
mamba env create -f prep.yaml
mamba activate prep
```


To preprocess the output, convert cofolding predictions to the openfold 3 format:
```
python3 Py_convert_boltz_to_of3_fmt.py -bd=examples/boltz-2_results/ -od=examples/boltz-2_results_reformatted
```

Next, extract sdf and pdb files for the ligands and receptors:
```
python3 util01_Py_extract_of3_ligand_sdfs.py -r=examples/boltz-2_results_reformatted/ -fd=examples/fragalysis_data/
```

## Code for running OST evaluations

Install the environment for ost scripts
```
mamba env create -f ost211.yaml
mamba activate ost211
```

Calculate posebusters checks on the cofolded models:
```
python3 util01b_Pt_calc_posebusters.py -r=examples/boltz-2_results_reformatted/
```
This script outputs `posebusters_data.json` in each results folder


Run the ost-ligand comparison script:
```
python3 util02_Py_score_of3_with_ost_v2.py -r=examples/boltz-2_results_reformatted/ -f=examples/fragalysis_data/ -o=examples/boltz-2_ost-ligand_results -m=boltz
```

Run the ost-recepotr comparison script:
```
python3 util03_Py_get_receptor_ost_metrics.py -r=examples/boltz-2_results_reformatted/ -f=examples/fragalysis_data/ -o=examples/
```

## Code for compiling OST metrics into a tsv file

Compile OST metrics
```
python3 Py_compile_performance_metrics.py -r=examples/boltz-2_results_reformatted/ -ol=examples/boltz-2_ost-ligand_results/ -or=examples/boltz-2_ost-receptor_results/ -st=examples/tsv_similarity_data_2023-06-01.tsv -m=boltz-2 -o=examples/metrics_boltz2_ost.tsv
```

Code for calculating the SuCOS-Pocket similarity data provided as the `--similarity-tsv` or `-st` can be installed via out [SuCOS_pocket code repository](https://github.com/OmeirK/pocket_sucos_code/tree/main)

If OST failed to score ligands, or there are any other issues, they will be stored in the in a error file that shares the name of the metrics csv. (i.e. `examples/metrics_boltz2_ost.tsv` has error file `examples/metrics_boltz2_ost.err`)

Calculate the `pocket_recall` metrics and append failure modes to the tsv file:
```
python3 util04_Py_get_pocket_recall.py -t=examples/metrics_boltz2_ost.tsv -r=examples/boltz-2_results_reformatted/ -f=examples/fragalysis_data/
```
If no output file is specified, the appended .tsv is stored in the same location as the input with a modified name: `examples/metrics_boltz2_ost_failure_modes.tsv`


