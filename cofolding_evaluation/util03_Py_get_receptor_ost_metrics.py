import os
import tqdm
import glob
import argparse
import subprocess

parser = argparse.ArgumentParser()

parser.add_argument('--of3_results', '-r', help='Directory containing of3 results')
parser.add_argument('--fragalysis_dir', '-f', help='Fragalysis aligned/ directory')
parser.add_argument('--out_dir', '-o', help='Path to output directory')

args = parser.parse_args()

def main():
    seeds = [1370180479, 1449838082, 1832854922, 1880307061, 2012026466]
    case_l = []
    for ff in os.listdir(args.fragalysis_dir):
        if os.path.isdir(os.path.join(args.fragalysis_dir,ff)):
            case_l.append(ff)
            print(ff)
    
    for case in case_l:
        fragalysis_rec = f'{args.fragalysis_dir}/{case}/{case}_apo.pdb'

        for seed in seeds:
            model_recs = glob.glob(f'{args.of3_results}/{case}/seed_{seed}/*_rec.pdb')

            outdir = f'{args.out_dir}/{case}/seed_{seed}'
            os.makedirs(outdir, exist_ok=True)
            for mr in model_recs:
                mr_name = os.path.basename(mr)[:-8]
                print(mr_name)
                outfile = f'{outdir}/ost-{mr_name}.json'

                if os.path.exists(outfile):
                    continue

                ost_cmd = f'ost compare-structures -m {mr} -r {fragalysis_rec} -rna --tm-score --lddt --bb-lddt --rigid-scores --local-lddt --aa-local-lddt -o {outfile} -v 0'
                subprocess.run(ost_cmd.split())

    #ost compare-structures -m test_results/A71EV2A-x6035a/seed_2012026466/A71EV2A-x6035a_seed_2012026466_sample_1_model_rec.pdb -r test_fragalysis/A71EV2A-x6035a/A71EV2A-x6035a_apo.pdb -rna --tm-score --lddt --bb-lddt --rigid-scores --local-lddt --aa-local-lddt -o test_rec_ost.json -v 1

if __name__=='__main__':
    main()

