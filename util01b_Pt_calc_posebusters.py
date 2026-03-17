
import os
import tqdm
import glob
import json
import argparse
import numpy as np
import pandas as pd
from posebusters import PoseBusters

parser = argparse.ArgumentParser()

parser.add_argument('--result_dir', '-r', help='Directory with cofolding results, in OF3 format', required=True)
#parser.add_argument('--fragalysis_dir', '-f', help='Path to aligned/ fragalysis folder', required=True)
#parser.add_argument('--outfile', '-o', help='Name of the output .tsv file', required=True)

args = parser.parse_args()
buster = PoseBusters(config="dock")

def check_posebusters(mdl_lig, mdl_rec):
    df = buster.bust([mdl_lig], None, mdl_rec)

    pb_valid = True
    for h in df.head():
        if df[h].iloc[0] == False:
            pb_valid = False

    return pb_valid

def main():
    for case in tqdm.tqdm(os.listdir(args.result_dir)):
        # Skip non directory items
        try:
            os.listdir(f'{args.result_dir}/{case}/')
        except:
            continue

        for seed in os.listdir(f'{args.result_dir}/{case}/'):
            print(case, seed)
            out_json = f'{args.result_dir}/{case}/{seed}/posebusters_data.json'
            pb_data = {}
            # A71EV2A-x7597a_seed_2012026466_sample_1_model.cif
            for model in glob.glob(f'{args.result_dir}/{case}/{seed}/*_model.cif'):
                model_n = os.path.basename(model)[:-4]
                sample = model_n.split('_')[-2]
                
                mdl_rec = f'{args.result_dir}/{case}/{seed}/{model_n}_rec.pdb'
                ligs = glob.glob(f'{args.result_dir}/{case}/{seed}/{model_n}_*-lig.sdf')

                for l in ligs:
                    pb_valid = check_posebusters(l, mdl_rec)
                    pb_data[os.path.basename(l)] = pb_valid

            with open(out_json, 'w') as fo:
                json.dump(pb_data, fo, indent=4)

if __name__=='__main__':
    main()
