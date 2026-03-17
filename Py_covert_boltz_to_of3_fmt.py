import os
import glob
import shutil
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--boltz_dir', '-bd', help='Boltz results directory')
parser.add_argument('--out_dir', '-od', help='Output directory with boltz results reformatted')

args = parser.parse_args()

def copy_files(file_l, case, seed, outdir):
    for f in file_l:
        fname = os.path.basename(f)
        print('\t', fname, fname[-4:])
        if fname[-4:] == '.cif':
            print('\t', f)
            sample = fname[:-4].split('_')[-1]
            new_fname = f'{case}_{seed}_sample_{sample}_model.cif'
            print('\t\t', new_fname)

            shutil.copy(f, f'{outdir}/{new_fname}')

        elif fname[-5:] == '.json':
            sample = fname[:-5].split('_')[-1]
            new_fname = f'{case}_{seed}_sample_{sample}_confidences_aggregated.json'
            print('\t\t', new_fname)
            
            shutil.copy(f, f'{outdir}/{new_fname}')

    pass

def main():
    case_list = os.listdir(args.boltz_dir)

    os.makedirs(args.out_dir, exist_ok=True)
    
    # Debug
    n = 0
    max_cases = 10 
    for case in case_list:
        seed_list = os.listdir(f'{args.boltz_dir}/{case}/')
        
        for seed in seed_list:
            print(case, seed)
            case_outdir = f'{args.out_dir}/{case}/{seed}/'
            os.makedirs(case_outdir, exist_ok=True)
            output_files = glob.glob(f'{args.boltz_dir}/{case}/{seed}/*/predictions/*/*')

            copy_files(output_files, case, seed, case_outdir)
        
        #Debug
        #n += 1
        #if n >= max_cases:
        #    return

if __name__=='__main__':
    main()

