import os
import glob
import shutil
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('--af3_dir', '-ad', help='Protenix results directory')
parser.add_argument('--out_dir', '-od', help='Output directory with boltz results reformatted')

args = parser.parse_args()

def copy_files(file_l, case, seed, outdir):
    for f in file_l:
        fname = os.path.basename(f)
        print('\t', fname, fname[-4:])
        if fname[-4:] == '.cif':
            print('\t', f)
            sample = fname[:-4].split('_')[-2].split('-')[1]
            new_fname = f'{case}_{seed}_sample_{sample}_model.cif'
            print('\t\t', new_fname)

            shutil.copy(f, f'{outdir}/{new_fname}')

        elif 'summary_confidences.json' in fname:
            sample = fname.split('_')[-3].split('-')[-1]
            new_fname = f'{case}_{seed}_sample_{sample}_confidences_aggregated.json'
            print('\t\t', new_fname)
            
            shutil.copy(f, f'{outdir}/{new_fname}')

    pass

def listdir_directories(query_dir):
    results = os.listdir(query_dir)
    
    dir_l = []
    for d in results:
        if os.path.isdir(f'{query_dir}/{d}'):
            dir_l.append(d)
            
    return dir_l

def main():
    case_list = os.listdir(args.af3_dir)

    os.makedirs(args.out_dir, exist_ok=True)
    
    # Debug
    n = 0
    max_cases = 1 
    for case in case_list:
        print(case)
        #seed_list = os.listdir(f'{args.af3_dir}/{case}/')
        seed_sample_list = listdir_directories(f'{args.af3_dir}/{case}/')
        print(seed_sample_list)


        for seed_sample in seed_sample_list:
            seed_n, sample_n = seed_sample.split('_')

            seed = '_'.join(seed_n.split('-'))
            print(case, seed)
            case_outdir = f'{args.out_dir}/{case}/{seed}/'
            os.makedirs(case_outdir, exist_ok=True)
            output_files = glob.glob(f'{args.af3_dir}/{case}/{seed_sample}/*')
            print(output_files)
            
            #print(output_files)
            copy_files(output_files, case, seed, case_outdir)
        
        #Debug
        #n += 1
        #if n >= max_cases:
        #    return

if __name__=='__main__':
    main()

