import os
import tqdm
import json
import glob
import argparse
import subprocess
import multiprocessing as mp
from rdkit import Chem

parser = argparse.ArgumentParser()

parser.add_argument('--of3_results', '-r', help='Directory with cofolding results and extracted ligand sdf files')
parser.add_argument('--fragalysis_dir', '-f', help='Path to aligned_files/ directory with fragalysis results')
parser.add_argument('--outdir', '-o', help='Directory to store ost comparison results')
parser.add_argument('--cofolding_model', '-m', help='Specify which cofolding model was used to predict structures. This parameter determines the ligand matching criteria', default='of3', choices=['of3', 'protenix', 'boltz', 'rf3', 'af3'])
parser.add_argument('--target_merge', '-tm', help='Enable this flag to compare to alternate fragalysis conformations of the target (i.e. compare x0152a to x0152a, x0152b, x0152c, etc)', default=False, action='store_true')

args = parser.parse_args()

# Identify ground-truth ligands that map back to of3 ligands
# Canonical RDKit smiles are used to determine matches between
# molecules
def check_lig_match(gt_ligs, gt_smis, gt_lines, of3_lig):
    fail_log = []

    of_mol = Chem.MolFromMolFile(of3_lig)
    of_smi = Chem.MolToSmiles(of_mol)
    
    tmp_lines = []
    match_ligs = []
    for i, smi in enumerate(gt_smis):
        if smi == of_smi:
            print('\t', i, smi, of_smi)
            match_ligs.append(gt_ligs[i])
            tmp_lines += gt_lines[gt_ligs[i]]

    if len(tmp_lines) == 0:
        fail_log.append(f'{of3_lig} failed!\n\tLig_SMI: {of_smi}\n\tGT_SMI_L: {" ".join(gt_smis)}\n')


    return match_ligs, tmp_lines, fail_log

def check_lig_match_resn(gt_ligs, gt_lines, of3_lig):
    of3_ligname = os.path.basename(of3_lig)
    of3_lig_resn = of3_ligname.split('_')[-1].split('-')[0]

    tmp_lines = []
    match_ligs = []
    fail_log = []
    
    #print('\t', of3_lig_resn)
    for gt_sdf in gt_lines:
        gt_lig_resn = gt_sdf.split('_')[-1].split('-')[0]
        
        # Custom ligand matching (assumes only one ligand of interest!!)
        if (args.cofolding_model == 'of3') or (args.cofolding_model == 'boltz') or (args.cofolding_model == 'af3'):
            if of3_lig_resn.startswith('LIG') and gt_lig_resn == 'LIG':
                match_ligs.append(gt_sdf)
                tmp_lines += gt_lines[gt_sdf]
        if args.cofolding_model == 'protenix':
            if of3_lig_resn.startswith('l0') and gt_lig_resn == 'LIG':
                match_ligs.append(gt_sdf)
                tmp_lines += gt_lines[gt_sdf]
        if args.cofolding_model == 'rf3':
            if of3_lig_resn.startswith('L:') and gt_lig_resn == 'LIG':
                match_ligs.append(gt_sdf)
                tmp_lines += gt_lines[gt_sdf]

        if of3_lig_resn == gt_lig_resn:
            match_ligs.append(gt_sdf)
            tmp_lines += gt_lines[gt_sdf]


    if len(tmp_lines) == 0:
        fail_log.append(f'{of3_lig} failed!\n\tLig_RESN: {of3_lig_resn}\n\tGT_LIG_RESN: {" ".join(gt_lines)}\n')
        
    return match_ligs, tmp_lines, fail_log


    '''
    # Compile all ground truth ligands into one sdf
    tmplines = []
    for fl in fragalysis_ligs:
        with open(fl) as f:
            fl_lines = f.readlines()
        tmplines += fl_lines
    
    with open('fl_tmp.sdf', 'w') as fo:
        fo.write(''.join(tmplines))
    '''

def mp_func(case):
    seeds = [1370180479, 1449838082, 1832854922, 1880307061, 2012026466]
    # Receptor file naming seems to be inconsistent between fragalysis systems :\
    fragalysis_rec = f'{args.fragalysis_dir}/{case}/{case}.pdb'

    if args.target_merge:
        case_prefix = case[:-1] # Remove suffix at the end
        fragalysis_ligs = glob.glob(f'{args.fragalysis_dir}/{case_prefix}*/ligand_sdfs/*.sdf')
    else:
        lig_path = f'{args.fragalysis_dir}/{case}/ligand_sdfs/'
        fragalysis_ligs = glob.glob(f'{lig_path}/*.sdf')
    
    #print(case)
    #print(fragalysis_ligs)

    fragalysis_lines = {}
    for sdf in fragalysis_ligs:
        with open(sdf) as f:
            fragalysis_lines[sdf] = f.readlines()
            

    for s in seeds:
        results_dir = f'{args.of3_results}/{case}/seed_{s}/'
        model_pdbs = glob.glob(f'{results_dir}/*_model_rec.pdb')
        case_outdir = f'{args.outdir}/{case}/seed_{s}'
        os.makedirs(case_outdir, exist_ok=True)

        for m in model_pdbs:
            m_name = os.path.basename(m).strip('_rec.pdb')

            lig_sdfs = glob.glob(f'{results_dir}/{m_name}*lig.sdf')
            #continue #Debug

            #print(m_name, lig_sdfs)

            for ml in lig_sdfs:
                ml_name = os.path.basename(ml).strip('.sdf')
                outfile = os.path.abspath(f'{case_outdir}/ost-{ml_name}.json')
                if os.path.exists(outfile):
                    continue

                matching_fragalysis, tmp_lines, fail_log = check_lig_match_resn(fragalysis_ligs, fragalysis_lines, ml)

                #print(ml)
                #print('\t', matching_fragalysis)

                if len(tmp_lines) == 0:
                    with open(f'{args.outdir}/ost_fails.log', 'a') as fo:
                        fo.write(''.join(fail_log))
                
                # Save temporary file with all matching ground truth ligands
                tmp_sdf = f'{case_outdir}/fl_tmp.sdf'
                with open(tmp_sdf, 'w') as f:
                    f.write(''.join(tmp_lines))

                ost_cmd = f'ost compare-ligand-structures -m {os.path.abspath(m)} -ml {os.path.abspath(ml)} -r {os.path.abspath(fragalysis_rec)} -rl {os.path.abspath(tmp_sdf)} --lddt-pli --rmsd -o {outfile} -v 0'.split()
                subprocess.run(ost_cmd)
                os.remove(tmp_sdf)
    
        
def main():
    os.makedirs(args.outdir, exist_ok=True)
        
    case_l = []
    
    for ff in os.listdir(args.of3_results):
        if os.path.isdir(os.path.join(args.of3_results,ff)):
            if os.path.exsts(fragalysis_rec = f'{args.fragalysis_dir}/{ff}/{ff}.pdb'):
                case_l.append(ff)
                print(ff)

    #with open(args.of3_json) as f:
    #    of3_inps = json.load(f)

    print(f'OST-lig eval for {len(case_l)} targets...')
    
    # Quick multiprocessing implementation
    with mp.Pool(mp.cpu_count()) as pool:
        r = list(tqdm.tqdm(pool.imap(mp_func, case_l, chunksize=1)))
    

if __name__=='__main__':
    main()
