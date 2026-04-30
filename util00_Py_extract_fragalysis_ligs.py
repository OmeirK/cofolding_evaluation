import os
import glob
import shutil
import argparse
from pymol import cmd, stored

parser = argparse.ArgumentParser()
parser.add_argument('--fragalysis_dir', '-r', help='aligned/ directory where fragalysis results are stored, and updated ligand .sdf files will be saved')

args = parser.parse_args()

def main():
    case_l = []
    
    for ff in os.listdir(args.fragalysis_dir):
        if os.path.isdir(os.path.join(args.fragalysis_dir,ff)):
            case_l.append(ff)
            print(ff)
    

    for case in case_l:
        rec_pdb = f'{args.fragalysis_dir}/{case}/{case}.pdb'
        print(rec_pdb)
        
        cmd.reinitialize()
        cmd.load(rec_pdb)
        cmd.remove('solvent')

        stored.lig_data = []
        cmd.iterate('hetatm', 'stored.lig_data.append("_".join([resn, resi, chain]))')

        print(len(stored.lig_data))
        lig_data = list(set(stored.lig_data))
        print(lig_data)
        
        lig_out = f'{args.fragalysis_dir}/{case}/ligand_sdfs/'
        os.makedirs(lig_out, exist_ok=True)

        for info in lig_data:
            lign, ligi, lc = info.split('_')
            lig_sdf = f'{lig_out}/{case}_{lign}-{ligi}-{lc}-lig.sdf'

            if lign == 'LIG':
                print(f'\tCopy: {lig_sdf}')
                shutil.copy(f'{args.fragalysis_dir}/{case}/{case}_ligand.sdf', lig_sdf)
            else:
                cmd.save(lig_sdf, f'chain {lc} and resi {ligi} and resn {lign}')


if __name__=='__main__':
    main()
