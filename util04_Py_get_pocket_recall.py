import os
import tqdm
import argparse
import pandas as pd
from scipy.spatial import cKDTree
import biotite.structure.io.mol as mol
import biotite.structure.io.pdbx as pdbx
import biotite.structure.io.pdb as pdb

parser = argparse.ArgumentParser(description='Calculate pocket recall from the fragalysis cofolding benchmark data. Then, assign failure modes to is_proper ligands and create an updated metrics_tsv. NOTE: This code assumes that residue numbering is the same between the model and ground truth')

parser.add_argument('--metrics_tsv', '-t', help='tsv file with ost metrics compiled for all systems')
parser.add_argument('--result_dir', '-r', help='Cofolding results in OST format')
parser.add_argument('--fragalysis_dir', '-f', help='Path to aligned/ directory with fragalysis structures')
parser.add_argument('--outfile', '-o', help='(optional) Name of the outputput file (default = {args.metrics_tsv}_failure_modes.tsv', default=None)

args = parser.parse_args()

def ch_map_as_dict(ch_map_str):
    map_dict = {}
    for cmap in ch_map_str.split(','):
        ref_ch, mdl_ch = cmap.split(':')
        map_dict[ref_ch] = mdl_ch

    return map_dict

def get_pocket(rec_pdb, lig_sdf, ch_map, ref=True, cutoff=6.0):
    pocket_atom_list = []
    

    rec = pdb.PDBFile.read(rec_pdb)
    lig = mol.SDFile.read(lig_sdf)

    rec_atoms = pdb.get_structure(rec, model=1)
    lig_atoms = mol.get_structure(lig)

    lig_tree = cKDTree(lig_atoms.coord)
    
    # For each relevant chain, get residues in the pocket
    pocket = []
    for ref_ch in ch_map:
        mdl_ch = ch_map[ref_ch]
        if ref == True:
            ch = ref_ch
        else:
            ch = mdl_ch
        

        rec_ch_atoms = rec_atoms[rec_atoms.chain_id == ch]
        dists, _ = lig_tree.query(rec_ch_atoms.coord)
        mask = dists <= cutoff
        
        pocket_atoms = rec_ch_atoms[mask]
        seen = set()
        for i in range(len(pocket_atoms)):
            rid, rname = pocket_atoms.res_id[i], pocket_atoms.res_name[i]

            r_code = f'{ref_ch}.{rid}.{rname}' # Always use the reference chain name
            if r_code not in seen:
                seen.add(r_code)
                pocket.append(r_code)
    
    #print(pocket)
    return sorted(pocket, key=lambda x: x[0])

# Calc pocket recall: The fraction of true positive
# pocket residues recovered in the model
def calc_pocket_recall(ref_pocket, mdl_pocket):
    true_positives = list(set(ref_pocket) & set(mdl_pocket))
    true_positives = len(true_positives)
    n_positives = len(ref_pocket)
    pocket_recall = true_positives/n_positives
    
    return pocket_recall

def get_failure_mode(rmsd, lddt_pli, lddt_lp, pocket_recall):
    conf_fail = False
    pocket_recall_fail = False
    pose_fail = False

    if lddt_lp < 0.8:
        conf_fail = True
    if pocket_recall < 0.65:
        pocket_recall_fail = True
    if (rmsd > 2.0) or (lddt_pli < 0.8):   
        pose_fail = True

    return conf_fail, pocket_recall_fail, pose_fail


def main():
    if args.outfile == None:
        outfile = f'{args.metrics_tsv[:-4]}_failure_modes.tsv'
    else:
        outfile = args.outfile

    metric_df = pd.read_csv(args.metrics_tsv, delimiter='\t')
    
    df_head = metric_df.head()

    header = ''
    for h in df_head:
        header += f'{h}\t'

    header = header[:-1] + '\tpocket_recall\tconf_fail\tpocket_recall_fail\tpose_fail'
    outlines = [header]


    case_list = os.listdir(f'{args.fragalysis_dir}/')

    for i, target in enumerate(tqdm.tqdm(metric_df['target'])):
        seed = metric_df['seed'].iloc[i]
        sample = metric_df['sample'].iloc[i]
        lig_id = metric_df['lig_id'].iloc[i]
        is_proper = metric_df['is_proper'].iloc[i]
        
        #Append empty lines if it is not a main ligand
        if not is_proper:
            outstr = ''
            for h in df_head:
                outstr += f'{metric_df[h].iloc[i]}\t'

            outstr = outstr[:-1] + f'\t{None}\t{None}\t{None}\t{None}'
            outlines.append(outstr)
            continue
        
        ch_map_str = metric_df['lig_rmsd_ch_mapping'].iloc[i]

        try:
            ch_map = ch_map_as_dict(ch_map_str)
        except:
            outstr = ''
            for h in df_head:
                outstr += f'{metric_df[h].iloc[i]}\t'

            outstr = outstr[:-1] + f'\t{None}\t{None}\t{None}\t{None}'
            outlines.append(outstr)
            continue
        
        result_path = f'{args.result_dir}/{target}/{seed}'
        #print(i, target, seed, sample, lig_id, is_proper)
        #print(ch_map)
        #print(result_path)

        mdl_rec = f'{result_path}/{target}_{seed}_sample_{sample}_model_rec.pdb'
        mdl_lig = f'{result_path}/{target}_{seed}_sample_{sample}_model_{lig_id}.sdf'

        ref_rec = f'{args.fragalysis_dir}/{target}/{target}.pdb'
        ref_lig = f'{args.fragalysis_dir}/{target}/{target}_ligand.sdf'

        ref_pocket = get_pocket(ref_rec, ref_lig, ch_map, ref=True, cutoff=6.0)
        mdl_pocket = get_pocket(mdl_rec, mdl_lig, ch_map, ref=False, cutoff=6.0)
        
        pocket_recall = calc_pocket_recall(ref_pocket, mdl_pocket)

        conf_fail, pocket_fail, pose_fail = get_failure_mode(metric_df['lig_rmsd'].iloc[i], metric_df['lddt_pli'].iloc[i], metric_df['lddt_lp'].iloc[i], pocket_recall)

        #print('\t', pocket_recall, metric_df['lig_rmsd'].iloc[i], metric_df['lddt_pli'].iloc[i], metric_df['lddt_lp'].iloc[i])
        #print('\t', pocket_fail, pose_fail,'\t', conf_fail)

        #Add metrics to the output
        outstr = ''
        for h in df_head:
            outstr += f'{metric_df[h].iloc[i]}\t'

        outstr = outstr[:-1] + f'\t{pocket_recall}\t{conf_fail}\t{pocket_fail}\t{pose_fail}'
        outlines.append(outstr)
        #'\tpocket_recall\tconf_fail\tpocket_recall_fail\tpose_fail'

    with open(outfile, 'w') as fo:
        fo.write('\n'.join(outlines))

if __name__=='__main__':
    main()
