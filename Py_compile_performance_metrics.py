import os
import tqdm
import glob
import json
import argparse
import numpy as np
import pandas as pd
from rdkit import Chem
from posebusters import PoseBusters

parser = argparse.ArgumentParser()

parser.add_argument('--result_dir', '-r', help='Directory with cofolding results, in OF3 format', required=True)
parser.add_argument('--ost_ligand', '-ol', help='Directory with ost ligand results', required=True)
parser.add_argument('--ost_receptor', '-or', help='Directory with ost receptor results', required=True)
#parser.add_argument('--fragalysis_dir', '-f', help='Path to aligned/ fragalysis folder', required=True)
parser.add_argument('--similarity_tsv', '-st', help='.tsv file with sucos_pocket_qcov similarity data', required=True)
parser.add_argument('--method', '-m', choices=['of3p', 'rf3', 'protenix', 'boltz-1', 'boltz-2', 'af3'], help='Name of the method that was benchmarked', required=True)
parser.add_argument('--outfile', '-o', help='Name of the output .tsv file', required=True)

args = parser.parse_args()

ALPHABET=list('ABCDEFGHIJKLMNOPQRSTUVWXYZ')

def parse_rec_ost(ost_f):
    with open(ost_f) as f:
        ost_data = json.load(f)
    
    if ost_data['status'] == 'FAILURE':
        return None, None, None, None

    rec_rmsd = ost_data['rmsd']
    rec_lddt = ost_data['lddt']
    tm_score = ost_data['tm_score']
    chain_mapping = ost_data['chain_mapping'] # ref_chain: model_chain
    
    # Convert chain mapping to a string
    chain_map_str = []
    for r_ch in chain_mapping:
        m_ch = chain_mapping[r_ch]

        chain_map_str.append(f'{r_ch}:{m_ch}')

    chain_map_str = ','.join(chain_map_str)
    
    return rec_rmsd, rec_lddt, tm_score, chain_map_str

# Read confidence metrics file
def parse_confidence_metrics(conf_data, method, lig_ch, lig_rmsd_chain_mapping):

    iptm = conf_data['iptm']
    pair_iptm = None # Ignore this for now
    
    # Convert chain mapping string to list
    try:
        rec_lig_ch_l = []
        for ch_pair in lig_rmsd_chain_mapping.split(','):
            mdl_ch = ch_pair.split(':')[1]
            rec_lig_ch_l.append(mdl_ch)
    except: #Fails if the ligand does not form a contact with the receptor
        pair_iptm = None
        return iptm, pair_iptm


    if method == 'of3p':
        pair_iptm_l = []
        for ch in rec_lig_ch_l:
            ch_pair = f'({ch}, {lig_ch})'
            pair_iptm_l.append(conf_data['chain_pair_iptm'][ch_pair])

        pair_iptm = np.average(pair_iptm_l)
        #print('\tiptm:', iptm, 'pair_iptm:', pair_iptm)

    if method == 'protenix':
        lig_ch_idx = ALPHABET.index(lig_ch)
        pair_iptm_l = []
        for ch in rec_lig_ch_l:
            ch_idx = ALPHABET.index(ch)
            pair_val = conf_data["chain_pair_iptm"][lig_ch_idx][ch_idx]
            pair_iptm_l.append(pair_val)
            #print(f'\t{lig_ch}:{lig_ch_idx}\t{ch}:{ch_idx}\t{pair_val}')

        pair_iptm = np.average(pair_iptm_l)
        #print(f'\t{lig_ch}:{lig_ch_idx}\t{pair_iptm}')

    if method == 'rf3':
        pair_iptm = iptm # RF3 does not report a pair iptm :(
        
    if method in ['boltz-1', 'boltz-2']:
        lig_ch_idx = str(ALPHABET.index(lig_ch))
        pair_iptm_l = []
        for ch in rec_lig_ch_l:
            ch_idx = str(ALPHABET.index(ch))

            # iptm scores are not symmetric in Boltz :(
            # Get the max value for a pair, and use that
            val1 = conf_data["pair_chains_iptm"][lig_ch_idx][ch_idx]
            val2 = conf_data["pair_chains_iptm"][ch_idx][lig_ch_idx]
            
            pair_val = max([val1, val2])
            pair_iptm_l.append(pair_val)

        pair_iptm = np.average(pair_iptm_l)
        

    #print('\tiptm:', iptm, 'pair_iptm:', pair_iptm)
    return iptm, pair_iptm


def parse_lig_ost(ost_f):
    with open(ost_f) as f:
        ost_data = json.load(f)
    
    if ost_data['status'] == 'FAILURE':
        return None, None, None, None, None 
    
    try:
        lddt_pli = ost_data['lddt_pli']['assigned_scores'][0]['score']
        lig_rmsd = ost_data['rmsd']['assigned_scores'][0]['score']
        rmsd_chain_mapping = ost_data['rmsd']['assigned_scores'][0]['chain_mapping'] #ref_chain: model_chain
        lddt_lp = ost_data['rmsd']['assigned_scores'][0]['lddt_lp']
        pocket_bb_rmsd = ost_data['rmsd']['assigned_scores'][0]['bb_rmsd']
    except:
        return None, None, None, None, None

    # Convert chain mapping to a string
    chain_map_str = []
    for r_ch in rmsd_chain_mapping:
        m_ch = rmsd_chain_mapping[r_ch]

        chain_map_str.append(f'{r_ch}:{m_ch}')

    chain_map_str = ','.join(chain_map_str)

    return lig_rmsd, lddt_pli, lddt_lp, pocket_bb_rmsd, chain_map_str

def check_posebusters(mdl_lig, mdl_rec):
    buster = PoseBusters(config="dock")
    df = buster.bust([mdl_lig], None, mdl_rec)

    pb_valid = True
    for h in df.head():
        if df[h].iloc[0] == False:
            pb_valid = False

    return pb_valid

def read_training_similarity(sim_tsv):
    df = pd.read_csv(sim_tsv, delimiter='\t')
    sim_data = {}

    for i, case in enumerate(df['query']):
        sim_data[case] = {}
        for metric in 'sucos_shape pocket_qcov sucos_shape_pocket_qcov'.split():
            sim_data[case][metric] = df[metric].iloc[i]
            
    return sim_data

def main():
    
    sim_data = read_training_similarity(args.similarity_tsv)

    case_data = {}
    #outlines = ['target\tseed\tsample\tmethod\tlig_id\tlig_ref_ch\tlig_mdl_ch\tlig_rmsd\tlddt_pli\tlddt_lp\tbb_rmsd_lp\t']
    #sucos_shape pocket_qcov sucos_shape_pocket_qcov
    outlines = ['target\tseed\tsample\tmethod\tlig_id\tpb_valid\tiptm\tpair_iptm\tis_succ\tis_proper\tlig_rmsd_ch_mapping\tlig_rmsd\tlddt_pli\tlddt_lp\tbb_rmsd_lp\trec_rmsd\trec_lddt\ttm_score\trec_ch_mapping\tsucos_shape\tpocket_qcov\tsucos_shape_pocket_qcov\tligand_smiles']
    err_log = []
    # Parse all cases
    for case in tqdm.tqdm(os.listdir(args.ost_receptor)):
        #print(case)
        if case not in case_data:
            case_data[case] = {}
        
        for seed in os.listdir(f'{args.ost_receptor}/{case}/'):
            #print('\t', seed)

            if seed not in case_data[case]:
                case_data[case][seed] = {}
            
            # Read confidence score jsons and store info, with sample number as the key
            conf_dict = {}
            conf_files = glob.glob(f'{args.result_dir}/{case}/{seed}/*confidences_aggregated.json')
            for cf in conf_files:
                sample_id = os.path.basename(cf).split('_')[-3]
                with open(cf) as f:
                    data = json.load(f)
                    conf_dict[sample_id] = data

            for ost_result in os.listdir(f'{args.ost_receptor}/{case}/{seed}/'):
                sample = ost_result.split('_')[-2]
                #print('\t\t', ost_result, sample)

                if sample not in case_data[case][seed]:
                    case_data[case][seed][sample] = {}

                ost_result_path = f'{args.ost_receptor}/{case}/{seed}/{ost_result}'


                lig_ost_results = glob.glob(f'{args.ost_ligand}/{case}/{seed}/ost-{case}_{seed}_sample_{sample}_model_*.json')
                lig_ost_results.sort()
                
                # Collect posebusters data for the system
                pb_file = f'{args.result_dir}/{case}/{seed}/posebusters_data.json'
                with open(pb_file) as f:
                    pb_data = json.load(f)

                rec_rmsd, rec_lddt, tm_score, rec_ch_mapping = parse_rec_ost(ost_result_path)

                if rec_rmsd == None:
                        err_log.append(f'ERR_OST-RECEPTOR: ost metric failure for {case} {seed} {sample}')
                
                for l_ost in lig_ost_results:
                    #print('\t', l_ost)
                    is_proper = False
                    lig_id = os.path.basename(l_ost).split('_')[-1][:-5]
                    
                    #mdl_rec = f'{args.result_dir}/{case}/{seed}/{case}_{seed}_sample_{sample}_model_rec.pdb'
                    #pb_valid = check_posebusters(mdl_lig, mdl_rec)
                    mdl_lig_name = f'{case}_{seed}_sample_{sample}_model_{lig_id}.sdf'
                    mdl_lig = f'{args.result_dir}/{case}/{seed}/{case}_{seed}_sample_{sample}_model_{lig_id}.sdf'
                    pb_valid = pb_data[mdl_lig_name]
                    lig_mol = Chem.MolFromMolFile(mdl_lig)
                    lig_smi = Chem.MolToSmiles(lig_mol)
                    

                    # Weird janky way to check if the ligand of interest is in this file, but
                    # it should work
                    if (lig_id.startswith('LIG')) or (lig_id.startswith('l0')) or (lig_id.startswith('L:')):
                        is_proper = True

                    is_succ = False
                    lig_rmsd, lddt_pli, lddt_lp, pocket_bb_rmsd, lig_rmsd_chain_mapping = parse_lig_ost(l_ost)
                    
                    if lig_rmsd != None:
                        if (lig_rmsd <= 2.0) and (lddt_pli >= 0.8):
                            is_succ = True
                    else:
                        err_log.append(f'ERR_OST-LIGAND: ost metric failure for {case} {seed} {sample} {lig_id}')
                    
                    
                    # Extract confidence metrics
                    lig_chain = lig_id.split('-')[2]
                    iptm, pair_iptm = parse_confidence_metrics(conf_dict[sample], args.method, lig_chain, lig_rmsd_chain_mapping)

                    # Save output data
                    outline = f'{case}\t{seed}\t{sample}\t{args.method}'
                    outline += f'\t{lig_id}\t{pb_valid}\t{iptm}\t{pair_iptm}\t{is_succ}\t{is_proper}\t{lig_rmsd_chain_mapping}'
                    outline += f'\t{lig_rmsd}\t{lddt_pli}\t{lddt_lp}\t{pocket_bb_rmsd}\t{rec_rmsd}\t{rec_lddt}\t{tm_score}\t{rec_ch_mapping}'
                    
                    if is_proper:
                        outline += f'\t{sim_data[case]["sucos_shape"]}\t{sim_data[case]["pocket_qcov"]}\t{sim_data[case]["sucos_shape_pocket_qcov"]}\t{lig_smi}'
                    else:
                        outline += f'\t{None}\t{None}\t{None}\t{lig_smi}'

                    #print('\t\t\t', lig_id, lig_rmsd, lddt_pli, lig_rmsd_chain_mapping, is_succ, is_proper)
                    outlines.append(outline)
        #break # Debug
    
    # Save output
    with open(args.outfile, 'w') as fo:
        fo.write('\n'.join(outlines))

    # Save error log
    with open(args.outfile[:-4] + '.err', 'w') as fo:
        fo.write('\n'.join(err_log))


if __name__=='__main__':
    main()
