from argparse import ArgumentParser
import gzip, pickle, os
import numpy as np
import pandas as pd
from tqdm import tqdm
from collections import defaultdict

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from reaction_center_extractor import extract_from_reaction

import gzip
import pickle

def get_scores(fp_cnts, data_num, common_ratio=0.001, min_cnt=1, max_score=3, min_score=-3):
    scores_dict = {}
    common_n = data_num*common_ratio
    max_value = np.log(max(fp_cnts.values())/common_n)
    min_value = np.log(1/common_n)
    for fp, cnt in fp_cnts.items():
        if cnt <= min_cnt:
            continue
        if cnt > common_n:
            score = np.log(cnt/common_n)*(max_score/max_value)
        else:
            score = np.log(cnt/common_n)*(min_score/min_value)
        scores_dict[fp] = score
    return scores_dict

def save_scores_to_pickle(scores_dict, output_path):
    with gzip.open(output_path, 'wb') as f:
        pickle.dump(scores_dict, f)
    return 

def get_RFrags(smi, changed_prod_idxs):
    mol = Chem.MolFromSmiles(smi)
    bi = {}
    RFrags = defaultdict(int)
    fps = rdMolDescriptors.GetMorganFingerprint(mol, 2, useChirality=True, bitInfo=bi)
    for fp, fs in bi.items():
        for (idx, radius) in fs:
            if radius != 2:
                break
            elif idx in changed_prod_idxs:
                RFrags[fp] += 2**(-changed_prod_idxs[idx])
    return RFrags

def get_BFrags(smi):
    mol = Chem.MolFromSmiles(smi)
    n_atoms = len(mol.GetAtoms())
    bi = {}
    BFrags = {}
    fps = rdMolDescriptors.GetMorganFingerprint(mol, 2, useChirality=True, bitInfo=bi)
    for fp, fs in bi.items():
        if fs[0][1] == 2:
            BFrags[fp] = len(fs)/n_atoms
    return BFrags

def prepare_RScores(args):
    df = pd.read_csv('%s' % args['Reaction_path'])
    rxns = df['reactants>reagents>production'].tolist()
    n_rxns = 0
    RFrags_all = defaultdict(int)
    for i, rxn in tqdm(enumerate(rxns), total=len(rxns)):
        changed_atom_tags, changed_prod_idxs = extract_from_reaction(rxn)
        if changed_prod_idxs:
            reactants, reagents, products = rxn.split('>')
            RFrags = get_RFrags(products, changed_prod_idxs)
            for fp, cnt in RFrags.items():
                RFrags_all[fp] += cnt
            n_rxns += 1
        if len(RFrags_all) >= 100:
            break
    RScores = get_scores(RFrags_all, n_rxns)
    save_scores_to_pickle(RScores, args['RScore_path']) 
    return RScores
    
def prepare_BScores(args):
    df = pd.read_csv('%s' % args['Buildingblock_path'])
    mols = df['SMILES'].tolist()
    n_mols = 0
    BFrags_all = defaultdict(int)
    for smi in tqdm(mols, total=len(mols)):
        if '[H]' not in smi and '*' not in smi:
            BFrags = get_BFrags(smi)
            for fp, cnt in BFrags.items():
                BFrags_all[fp] += cnt
            n_mols += 1
        if len(BFrags_all) >= 100:
            break
    BScores = get_scores(BFrags_all, n_mols)
    save_scores_to_pickle(BScores, args['BScore_path']) 
    return BScores

def get_max_value(dict1, dict2):
    for k, v in dict2.items():
        if k in dict1:
            dict1[k] = max([dict1[k], v])
        else:
            dict1[k] = v
    return dict1

def main(args):
    if os.path.exists(args['RScore_path']):
        print('Loading previously saved RScores...')
        RScores = pickle.load(gzip.open(args['RScore_path']))
    else:
        print('Preparing RScores from scratch...')
        RScores = prepare_RScores(args)
        
    if os.path.exists(args['BScore_path']):
        print('Loading previously saved BScores...')
        BScores = pickle.load(gzip.open(args['BScore_path']))
    else:
        print('Preparing BScores from scratch...')
        BScores = prepare_BScores(args)
        
    R2Scores = get_max_value(RScores, BScores)    
    save_scores_to_pickle(R2Scores, args['R2Score_path']) 
    return

if __name__ == '__main__':  
    parser = ArgumentParser()
    parser.add_argument('-r', '--Reaction-name', default='uspto', help='Path of reaction knowledge')
    parser.add_argument('-b', '--Buildingblock-name', default='emolecules',  help='Path of resource knowledge')
    args = parser.parse_args().__dict__
    
    for data in ['Reaction', 'Buildingblock']:
        args['%s_path' % data] = 'data/%s.csv' % args['%s_name' % data]
        args['%sScore_path' % data[0]] = 'R2SAScore/pickle/%sScores_%s.pkl.gz' % (data[0], args['%s_name' % data])
    args['R2Score_path'] = 'R2SAScore/pickle/R2SAScores_%s_%s.pkl.gz' % (args['Reaction_name'], args['Buildingblock_name'])
        
    main(args)
