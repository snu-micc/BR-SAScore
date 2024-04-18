#
# This python script is modified from rdchiral template extractor 
# https://github.com/connorcoley/rdchiral/blob/master/rdchiral/template_extractor.py
#

import re
from rdkit import Chem
from rdkit.Chem import AllChem

def clean_map_and_sort(smiles_list, no_clean_numbers = [], return_mols = False):
    mols = []
    for smiles in smiles_list:
        if not smiles: continue
        mol = Chem.MolFromSmiles(smiles)
        [atom.SetAtomMapNum(0) for atom in mol.GetAtoms() if atom.GetAtomMapNum() not in no_clean_numbers]
        mols.append(Chem.MolFromSmiles(Chem.MolToSmiles(mol)))
    mols = sorted(mols, key= lambda m: m.GetNumAtoms(), reverse=True)
    if return_mols:
        return mols                                                                                
    else:
        return [Chem.MolToSmiles(mol) for mol in mols]

def replace_deuterated(smi):
    return re.sub('\[2H\]', r'[H]', smi)

def clear_mapnum(mol):
    [a.ClearProp('molAtomMapNumber') for a in mol.GetAtoms() if a.HasProp('molAtomMapNumber')]
    return mol

def get_tagged_atoms_from_mols(mols):
    '''Takes a list of RDKit molecules and returns total list of
    atoms and their tags'''
    atoms = []
    atom_tags = []
    for mol in mols:
        new_atoms, new_atom_tags = get_tagged_atoms_from_mol(mol)
        atoms += new_atoms 
        atom_tags += new_atom_tags
    return atoms, atom_tags

def get_tagged_atoms_from_mol(mol):
    '''Takes an RDKit molecule and returns list of tagged atoms and their
    corresponding numbers'''
    atoms = []
    atom_tags = []
    for atom in mol.GetAtoms():
        if atom.HasProp('molAtomMapNumber'):
            atoms.append(atom)
            atom_tags.append(str(atom.GetProp('molAtomMapNumber')))
    return atoms, atom_tags

def atoms_are_different(atom1, atom2): 
    '''Compares two RDKit atoms based on basic properties'''

    if atom1.GetAtomicNum() != atom2.GetAtomicNum(): return True # must be true for atom mapping
    
    if atom1.GetNumRadicalElectrons() != atom2.GetNumRadicalElectrons(): return True
    if atom1.GetFormalCharge() != atom2.GetFormalCharge(): return True
    if atom1.GetTotalNumHs() != atom2.GetTotalNumHs(): return True
    if atom_neighbors(atom1) != atom_neighbors(atom2): return True 
    
    # change bonds
    bonds1 = sorted([bond_to_smarts(bond) for bond in atom1.GetBonds()]) 
    bonds2 = sorted([bond_to_smarts(bond) for bond in atom2.GetBonds()]) 
    if bonds1 != bonds2: return True

    return False

def find_map_num(mol, mapnum):
    return [(a.GetIdx(), a) for a in mol.GetAtoms() if a.HasProp('molAtomMapNumber') 
         and a.GetProp('molAtomMapNumber') == str(mapnum)][0]

def atom_neighbors(atom):
    neighbor = []
    for n in atom.GetNeighbors():
        neighbor.append(n.GetAtomMapNum())
    return sorted(neighbor)

def bond_to_smarts(bond):
    '''This function takes an RDKit bond and creates a label describing
    the most important attributes'''
    a1_label = str(bond.GetBeginAtom().GetAtomicNum())
    a2_label = str(bond.GetEndAtom().GetAtomicNum())
    if bond.GetBeginAtom().HasProp('molAtomMapNumber'):
        a1_label += bond.GetBeginAtom().GetProp('molAtomMapNumber')
    if bond.GetEndAtom().HasProp('molAtomMapNumber'):
        a2_label += bond.GetEndAtom().GetProp('molAtomMapNumber')
    atoms = sorted([a1_label, a2_label])
    bond_smarts = bond.GetSmarts()
    if bond_smarts == '':
        bond_smarts = '-'
    
    return '{}{}{}'.format(atoms[0], bond_smarts, atoms[1])

def get_changed_atoms(reactants, products, radius=2):
    '''Looks at mapped atoms in a reaction and determines which ones changed'''

    err = 0
    prod_atoms, prod_atom_tags = get_tagged_atoms_from_mols(products)
    reac_atoms, reac_atom_tags = get_tagged_atoms_from_mols(reactants)

    # Find differences 
    changed_atoms = [] # actual reactant atom species
    changed_atom_tags = [] # atom map numbers of those atoms
    changed_atom_idxs = []

    # Product atoms that are different from reactant atom equivalent
    for i, prod_tag in enumerate(prod_atom_tags):
        for j, reac_tag in enumerate(reac_atom_tags):
            if reac_tag != prod_tag: continue
            if reac_tag not in changed_atom_tags: # don't bother comparing if we know this atom changes
                # If atom changed, add
                prod_atom, reac_atom = prod_atoms[i], reac_atoms[j]
                if atoms_are_different(prod_atom, reac_atom):
                    changed_atoms.append(prod_atom)
                    changed_atom_tags.append(reac_tag)
                    changed_atom_idxs.append(prod_atom.GetIdx())
                    break
            
    radius_atoms_idx = {idx:0 for idx in changed_atom_idxs}
    for k in range(radius):
        changed_atoms, changed_atom_idxs = expand_neighbors(changed_atoms, changed_atom_idxs)
        for new_idx in changed_atom_idxs:
            if new_idx not in radius_atoms_idx:
                radius_atoms_idx[new_idx] = k+1
    
    return changed_atom_tags, radius_atoms_idx, err

def expand_neighbors(changed_atoms, changed_atom_idxs):
    new_changed_atoms = changed_atoms[:]
    new_changed_atom_idxs = changed_atom_idxs[:]
    for atom in changed_atoms:
        for natom in atom.GetNeighbors():
            natom_idx = natom.GetIdx()
            if natom_idx not in changed_atom_idxs:
                new_changed_atoms.append(natom)
                new_changed_atom_idxs.append(natom_idx)
    return new_changed_atoms, new_changed_atom_idxs
    
def expand_atoms_to_use(mol, atoms_to_use, groups=[], symbol_replacements=[]):
    '''Given an RDKit molecule and a list of AtomIdX which should be included
    in the reaction, this function expands the list of AtomIdXs to include one 
    nearest neighbor with special consideration of (a) unimportant neighbors and
    (b) important functional groupings'''
    
    new_atoms_to_use = atoms_to_use[:]
    
    radius_atoms = {idx:0 for idx in atoms_to_use}
    # Look for all atoms in the current list of atoms to use
    for atom in mol.GetAtoms():
        if atom.GetIdx() not in atoms_to_use: continue
        # Ensure membership of changed atom is checked against group
        for group in groups:
            if int(atom.GetIdx()) in group[0]:
                if VERBOSE: 
                    print('adding group due to match')
                    try:
                        print('Match from molAtomMapNum {}'.format(
                            atom.GetProp('molAtomMapNumber'),
                        ))
                    except KeyError:
                        pass
                for idx in group[1]:
                    if idx not in atoms_to_use:
                        new_atoms_to_use.append(idx)
                        symbol_replacements.append((idx, convert_atom_to_wildcard(mol.GetAtomWithIdx(idx))))
        # Look for all nearest neighbors of the currently-included atoms
        for r in range(radius):
            for neighbor in atom.GetNeighbors():
                # Evaluate nearest neighbor atom to determine what should be included
                new_atoms_to_use, symbol_replacements = \
                        expand_atoms_to_use_atom(mol, new_atoms_to_use, neighbor.GetIdx(), 
                            groups=groups, symbol_replacements=symbol_replacements)
                
            for new_atom in new_atoms_to_use:
                if new_atom not in radius_atoms or r < radius_atoms[new_atom]:
                    radius_atoms[new_atom] = r
    return radius_atoms, symbol_replacements
    
def split_reagents(reaction):
    rs, ps = replace_deuterated(reaction['reactants']).split('.'), replace_deuterated(reaction['products']).split('.')
    reagents = [smiles for smiles in rs if smiles in ps]
    return [r for r in rs if r not in reagents], [p for p in ps if p not in reagents], reagents                                                     
                                                                                        
def extract_from_reaction(reaction):
    if type(reaction) == type('string'):
        reaction = {'reactants': reaction.split('>>')[0], 'products': reaction.split('>>')[1], '_id' : 0}
    reactants_list, products_list, reagents_list = split_reagents(reaction)
    product_maps = [atom.GetAtomMapNum() for products in products_list for atom in Chem.MolFromSmiles(products).GetAtoms()]
    products = clean_map_and_sort(products_list, product_maps, return_mols = True)
    reactants = clean_map_and_sort(reactants_list, product_maps, return_mols = True)
    changed_atom_tags, radius_atoms_idx, err = get_changed_atoms(reactants, products)
    
    if err: 
        return None, None
    
    else:
        return changed_atom_tags, radius_atoms_idx