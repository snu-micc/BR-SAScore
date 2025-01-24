#
# calculation of synthetic accessibility score as described in:
#
# Estimation of Synthetic Accessibility Score of Drug-like Molecules based on Molecular Complexity and Fragment Contributions
# Peter Ertl and Ansgar Schuffenhauer
# Journal of Cheminformatics 1:8 (2009)
# http://www.jcheminf.com/content/1/1/8
#
# several small modifications to the original paper are included
# particularly slightly different formula for marocyclic penalty
# and taking into account also molecule symmetry (fingerprint density)
#
# for a set of 10k diverse molecules the agreement between the original method
# as implemented in PipelinePilot and this implementation is r2 = 0.97
#
# peter ertl & greg landrum, september 2013
#

import math
import os.path as op
import matplotlib.cm as cm
import pickle, gzip
from collections import defaultdict
from pkg_resources import resource_filename

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

def numBridgeheadsAndSpiro(mol):
    nSpiro = rdMolDescriptors.CalcNumSpiroAtoms(mol)
    nBridgehead = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
    return nBridgehead, nSpiro

def numMacroAndMulticycle(mol, nAtoms):
    ri = mol.GetRingInfo()
    nMacrocycles = 0
    multi_ring_atoms = {i:0 for i in range(nAtoms)}
    for ring_atoms in ri.AtomRings():
        if len(ring_atoms) > 6:
            nMacrocycles += 1
        for atom in ring_atoms:
            multi_ring_atoms[atom] += 1
    nMultiRingAtoms = sum([v-1 for k, v in multi_ring_atoms.items() if v > 1])
    return nMacrocycles, nMultiRingAtoms

class SAScorer():
    def __init__(self, reaction_from='uspto', buildingblock_from='emolecules', frag_penalty=-6.0, complexity_buffer=1.0):
        if reaction_from == 'uspto' and buildingblock_from == 'emolecules':
            pickle_path = resource_filename('BRSAScore', 'pickle/BRScores_%s_%s.pkl.gz' % (reaction_from, buildingblock_from))
        else:
            pickle_path = 'BRSAScore/pickle/BRScores_%s_%s.pkl.gz' % (reaction_from, buildingblock_from)
        self._fscores = pickle.load(gzip.open(pickle_path))
        self.frag_penalty = frag_penalty
        self.max_score = 0
        self.min_score = frag_penalty-complexity_buffer
        
    def calculateScore(self, smi):
        sascore = 0        
        m = Chem.MolFromSmiles(smi)
        contribution = {}
        
        # fragment score
        bi = {}
        fp = rdMolDescriptors.GetMorganFingerprint(m, 2, useChirality=True, bitInfo=bi)

        fps = fp.GetNonzeroElements()
        score1 = 0.
        nf, nn = 0, 0
        for bitId, vs in bi.items():
            if vs[0][1] != 2:
                continue
            fscore = self._fscores.get(bitId, self.frag_penalty)
            if fscore < 0:
                nf += 1
                score1 += fscore
                for v in vs:
                    contribution[v[0]] = fscore
            if fscore == self.frag_penalty:
                nn += len(vs)
        if nf != 0:
            score1 /= nf
        sascore += score1

        # features score
        nAtoms = m.GetNumAtoms()
        nChiralCenters = len(Chem.FindMolChiralCenters(m, includeUnassigned=True))
        nBridgeheads, nSpiro = numBridgeheadsAndSpiro(m)
        nMacrocycles, nMulticycleAtoms = numMacroAndMulticycle(m, nAtoms)
            
        sizePenalty = nAtoms**1.005 - nAtoms
        stereoPenalty = math.log10(nChiralCenters + 1)
        spiroPenalty = math.log10(nSpiro + 1)
        bridgePenalty = math.log10(nBridgeheads + 1)
        macrocyclePenalty = math.log10(2) if nMacrocycles > 0 else 0
        multicyclePenalty = math.log10(nMulticycleAtoms + 1)

        score2 = 0. - sizePenalty - stereoPenalty - spiroPenalty - bridgePenalty - macrocyclePenalty - multicyclePenalty
        sascore += score2
           
        # correction for the fingerprint density
        # not in the original publication, added in version 1.1
        # to make highly symmetrical molecules easier to synthetise
        score3 = 0.
        if nAtoms > len(fps):
            score3 = math.log(float(nAtoms) / len(fps)) * .5
        sascore += score3

        # need to transform "raw" value into scale between 0 and 1
        sascore = (sascore - self.min_score) / (self.max_score - self.min_score)            
        if sascore > 1:
            sascore = 1
        elif sascore < 0.:
            sascore = 0
        sascore = 10 - sascore*9
        return sascore, contribution
    
    def contribution_to_svg(self, smiles, contribution):
        mol = Chem.MolFromSmiles(smiles)
        norm = cm.colors.Normalize(vmin=0, vmax=1)
        cmap = cm.get_cmap('OrRd') 
        plt_colors = cm.ScalarMappable(norm=norm, cmap=cmap)
        n_atoms = len(mol.GetAtoms())
        weights = [(-contribution.get(i, 0)/6) for i in range(n_atoms)]
        atom_colors = {i: plt_colors.to_rgba(w) for i, w in enumerate(weights)}
        rdDepictor.Compute2DCoords(mol)
        dr = rdMolDraw2D.MolDraw2DSVG(400, 370)
        do = rdMolDraw2D.MolDrawOptions()
        do.bondLineWidth = 4
        do.fixedBondLength = 30
        do.highlightRadius = 4
        mol = rdMolDraw2D.PrepareMolForDrawing(mol)
        dr.DrawMolecule(mol, highlightAtoms=range(n_atoms),
                            highlightBonds=[],
                            highlightAtomColors=atom_colors)
        dr.FinishDrawing()
        svg = dr.GetDrawingText()
        svg = svg.replace('svg:', '')
        return svg
    
#
#  Copyright (c) 2013, Novartis Institutes for BioMedical Research Inc.
#  All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#     * Neither the name of Novartis Institutes for BioMedical Research Inc.
#       nor the names of its contributors may be used to endorse or promote
#       products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#