#!/usr/bin/env python

# Author: Hongjun Bai @ VGS (Rolland Lab), MHRP, HJF (hbai@hivresearch.org)

import sys
import argparse
import collections
import math
import numpy
import itertools as it
import scipy.cluster.hierarchy as hac
from scipy.spatial.distance import squareform

import Bio.PDB
import warnings
from Bio import BiopythonWarning
warnings.simplefilter('ignore', BiopythonWarning)

# Number of atoms and number of neigbhors are censored at 98% percentile, if a 
# raw count >= 98%_percentile. weights_over_98%_percentile is set to threshold
th_natoms = 12.0
th_nnbs = 8.0

# The number of contact and number of neighbor residues for glycans on a 
# glycosylation site can be overwhelming when comparing with that for amino 
# acids, scale it down so that they are comparable to the numbers of a given
# amino acid. The empirical value is calculated based on a list of antibodies
# in  Bai, et. al, PLoS Comp. Bio. 2019.  The scale factor will make
# 98%_of_glycan_terms * scale == 98%_of_aa_terms
s_natoms = 4.7
s_nnbs = 2.5
roundit = lambda x: round(x*10)/10.0
#

AA3to1 = {'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
          'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
          'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
          'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
          'MSE': 'M', 'UNK': 'X'}
Glycans = set('NAG MAN BMA FUC BGC NDG FUL GAL SIA'.split())
AAs = set('ALA CYS ASP GLU PHE GLY HIS ILE LYS LEU MET ASN PRO GLN ARG SER THR VAL TRP TYR'.split())

def parse():
    parser = argparse.ArgumentParser(description="Get weight of epitope sites according to specified Ab:Ag structures")
    parser.add_argument('pdbfile', help='pdbfile of complex structure (.pdb, .pdb.gz, .cif, .cif.gz)')
    parser.add_argument('ag_chains', help='chain ids of antigen in the complex structure')
    parser.add_argument('ab_chains', help='chain ids of antibody in the complex structure')
    parser.add_argument('-range', help='start-end of residue id in the chimeric antigen chain (e.g. 123-194, for HIV Env in chain G of 5ESZ) (default: None)', default=None)
    parser.add_argument('--epsites', help='Print only epitope sites, instead of all sites (default: False)', default=False, action='store_true')
    parser.add_argument('-o', '--outfile', help='Output file (default: stdout)', default=sys.stdout, type=argparse.FileType('wt'))
    return parser.parse_args()

def main():
    para = parse()

    if para.range:
        pair = [int(x) for x in para.range.strip().split('-')]
        assert(len(pair) == 2)
        segment = pair
    else:
        segment = None

    parser = Bio.PDB.PDBParser()  # default: PDB format
    if '.cif' in para.pdbfile:  # Not tested yet
        parser = Bio.PDB.MMCIFParser()

    with zopen(para.pdbfile) as infile:
        structure = parser.get_structure('', infile)[0]

    w_natoms = extr_epitope_natoms(structure, para.ag_chains, para.ab_chains)
    w_nnbs = get_nnbs(structure, para.ag_chains, para.ab_chains, segment, 8.5)

    w_natoms = dict((rec[0], rec[1:]) for rec in w_natoms)
    weights = {}
    for resid, nnbs_aa, nnbs_glycan, nb_ag_sites in w_nnbs:
        res = get_by_full_id(structure, resid)
        if res.resname not in AA3to1:
            continue
        #
        if resid in w_natoms:
            natoms_aa, natoms_glycan = w_natoms[resid]
            #print(resid, nnbs_aa, nnbs_glycan, natoms_aa, natoms_glycan, nb_ag_sites)
        else:
            natoms_aa, natoms_glycan = 0, 0
        #
        w_natom = natoms_aa + natoms_glycan/s_natoms
        if w_natom > th_natoms:
            w_natom = th_natoms
        w_nnb = nnbs_aa + nnbs_glycan/s_nnbs
        if w_nnb > th_nnbs:
            w_nnb = th_nnbs
        #
        weights[resid] = (w_natom, w_nnb)
    #
    epsites = [resid for resid, (w_natom, w_nnb) in weights.items() if w_natom+w_nnb > 0]
    scale_w_natoms = sum(weights[rid][0] for rid in epsites) / len(epsites) * 2
    scale_w_nnbs = sum(weights[rid][1] for rid in epsites) / len(epsites) * 2
    # Output
    outs = para.outfile
    outs.write("#$ {}\n".format(' '.join(sys.argv)))
    #outs.write("site\tAA\tweight\tw_natom\tw_nnb\tneighbor_sites\n")
    outs.write("site\tAA\tweight\tw_natom\tw_nnb\n")
    for resid, nnbs_aa, nnbs_glycan, nb_sites in w_nnbs:
        res = get_by_full_id(structure, resid)
        if res.resname not in AA3to1:
            continue
        if para.epsites and (nnbs_aa == 0 and nnbs_glycan == 0 and resid not in w_natoms):
            continue
        #
        site = to_site(resid)
        aa = AA3to1[res.resname]
        w_natom, w_nnb = weights[resid]
        weight = w_natom/scale_w_natoms + w_nnb/scale_w_nnbs
        #outs.write(f'{site}\t{aa}\t{weight:.3f}\t{w_natom:.1f}\t{w_nnb:.1f}\t{nb_sites}\n')
        outs.write(f'{site}\t{aa}\t{weight:.3f}\t{w_natom:.1f}\t{w_nnb:.1f}\n')

#
EpNote = collections.namedtuple('EpNote', 'res_id w_aa w_glycan')
def extr_epitope_natoms(model, Ag_chains, Ab_chains, threshold=4.0, threshold_interactions=4.0, sugar_set=Glycans, whitelist=Glycans|AAs):
    """ Extract epitope sites from Ag_chains (#ab_atoms to aa and glycan are annotated)"""
    # Clean up the structure: remove H atoms + remove non-(AA/Glycan) residues
    clean_stru(model, whitelist)
    ag_atoms = [atom for chain in Ag_chains for atom in model[chain].get_atoms()]
    ab_atoms = [atom for chain in Ab_chains for atom in model[chain].get_atoms()]
    glycan_groups = group_glycan(model)
    glycan2ngs = dict((glycan, ngs) for ngs, glycans in glycan_groups.items() for glycan in glycans)
    ep_sites = get_ep_sites(ag_atoms, ab_atoms, glycan2ngs, threshold, sugar_set)
    # Calculate and assign natoms
    epitope_atoms = get_epitope_atoms(model, ep_sites, glycan_groups)
    nb_ab_atoms = get_nb_atoms(epitope_atoms, ab_atoms, threshold_interactions)
    # assign contacted atoms to site
    w_aa_dict = collections.defaultdict(set)
    w_glycan_dict = collections.defaultdict(set)
    for atom in epitope_atoms:
        aid = atom.get_full_id()
        if aid not in nb_ab_atoms: continue
        if atom.parent.resname in sugar_set:  # Sugar residues
            try:
                ngs_id = glycan2ngs[atom.parent.get_full_id()]
                w_glycan_dict[ngs_id].update(nb_ab_atoms[aid])
            except KeyError:
                warnings.warn('No bonded ASN found for glycan {}.{}!'.format(encode_resid(atom.parent.get_full_id()), resname))
                pass
        else:
            res_id = atom.parent.get_full_id()
            w_aa_dict[res_id].update(nb_ab_atoms[aid])
    epinfo_list = [EpNote(x, len(w_aa_dict[x]), len(w_glycan_dict[x])) for x in sorted(ep_sites, key=lambda x: (x[3][0], x[3][1]))]
    return epinfo_list

def get_epitope_atoms(model, ep_sites, glycan_groups):
    epitope_atoms = []
    for res_id in ep_sites: 
        epitope_atoms += get_by_full_id(model, res_id).child_list
        if res_id in glycan_groups:
            for g_id in glycan_groups[res_id]:
                epitope_atoms += get_by_full_id(model, g_id).child_list
    return epitope_atoms

def get_ep_sites(ag_atoms, ab_atoms, glycan2ngs, threshold=4.0, sugar_set=Glycans):
    epitope_atoms = select_by_dist(ag_atoms, ab_atoms, threshold)
    ep_sites = set()
    uniq_sites = set((atom.parent.get_full_id(), atom.parent.resname) for atom in epitope_atoms)
    for res_id, resname in uniq_sites:
        if resname in sugar_set:  # Sugar residues
            try:
                ep_sites.add(glycan2ngs[res_id])
            except KeyError:
                #warnings.warn('No bonded ASN found for glycan {}.{}!'.format(encode_resid(res_id), resname))
                pass
        else:
            ep_sites.add(res_id)
    return ep_sites

#
to_site = lambda resid: resid[2]+'.'+''.join(str(x) for x in resid[3]).replace(' ','')

def get_nnbs(model, Ag_chains, Ab_chains, segment=None, threshold=8.5, whitelist=Glycans|AAs):
    """ Get Ab nb sites for all ag sites """
    # Clean up the structure: remove H atoms + remove non-(AA/Glycan) residues
    clean_stru(model, whitelist)
    glycan_groups = group_glycan(model)
    glycan2ngs = dict((glycan, ngs) for ngs, glycans in glycan_groups.items() for glycan in glycans)
    glycan_chains = set(glycan[2] for ngs, glycans in glycan_groups.items() for glycan in glycans if ngs[2] in Ag_chains)
    # Get nb_dict
    all_chains = set(Ag_chains+Ab_chains) | glycan_chains
    residues = [r for r in model.get_residues() if r.parent.id in all_chains]
    #print(residues)
    nb_dict = get_delaunay_nbs(residues, threshold, scc=True)  # {res_full_id: [nb_full_id1, nb_full_id2, ...]}
    # Get result
    if segment is not None:
        start, end = segment
        within_range = lambda i: i >= start and i <= end
    note_nnbs = []
    for res in residues:
        if res.parent.id in Ab_chains: continue
        if (segment is not None) and (not within_range(res.id[1])): continue
        resid = res.get_full_id()
        # if glycosylated, count contacts with glycans
        nbs_glycan = set()
        if resid in glycan_groups:
            for glycan_id in glycan_groups[resid]:
                for nb_resid in nb_dict[glycan_id]:
                    if nb_resid[2] in Ab_chains:
                        nbs_glycan.add(nb_resid)
        nnbs_glycan = len(nbs_glycan)
        # get valid Ab+Ag neighbor sites, if segment is specified
        if segment is not None:
            all_current_nbs = []
            for rid in nb_dict[resid]:
                if rid[2] in Ab_chains:
                    all_current_nbs.append(rid)
                else:
                    if within_range(rid[3][1]):
                        all_current_nbs.append(rid)
        else:
            all_current_nbs = nb_dict[resid]
        # 
        nnbs_aa = len([rid for rid in all_current_nbs if rid[2] in Ab_chains])
        #nb_ag_sites = ','.join(to_site(rid) for rid in all_current_nbs if rid[2] in Ag_chains)
        nb_ag_sites = ','.join(to_site(rid) for rid in all_current_nbs)
        note_nnbs.append((resid, nnbs_aa, nnbs_glycan, nb_ag_sites))
    return note_nnbs

def get_delaunay_nbs(residues, threshold=8.5, scc=False):
    if scc:
        xyzs = get_sc_centers(residues)
    else:
        xyzs = get_cb_xyzs(residues)
    qh_out = DelaunayTHD(xyzs, threshold)
    edges = sorted(qh_out.edges)
    nb_dict = collections.OrderedDict()
    nb_dict.update((r.get_full_id(), []) for r in residues)
    for (i, j), d in edges:
        ri = residues[i]
        rj = residues[j]
        nb_dict[ri.get_full_id()].append(rj.get_full_id())
        nb_dict[rj.get_full_id()].append(ri.get_full_id())
    return nb_dict

BondThreshold = 2.0
def group_glycan(prot, sugar_set=Glycans, bond_length_th=BondThreshold):
    """ Group glycans in protein and return {NGS: glycan_list, ...}.

    Input:
    prot -- Bio.PDB stru/model/chain instance
    sugar_set -- list of sugar names
    bond_length_th -- lenght of bond to link two glycans together

    Output: {NGS_full_id: glycan_full_id_list, ...}, ordered Dictionary
    NGS_full_id -- N-linked glycosylation site (residue N), full_id by Bio.PDB.X.get_full_id()
    glycan_full_id_list -- A list of glycans attaced to NGS, full_id by Bio.PDB.X.get_full_id()
    """
    def find_ngs(atoms_prot, atoms_sugar, bond_length_th):
        """ Tested: get_nb_atoms() is slower than select_by_dist() """
        nd2_candidates = select_by_dist(atoms_prot, atoms_sugar, threshold=bond_length_th)
        ngs = None
        for atom in nd2_candidates:
            if atom.name == 'ND2':
                ngs = atom.parent
                break
        if ngs is None:
            glycans = set(a.parent.get_full_id() for a in atoms_sugar)
            glycan_names = ','.join(gid[-2] + ''.join(str(x).replace(' ', '') for x in gid[-1])
                                        for gid in glycans)
            #warnings.warn('No bonded ASN found for glycan(s): {}'.format(glycan_names))
        return ngs
    sugar_res = [res for res in prot.get_residues() if res.resname in sugar_set]
    groups = []
    if len(sugar_res) < 1:
        pass
    elif len(sugar_res) == 1:
        atoms_prot = [atom for atom in prot.get_atoms() if atom.parent.resname not in sugar_set]
        atoms_sugar = [atom for residue in sugar_res for atom in residue]
        ngs = find_ngs(atoms_prot, atoms_sugar, bond_length_th)
        if ngs is not None:
            groups.append((ngs, sugar_res))
    else:
        # Build distance matrix between glycan residues
        min_sqdist = lambda r_a, r_b: numpy.min(get_sqdist_mat(r_a.child_list, r_b.child_list))
        sqdist_mat = numpy.zeros((len(sugar_res), len(sugar_res)))
        for i in range(len(sugar_res)):
            for j in range(i+1, len(sugar_res)):
                sqdist_mat[i, j] = sqdist_mat[j, i] = min_sqdist(sugar_res[i], sugar_res[j])
        # Cluster glycans according to distance matrix
        pdist = squareform(sqdist_mat)  # transform distance matrix to condensed form
        z = hac.linkage(pdist, 'single')  # single linkage method to cluster all things together
        # extract cluster by a distance cutoff 2.0*2.0
        cls_ids = hac.fcluster(z, bond_length_th*bond_length_th, 'distance') 
        glycan_clusters = collections.defaultdict(list)
        for cls_id, sugar in zip(cls_ids, sugar_res):
            glycan_clusters[cls_id].append(sugar)
        # Get the NGS of each cluster
        atoms_prot = [atom for atom in prot.get_atoms() if atom.parent.resname not in sugar_set]
        for cls_id, glycans in glycan_clusters.items():
            atoms_sugar = [atom for residue in glycans  for atom in residue]
            ngs = find_ngs(atoms_prot, atoms_sugar, bond_length_th)
            if ngs is not None:
                groups.append((ngs, glycans))
    groups.sort(key=lambda x: x[0].id)
    result = collections.OrderedDict((ngs.get_full_id(), [glycan.get_full_id() for glycan in glycans]) for ngs, glycans in groups)
    return result  

# Utility functions
def zopen(fname, specifications='rt', msg=None):
    " open compressed files "
    if fname.endswith('.gz'):
        import gzip
        open_fun = gzip.open
    elif fname.endswith('.bz2'):
        import bz2
        open_fun = bz2.open
    else:
        open_fun = open
    try:
        return open_fun(fname, specifications)
    except IOError:
        print('File open error. %s cannot be opened' %(fname))
        if msg is not None: print(msg)
        sys.exit()

QhOut =  collections.namedtuple('QhOut', 'tetrahedra edges')
def DelaunayTHD(xyzs, threshold):
    """ Wrapper to use scipy.spatial.Delaunay

    Keyword arguments:
    xyzs -- [[x0, y0, z0], ...]
    threshold -- distance threshold to exclude tetrahedra

    Output:
    1, qhOut.tetraphedrons: a list of tetrahedrons, in which tetrahedron is recorded as a set of
       four ordered indices of vertexes, e.g. (Vidx0, Vidx1, Vidx2, Vidx3).  Index referred to
       coords array (zero-based, not residue ID).
    2, qhOut.edges: a list of edges, in which edges is recorded as a pair of vertexes,
       e.g. ((Vidx0, Vidx1), dist).  Indices are referred to coords array.

    """
    try:
        from scipy.spatial import Delaunay
    except ModuleNotFoundError:
        print("scipy.spatial.Delaunay cann't be imported.  Please check if scipy is installed.")
        sys.exit()
    # check if it is valid input
    if isinstance(xyzs, numpy.ndarray):
        working_xyzs = xyzs
    else:
        working_xyzs = numpy.array(xyzs)
    assert(len(working_xyzs.shape) == 2 and working_xyzs.shape[0] >= 4 and working_xyzs.shape[1] == 3)
    # 
    thd = Delaunay(working_xyzs)
    #
    # All uniq paris
    all_pairs = set()
    for ijkl in thd.simplices:
        for a, b in it.combinations(ijkl, 2):
            all_pairs.add((a, b) if a < b else (b, a))
    all_pairs = list(all_pairs)
    # Filter edges by distance
    xa = numpy.array([working_xyzs[a] for a, b in all_pairs])
    xb = numpy.array([working_xyzs[b] for a, b in all_pairs])
    dx = xa - xb
    distances = numpy.sqrt(numpy.sum(dx*dx, axis=1))
    valid_pairs = dict((pair, dist) for pair, dist in zip(all_pairs, distances) if dist < threshold)
    edges = sorted(valid_pairs.items(), key=lambda x: x[0])
    # Filter tetrahedra by edges
    tetrahedra = []
    for ijkl in thd.simplices:
        if all(pair in valid_pairs for pair in it.combinations(ijkl, 2)):
            tetrahedra.append(ijkl)
    return QhOut(tetrahedra, edges)

EntLevels='SMCRA'  # Bio.PDB organizes protein structure as SMCRA (Structure-Model-Chain-Residue-Atom)

def get_by_full_id(entity, full_id):
    """ Get elment in entity by full_id (entity could S/M/C/R/A) """
    def _get_by_corresponding_level_id(entity, entry_id):
        if len(entry_id) > 1:
            return _get_by_corresponding_level_id(entity[entry_id[0]], entry_id[1:])
        else:
            if entity.level == 'R':
                return entity[entry_id[0][0]]
            else:
                return entity[entry_id[0]]

    corresponding_level_id = full_id[EntLevels.index(entity.level)+1:]
    return _get_by_corresponding_level_id(entity, corresponding_level_id)

def encode_id(x_id):
    encoded = ';'.join(','.join(str(y) for y in x) if type(x) == tuple else str(x) for x in x_id)
    return f'({encoded})'

def encode_resid(resid):
    return encode_id(resid)

def get_nb_atoms(atoms, atoms_ref, threshold=4.0):
    """ Retrive IDs of atoms in "atoms_ref" within "threshold" of any atom in "atoms"

    Input:
    atoms -- atom whoes neighbor should be retrived, list of Bio.PDB atom instances
    atoms_ref -- neighbor atoms are subject to extraction, list of Bio.PDB atom instances
    threshold -- square of distance threshold to define neighbor atoms

    Output: {atom_id: id_nb_atom_in_ref, ...}
    atom_id  -- id of atoms have neighbors in atoms_ref
    id_nb_atom_in_ref -- list of ids of neighbor atoms in atoms_ref
    """
    id2coord = dict((atom.get_full_id(), atom.coord) for atom in atoms_ref)
    ref_frame = CubicIdx(atoms_ref, threshold)
    result = collections.OrderedDict()
    threshold2 = threshold * threshold
    result = {}
    for atom in atoms:
        potential_nb_ids = ref_frame.query(atom.coord, threshold)
        if len(potential_nb_ids) > 0:
            potential_nb_coords = numpy.array([id2coord[atom_id] for atom_id in potential_nb_ids])
            d_diff = potential_nb_coords - atom.coord
            dist2 = numpy.sum(d_diff * d_diff, axis=1)
            nb_ids = [atom_id for d2, atom_id in  zip(dist2, potential_nb_ids) if d2 < threshold2]
            if nb_ids:
                result[atom.get_full_id()] = nb_ids
    return result

class CubicIdx:
    """ Simple yet powerful spacial index """
    def __init__(self, atoms, grid):
        assert(grid > 0.0)
        self._grid = grid
        self._idx = collections.defaultdict(list)
        for atom in atoms:
            self._idx[self._index(atom.coord)].append(atom.get_full_id())

    def query(self, q_xyz, q_range):
        n = int(math.ceil(q_range/self._grid))
        search_range = (range(i-n, i+n+1) for i in self._index(q_xyz))
        search_grids = it.product(*search_range)
        result = sum((self._idx[grid] for grid in search_grids if grid in self._idx), [])
        return result

    def _index(self, xyz):
        return tuple(int(math.floor(x/self._grid)) for x in xyz)

def get_sc_centers(residues):
    """ Get center of side chain (CA is considered as side chain atom) """
    xyzs = []
    mc_atoms = set('C N O'.split())
    for r in residues:
        sc_coords = [atom.get_coord() for atom in r if atom.name not in mc_atoms]
        sc_enter = numpy.mean(sc_coords, axis=0)
        xyzs.append(sc_enter)
    return xyzs

def get_cb_xyzs(residues):
    xyzs = []
    for r in residues:
        if 'CB' in r:  # Use CB atoms when possible.
            xyzs.append(r['CB'].get_coord())
        elif 'CA' in r:  # If it is not possible (e.g. GLY, missing CB in resolved structure, etc), use CA atoms
            xyzs.append(r['CA'].get_coord())
        else:  # If both 'CA' and 'CB' are missing, use the center 
            coor = [0.0, 0.0, 0.0]
            n = len(r)
            for atom in r:
                current_coor = atom.get_coord()
                coor = [x+x_i/n for x, x_i in zip(coor, current_coor)]
            xyzs.append(coor)
            #xyzs.append(numpy.mean([atom.get_coord() for atom in r], axis=0))
    return xyzs

def clean_stru(struture, whitelist=Glycans|AAs):
    filter_res(struture, whitelist)
    remove_h(struture)

def filter_res(struture, whitelist=Glycans|AAs):
    for chain in Bio.PDB.Selection.unfold_entities(struture, 'C'):
        residues_to_be_removed = [res.id for res in chain if res.resname not in whitelist]
        for res_id in residues_to_be_removed:
            chain.detach_child(res_id)

def remove_h(struture):
    for residue in Bio.PDB.Selection.unfold_entities(struture, 'R'):
        atoms_to_be_removed = [atom.id for atom in residue if atom.element == 'H']
        for atom_id in atoms_to_be_removed:
            residue.detach_child(atom_id)  

def select_by_dist(atoms, atoms_ref, threshold=4.0):
    """ select atom in "atoms" within "threshold" of "atoms_ref" """
    sqdist_mat = get_sqdist_mat(atoms, atoms_ref)
    sqdist_min = numpy.min(sqdist_mat, axis=1)
    threshold2 = threshold * threshold
    result = [(atom, min_dist) for atom, min_dist in zip(atoms, sqdist_min) if min_dist < threshold2]
    result.sort(key=lambda atom_dist_pair:atom_dist_pair[1])
    return [atom_dist_pair[0] for atom_dist_pair in result]

def get_sqdist_mat(atoms_a, atoms_b):
    coordsa = numpy.array([atom.coord for atom in atoms_a])
    coordsb = numpy.array([atom.coord for atom in atoms_b])
    return sqdist(coordsa, coordsb)

def sqdist(xyza, xyzb):
    ''' Get the distance matrix between coords array xyza and xyzb.

    Input: 
        xyza: [[xa1, ya1, za1], [xa2, ya2, za2], ...]
        xyzb: [[xb1, yb1, zb1], [xb2, yb2, zb2], ...]

    Output:
        distmatrix: (an x bn)
        [[D_a1_b1, D_a1_b2, D_a1_b3, ..., D_a1_bn], 
         [D_a2_b1, D_a2_b2, D_a2_b3, ..., D_a2_bn], 
         .
         .
         .
         [D_an_b1, D_an_b2, D_an_b3, ..., D_an_bn], 
    '''
    sizea = xyza.shape[0]
    sizeb = xyzb.shape[0]
    mat_a = xyza.reshape(sizea, 1, 3)
    mat_a = mat_a.repeat(sizeb, axis=1)
    # mat_a:
    # [[[xa1, ya1, za1], [[xa1, ya1, za1], ...],
    #  [[xa2, ya2, za2], [[xa2, ya2, za2], ...], 
    #  .
    #  .
    #  .
    #  [[xan, yan, zan], [[xan, yan, zan], ...]]
    mat_b = xyzb.reshape(1, sizeb, 3)
    mat_b = mat_b.repeat(sizea, axis=0)
    # mat_b:
    # [[[xb1, yb1, zb1], [xb2, yb2, zb2], ...],
    #  [[xb1, yb1, zb1], [xb2, yb2, zb2], ...],
    #  .
    #  .
    #  .
    #  [[xb1, yb1, zb1], [xb2, yb2, zb2], ...]]
    dist = mat_a - mat_b
    dist = numpy.sum(dist * dist, axis=2)
    return dist

if __name__ == '__main__':
    main()
