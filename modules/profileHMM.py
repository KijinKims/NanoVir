try:
    import numpy as np
    import networkx as nx
    from hmm_profile import reader
    from hmm_profile.models import HMM
except ImportError:
    print('[Error] Seems you do not have the required python packages. Please check it.')

# python modules
from math import log
from typing import List, Dict, Tuple, Any
from collections.abc import Mapping
from typing import NewType

# NanoVir modules
from codon_table import codon_read_only
from correct import DAG, Idx, NodeId, Base

ProteinCode = NewType('ProteinCode', str)

class PHMM:
    """The class encapsulating HMM_profile HMM object.

    The overall structure of this class is adapted from https://github.com/janmax/Profile-HMM.

    Attributes:
        _phmm: A HMM_profile HMM object.
        _alphabet: A list of alphabets used in HMM. In this module, protein code.
        _alphabet_to_index: A list of index coverted from each protein code.
        _len: An integer count of the residues in HMM.
        _transmissions: A numpy array containing the transmission probabilities from each residue to residue.
        _emissions_from_M: A numpy array containing the match state emission probabilities at each residue.
        _emissions_from_I: A numpy array containing the insertion state emission probabilities at each residue.
    """

    def __init__(self, phmm_ : HMM):
        """Inits PHMM with a given HMM_profile HMM object."""
        # Getting the length
        self._phmm : HMM = phmm_
        self._alphabet : List[ProteinCode] = [ProteinCode(x) for x in self._phmm.metadata.alphabet]
        self._alphabet_to_index : Dict[ProteinCode, Idx] = {value: Idx(index) for index, value in enumerate(self._alphabet)}
        self._len : int = self._phmm.metadata.length

        # Transferting the probabilities
        self._transmissions : np.ndarray = self.transfer_transmissions()
        self._emissions_from_M : np.ndarray
        self._emissions_from_I : np.ndarray
        self._emissions_from_M, self._emissions_from_I = self.transfer_emissons()

    def __len__(self) -> int:
        """Return the number of residues in the HMM."""
        return self._len

    def modified_viterbi(self, dag_ : DAG) -> List[NodeId]:
        """Return the path corrected with viterbi algorithm.
           Generate data objects to store predecessors, ancestors and base of each node where every node id is converted into index in ordering for numpy operation.    
        """
        predecessors : List[List[Idx]] = [[]] * len(dag_)
        for i in range(len(dag_)):
            predecessors[i] = [ dag_.node_id_to_index(predecessor_id) for predecessor_id in dag_.predecessors(dag_.index_to_node_id(i)) ]

        ancestors : List[List[List[Idx]]] = [[]]* len(dag_)
        for i in range(len(dag_)):
            for ancestor_list in dag_.ancestors(dag_.index_to_node_id(i)):
                index_converted_ancestor_list : List[Idx] = [ dag_.node_id_to_index(ancestor_id) for ancestor_id in ancestor_list ]

                if not ancestors[i]:
                    ancestors[i] = [ index_converted_ancestor_list ]
                else:
                    ancestors[i].append(index_converted_ancestor_list)

        bases : List[Base]= [''] * len(dag_)
        for i in range(len(dag_)):
            bases[i] = Base(dag_.base(dag_.index_to_node_id(i)))

        tr, max_tr_idx = self._modified_viterbi(
            predecessors,
            ancestors,
            bases,
            self._emissions_from_M,
            self._emissions_from_I,
            self._transmissions,
            codon_read_only,
            self._alphabet_to_index,
            len(dag_),
            len(self)
        )

        corrected_path : List[NodeId] = [ dag_.index_to_node_id(x) for x in self.traceback(tr, max_tr_idx) ]

        return corrected_path

    def _modified_viterbi(self, predecessors_ : List[List[Idx]], ancestors_ : List[List[List[Idx]]], bases_ : List[Base], e_M_ : np.ndarray, e_I_ : np.ndarray, a_ : np.ndarray, codon_dict_ : Mapping, alphabet_to_index_ : Dict[ProteinCode, Idx], N_ : int, L_ : int):
        """Inner function for Viterbi algorithm.
           TO DO: optimize the performance with numba compile strategy.
        """
        V_M : np.ndarray = np.array(np.ones((N_, L_)) * -np.inf)
        V_I : np.ndarray = np.array(np.ones((N_, L_)) * -np.inf)
        V_D : np.ndarray = np.array(np.ones((N_, L_)) * -np.inf)
        V_N : np.ndarray = np.ones(N_) * -np.inf
        V_C : np.ndarray = np.ones(N_) * -np.inf

        V_M_tr : np.ndarray = np.array(np.zeros((N_, L_)), dtype=[('alignment_type', 'i'), ('dag_node', 'i'), ('hmm_residue', 'i'), ('parent', 'i'), ('grandparent', 'i')])   # first: type of alignment 0 - M, 1 - I, 2 - D, 3 - N, 4 - C
        V_I_tr : np.ndarray = np.array(np.zeros((N_, L_)), dtype=[('alignment_type', 'i'), ('dag_node', 'i'), ('hmm_residue', 'i'), ('parent', 'i'), ('grandparent', 'i')])   # second: DAG node index
        V_D_tr : np.ndarray = np.array(np.zeros((N_, L_)), dtype=[('alignment_type', 'i'), ('dag_node', 'i'), ('hmm_residue', 'i'), ('parent', 'i'), ('grandparent', 'i')])       # third: hmm residue index
        V_N_tr : np.ndarray = np.zeros(N_, dtype=[('alignment_type', 'i'), ('dag_node', 'i'), ('hmm_residue', 'i')])                         # fourth: parent node ordering index
        V_C_tr : np.ndarray = np.zeros(N_, dtype=[('alignment_type', 'i'), ('dag_node', 'i'), ('hmm_residue', 'i')])                         # fifth: grandparent node ordering index
        tr : List[np.ndarray] = [V_M_tr, V_I_tr, V_D_tr, V_N_tr, V_C_tr]

        V_N[0] = 0

        for i in range(1, N_): # Node index in topological order
            for p in predecessors_[i]: # x_i^(1)
                assert p < i

                n_to_n = V_N[p]    # N->N
                if n_to_n > V_N[i]:
                    V_N[i] = n_to_n
                    tr[3][i] = (3,p,0)

                max_idx = np.argmax(V_M[p, :])
                m_to_c = V_M[p, max_idx]   # M->C
                if m_to_c > V_C[i]:
                    V_C[i] = m_to_c
                    tr[4][i] = (0,p,max_idx)

                c_to_c = V_C[p]   # C->C
                if c_to_c > V_C[i]:
                    V_C[i] = c_to_c
                    tr[4][i] = (4,p,0)

            if i == N_-1 or bases_[i] =='^' or bases_[i] == '$':
                continue

            if ancestors_[i]:

                for ancestor_list in ancestors_[i]:
                    gg = ancestor_list[2] # x_i^(3) grandgrandparent
                    g = ancestor_list[1] # x_i^(2) grandparent
                    p = ancestor_list[0] # x_i^(1) parent
                    assert gg < i
                    assert g < i
                    assert p < i

                    codon = bases_[g] + bases_[p] + bases_[i]
                    T = codon_dict_[codon]

                    # skip if stop codon
                    if T == '*':
                        continue

                    x = alphabet_to_index_[T]

                    for j in range(L_): # HMM residue index

                        if j != 0: # skip first residue
                            m_to_m = log(e_M_[x][j+1]) - log(e_M_[x][0]) + V_M[gg, j-1] + log(a_[0][j]) # M->M
                            if m_to_m > V_M[i, j]: 
                                V_M[i, j] = m_to_m
                                tr[0][i,j] = (0, gg, j-1, p, g)
                            
                            i_to_m = log(e_M_[x][j+1]) - log(e_M_[x][0]) + V_I[gg, j-1] + log(a_[3][j]) # I->M
                            if i_to_m > V_M[i, j]: 
                                V_M[i, j] = i_to_m
                                tr[0][i,j] = (1, gg, j-1, p, g)
                            
                            d_to_m = log(e_M_[x][j+1]) - log(e_M_[x][0]) + V_D[gg, j-1] + log(a_[5][j]) # D->M
                            if d_to_m > V_M[i, j]: 
                                V_M[i, j] = d_to_m
                                tr[0][i,j] = (2, gg, j-1, p, g)

                        n_to_m = log(e_M_[x][j+1]) - log(e_M_[x][0]) + V_N[gg] # N->M
                        if n_to_m > V_M[i, j]: 
                            V_M[i, j] = n_to_m
                            tr[0][i,j] = (3, gg, 0, p, g)

                        m_to_i = log(e_I_[x][j+1]) - log(e_I_[x][0]) + V_M[gg, j] + log(a_[1][j+1]) # M->I
                        if m_to_i > V_I[i, j]:
                            V_I[i, j] = m_to_i
                            tr[1][i,j] = (0, gg, j, p, g)
                        
                        i_to_i = log(e_I_[x][j+1]) - log(e_I_[x][0]) + V_I[gg, j] + log(a_[4][j+1]) # I->I
                        if i_to_i > V_I[i, j]:
                            V_I[i, j] = i_to_i
                            tr[1][i,j] = (1, gg, j, p, g)

                        if j != 0 and j != L_-1: # skip first and last residues
                            m_to_d = V_M[i, j-1] + log(a_[2][j+1]) # M->D
                            if m_to_d > V_D[i, j]:
                                V_D[i, j] = m_to_d
                                tr[2][i,j] = (0, i, j-1)
                            
                            d_to_d = V_D[i, j-1] + log(a_[6][j+1]) # D->D
                            if d_to_d > V_D[i, j]:
                                V_D[i, j] = d_to_d
                                tr[2][i,j] = (2, i, j-1)

        max_tr_idx = (4,N_-1,0)

        return tr, max_tr_idx

    def traceback(self, tr_, tr_start_idx_):
        """Trace back the traceback matrix so that we could identify the path with best score."""
        
        traceback_index_list = []

        t, i, j = tr_start_idx_

        while i != 0:
            traceback_index_list.append(i)

            if t == 0 or t == 1:
                traceback_index_list.append(tr_[t][i,j]['parent'])
                traceback_index_list.append(tr_[t][i,j]['grandparent'])
                t, i, j = (tr_[t][i,j]['alignment_type'], tr_[t][i,j]['dag_node'], tr_[t][i,j]['hmm_residue'])
            
            elif t == 2:
                t, i, j = (tr_[t][i,j]['alignment_type'], tr_[t][i,j]['dag_node'], tr_[t][i,j]['hmm_residue'])

            else: # t == 3 or t == 4:
                t, i, j = (tr_[t][i]['alignment_type'], tr_[t][i]['dag_node'], tr_[t][i]['hmm_residue'])

        # begin node
        traceback_index_list.append(i)
        
        traceback_index_list.reverse()
        return traceback_index_list

    def transfer_emissons(self) -> np.ndarray:
        """Transfer the emission probabilites into numpy arrays from HMM_profile HMM object."""
        emissions_from_M : Dict[ProteinCode, np.ndarray] = {ProteinCode(char): np.zeros(len(self)+1) for char in self._alphabet}
        emissions_from_I : Dict[ProteinCode, np.ndarray] = {ProteinCode(char): np.zeros(len(self)+1) for char in self._alphabet}

        for i, alphabet in enumerate(self._alphabet):
            emissions_from_M[alphabet][0] = self._phmm.start_step.p_emission_char[i]
            emissions_from_I[alphabet][0] = self._phmm.start_step.p_insertion_char[i]

        for i, alphabet in enumerate(self._alphabet):
            for j in range(1, len(self)+1):
                emissions_from_M[alphabet][j] = self._phmm.steps[j-1].p_emission_char[i]
                emissions_from_I[alphabet][j] = self._phmm.steps[j-1].p_insertion_char[i]

        # return 2D arrays for performance
        return \
            np.vstack([emissions_from_M[c] for c in self._alphabet]), \
            np.vstack([emissions_from_I[c] for c in self._alphabet])
    
    def transfer_transmissions(self) -> np.ndarray:
        """Transfer the transmission probabilites into numpy arrays from HMM_profile HMM object."""
        # these are all the transmissions we want to observe
        transmission_list = [
            'm->m', 'm->i', 'm->d', 'i->m', 'i->i', 'd->m', 'd->d'
            ]
        transmissions : Dict[str, np.ndarray]= {t: np.zeros(len(self)+1) for t in transmission_list}

        for i in range(1, len(self)+1):
            transmissions['m->m'][i] = self._phmm.steps[i-1].p_emission_to_emission
            transmissions['m->i'][i] = self._phmm.steps[i-1].p_emission_to_insertion
            transmissions['m->d'][i] = self._phmm.steps[i-1].p_emission_to_deletion
            transmissions['i->m'][i] = self._phmm.steps[i-1].p_insertion_to_emission
            transmissions['i->i'][i] = self._phmm.steps[i-1].p_insertion_to_insertion
            transmissions['d->m'][i] = self._phmm.steps[i-1].p_deletion_to_emission
            transmissions['d->d'][i] = self._phmm.steps[i-1].p_deletion_to_deletion

        # return everything as a 2D array for performance
        return np.vstack([transmissions[t] for t in transmission_list])

#################################### just for testing
def locate_grand_parents(DAG_):
    ordering = nx.algorithms.dag.topological_sort(DAG_)

    for node_id in ordering:
        for parent in DAG_.predecessors(node_id):
            for grandparent in DAG_.predecessors(parent):
                for grandgrandparent in DAG_.predecessors(grandparent):
                    ancestors = [parent, grandparent, grandgrandparent]

                    if 'ancestors' not in DAG_.nodes[node_id]:
                        DAG_.nodes[node_id]['ancestors'] = [ancestors]
                    else:
                        DAG_.nodes[node_id]['ancestors'].append(ancestors)
#############################################

if __name__ == "__main__":
    G = nx.DiGraph()

    G.add_node("^_0", base="^")
    G.add_node("A_1", base="A")
    G.add_node("T_2", base="T")
    G.add_node("G_3", base="G")
    G.add_node("C_8", base="C")
    G.add_node("G_4", base="G")
    G.add_node("T_5", base="T")
    G.add_node("T_6", base="T")
    G.add_node("$_7", base="$")
    G.add_edge("^_0", "A_1")
    G.add_edge("A_1", "T_2")
    G.add_edge("T_2", "G_3")
    G.add_edge("G_3", "G_4")
    G.add_edge("G_3", "C_8")
    G.add_edge("C_8", "T_5")
    G.add_edge("G_4", "T_5")
    G.add_edge("T_5", "T_6")
    G.add_edge("T_6", "$_7")
    locate_grand_parents(G)

    f = open("test.hmm")
    model_generator = reader.read_all(f)
    phmms = list(model_generator)
    f.close()
    phmm = phmms[0]

    hmm = PHMM(phmm)
    print(hmm.modified_viterbi(G))