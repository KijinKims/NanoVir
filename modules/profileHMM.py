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
from correct import DAG

NodeId = NewType('NodeId', str)
ProteinCode = NewType('ProteinCode', str)
Idx = NewType('Idx', int)
Base = NewType('Base', str)

class PHMM:
    """docstring for PHMM"""
    def __init__(self, phmm_ : HMM):

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
        return self._len

    def modified_viterbi(self, DAG_ : DAG) -> Tuple[DAG, List[NodeId]]:

        ordering = DAG_.ordering

        node_id_to_index : Dict[NodeId, Idx] = { NodeId(val) : Idx(idx) for idx, val in enumerate(ordering) }

        predecessors : List[List[Idx]] = [[]] * len(DAG_)
        for i in range(len(DAG_)):
            index_to_node_id = ordering[i]
            predecessors[i] = [ node_id_to_index[predecessor_id] for predecessor_id in DAG_.predecessors(index_to_node_id) ]

        ancestors : List[List[List[Idx]]] = [[]]* len(DAG_)
        for i in range(len(DAG_)):
            index_to_node_id = ordering[i]
            for ancestor_list in DAG_.ancestors(index_to_node_id):
                index_converted_ancestor_list : List[Idx] = [ node_id_to_index[ancestor_id] for ancestor_id in ancestor_list ]

                if not ancestors[i]:
                    ancestors[i] = [ index_converted_ancestor_list ]
                else:
                    ancestors[i].append(index_converted_ancestor_list)

        bases : List[Base]= [''] * len(DAG_)
        for i in range(len(DAG_)):
            index_to_node_id = ordering[i]
            bases[i] = Base(DAG_.base(ordering[i]))

        tr, max_tr_idx = self._modified_viterbi(
            predecessors,
            ancestors,
            bases,
            self._emissions_from_M,
            self._emissions_from_I,
            self._transmissions,
            codon_read_only,
            self._alphabet_to_index,
            len(DAG_),
            len(self)
        )

        corrected_path = [ ordering[x] for x in self.traceback(tr, max_tr_idx) ] # list of node id

        return DAG_, corrected_path

    def _modified_viterbi(self, predecessors : List[List[Idx]], ancestors : List[List[List[Idx]]], bases : List[Base], e_M : np.ndarray, e_I : np.ndarray, a : np.ndarray, codon_dict : Mapping, alphabet_to_index : Dict[ProteinCode, Idx], N : int, L : int):

        V_M : np.ndarray = np.array(np.ones((N, L)) * -np.inf)
        V_I : np.ndarray = np.array(np.ones((N, L)) * -np.inf)
        V_D : np.ndarray = np.array(np.ones((N, L)) * -np.inf)
        V_N : np.ndarray = np.ones(N) * -np.inf
        V_C : np.ndarray = np.ones(N) * -np.inf

        V_M_tr : np.ndarray = np.array(np.zeros((N, L)), dtype='i,i,i,i,i')   # first: type of alignment 0 - M, 1 - I, 2 - D, 3 - N, 4 - C
        V_I_tr : np.ndarray = np.array(np.zeros((N, L)), dtype='i,i,i,i,i')   # second: DAG node index
        V_D_tr : np.ndarray = np.array(np.zeros((N, L)), dtype='i,i,i')       # third: hmm residue index
        V_N_tr : np.ndarray = np.zeros(N, dtype='i,i,i')                         # fourth: parent node ordering index
        V_C_tr : np.ndarray = np.zeros(N, dtype='i,i,i')                         # fifth: grandparent node ordering index
        tr : List[np.ndarray] = [V_M_tr, V_I_tr, V_D_tr, V_N_tr, V_C_tr]

        V_N[0] = 0

        for i in range(1, N): # Node index in topological order
            for p in predecessors[i]: # x_i^(1)
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

            if i == N-1 or bases[i] =='^' or bases[i] == '$':
                continue

            if ancestors[i]:

                for ancestor_list in ancestors[i]:
                    gg = ancestor_list[2] # x_i^(3) grandgrandparent
                    g = ancestor_list[1] # x_i^(2) grandparent
                    p = ancestor_list[0] # x_i^(1) parent
                    assert gg < i
                    assert g < i
                    assert p < i

                    codon = bases[g] + bases[p] + bases[i]
                    T = codon_dict[codon]

                    # skip if stop codon
                    if T == '*':
                        continue

                    x = alphabet_to_index[T]

                    for j in range(L): # HMM residue index

                        if j != 0: # skip first residue
                            m_to_m = log(e_M[x][j+1]) - log(e_M[x][0]) + V_M[gg, j-1] + log(a[0][j]) # M->M
                            if m_to_m > V_M[i, j]: 
                                V_M[i, j] = m_to_m
                                tr[0][i,j] = (0, gg, j-1, p, g)
                            
                            i_to_m = log(e_M[x][j+1]) - log(e_M[x][0]) + V_I[gg, j-1] + log(a[3][j]) # I->M
                            if i_to_m > V_M[i, j]: 
                                V_M[i, j] = i_to_m
                                tr[0][i,j] = (1, gg, j-1, p, g)
                            
                            d_to_m = log(e_M[x][j+1]) - log(e_M[x][0]) + V_D[gg, j-1] + log(a[5][j]) # D->M
                            if d_to_m > V_M[i, j]: 
                                V_M[i, j] = d_to_m
                                tr[0][i,j] = (2, gg, j-1, p, g)

                        n_to_m = log(e_M[x][j+1]) - log(e_M[x][0]) + V_N[gg] # N->M
                        if n_to_m > V_M[i, j]: 
                            V_M[i, j] = n_to_m
                            tr[0][i,j] = (3, gg, 0, p, g)

                        m_to_i = log(e_I[x][j+1]) - log(e_I[x][0]) + V_M[gg, j] + log(a[1][j+1]) # M->I
                        if m_to_i > V_I[i, j]:
                            V_I[i, j] = m_to_i
                            tr[1][i,j] = (0, gg, j, p, g)
                        
                        i_to_i = log(e_I[x][j+1]) - log(e_I[x][0]) + V_I[gg, j] + log(a[4][j+1]) # I->I
                        if i_to_i > V_I[i, j]:
                            V_I[i, j] = i_to_i
                            tr[1][i,j] = (1, gg, j, p, g)

                        if j != 0 and j != L-1: # skip first and last residues
                            m_to_d = V_M[i, j-1] + log(a[2][j+1]) # M->D
                            if m_to_d > V_D[i, j]:
                                V_D[i, j] = m_to_d
                                tr[2][i,j] = (0, i, j-1)
                            
                            d_to_d = V_D[i, j-1] + log(a[6][j+1]) # D->D
                            if d_to_d > V_D[i, j]:
                                V_D[i, j] = d_to_d
                                tr[2][i,j] = (2, i, j-1)

        print("final:")
        print(V_M)
        print(V_I)
        print(V_D)
        print(V_N)
        print(V_C)

        max_tr_idx = (4,N-1,0)

        return tr, max_tr_idx

    def traceback(self, tr, tr_start_idx):
        
        node_index_list = []

        t, i, j = tr_start_idx

        while i != 0:
            node_index_list.append(i)

            if t == 0 or t == 1:
                node_index_list.append(tr[t][i,j][3])
                node_index_list.append(tr[t][i,j][4])
                t, i, j = (tr[t][i,j][0], tr[t][i,j][1], tr[t][i,j][2])
            
            elif t == 2:
                t, i, j = (tr[t][i,j][0], tr[t][i,j][1], tr[t][i,j][2])

            else: # t == 3 or t == 4:
                t, i, j = (tr[t][i][0], tr[t][i][1], tr[t][i][2])

        # begin node
        node_index_list.append(i)
        
        node_index_list.reverse()
        return node_index_list

    def transfer_emissons(self) -> np.ndarray:
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