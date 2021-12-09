# external Libraries
try:
    import numpy as np
    import networkx as nx
except ImportError:
    print('[Error] Seems you do not have numpy, numba, networkx. Too bad - Exiting...')

# python modules
from math import log
from hmm_profile import reader

# Standard genetic code:
codon_dict = {}
# nTn
codon_dict["TTT"] = codon_dict["TTC"] = "F"
codon_dict["TTA"] = codon_dict["TTG"] = codon_dict["CTT"] = codon_dict["CTC"] = codon_dict["CTA"] = codon_dict["CTG"] = "L"
codon_dict["ATT"] = codon_dict["ATC"] = codon_dict["ATA"] = "I"
codon_dict["ATG"] = "M"
codon_dict["GTT"] = codon_dict["GTC"] = codon_dict["GTA"] = codon_dict["GTG"] = "V"

# nCn
codon_dict["TCT"] = codon_dict["TCC"] = codon_dict["TCA"] = codon_dict["TCG"] = "S"
codon_dict["CCT"] = codon_dict["CCC"] = codon_dict["CCA"] = codon_dict["CCG"] = "P"
codon_dict["ACT"] = codon_dict["ACC"] = codon_dict["ACA"] = codon_dict["ACG"] = "T"
codon_dict["GCT"] = codon_dict["GCC"] = codon_dict["GCA"] = codon_dict["GCG"] = "A"

# nAn
codon_dict["TAT"] = codon_dict["TAC"] = "Y"
codon_dict["TAA"] = codon_dict["TAG"] = "*" 	#Stop
codon_dict["CAT"] = codon_dict["CAC"] = "H"
codon_dict["CAA"] = codon_dict["CAG"] = "Q"
codon_dict["AAT"] = codon_dict["AAC"] = "N"
codon_dict["AAA"] = codon_dict["AAG"] = "K"
codon_dict["GAT"] = codon_dict["GAC"] = "D"
codon_dict["GAA"] = codon_dict["GAG"] = "E"

# nGn
codon_dict["TGT"] = codon_dict["TGC"] = "C"
codon_dict["TGA"] = "*"  			#Stop
codon_dict["TGG"] = "W"
codon_dict["CGT"] = codon_dict["CGC"] = codon_dict["CGA"] = codon_dict["CGG"] = "R"
codon_dict["AGT"] = codon_dict["AGC"] = "S"
codon_dict["AGA"] = codon_dict["AGG"] = "R"
codon_dict["GGT"] = codon_dict["GGC"] = codon_dict["GGA"] = codon_dict["GGG"] = "G"

class HMM:
    """docstring for HMM"""
    def __init__(self, phmm_):

        # Getting the length
        self.phmm = phmm_
        self.alphabet = self.phmm.metadata.alphabet
        self.alphabet_dict = {value: index for index, value in enumerate(self.alphabet)}
        self.n = self.phmm.metadata.length

        # Transferting the probabilities
        self.transmissions = self.transfer_transmissions()
        self.emissions_from_M, self.emissions_from_I = self.transfer_emissons()

    def modified_viterbi(self, DAG_):

        ordering = list(nx.algorithms.dag.topological_sort(DAG_))
        node_id_to_sorted_order_dict = { val : idx for idx, val in enumerate(ordering) }

        predecessors = [[]] * len(ordering) # list of list
        for i in range(len(ordering)):
            predecessors[i] = [ node_id_to_sorted_order_dict[predecessor_id] for predecessor_id in DAG_.predecessors(ordering[i]) ]

        ancestors = [[]] * len(ordering) # list of list of list
        for i in range(len(ordering)):
            if 'ancestors' in DAG_.nodes[ordering[i]]:
                for ancestor_set in DAG_.nodes[ordering[i]]['ancestors']:
                    sorted_order_converted = [ node_id_to_sorted_order_dict[ancestor_id] for ancestor_id in ancestor_set ]
                    
                    if not ancestors[i]:
                        ancestors[i] = [ sorted_order_converted ]
                    else:
                        ancestors[i].append(sorted_order_converted)

        bases = [''] * len(ordering)
        for i in range(len(ordering)):
            bases[i] = DAG_.nodes[ordering[i]]['base']

        tr, max_tr_idx = self._modified_viterbi(
            predecessors,
            ancestors,
            bases,
            self.emissions_from_M,
            self.emissions_from_I,
            self.transmissions,
            codon_dict,
            self.alphabet_dict,
            len(ordering),
            self.n
        )

        corrected_path = [ ordering[x] for x in self.traceback(tr, max_tr_idx) ] # list of node id
        attrs = {}
        for node_id in corrected_path:
            attrs[node_id] = { "corrected_consensus": 'true' }
        for node_id in DAG_.nodes():
            if node_id not in attrs:
                attrs[node_id] = { "corrected_consensus": 'false' }
        nx.set_node_attributes(DAG_, attrs)
        return DAG_, corrected_path

    def _modified_viterbi(self, predecessors, ancestors, bases, e_M, e_I, a, codon_dict, alphabet_dict, N, L):

        V_M = np.matrix(np.ones((N, L)) * -np.inf)
        V_I = np.matrix(np.ones((N, L)) * -np.inf)
        V_D = np.matrix(np.ones((N, L)) * -np.inf)
        V_N = np.ones(N) * -np.inf
        V_C = np.ones(N) * -np.inf

        V_M_tr = np.matrix(np.zeros((N, L)), dtype='i,i,i,i,i')   # first: type of alignment 0 - M, 1 - I, 2 - D, 3 - N, 4 - C
        V_I_tr = np.matrix(np.zeros((N, L)), dtype='i,i,i,i,i')   # second: DAG node index
        V_D_tr = np.matrix(np.zeros((N, L)), dtype='i,i,i')       # third: hmm residue index
        V_N_tr = np.zeros(N, dtype='i,i,i')                         # fourth: parent node ordering index
        V_C_tr = np.zeros(N, dtype='i,i,i')                         # fifth: grandparent node ordering index
        tr = [V_M_tr, V_I_tr, V_D_tr, V_N_tr, V_C_tr]

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

            if i == N-1:
                break

            if ancestors[i]:

                for ancestor_set in ancestors[i]:
                    gg = ancestor_set[2] # x_i^(3) grandgrandparent
                    g = ancestor_set[1] # x_i^(2) grandparent
                    p = ancestor_set[0] # x_i^(1) parent

                    assert gg < i
                    assert g < i
                    assert p < i

                    codon = bases[g] + bases[p] + bases[i]
                    T = codon_dict[codon]

                    # skip if stop codon
                    if T == '*':
                        continue

                    x = alphabet_dict[T]

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

    def transfer_emissons(self):
        emissions_from_M = {char: np.zeros(self.n+1) for char in self.alphabet}
        emissions_from_I = {char: np.zeros(self.n+1) for char in self.alphabet}

        for alphabet in self.alphabet:
            emissions_from_M[alphabet][0] = self.phmm.start_step.p_emission_char[alphabet]
            emissions_from_I[alphabet][0] = self.phmm.start_step.p_insertion_char[alphabet]

        for alphabet in self.alphabet:
            for i in range(1, self.n+1):
                emissions_from_M[alphabet][i] = self.phmm.steps[i-1].p_emission_char[alphabet]
                emissions_from_I[alphabet][i] = self.phmm.steps[i-1].p_insertion_char[alphabet]

        # return 2D arrays for performance
        return \
            np.vstack([emissions_from_M[c] for c in self.alphabet]), \
            np.vstack([emissions_from_I[c] for c in self.alphabet])
    
    def transfer_transmissions(self):
        # these are all the transmissions we want to observe
        transmission_list = [
            'm->m', 'm->i', 'm->d', 'i->m', 'i->i', 'd->m', 'd->d'
            ]
        transmissions = {t: np.zeros(self.n+1) for t in transmission_list}

        for i in range(1, self.n+1):
            transmissions['m->m'][i] = self.phmm.steps[i-1].p_emission_to_emission
            transmissions['m->i'][i] = self.phmm.steps[i-1].p_emission_to_insertion
            transmissions['m->d'][i] = self.phmm.steps[i-1].p_emission_to_deletion
            transmissions['i->m'][i] = self.phmm.steps[i-1].p_insertion_to_emission
            transmissions['i->i'][i] = self.phmm.steps[i-1].p_insertion_to_insertion
            transmissions['d->m'][i] = self.phmm.steps[i-1].p_deletion_to_emission
            transmissions['d->d'][i] = self.phmm.steps[i-1].p_deletion_to_deletion

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

    hmm = HMM(phmm)
    print(hmm.modified_viterbi(G))