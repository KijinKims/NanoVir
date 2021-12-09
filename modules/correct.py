from Bio import SeqIO
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from hmm_profile import reader
import pydot
import networkx as nx
import modules.profileHMM as profileHMM

def read_contigs(contigs_path_):
    return SeqIO.parse(contigs_path_, "fasta")

def read_phmmDB(phmmDB_path_):
    f = open(phmmDB_path_)
    model_generator = reader.read_all(f)

    phmms = list(model_generator)
    f.close()
    
    return  {x.metadata.model_name: x for x in phmms}

def read_dot(dot_path_):
    graphs = pydot.graph_from_dot_file(dot_path_)
    
    graphs = map(nx.nx_pydot.from_pydot, graphs)
    
    return graphs

def prun_DAG(DAGS_):
    pass

def parse_hmmscanresult(hmmscan_domtbl_path_):
    domtbl = SearchIO.parse(hmmscan_domtbl_path_, "hmmscan3-domtab")
    
    uniq_contig_HMM_dict = {}
    # for each pair of contig and HMM,
    # only the longest hsp is left
    for query in domtbl:
        tig_id = query.id.split("_rframe")[0]
        for hit in query:
            hmm_id = hit.id
            for hsp in hit:
                if (tig_id, hmm_id) not in uniq_contig_HMM_dict:
                    uniq_contig_HMM_dict[(tig_id, hmm_id)] = hsp
                else:
                    curr_best_hsp = uniq_contig_HMM_dict[(tig_id, hmm_id)]
                    if hsp.bitscore > curr_best_hsp.bitscore:
                        uniq_contig_HMM_dict[(tig_id, hmm_id)] = hsp
    
    # extract reading frame and direction for later use
    for hsp in uniq_contig_HMM_dict.values():
        hsp.tig_id = hsp.query_id.split("_rframe")[0]
        hsp.reading_frame = hsp.query_id.split("_rframe")[1]
        hsp.direction = "-" if hsp.reading_frame[0]  == "-" else "+"
        hsp.reading_frame = abs(int(hsp.query_id.split("_rframe")[1]))
    
    return list(uniq_contig_HMM_dict.values())

def pair_contig_DAG(contigs_, DAGs_):
    
    contigs_dict = {x.id: x for x in contigs_}
    DAGs_dict = {x.graph['name']: x for x in DAGs_}
    pair_dict = {}

    for key in contigs_dict:
        pair_dict[key] = (contigs_dict[key], DAGs_dict[key])

    return pair_dict

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

def locate_node_in_DAG(DAG_, cns_pos_):
    for node in DAG_.nodes(data="cns_pos"):
        if(int(node[1]) == cns_pos_):   # node[1] : position in the context of contig
            return node[0]              # node[0] : node ID
    return None

def base_quote_strip(DAG_):
    attrs = {}
    for node, original_base in nx.get_node_attributes(DAG_, 'base').items():
        stripped_base = original_base.strip('"')

        attrs[node] = {}
        attrs[node]['base'] = stripped_base
        nx.set_node_attributes(DAG_, attrs)

    return DAG_

def complement(DAG_):
    
    complement_dict = { "A": "T",
                        "T": "A",
                        "G": "C",
                        "C": "G",
                        } 
    
    attrs = {}
    for node, original_base in nx.get_node_attributes(DAG_, 'base').items():
        if original_base in ["^", "$"]:
            continue
        
        attrs[node] = {}
        attrs[node]['base'] = complement_dict[original_base]
        nx.set_node_attributes(DAG_, attrs)
    
    return DAG_

def correct_contig(DAG_, phmm_):
    locate_grand_parents(DAG_)
    PHMM = profileHMM.HMM(phmm_)
    return PHMM.modified_viterbi(DAG_)
    

def export_to_DOT(DAG_):

    # set label to be displayed in graphviz
    attrs = {}
    for node in DAG_.nodes(data="base"):
        attrs[node[0]] = {}
        attrs[node[0]]['label'] = node[1]
    
    for node_id in DAG_.nodes():
        if DAG_.nodes[node_id]['initial_consensus'] == 'true' or DAG_.nodes[node_id]['corrected_consensus'] == 'true':
            attrs[node_id]['style'] = 'filled'

            if DAG_.nodes[node_id]['initial_consensus'] == 'true' :
                if DAG_.nodes[node_id]['corrected_consensus'] == 'true' :
                    attrs[node_id]['color'] = 'purple'
                else:
                    attrs[node_id]['color'] = 'blue'
            else:
                attrs[node_id]['color'] = 'red'
    
    nx.set_node_attributes(DAG_, attrs)

    nx.drawing.nx_pydot.write_dot(DAG_, "example_corrected_graph.dot")


def write_corrected_contig(corrected_DAG_, corrected_path_, contig_id_, phmm_id_, file_handle_):
    
    correcte_sequence = ""
    for node_id in corrected_path_:
        node_bases = corrected_DAG_.nodes(data='base')
        correcte_sequence += node_bases[node_id]
    
    seqid = contig_id_ + "_corrected_" + phmm_id_ + '\n'
    record = SeqRecord(
        Seq(correcte_sequence),
        id=seqid)
    SeqIO.write(record, output_handle, "fasta")

    

if __name__ == '__main__':
    contigs_path        = "contigs.fasta"
    phmmDB_path         = "FAM173.hmm"
    dot_path            = "graph.dot"
    hmmscan_domtbl_path = "hmmscan.domtbl"
    output_path         = "corrected_contigs.fasta"
    output_dot          = "output.dot"
    minimum_edge_weight = 1
    
    phmms           = read_phmmDB(phmmDB_path) # dictionary
                                               # key: HMM name
                                               # value: HMM obj
    hmmscan_domtbl  = parse_hmmscanresult(hmmscan_domtbl_path) # list of hsp obj

    contigs         = read_contigs(contigs_path) # list of SeqRecord obj
    DAGs            = read_dot(dot_path) # list of networkx obj
    pruned_DAGs     = prun_DAG(DAGs, minimum_edge_weight) # prun the edges whose weights are smaller than threshold
    contig_DAG      = pair_contig_DAG(contigs, pruned_DAGs) # dictionary 
                                                     # key: contig ID
                                                     # value: (contig's SeqRecord obj, networkx obj of its DAG)
    with open(output_path, "w") as output_handle:

        for hsp in hmmscan_domtbl:
            contig_id       = hsp.query_id.split("_rframe")[0]
            matched_phmm_id = hsp.hit_id

            contig          = contig_DAG[contig_id][0]
            DAG             = contig_DAG[contig_id][1]
            matched_phmm    = phmms[matched_phmm_id]

            DAG = base_quote_strip(DAG)

            if hsp.direction == "-":
                DAG = DAG.reverse()
                DAG = complement(DAG)

            corrected_DAG, corrected_path = correct_contig(DAG, matched_phmm)
            export_to_DOT(corrected_DAG)

            write_corrected_contig(corrected_DAG, corrected_path, contig_id, matched_phmm_id, output_handle)