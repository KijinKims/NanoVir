from typing import DefaultDict, List, Dict, Tuple
from typing import NewType
from collections import defaultdict
from pathlib import Path
import copy
from Bio import SeqIO
from Bio import SearchIO
from Bio.SearchIO._model import HSP
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from hmm_profile import reader
from hmm_profile.models import HMM
import pydot
import networkx as nx
import profileHMM

GraphId = NewType('GraphId', str)
NodeId = NewType('NodeId', str)
EdgeId = NewType('EdgeId', str)
TigId = NewType('TigId', str)
HMMId = NewType('HMMId', str)
Base = NewType('Base', str)

class DAG:
    def __init__(self, digraph : nx.DiGraph, edge_threshold = 1):
        self._name : GraphId = GraphId(digraph.graph['name'])
        self._graph : nx.DiGraph = nx.DiGraph(digraph)
        self._ordering : List[NodeId] = [NodeId(node_id) for node_id in nx.algorithms.dag.topological_sort(self._graph)]
        self._len : int = len(self._graph)
        self._source : NodeId = NodeId(self._ordering[0])
        self._sink : NodeId = NodeId(self._ordering[len(self)-1])
        self._edge_threshold : int = edge_threshold

    def __len__(self) -> int :
        return self._len

    def len(self, new_len : int) :
        assert new_len > 0
        self._len = new_len

    @property
    def name(self) -> GraphId :
        return self._name

    @property
    def ordering(self) -> List[NodeId] :
        return self._ordering

    @ordering.setter
    def ordering(self, new_ordering: List[NodeId]):
        self._ordering = new_ordering

    @property
    def source(self) -> NodeId :
        return self._source

    @source.setter
    def source(self, new_source : NodeId) :
        assert self._graph.has_node(new_source)
        self._source = new_source

    @property
    def sink(self) -> NodeId :
        return self._sink

    @sink.setter
    def sink(self, new_sink : NodeId) :
        assert self._graph.has_node(new_sink)
        self._sink = new_sink

    def nodes(self, *args, **argv) :
        return self._graph.nodes(*args, **argv)

    def edges(self, *args, **argv) :
        return self._graph.edges(*args, **argv)

    def in_degree(self, node_id : NodeId):
        return self._graph.in_degree(node_id)

    def out_degree(self, node_id : NodeId):
        return self._graph.out_degree(node_id)

    def prun(self, edge_threshold : int) :

        self._edge_threshold = edge_threshold
        self.filter_by_edge_weight(self._edge_threshold)
        self.clean_obsolete_nodes()
        self.len(len(self._graph))
        self.leave_only_one_source()
        self.leave_only_one_sink()

        return self

    def filter_by_edge_weight(self, threshold: int) :
        """ Remove all edges with weight<threshold from DAG. """
        edges = [(n1, n2) for n1, n2, w in self.edges(data="weight") if int(w) < threshold]
        self._graph.remove_edges_from(edges)

    def clean_obsolete_nodes(self):
        """ Remove disconnected nodes from DAG. """
        obsolete_nodes = []
        for node_id in self.nodes():
            if node_id == self.source or node_id == self.sink:
                continue
            if self.in_degree(node_id) == 0 and self.out_degree(node_id) == 0:
                obsolete_nodes.append(node_id)

        self._graph.remove_nodes_from(obsolete_nodes)
    
    def leave_only_one_source(self):
        """Connect all sources to original source"""
        sources = []
        for node_id in self.nodes():
            if node_id == self.source or node_id == self.sink:
                continue
            if self.in_degree(node_id) == 0:
                sources.append(node_id)
        
        self._graph.add_edges_from([(self.source, x) for x in sources])

    def leave_only_one_sink(self):
        """Connect all sinks to original sink"""
        sinks = []
        for node_id in self.nodes():
            if node_id == self.source or node_id == self.sink:
                continue
            if self.out_degree(node_id) == 0:
                sinks.append(node_id)
        
        self._graph.add_edges_from([(x, self.sink) for x in sinks])

    def set_node_attributes(self, *args, **argv) :
        nx.set_node_attributes(self._graph, *args, **argv)

    def get_node_attributes(self, *args, **argv) :
        return nx.get_node_attributes(self._graph, *args, **argv)

    def set_edge_attributes(self, *args, **argv) :
        nx.set_edge_attributes(self._graph, *args, **argv)

    def get_edge_attributes(self, *args, **argv) :
        return nx.get_edge_attributes(self._graph, *args, **argv)

    def predecessors(self, node_id : NodeId) -> List[NodeId] :
        return self._graph.predecessors(node_id)

    def set_ancestors(self):
        """find parent, grandparent, grandgrandparent of each node
           and set them as member variable"""
        attrs = defaultdict(lambda:{'ancestors': []})

        for node_id in self.ordering:
            attrs[node_id]['ancestors'] = []
            for parent in self.predecessors(node_id):
                for grandparent in self.predecessors(parent):
                    for grandgrandparent in self.predecessors(grandparent):
                        ancestors = [parent, grandparent, grandgrandparent]
                        attrs[node_id]['ancestors'].append(ancestors)
        
        self.set_node_attributes(attrs)

    def ancestors(self, node_id : NodeId) -> List[NodeId] :
        return self.get_node_attributes('ancestors')[node_id]

    def bases(self) -> Dict[NodeId, Base] :
        """return dictionary of base of node with node id as a key"""
        return self.nodes(data="base")

    def base(self, node_id : NodeId) -> Base :
        """return base of the node with given node id"""
        return self.bases()[node_id]

    def export_dot(self, dot_path : Path) :
        """export corrected path in DAG as graphviz DOT format"""
        # set node label to be displayed in graphviz
        attrs = defaultdict(dict)
        for u, base in self.bases():
            attrs[u]['label'] = base

        # set node color to be displayed in graphviz
        for node_id in self.nodes():
            if self.nodes(data='initial_consensus')[node_id] == 'true' or self.nodes(data='corrected_consensus')[node_id] == 'true':
                attrs[node_id]['style'] = 'filled'

                if self.nodes(data='initial_consensus')[node_id] == 'true' :
                    if self.nodes(data='corrected_consensus')[node_id] == 'true' :
                        attrs[node_id]['color'] = 'purple'
                    else:
                        attrs[node_id]['color'] = 'blue'
                else:
                    attrs[node_id]['color'] = 'red'
        
        self.set_node_attributes(attrs)

        # set edge label to be displayed in graphviz
        attrs = defaultdict(dict)
        for u,v,w in self.edges(data="weight"):
            attrs[(u, v)]['label'] = w

        self.set_edge_attributes(attrs)

        nx.drawing.nx_pydot.write_dot(self._graph, str(dot_path))

    def base_quote_strip(self):
        """strip the double quotes on both ends of base attribute of each node"""
        attrs = defaultdict(dict)
        for node_id, original_base in self.bases():
            stripped_base = original_base.strip('"')
            attrs[node_id]['base'] = stripped_base
            self.set_node_attributes(attrs)

    def reverse(self):
        """reverse the graph"""
        self._graph = self._graph.reverse()

    def complement(self):
        """change each base as complimentary one"""
        complement_dict = { "A": "T",
                            "T": "A",
                            "G": "C",
                            "C": "G",
                            } 
        
        attrs = defaultdict(dict)
        for node_id, original_base in self.bases():
            if original_base in ["^", "$"]:
                continue
            
            attrs[node_id]['base'] = complement_dict[original_base]
            self.set_node_attributes(attrs)

    def swap_source_sink(self):
        """swap source node and sink node"""
        tmp : NodeId = self.sink
        self.sink = copy.copy(self.source)
        self.source = tmp

    def reverse_complement(self):
        self.reverse()
        self.complement()
        self.swap_source_sink()

    def leave_only_indels(self, path_ : List[NodeId]):
        pass

    def trim_non_consensus(self, path_ : List[NodeId]):
        """corrected consensus sequence should be shorter than initial consensus sequence
           so trim the both ends of consensus sequence until any initial consensus node is encountered"""
        initial_consensus_tf = dict(self.nodes(data='initial_consensus', default=False))

        start_index = next(index for index, node_id in enumerate(path_) if initial_consensus_tf[node_id])

        end_index = next(index for index, node_id in reversed(list(enumerate(path_))) if initial_consensus_tf[node_id])

        consensus_path = path_[start_index:end_index+1]

        attrs = defaultdict(dict)
        for node_id in consensus_path:
            attrs[node_id]["corrected_consensus"] = 'true'
        for node_id in self.nodes():
            if node_id not in attrs:
                attrs[node_id]["corrected_consensus"] = 'false'
        dag.set_node_attributes(attrs)

        return consensus_path

    def topological_sort(self) -> List[NodeId] :
        return [NodeId(node_id) for node_id in nx.algorithms.dag.topological_sort(self._graph)]

def read_contigs(contigs_path_) -> List[SeqRecord]:
    """read the contigs from fasta
       store them as list of SeqRecord objects in Biopython"""
    return SeqIO.parse(contigs_path_, "fasta")

def read_phmmDB(phmmDB_path_ : Path) -> Dict[HMMId, HMM]:
    """read the profile HMMs
       store them as dictionary of HMM objects in profileHMM module with each HMM's id as a key"""
    f = phmmDB_path_.open()
    model_generator = reader.read_all(f)
    
    return  {HMMId(x.metadata.model_name): x for x in model_generator}

def read_dot(dot_path_) -> List[DAG] :
    """read the DAGs from dot
       store them as list of DAG objects"""
    graphs = pydot.graph_from_dot_file(dot_path_)
    
    graphs = map(nx.nx_pydot.from_pydot, graphs)
    
    dags = [DAG(graph) for graph in graphs]

    return dags

def parse_hmmscanresult(hmmscan_domtbl_path_ : Path) -> List[HSP]:
    """read the hmmscan domain hits from domtbl
       store them as list of HSPs object in Biopython SearchIO module
       for unique contig and hmm pair, only a hit with the best score is stored for reducing redundancies"""
    domtbl = SearchIO.parse(hmmscan_domtbl_path_, "hmmscan3-domtab")
    
    uniq_contig_HMM_dict = {}
    # for each pair of contig and HMM,
    # only the hsp with the best score is left
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

    def extract_rframe(hsp : HSP):
        """extract reading frame and direction for later use"""
        hsp.tig_id = hsp.query_id.split("_rframe")[0]
        hsp.reading_frame = hsp.query_id.split("_rframe")[1]
        hsp.direction = "-" if hsp.reading_frame[0]  == "-" else "+"
        hsp.reading_frame = abs(int(hsp.query_id.split("_rframe")[1]))
    
    for hsp in uniq_contig_HMM_dict.values():
        extract_rframe(hsp)
    
    return list(uniq_contig_HMM_dict.values())

def pair_contig_DAG(contigs_ : List[SeqRecord], dags_ : List[DAG]) -> DefaultDict[TigId, Tuple[SeqRecord, DAG]]:
    """pair a contig and DAG with a matching contig id"""
    contigs_dict = {x.id: x for x in contigs_}
    dags_dict = {x.name: x for x in dags_}
    pair_dict : DefaultDict[TigId, Tuple[SeqRecord, DAG]] = defaultdict(tuple)

    for key in contigs_dict:
        pair_dict[TigId(key)] = (contigs_dict[key], dags_dict[key])

    return pair_dict

def correct_consensus(dag_ : DAG, phmm_ : HMM) -> Tuple[DAG, List[NodeId]] :
    '''modify DAG to correct indels more agressively and return corrected consensus path'''
    PHMM = profileHMM.PHMM(phmm_)
    return PHMM.modified_viterbi(dag_)

def write_corrected_contig(corrected_dag_ : DAG, corrected_path_ : List[NodeId], contig_id_ : TigId, phmm_id_ : HMMId, f_ : Path):
    """write the corrected sequence as a FASTA format
       with the header including contig ID and HMM ID used for correction"""
    corrected_sequence = ""
    for node_id in corrected_path_:
        corrected_sequence += corrected_dag_.base(node_id)
    
    seqid = contig_id_ + "_corrected_" + phmm_id_ + '\n'
    record = SeqRecord(
        Seq(corrected_sequence),
        id=seqid)

    SeqIO.write(record, f_, "fasta")

if __name__ == '__main__':
    contigs_path : Path = Path("contigs.fasta")
    phmmDB_path : Path = Path("FAM173.hmm")
    dot_path : Path = Path("graph.dot")
    hmmscan_domtbl_path : Path = Path("hmmscan.domtbl")
    output_path : Path = Path("corrected_contigs.fasta")
    output_dot : Path = Path("output.dot")
    minimum_edge_weight : int = 2
    
    phmms = read_phmmDB(phmmDB_path)
    hmmscan_domtbl = parse_hmmscanresult(hmmscan_domtbl_path)

    contigs = read_contigs(contigs_path)
    dags = read_dot(dot_path)
    contig_dag      = pair_contig_DAG(contigs, dags)

    for hsp in hmmscan_domtbl:
        contig_id       = hsp.query_id.split("_rframe")[0]
        matched_phmm_id = hsp.hit_id

        contig          = contig_dag[contig_id][0]
        dag             = contig_dag[contig_id][1]
        matched_phmm    = phmms[matched_phmm_id]

        dag.base_quote_strip()

        if hsp.direction == "-":
            dag.reverse_complement()

        dag.prun(minimum_edge_weight)
        # after prunning (and revserse complement, if happend)
        # new ordering and ancestor calcuation is needed
        dag.ordering = dag.topological_sort()
        dag.set_ancestors()

        corrected_dag, corrected_path = correct_consensus(dag, matched_phmm)
        consensus_path = corrected_dag.trim_non_consensus(corrected_path)
        corrected_dag.export_dot(output_dot)

        write_corrected_contig(corrected_dag, consensus_path, contig_id, matched_phmm_id, output_path)