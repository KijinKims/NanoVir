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
Idx = NewType('Idx', int)

class DAG:
    """The class encapsulating Networkx Digraph object.

    Attributes:
        _name: A string indicating the name of the graph.
        _graph: A Networkx Digraph object.
        _ordering: A list of strings that are node ids of the graph.
        _len: An integer count of the nodes in the graph.
        _source: A string indicating the node id of source node in the graph.
        _sink: A string indicating the node id of sink node in the graph.
        _edge_threshold: An integer indicating the threshold of weight for filtering edges in the graph.
        _node_id_to_index_dict: A dictionary matching every node id to index in self._ordering list. It is updated every time ordering is updated.
    """

    def __init__(self, digraph_ : nx.DiGraph):
        """Inits DAG with a given Networkx Digraph object."""
        self._name : GraphId = GraphId(digraph_.graph['name'])
        self._graph : nx.DiGraph = nx.DiGraph(digraph_)
        self._ordering : List[NodeId] = [NodeId(node_id) for node_id in nx.algorithms.dag.topological_sort(self._graph)]
        self._len : int = len(self._graph)
        self._source : NodeId = NodeId(self._ordering[0])
        self._sink : NodeId = NodeId(self._ordering[len(self)-1])

    def __len__(self) -> int :
        """Return the number of nodes in the graph."""
        return self._len

    def len(self, len_ : int) :
        """Set the number of nodes in the graph."""
        assert len_ > 0
        self._len = len_

    @property
    def name(self) -> GraphId :
        """Return the name of the graph."""
        return self._name

    @property
    def ordering(self) -> List[NodeId] :
        """Return the sorted order of nodes in the graph. Typically, they are topologically sorted."""
        return self._ordering

    @ordering.setter
    def ordering(self, ordering_: List[NodeId]):
        """Set orted order of nodes in the graph. Typically, they are topologically sorted. Afterwards, update self._node_id_to_index_dict."""
        self._ordering = ordering_
        self.set_node_id_to_index()

    @property
    def source(self) -> NodeId :
        """Return the source node in the graph."""
        return self._source

    @source.setter
    def source(self, source_ : NodeId) :
        """Set the source node in the graph. The source node should have no indegree edges."""
        assert self._graph.has_node(source_)
        assert self._graph.in_degree(source_) == 0
        self._source = source_

    @property
    def sink(self) -> NodeId :
        """Return the sink node in the graph."""
        return self._sink

    @sink.setter
    def sink(self, sink_ : NodeId) :
        """Set the sink node in the graph. The sink node should have no outdegree edges."""
        assert self._graph.has_node(sink_)
        assert self._graph.out_degree(sink_) == 0
        self._sink = sink_

    def nodes(self, *args, **argv) :
        """Return the list of nodes with or without their attributes. It depends on the given arguments. Actually, just forward the arguments to the nodes() method of Networkx Digraph class."""
        return self._graph.nodes(*args, **argv)

    def edges(self, *args, **argv) :
        """Return the list of edges with or without their attributes. It depends on the given arguments. Actually, just forward the arguments to the edges() method of Networkx Digraph class."""
        return self._graph.edges(*args, **argv)

    def in_degree(self, node_id_ : NodeId):
        """Return the indegree of the given node."""
        return self._graph.in_degree(node_id_)

    def out_degree(self, node_id_ : NodeId):
        """Return the outdegree of the given node."""
        return self._graph.out_degree(node_id_)

    def prun(self, edge_threshold_ : int) :
        """Remove the edges with weights lower than given threshold weight and update the graph and its attribute to make it plausible.
        
        Initially, edges are filtered. It make some obsolete nodes, which has no zero in- and out-degree other than source and sink nodes. Such nodes are excluded from the graph.
        Update the self._len. Then, connect all zero indegree nodes other than source to source node by generating an edge between them. Likewise, connect all zero outdegree nodes other than sink to sink node.

        Args:
            edge_threshold_:
                An integer indicating the threshold of weight for filtering edges in the graph.
        Returns:
            Updated object itself.
        """
        assert edge_threshold_ > 0
        self._edge_threshold = edge_threshold_
        self.filter_by_edge_weight(self._edge_threshold)
        self.clean_obsolete_nodes()
        self.len(len(self._graph))
        self.leave_only_one_source()
        self.leave_only_one_sink()

        return self

    def filter_by_edge_weight(self, threshold_: int) :
        """Remove all edges with weight<threshold from the graph."""
        edges = [(n1, n2) for n1, n2, w in self.edges(data="weight") if int(w) < threshold_]
        self._graph.remove_edges_from(edges)

    def clean_obsolete_nodes(self):
        """Remove zero in- and out-degree nodes other than source and sink from the graph."""
        obsolete_nodes = []
        for node_id in self.nodes():
            if node_id == self.source or node_id == self.sink:
                continue
            if self.in_degree(node_id) == 0 and self.out_degree(node_id) == 0:
                obsolete_nodes.append(node_id)

        self._graph.remove_nodes_from(obsolete_nodes)
    
    def leave_only_one_source(self):
        """Connect all zero indegree nodes other than source to source"""
        sources = []
        for node_id in self.nodes():
            if node_id == self.source or node_id == self.sink:
                continue
            if self.in_degree(node_id) == 0:
                sources.append(node_id)
        
        self._graph.add_edges_from([(self.source, x) for x in sources])

    def leave_only_one_sink(self):
        """Connect all zero outdegree nodes other than sink to sink"""
        sinks = []
        for node_id in self.nodes():
            if node_id == self.source or node_id == self.sink:
                continue
            if self.out_degree(node_id) == 0:
                sinks.append(node_id)
        
        self._graph.add_edges_from([(x, self.sink) for x in sinks])

    def set_node_attributes(self, *args, **argv) :
        """Set the node attributes. Actually, just forward the arguments to the set_node_attributes() function of Networkx library."""
        nx.set_node_attributes(self._graph, *args, **argv)

    def get_node_attributes(self, *args, **argv) :
        """Return the node attributes. Actually, just forward the arguments to the get_node_attributes() function of Networkx library."""
        return nx.get_node_attributes(self._graph, *args, **argv)

    def set_edge_attributes(self, *args, **argv) :
        """Set the edge attributes. Actually, just forward the arguments to the set_edge_attributes() function of Networkx library."""
        nx.set_edge_attributes(self._graph, *args, **argv)

    def get_edge_attributes(self, *args, **argv) :
        """Return the edge attributes. Actually, just forward the arguments to the get_edge_attributes() function of Networkx library."""
        return nx.get_edge_attributes(self._graph, *args, **argv)

    def predecessors(self, node_id_ : NodeId) -> List[NodeId] :
        """Return the predecessor of the node with given node id. A predecessor of n is a node m such that there exists a directed edge from m to n."""
        return self._graph.predecessors(node_id_)

    def set_ancestors(self):
        """Find parent, grandparent, grandgrandparent of each node and set the list of them as attribute called 'ancestors' to each node in the graph.
           A parent is equivalent to a predecessor of a node.
           A grandparent is equivalent to a predecessor of parent.
           A grandgrandparent is equivalent to a predecessor of a grandparent.
        """
        attrs = defaultdict(lambda:{'ancestors': []})

        for node_id in self.ordering:
            attrs[node_id]['ancestors'] = []
            for parent in self.predecessors(node_id):
                for grandparent in self.predecessors(parent):
                    for grandgrandparent in self.predecessors(grandparent):
                        ancestors = [parent, grandparent, grandgrandparent]
                        attrs[node_id]['ancestors'].append(ancestors)
        
        self.set_node_attributes(attrs)

    def ancestors(self, node_id_ : NodeId) -> List[NodeId] :
        """Return the ancestors of the node with given node id."""
        return self.get_node_attributes('ancestors')[node_id_]

    def bases(self) -> Dict[NodeId, Base] :
        """Return dictionary of bases of nodes with their node id as keys."""
        return self.nodes(data="base")

    def base(self, node_id_ : NodeId) -> Base :
        """Return the base of the node with given node id."""
        return self.bases()[node_id_]

    def export_dot(self, dot_path_ : Path) :
        """Export the DAG with updated attributes as graphviz DOT format.
           Set node's base as node's label and edge's weight as edge's label.
           Fill the node that is in consensus sequene only before correction with blue color.
           Fill the node that is in consensus sequene only after correction with red color.
           Fill the node that is in consensus sequene both before and after correction with purple color.

           Args:
              dot_path_: A path obeject indicating the path of the file to which graph notation in DOT format is written.
        """
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

        nx.drawing.nx_pydot.write_dot(self._graph, str(dot_path_))

    def base_quote_strip(self):
        """Strip the double quotes on both ends of base attribute of each node in the graph."""
        attrs = defaultdict(dict)
        for node_id, original_base in self.bases():
            stripped_base = original_base.strip('"')
            attrs[node_id]['base'] = stripped_base
            self.set_node_attributes(attrs)

    def reverse(self):
        """Reverse the graph."""
        self._graph = self._graph.reverse()

    def complement(self):
        """Change each base of nodes in the graph as complimentary base."""
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
        """Swap source and sink node."""
        tmp : NodeId = self.sink
        self.sink = copy.copy(self.source)
        self.source = tmp

    def reverse_complement(self):
        """Reverse and complement the graph and swap source and sink to make it plausible."""
        self.reverse()
        self.complement()
        self.swap_source_sink()

    def leave_only_indels(self, path_ : List[NodeId]):
        pass

    def trim_path(self, path_ : List[NodeId]):
        """Trim the both ends of the given path until any initial consensus node is encountered.
           Corrected consensus sequence should be shorter than initial consensus sequence.
           Exclude source and sink nodes, if included.
        """
        initial_consensus_tf = dict(self.nodes(data='initial_consensus', default=False))
        start_index = next(index for index, node_id in enumerate(path_) if initial_consensus_tf[node_id])
        end_index = next(index for index, node_id in reversed(list(enumerate(path_))) if initial_consensus_tf[node_id])

        if path_[start_index] == self.source:
            start_index += 1
        if path_[end_index] == self.sink:
            end_index -= 1

        return path_[start_index:end_index+1]

    def set_corrected_consensus_attribute(self, path_ : List[NodeId]):
        """Set the attribute called 'corrected_consensus' as true if a node is included in the given path."""
        attrs = defaultdict(dict)
        for node_id in path_:
            attrs[node_id]["corrected_consensus"] = 'true'
        for node_id in self.nodes():
            if node_id not in attrs:
                attrs[node_id]["corrected_consensus"] = 'false'
        dag.set_node_attributes(attrs)

    def topological_sort(self) -> List[NodeId] :
        """Sort the nodes in topogical order."""
        return [NodeId(node_id) for node_id in nx.algorithms.dag.topological_sort(self._graph)]

    def extract_base_with_node_id_list(self, node_id_list_ : List[NodeId]) -> str:
        """Return the sequence concatenating bases extracted from each node with the node id in given list."""
        sequence : str = ""
        for node_id in node_id_list_:
            sequence += self.base(node_id)
        return sequence

    def set_node_id_to_index(self) :
        """Indexize node ids to index in self._ordering in order to make them fit to numpy operations."""
        self._node_id_to_index_dict : Dict[NodeId, Idx] = { NodeId(val) : Idx(idx) for idx, val in enumerate(self.ordering) }

    def node_id_to_index(self, node_id_ : NodeId) -> Idx :
        """Return index in self._ordering of the node with given node id."""
        return self._node_id_to_index_dict[node_id_]

    def index_to_node_id(self, idx_ : Idx) -> NodeId :
        """Return the node id corresponding given index in self._ordering."""
        return self.ordering[idx_]
        

def read_consensuses(consensuses_path_) -> List[SeqRecord]:
    """Read the consensuses from fasta.
       Store them as list of SeqRecord objects in Biopython.
    """
    return SeqIO.parse(consensuses_path_, "fasta")

def read_phmmDB(phmmDB_path_ : Path) -> Dict[HMMId, HMM]:
    """Read the profile HMMs.
       Store them as dictionary of HMM objects in profileHMM module with each HMM's id as a key.
    """
    f = phmmDB_path_.open()
    model_generator = reader.read_all(f)
    
    return  {HMMId(x.metadata.model_name): x for x in model_generator}

def read_dot(dot_path_ : Path) -> List[DAG] :
    """Read the DAGs from given DOT file.
       store them as list of DAG objects."""
    graphs = pydot.graph_from_dot_file(dot_path_)
    
    graphs = map(nx.nx_pydot.from_pydot, graphs)
    
    dags = [DAG(graph) for graph in graphs]

    return dags

def parse_hmmscanresult(hmmscan_domtbl_path_ : Path) -> List[HSP]:
    """Read the hmmscan domain hits from domtbl file.
       Store them as list of HSPs object in Biopython SearchIO module.
       For unique consensus and hmm pair, only a hit with the best score is stored for reducing redundancies."""
    domtbl = SearchIO.parse(hmmscan_domtbl_path_, "hmmscan3-domtab")
    
    uniq_consensus_HMM_dict = {}
    # for each pair of consensus and HMM,
    # only the hsp with the best score is left
    for query in domtbl:
        tig_id = query.id.split("_rframe")[0]
        for hit in query:
            hmm_id = hit.id
            for hsp in hit:
                if (tig_id, hmm_id) not in uniq_consensus_HMM_dict:
                    uniq_consensus_HMM_dict[(tig_id, hmm_id)] = hsp
                else:
                    curr_best_hsp = uniq_consensus_HMM_dict[(tig_id, hmm_id)]
                    if hsp.bitscore > curr_best_hsp.bitscore:
                        uniq_consensus_HMM_dict[(tig_id, hmm_id)] = hsp

    def extract_rframe(hsp_ : HSP):
        """Extract reading frame and direction of hmmscan domain hit."""
        hsp_.tig_id = hsp_.query_id.split("_rframe")[0]
        hsp_.reading_frame = hsp_.query_id.split("_rframe")[1]
        hsp_.direction = "-" if hsp_.reading_frame[0]  == "-" else "+"
        hsp_.reading_frame = abs(int(hsp_.query_id.split("_rframe")[1]))
    
    for hsp in uniq_consensus_HMM_dict.values():
        extract_rframe(hsp)
    
    return list(uniq_consensus_HMM_dict.values())

def pair_consensus_DAG(consensuses_ : List[SeqRecord], dags_ : List[DAG]) -> DefaultDict[TigId, Tuple[SeqRecord, DAG]]:
    """Pair a consensus and DAG with a matching consensus id."""
    consensuses_dict = {x.id: x for x in consensuses_}
    dags_dict = {x.name: x for x in dags_}
    pair_dict : DefaultDict[TigId, Tuple[SeqRecord, DAG]] = defaultdict(tuple)

    for key in consensuses_dict:
        pair_dict[TigId(key)] = (consensuses_dict[key], dags_dict[key])

    return pair_dict

def correct_consensus(dag_ : DAG, phmm_ : HMM) -> List[NodeId] :
    '''Return corrected consensus path in DAG using Viterbi algorithm.'''
    PHMM = profileHMM.PHMM(phmm_)
    return PHMM.modified_viterbi(dag_)

def write_corrected_consensus(corrected_sequence_ : str, consensus_id_ : TigId, phmm_id_ : HMMId, f_ : Path):
    """Write the corrected consensus as a FASTA format with the header including consensus ID and HMM ID used for correction"""
    seqid = f'{consensus_id_}_corrected_with_{phmm_id_}'
    record = SeqRecord(
        Seq(corrected_sequence_),
        id=seqid,
        description=f'{consensus_id_} len={len(corrected_sequence_)}')

    SeqIO.write(record, f_, "fasta")

if __name__ == '__main__':
    consensuses_path : Path = Path("consensuses.fasta")
    phmmDB_path : Path = Path("FAM173.hmm")
    dot_path : Path = Path("graph.dot")
    hmmscan_domtbl_path : Path = Path("hmmscan.domtbl")
    output_path : Path = Path("corrected_consensuses.fasta")
    output_dot : Path = Path("output.dot")
    minimum_edge_weight : int = 2
    
    phmms = read_phmmDB(phmmDB_path)
    hmmscan_domtbl = parse_hmmscanresult(hmmscan_domtbl_path)

    consensuses = read_consensuses(consensuses_path)
    dags = read_dot(dot_path)
    consensus_dag = pair_consensus_DAG(consensuses, dags)

    for hsp in hmmscan_domtbl:
        consensus_id       = hsp.query_id.split("_rframe")[0]
        matched_phmm_id = hsp.hit_id

        consensus          = consensus_dag[consensus_id][0]
        dag             = consensus_dag[consensus_id][1]
        matched_phmm    = phmms[matched_phmm_id]

        # initial base attribute is wrapped with double quotes
        dag.base_quote_strip()

        # if hit appeared in reverse direction, dag also reversed
        if hsp.direction == "-":
            dag.reverse_complement()

        # filter edges with low weight
        dag.prun(minimum_edge_weight)

        # update dag's properties
        dag.ordering = dag.topological_sort()
        dag.set_ancestors()

        # correct initial consensus
        corrected_path = correct_consensus(dag, matched_phmm)
        corrected_path = dag.trim_path(corrected_path)
        dag.set_corrected_consensus_attribute(corrected_path)

        # export results
        dag.export_dot(output_dot)
        corrected_sequence = dag.extract_base_with_node_id_list(corrected_path)
        write_corrected_consensus(corrected_sequence, consensus_id, matched_phmm_id, output_path)