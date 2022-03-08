from re import L
from typing import DefaultDict, List, Dict, Set, Tuple
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
import time

GraphId = NewType('GraphId', str)
NodeId = NewType('NodeId', str)
EdgeId = NewType('EdgeId', str)
TigId = NewType('TigId', str)
HMMId = NewType('HMMId', str)
Base = NewType('Base', str)
Idx = NewType('Idx', int)

def resolve_cycle(digraph_: nx.DiGraph) -> nx.DiGraph :

    # make dictionary of node_id and node attributes
    node_attr_dict = {}
    for node_id, attr_dict in digraph_.nodes(data=True):
        node_attr_dict[node_id] = attr_dict

    # make dictionary of edge ids and edge attributes
    edge_attr_dict = {}
    for out_node_id, in_node_id, attr_dict in digraph_.edges(data=True):
        edge_attr_dict[(out_node_id, in_node_id)] = attr_dict

    print(edge_attr_dict)
    
    for e in list(digraph_.edges(data=False)):
        out_node, in_node = e

        # if cycle length of one or two
        if out_node == in_node or digraph_.has_edge(in_node, out_node):
            if not digraph_.has_edge(out_node, in_node):
                continue

            # create duplicate node and copy the node attributes
            duplicate_in_node_id = f"{in_node}_1"
            digraph_.add_node(duplicate_in_node_id)
            nx.set_node_attributes(digraph_, {duplicate_in_node_id: node_attr_dict[in_node]})

            for predecessor in digraph_.predecessors(in_node):
                if predecessor == in_node or predecessor == out_node:
                    continue
                else:
                    digraph_.add_edge(predecessor, in_node)

            for successor in digraph_.successors(in_node):
                if successor == in_node or successor == out_node:
                    continue
                else:
                    digraph_.add_edge(in_node, successor)

            if out_node == in_node:
                digraph_.add_edge(in_node, duplicate_in_node_id)
                nx.set_edge_attributes(digraph_, {(in_node,duplicate_in_node_id): edge_attr_dict[(in_node, in_node)]})
                digraph_.remove_edge(in_node, in_node)
                
            if digraph_.has_edge(in_node, out_node):
                duplicate_out_node_id = f"{out_node}_1"
                digraph_.add_node(duplicate_out_node_id)
                nx.set_node_attributes(digraph_, {duplicate_out_node_id: node_attr_dict[out_node]})

                for predecessor in digraph_.predecessors(out_node):
                    if predecessor == in_node:
                        continue
                    else:
                        digraph_.add_edge(predecessor, out_node)
                        nx.set_edge_attributes(digraph_, {(predecessor, out_node): edge_attr_dict[(predecessor, out_node)]})

                for successor in digraph_.successors(out_node):
                    if successor == in_node:
                        continue
                    else:
                        digraph_.add_edge(out_node, successor)
                        nx.set_edge_attributes(digraph_, {(out_node, successor): edge_attr_dict[(out_node, successor)]})

                # transfer the incoming edges to duplicate node
                digraph_.add_edge(in_node, duplicate_out_node_id)
                nx.set_edge_attributes(digraph_, {(in_node, duplicate_out_node_id): edge_attr_dict[(in_node, out_node)]})
                digraph_.add_edge(out_node, duplicate_in_node_id)
                nx.set_edge_attributes(digraph_, {(out_node, duplicate_in_node_id): edge_attr_dict[(out_node, in_node)]})
                digraph_.remove_edge(in_node, out_node)
                digraph_.remove_edge(out_node, in_node)

    return digraph_
            
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
        self._ordering = list(nx.algorithms.dag.topological_sort(self._graph))
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

    def update_ordering(self):
        """Update ordering"""
        self.ordering = list(nx.algorithms.dag.topological_sort(self._graph))

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
        """
        assert edge_threshold_ > 0
        self._edge_threshold = edge_threshold_
        self.filter_by_edge_weight(self._edge_threshold)
        self.clean_obsolete_nodes()
        self.len(len(self._graph))
        self.leave_only_one_source()
        self.leave_only_one_sink()

    def filter_by_edge_weight(self, threshold_: int) :
        """Remove all edges with weight<threshold from the graph."""
        edges = []
        for (n1, n2, w) in self.edges(data="weight"):
            if w == None or int(w) < threshold_:
                edges.append((n1, n2))
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
        
        self.add_edges_from([(self.source, x) for x in sources])

    def leave_only_one_sink(self):
        """Connect all zero outdegree nodes other than sink to sink"""
        sinks = []
        for node_id in self.nodes():
            if node_id == self.source or node_id == self.sink:
                continue
            if self.out_degree(node_id) == 0:
                sinks.append(node_id)
        
        self.add_edges_from([(x, self.sink) for x in sinks])

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

    def add_node(self, *args, **argv) :
        """Add a node. Actually, just forward the arguments to the add_node() function of Networkx library."""
        return self._graph.add_node(*args, **argv)

    def add_nodes_from(self, *args, **argv) :
        """Add a node. Actually, just forward the arguments to the add_nodes_from() function of Networkx library."""
        return self._graph.add_nodes_from(*args, **argv)

    def add_edge(self, *args, **argv) :
        """Add an edge. Actually, just forward the arguments to the add_edge() function of Networkx library."""
        return self._graph.add_edge(*args, **argv)

    def add_edges_from(self, *args, **argv) :
        """Add an edge. Actually, just forward the arguments to the add_edges_from() function of Networkx library."""
        return self._graph.add_edges_from(*args, **argv)

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
    
    def update_base_dict(self) :
        self._base_dict = dict(self.bases())

    def base(self, node_id_ : NodeId) -> Base :
        """Return the base of the node with given node id."""
        return self._base_dict[node_id_]

    def neighbors(self, node_id_ : NodeId) :
        """Return an iterator over successor nodes of node_id."""
        return self._graph.neighbors(node_id_)

    def export_dot(self, dot_path_ : Path, horizontal_ = True) :
        """Export the DAG with updated attributes as graphviz DOT format.
           Set node's base as node's label and edge's weight as edge's label.
           Fill the node that is in consensus sequene only before correction with blue color.
           Fill the node that is in consensus sequene only after correction with red color.
           Fill the node that is in consensus sequene both before and after correction with purple color.

           Args:
              dot_path_: A path obeject indicating the path of the file to which graph notation in DOT format is written.
              
              Keyword:
                horizontal: A boolean indicating the orientation of graph when displayed
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

        # determine the orientation of graph when displayed
        if horizontal_:
            self._graph.graph['graph']={'rankdir':'LR'}

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
        self.set_node_attributes(attrs)

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
        self._node_id_to_index_dict : Dict[NodeId, Idx] = { NodeId(val) : Idx(idx) for idx, val in enumerate(self._ordering) }

    def node_id_to_index(self, node_id_ : NodeId) -> Idx :
        """Return index in self._ordering of the node with given node id."""
        return self._node_id_to_index_dict[node_id_]

    def index_to_node_id(self, idx_ : Idx) -> NodeId :
        """Return the node id corresponding given index in self._ordering."""
        return self._ordering[idx_]

    def aggressive_mode(self) :
        self.aggressive_deletion()
        self.aggressive_insertion()

    def aggressive_deletion(self) :
        """Add synthetic nodes and edges to cope with deletion errors in homopolymer regions"""
        new_nodes = []
        new_edges = []
        for node_id in self.nodes():
            self.set_node_attributes({node_id: {"synthetic":False}})
            if node_id == self.source or node_id == self.sink:
                continue
            else:
                duplicate_node_id = "dupicate_" + node_id
                new_nodes.append((duplicate_node_id, dict(base=self.base(node_id), synthetic=True, initial_consensus=False)))
                new_edges.append((node_id, duplicate_node_id))
                for neighbor_id in self.neighbors(node_id):
                    new_edges.append((duplicate_node_id, neighbor_id))
        
        self.add_nodes_from(new_nodes)
        self.add_edges_from(new_edges)

    def aggressive_insertion(self) :
        """Add synthetic nodes to cope with insertion errors in homopolymer regions"""
        for node_id in self.ordering:
            curr_base = self.base(node_id)

            curr_node_id = node_id
            homopol = False
            homopol_start_node_id = curr_node_id
            
            while True:
                for predecessor in self.predecessors(curr_node_id):
                    if self.get_node_attributes('synthetic')[predecessor]:
                        continue
                    if self.base(predecessor) == curr_base:
                        homopol_start_node_id = predecessor
                        homopol = True
                if homopol_start_node_id != curr_node_id:
                    curr_node_id = homopol_start_node_id
                else:
                    break

            if homopol:
                for predecessor in self.predecessors(homopol_start_node_id):
                    self.add_edge(predecessor, node_id)

    def leave_only_indel_errors(self, consensus_path : List[NodeId]) :
        """If corrected subpath has the same length with original subpath, correction is cancelled. 
           This is for excluding correction of substitution erorrs that could hide true variants.
        """
        for idx, node_id in enumerate(consensus_path):

            if self.in_degree(node_id) > 1:
                print("node_id:", node_id)
                init_consensus_predecessor_id=""
                corr_consensus_predecessor_id=""
                for predecessor_id in self.predecessors(node_id):
                    if self.get_node_attributes('initial_consensus')[predecessor_id] == 'true' and self.get_node_attributes('corrected_consensus')[predecessor_id] == 'false':
                        print(predecessor_id, " init:", self.get_node_attributes('initial_consensus')[predecessor_id])
                        init_consensus_predecessor_id=predecessor_id

                    if self.get_node_attributes('initial_consensus')[predecessor_id] == 'false' and self.get_node_attributes('corrected_consensus')[predecessor_id] == 'true':
                        print(predecessor_id, " corr:", self.get_node_attributes('corrected_consensus')[predecessor_id])
                        corr_consensus_predecessor_id=predecessor_id

                if init_consensus_predecessor_id != "" and corr_consensus_predecessor_id != "":
                    print(init_consensus_predecessor_id, corr_consensus_predecessor_id)
                    init_consensus_path = [init_consensus_predecessor_id]
                    curr_id = init_consensus_predecessor_id
                    continued = True
                    l = 1

                    while continued :
                        continued = False
                        for predecessor_id in self.predecessors(curr_id):
                            if self.get_node_attributes('initial_consensus')[predecessor_id] == 'true' and self.get_node_attributes('corrected_consensus')[predecessor_id] == 'false':
                                init_consensus_path.append(predecessor_id)
                                continued = True
                                curr_id = predecessor_id
                                l += 1

                    corr_consensus_path = [corr_consensus_predecessor_id]
                    curr_id = corr_consensus_predecessor_id
                    continued = True

                    while continued :
                        continued = False
                        for predecessor_id in self.predecessors(curr_id):
                            if self.get_node_attributes('initial_consensus')[predecessor_id] == 'false' and self.get_node_attributes('corrected_consensus')[predecessor_id] == 'true':
                                init_consensus_path.append(predecessor_id)
                                continued = True
                                curr_id = predecessor_id
                else:
                    continue

                print(init_consensus_path)
                print(corr_consensus_path)
                if len(init_consensus_path) == len(corr_consensus_path):
                    for node_id in corr_consensus_path:
                        self.set_node_attributes({node_id: {'corrected_consensus': 'false'}})

                    for node_id in init_consensus_path:
                        self.set_node_attributes({node_id: {'corrected_consensus': 'true'}})
                    
                    print(consensus_path[idx-l:idx])
                    consensus_path[idx-l:idx] = init_consensus_path[::-1]
                    print(consensus_path[idx-l:idx])

        return consensus_path
        

def read_consensuses(consensuses_path_, consensus_ids) -> Dict[TigId, SeqRecord]:
    """Read the consensuses from fasta.
       Store them as list of SeqRecord objects in Biopython.
    """
    seqrecords = {}
    for record in SeqIO.parse(consensuses_path_, "fasta"):
        if record.id in consensus_ids:
            seqrecords[TigId(record.id)] = record
    return seqrecords

def read_dot(dot_path_ : Path, consensus_ids) -> Dict[TigId, DAG] :
    """Read the DAGs from given DOT file.
       store them as list of DAG objects.
    """
    pydot_graphs = pydot.graph_from_dot_file(dot_path_)
    
    graphs = map(nx.nx_pydot.from_pydot, pydot_graphs)
    
    dags = {}
    for graph in graphs:
        dag = DAG(graph)
        if dag._name in consensus_ids:
            dags[TigId(dag._name)] = dag

    return dags

def read_phmmDB(phmmDB_path_ : Path, phmm_ids) -> Dict[HMMId, HMM]:
    """Read the profile HMMs.
       Store them as dictionary of HMM objects in profileHMM module with each HMM's id as a key.
    """
    f = phmmDB_path_.open()
    model_generator = reader.read_all(f)
    
    models = {}
    for x in model_generator:
        if x.metadata.model_name in phmm_ids:
            models[HMMId(x.metadata.model_name)] = x

    return models


def parse_hmmscanresult(hmmscan_domtbl_path_ : Path) -> List[HSP]:
    """Read the hmmscan domain hits from domtbl file.
       Store them as list of HSPs object in Biopython SearchIO module.
       For unique consensus and hmm pair, only a hit with the best score is stored for reducing redundancies.
    """
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

def correct_consensus(dag_ : DAG, phmm_ : HMM) -> List[NodeId] :
    """Return corrected consensus path in DAG using Viterbi algorithm."""
    PHMM = profileHMM.PHMM(phmm_)
    return PHMM.modified_viterbi(dag_)

def make_corrected_consensus_as_seqrecord(corrected_sequence_ : str, consensus_id_ : TigId, phmm_id_ : HMMId):
    """Write the corrected consensus as a FASTA format with the header including consensus ID and HMM ID used for correction."""
    seqid = f'{consensus_id_}_corrected_with_{phmm_id_}'
    record = SeqRecord(
        Seq(corrected_sequence_),
        id=seqid,
        description=f'{consensus_id_} len={len(corrected_sequence_)}')

    return record

if __name__ == '__main__':
    start = time.time()

    consensuses_path : Path = Path("PR8_H1N1.contigs.fasta")
    dot_path : Path = Path("PR8_H1N1.graph.dot")
    phmmDB_path : Path = Path("HMM/U-RVDBv23.0-prot.hmm")
    hmmscan_domtbl_path : Path = Path("U_RVDB_e1000.domtbl")
    output_path : Path = Path("corrected_consensuses.fasta")
    outdir : Path = Path("outdir_only_indel")
    outdir.mkdir(parents=True, exist_ok=True)
    minimum_edge_weight : int = 2
    aggr_mode : bool = False
    only_indel_mode : bool = True

    hmmscan_domtbl = parse_hmmscanresult(hmmscan_domtbl_path)

    hit_consensus_id_set : Set[TigId] = set()
    hit_phmm_id_set : Set[TigId] = set()

    for hsp in hmmscan_domtbl:
        hit_consensus_id_set.add(TigId(hsp.query_id.split("_rframe")[0]))
        hit_phmm_id_set.add(TigId(hsp.hit_id))

    consensuses = read_consensuses(consensuses_path, hit_consensus_id_set)
    dags = read_dot(dot_path, hit_consensus_id_set)
    phmms = read_phmmDB(phmmDB_path, hit_phmm_id_set)

    corrected_seqrecords = []

    for hsp in hmmscan_domtbl:
        consensus_id    = hsp.query_id.split("_rframe")[0]
        matched_phmm_id = hsp.hit_id

        consensus       = consensuses[consensus_id]
        dag             = dags[consensus_id]
        
        matched_phmm    = phmms[matched_phmm_id]
        print(f"{consensus_id}_{matched_phmm_id}")
        # initial base attribute is wrapped with double quotes
        dag.base_quote_strip()

        # if hit appeared in reverse direction, dag also reversed
        if hsp.direction == "-":
            dag.reverse_complement()
            
        # filter edges with low weight
        dag.prun(minimum_edge_weight)
        dag.update_ordering()
        dag.update_base_dict()
        
        # aggressive mode
        if aggr_mode:
            dag.aggressive_mode()
            dag.update_ordering()
            dag.update_base_dict()

        dag.set_ancestors()

        # correct initial consensus
        corrected_path = correct_consensus(dag, matched_phmm)
        dag.set_corrected_consensus_attribute(corrected_path)

        if only_indel_mode:
            corrected_path = dag.leave_only_indel_errors(corrected_path)


        # export results
        dag.export_dot(Path(outdir, f"{consensus_id}_{matched_phmm_id}.dot"), horizontal_=True)
        corrected_sequence = dag.extract_base_with_node_id_list(corrected_path)
        corrected_seqrecords.append(make_corrected_consensus_as_seqrecord(corrected_sequence, consensus_id, matched_phmm_id))

    
    SeqIO.write(corrected_seqrecords, Path(outdir, output_path), "fasta")

    end = time.time()
    print(end - start)