"""
Basically the debais prep code, made more readable with slight changes
USE ONTOLOGY.IO to do this?
"""

import networkx as nx
import pickle as cp
################################ GRAPH HELPER TOOLS ###############################
# Some constants
ROOT_BPO = 'GO:0008150'
ROOT_CCO = 'GO:0005575'
ROOT_MFO = 'GO:0003674'


def parseGOTerms(obo_path):
    '''
    Create graphs from an OBO file.   
    
    Input:
    obo_path : String      File location of OBO
    
    Output:
    [0]      : MFO Graph
    [1]      : BPO Graph
    [2]      : CCO Graph
    [3]      : AltID dictionary         
    '''
    
    obo_file = open(obo_path,"r") 
    # Generate graph objects
    mf_g = nx.DiGraph()
    cc_g = nx.DiGraph()
    bp_g = nx.DiGraph()
    # Read in all Terms
    allGOterms = obo_file.read().split("[Term]")
    # Intialize variables
    GO_term   = ""
    namespace = ""
    relation  = ""
    # Alt name dictionary 
    alt_id_to_id = dict()
    # Lists for each namespace
    mf = []
    bp = []
    cc = []
    # Add terms to respective namspace lists
    for term in allGOterms[1:]:
        split_term = term.split("\n")
        for line in split_term:
            if "id:" in line and "GO:" in line and "alt_id" not in line:
                GO_term = "GO:"+line.split("GO:")[-1].strip()
            if   "namespace: biological_process" in line:
                namespace = "bp"
                bp.append(GO_term)
            elif "namespace: cellular_component" in line:
                namespace = "cc"   
                cc.append(GO_term)
            elif "namespace: molecular_function" in line:
                namespace = "mf"
                mf.append(GO_term)
    # Populate graphs -> Have to do this before adding edges!
    mf_g.add_nodes_from(mf)
    bp_g.add_nodes_from(bp)
    cc_g.add_nodes_from(cc)
    # For each Term text object gotten form the OBO
    for term in allGOterms[1:]:
        split_term = term.split("\n")
        # Read term info in by line
        for line in split_term:
            # Get the Go Term
            if "id:" in line and "GO:" in line and "alt_id" not in line:
                GO_term = "GO:" + line.split("GO:")[-1].strip()
            # Add to ALT_ID if needed
            if "GO:" in line and "alt_id" in line:
                alt_id = "GO:" + line.split("GO:")[-1].strip()
                alt_id_to_id[alt_id] = GO_term
            # Get the Namespace for this entry
            if   "namespace: biological_process" in line:
                namespace = "bp"
            elif "namespace: cellular_component" in line:
                namespace = "cc"   
            elif "namespace: molecular_function" in line:
                namespace = "mf"
            # Add edges based on the relationships in the entry
            if(("is_a:" in line or "relationship: part_of" in line) and "GO" in line):
                
                if "is_a" in line:
                    relation = "is_a"
                else:
                    relation = "part_of"
                try:
                    # Parent GO Term is listed if in a relationship
                    parent_GO_term = "GO:" + line.split("GO:")[1][:7]
                except IndexError:
                    print(line)
                    exit()
                
                #print(GO_term,line)
                if   namespace == "bp" and parent_GO_term in bp:
                    bp_g.add_edge(GO_term, parent_GO_term, weight = relation)
                elif namespace == "cc" and parent_GO_term in cc:
                    cc_g.add_edge(GO_term, parent_GO_term, weight = relation)
                elif namespace == "mf" and parent_GO_term in mf:
                    mf_g.add_edge(GO_term, parent_GO_term, weight = relation)

    return mf_g, bp_g, cc_g, alt_id_to_id



def findAllAncestorsTerm(graph, node):
    '''
    Find all Ancestors for a given Term
    Recursive function  
    
    Input:
    graph : The graph to look on
    node  : The current Node
    
    Output:
    [0]   : List[Set of Ancestor Terms]
    '''
    # The Term is an ancestor of itself
    ancestors = [node]
    # Recursively extend ancestors
    for immediate_ancestor in graph.successors(node):
        # Add all immediate ancestors
        ancestors.append(immediate_ancestor)
        ancestor_temp = findAllAncestorsTerm(graph, immediate_ancestor)
        # Add all ancestors of that ancestor
        ancestors.extend(ancestor_temp)
    return list(set(ancestors))
    

def findAllAncestors(filename):
    '''
    Find Ancestors for all Terms in the Ontology
    
    Input:
    filename : The name of the graph
    
    Output:
    [0]      : Dictionary KEY: Nodes VALUE: List[Set of Ancestors]
    '''
    
    # Load graph
    graph = cp.load( open(filename, "rb" ) )
    ancestors = dict()
    # For all Nodes, find their ancestors
    for term in graph.nodes():
        ancestors[term] = findAllAncestorsTerm(graph, term)
        
    return ancestors


def setupGraphs(obo_path):
    '''
    Code used to create the ontology graphs based on a given OBO
    
    Input:
    obo_path : String      File path for the OBO
    '''
    
    # Get needed information
    mf_g, bp_g, cc_g, alt_id_to_id = parseGOTerms(obo_path)
    # Save to file
    cp.dump(mf_g, open("ICdata/MFO.graph", "wb"))
    cp.dump(bp_g, open("ICdata/BPO.graph", "wb"))
    cp.dump(cc_g, open("ICdata/CCO.graph", "wb"))
    cp.dump(alt_id_to_id, open("ICdata/ALTERNATE_ID_TO_ID.map", "wb"))
    # Load from file
    # Find ancestors
    # Save to File
    # Clear Memory
    mf_g = cp.load(       open("ICdata/MFO.graph", "rb" ) )
    mf_ancestors = findAllAncestors("ICdata/MFO.graph")
    cp.dump(mf_ancestors, open("ICdata/MFO_ANCESTORS.graph", "wb"))
    del mf_g, mf_ancestors
    
    bp_g = cp.load(       open("ICdata/BPO.graph", "rb" ) )
    bp_ancestors = findAllAncestors("ICdata/BPO.graph")
    cp.dump(bp_ancestors, open("ICdata/BPO_ANCESTORS.graph", "wb"))
    del bp_g, bp_ancestors
    
    cc_g = cp.load(       open("ICdata/CCO.graph", "rb" ) )
    cc_ancestors = findAllAncestors("ICdata/CCO.graph")
    cp.dump(cc_ancestors, open("ICdata/CCO_ANCESTORS.graph", "wb"))
    del cc_g, cc_ancestors