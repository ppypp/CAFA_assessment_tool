# -*- coding: utf-8 -*-
"""
Basically the debais prep code, made more readable with slight changes
USE ONTOLOGY.IO to do this?
"""

import networkx as nx
import pickle as cp
################################ GRAPH HELPER TOOLS ###############################
# Some filenames
ROOT_BPO='GO:0008150'
ROOT_CCO='GO:0005575'
ROOT_MFO='GO:0003674'


def parseGOTerms(obo_path):
    '''
    Create graphs from an OBO file.   
    
    Input:
    obo_path : String      File location of OBO
    
    Output:
    [0]     : MFO Graph
    [1]     : BPO Graph
    [2]     : CCO Graph
    [3]     : AltID map
                
    '''
    
    # USe same OBO as the assessment tool ##################Make this user configurable
    fhr = open(obo_path,"r") 
    
    mf_g = nx.DiGraph()
    cc_g = nx.DiGraph()
    bp_g = nx.DiGraph()

    allGOterms = fhr.read().split("[Term]")

    GO_term = ""
    namespace = ""
    relation = ""

    alt_id_to_id = dict()

    mf = []
    bp = []
    cc = []

    for term in allGOterms[1:]:
        split_term = term.split("\n")
        for line in split_term:
            if "id:" in line and "GO:" in line and "alt_id" not in line:
                GO_term = "GO:"+line.split("GO:")[-1].strip()
            if "namespace: biological_process" in line:
                namespace = "bp"
                bp.append(GO_term)
            elif "namespace: cellular_component" in line:
                namespace = "cc"   
                cc.append(GO_term)
            elif "namespace: molecular_function" in line:
                namespace = "mf"
                mf.append(GO_term)
    
    mf_g.add_nodes_from(mf)
    bp_g.add_nodes_from(bp)
    cc_g.add_nodes_from(cc)
    
    for term in allGOterms[1:]:
        split_term=term.split("\n")
        #alt_ids=[]
        for line in split_term:
            if "id:" in line and "GO:" in line and "alt_id" not in line:
                GO_term="GO:"+line.split("GO:")[-1].strip()
            if "GO:" in line and "alt_id" in line:
                alt_id="GO:"+line.split("GO:")[-1].strip()
                alt_id_to_id[alt_id] = GO_term
            if "namespace: biological_process" in line:
                namespace="bp"
            elif "namespace: cellular_component" in line:
                namespace = "cc"   
            elif "namespace: molecular_function" in line:
                namespace = "mf"
            
            if(("is_a:" in line or "relationship: part_of" in line) and "GO" in line):
                
                if "is_a" in line:
                    relation = "is_a"
                else:
                    relation = "part_of"
                try:
                    parent_GO_term = "GO:" + line.split("GO:")[1][:7]
                except IndexError:
                    print(line)
                    exit()
                
                #print(GO_term,line)
                if namespace == "bp" and parent_GO_term in bp:
                    bp_g.add_edge(GO_term, parent_GO_term, weight = relation)
                elif namespace == "cc" and parent_GO_term in cc:
                    cc_g.add_edge(GO_term, parent_GO_term, weight = relation)
                elif namespace == "mf" and parent_GO_term in mf:
                    mf_g.add_edge(GO_term, parent_GO_term, weight = relation)

    return mf_g, bp_g, cc_g, alt_id_to_id


def findAllCommonAncestorsAndDisjointCommonAncestors(mf_g, bp_g, cc_g, go1, go2): ################NEVER USED ###########################
    '''
    Description    
    
    Input:

    
    Output:
    [0]   :
    [1]   :
    '''
    
    mf, bp, cc = 0, 0, 0
    if mf_g.has_node(go1):
        all_paths_go1 = nx.all_simple_paths(mf_g,source = go1, target = ROOT_MFO)
        mf += 1
    elif bp_g.has_node(go1):
        all_paths_go1 = nx.all_simple_paths(bp_g,source = go1, target = ROOT_BPO)
        bp += 1
    elif cc_g.has_node(go1):
        all_paths_go1 = nx.all_simple_paths(cc_g,source = go1, target = ROOT_CCO)
        cc += 1

    
    if mf_g.has_node(go2):
        all_paths_go2 = nx.all_simple_paths(mf_g,source = go2,target = ROOT_MFO)
        mf += 1
    elif bp_g.has_node(go2):
        all_paths_go2 = nx.all_simple_paths(bp_g,source = go2,target = ROOT_BPO)
        bp += 1
    elif cc_g.has_node(go2):
        all_paths_go2 = nx.all_simple_paths(cc_g,source = go2,target = ROOT_CCO)
        cc += 1

    all_paths_go1 = list(all_paths_go1)
    all_paths_go2 = list(all_paths_go2)
    
    #print(mf, bp, cc)
    if mf == 2:
        current_graph = mf_g
        root_node = ROOT_MFO
    elif bp == 2:
        current_graph = bp_g
        root_node = ROOT_BPO
    elif cc == 2:
        current_graph = cc_g
        root_node = ROOT_CCO
    
        
    common_ancestors = []
    for path1 in all_paths_go1:
        for path2 in all_paths_go2:
            print(path1)
            print(path2)
            print(list(set(path1) & set(path2)))
            print("*"*50)
            if list(set(path1) & set(path2)) not in common_ancestors:
                common_ancestors.append(list(set(path1) & set(path2)))
            #common_ancestors=list(set(common_ancestors))
    dca=[]
    for each_path in common_ancestors:
        node_to_root_distance=dict()
        for node in each_path:
            #if node is not root_node:
            node_to_root_distance[node] = nx.shortest_path_length(current_graph,source = node,target = root_node)
        print(node_to_root_distance)
        if len(node_to_root_distance) != 0:
            key, _ = max(node_to_root_distance.iteritems(), key = lambda x:x[1])
            if key not in dca:
                dca.append(key)
    
    return common_ancestors,dca

def findAllAncestors(graph, node):
    '''
    Recursive function  
    
    Input:
    graph : The graph to look on
    node  : The current Node
    
    Output:
    [0]   : List[Set of Ancestors]
    '''
    
    ancestors = [node]
    for immediate_ancestor in graph.successors(node):
        ancestors.append(immediate_ancestor)
        ancestor_temp = findAllAncestors(graph, immediate_ancestor)
        ancestors.extend(ancestor_temp)
    return list(set(ancestors))
    

def findAllAncestorsForAllNodesForOntology(filename):
    '''
    Description
    
    Input:
    filename : The name of the graph
    
    Output:
    [0]      : Dictionary KEY: Nodes VALUE: List[Set of Ancestors]
    '''
    
    # Load graph
    graph = cp.load( open(filename, "rb" ) )
    graph_ancestors = dict()
    # For all Nodes, find their ancestors
    for nodes in graph.nodes():
        graph_ancestors[nodes] = findAllAncestors(graph, nodes)
        
    return graph_ancestors


def setupGraphs(obo_path):
    '''
    Code used to create the ontology graphs based on a given OBO
    
    Input:
    obo_path : String      File path for the OBO
    '''
    
    # Get needed information
    mf_g, bp_g, cc_g, alt_id_to_id = parseGOTerms("../" + obo_path)
    # Save to file
    cp.dump(mf_g, open("MFO.graph", "wb"))
    cp.dump(bp_g, open("BPO.graph", "wb"))
    cp.dump(cc_g, open("CCO.graph", "wb"))
    cp.dump(alt_id_to_id, open("ALTERNATE_ID_TO_ID.map", "wb"))
    # Load from file
    mf_g = cp.load( open("MFO.graph", "rb" ) )
    bp_g = cp.load( open("BPO.graph", "rb" ) )
    cc_g = cp.load( open("CCO.graph", "rb" ) )
    # Find ancestors
    mf_ancestors = findAllAncestorsForAllNodesForOntology("MFO.graph")
    bp_ancestors = findAllAncestorsForAllNodesForOntology("BPO.graph")
    cc_ancestors = findAllAncestorsForAllNodesForOntology("CCO.graph")
    # Save to File
    cp.dump(mf_ancestors, open("MFO_ANCESTORS.graph", "wb"))
    cp.dump(bp_ancestors, open("BPO_ANCESTORS.graph", "wb"))
    cp.dump(cc_ancestors, open("CCO_ANCESTORS.graph", "wb"))
    