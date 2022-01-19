#!/usr/bin/env python


import sys
import Bio.PDB
import commands
import numpy as np
from copy import deepcopy


class Weight:
    def __init__(self, vertex1, vertex2):
        self.residues = [vertex1.resname, vertex2.resname]
        # eudlidean distances between CA - CA and CA - CB
        self.S1S2 = np.linalg.norm(vertex1.S.coord - vertex2.S.coord)
        self.S1E2 = np.linalg.norm(vertex1.S.coord - vertex2.E.coord)
        self.S2E1 = np.linalg.norm(vertex1.E.coord - vertex2.S.coord)
        self.q = self.get_dihedral(vertex1.E.coord, vertex1.S.coord,
                                   vertex2.S.coord, vertex2.E.coord)
        #self.q = Bio.PDB.calc_dihedral()

    @classmethod
    def get_dihedral(cls, a,b,c,d):
        u = np.cross(a-b, b-c)
        v = np.cross(b-c, c-d)
        return np.rad2deg(np.arccos(np.dot(u/np.linalg.norm(u),
                                           v/np.linalg.norm(v))))

    def __repr__(self):
        return "%s, %s, %s, %s" % (self.S1S2, self.S1E2, self.S2E1, self.q)

    def reverse(self):
        self.residues.reverse()
        self.S1E2, self.S2E1 = self.S2E1, self.S1E2
        # self.q stays the same !


class Edge:
    def __init__(self, vertex1, vertex2):
        self.v1 = vertex1
        self.v2 = vertex2
        self.w = Weight(vertex1, vertex2)

    def __getitem__(self, index):
        if index == 0:
            return self.v1
        elif index == 1:
            return self.v2
        else:
            raise IndexError("index out of range")

    def __repr__(self):
        return "%s <--> %s [%s]" %(self.v1, self.v2, self.w)
    
    def reverse(self):
        self.v1, self.v2 = self.v2, self.v1
        self.w.reverse()
    
    @property
    def tabrepr(self):
        """ Vertex1 x1 y1 z1 Vertex2 x2 y2 z2 """
        return "\t".join([self.v1.shortname,
                          str(self.v1.x), str(self.v1.y), str(self.v1.z),
                          self.v2.shortname,
                          str(self.v2.x), str(self.v2.y), str(self.v2.z)]) + "\n"
    
        
class Vertex:
    def __init__(self, S, E):
        self.S = S
        self.E = E
        self.coord = S.coord
        self.x = S.coord[0]
        self.y = S.coord[1]
        self.z = S.coord[2]
        self.resname = S.parent.resname
        self.id = S.parent.id[1]

        
    def __getitem__(self, index):
        if index == 0:
            return self.S
        elif index == 1:
            return self.E
        elif index == 2:
            return self.id
        elif index == 3:
            return self.resname
        else:
            raise IndexError("index out of range")
    
    def __repr__(self):
        #return "<Vertex %s %s>" %(self.resname, self.id)
        return "%s %s" %(self.resname, self.id)
    
    @property
    def shortname(self):
        return "%s%s" % (Bio.PDB.protein_letters_3to1[self.resname], self.id)

        
def debug_vertex_type(edges_dict):
    for e in edges_dict.keys():
        if not isinstance(e,Vertex):
            print "--->", e
            break
        for ep in edges_dict[e] :
            if not isinstance(ep,Vertex):
                print "===>", ep
                exit()

                        
#### Part 1 : construction of graph
# 2CPK.pdb
# 3LCK.pdb
####


def load_PDB(pdb_fname):
    parser = Bio.PDB.PDBParser()
    structure = parser.get_structure(pdb_fname.replace(".pdb",""),pdb_fname)
    return structure


def select_aa_acessible(pdb_fname, chain_name=None, Wacc=1):

    dssp_fname = pdb_fname.replace(".pdb",".dssp")

    if commands.getstatusoutput("mkdssp -i %s -o %s"%(pdb_fname,dssp_fname))[0]:
        raise RuntimeError("mkdssp failed !")

    dssp_dict, key_list = Bio.PDB.make_dssp_dict(dssp_fname)

    if chain_name :
        #return (dssp_dict,[ key for key in key_list if chain_name in key and \
        #dssp_dict[key][2] >= Wacc ])
        return [k[1][1] for k in key_list if chain_name == k[0] and dssp_dict[k][2] >= Wacc ]
    else :
        #return (dssp_dict,[ key for key in key_list if dssp_dict[key][2] >= Wacc ])
        return [k[1][1] for k in key_list if dssp_dict[k][2] >= Wacc ]
    


def select_vertices(structure, accessible_aa = None):

    vertices_list = []

    if accessible_aa is None:
        is_accessible = lambda residu: True
    else:
        is_accessible = lambda residu: residu.id[1] in accessible_aa

    for residu in structure[0].get_residues() :

        if Bio.PDB.is_aa(residu.resname) and is_accessible(residu):
            atom_pair_dict = {}

            for a in residu :
                if a.name in ['CA', 'CB', 'N'] :
                    atom_pair_dict[a.name] = a 
            
            if atom_pair_dict.get("CA") != None and atom_pair_dict.get("CB") != None:
                vertices_list.append(Vertex(atom_pair_dict["CA"], atom_pair_dict["CB"]))

            elif atom_pair_dict.get("CA")!= None and atom_pair_dict.get("N")!= None:
                vertices_list.append(Vertex(atom_pair_dict["CA"], atom_pair_dict["N"]))

            else :
                print atom_pair_dict
                print residu.id[1]
                raise ValueError("Missing value in vertex.")
    return vertices_list
    

#def build_mtx_dist(vertices_list):
#    N = len(vertices_list)
#    mtx = np.zeros((N,N))
#    for i in range(N):
#        for j in range(i+1, N):
#            a = vertices_list[i][0]
#            b = vertices_list[j][0]
#            dist = np.linalg.norm(a-b)
#            mtx[i,j] = dist
#            mtx[j,i] = dist
#    return mtx


def build_graph(vertices_list, Dmax=8): 
    graph = []
    N = len(vertices_list)
    #mtx = build_mtx_dist(vertices_list)
    for i in range(N):
        for j in range(i+1, N):
            a = vertices_list[i].coord
            b = vertices_list[j].coord
            dist = np.linalg.norm(a-b)
            if dist <= Dmax:
                graph.append(Edge(vertices_list[i], vertices_list[j]))
    return graph


#### Part 2 : comparison of graph

# from: http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
_blosum62_aa =  ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K",
                 "M", "F", "P", "S", "T", "W", "Y", "V", "B", "Z", "X", "*" ]

    
_blosum62 = np.matrix(
""" 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4 ; 
   -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4 ;
   -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4 ;
   -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4 ;
    0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4 ;
   -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4 ;
   -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 ;
    0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4 ;
   -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4 ;
   -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4 ;
   -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4 ;
   -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4 ;
   -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4 ;
   -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4 ;
   -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4 ;
    1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4 ;
    0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4 ;
   -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4 ;
   -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4 ;
    0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4 ;
   -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4 ;
   -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4 ;
    0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4 ;
   -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1 """)


_compatible_aa = _blosum62 >= 2

def aa_identical(aa1, aa2, compatible_aa = _compatible_aa):
    """Using 3 letters amino acid names.
    compatible_aa: matrix of similar amino-acids. ex: BLOSUM62 >= 2"""
    a1 = Bio.PDB.protein_letters_3to1[aa1]
    a2 = Bio.PDB.protein_letters_3to1[aa2]
    return compatible_aa[_blosum62_aa.index(a1), _blosum62_aa.index(a2)]


#def align_edges(edge1, edge2, compatible_aa):
def edges_aa_identical(edge1, edge2):
    """compatible_aa: matrix of similar amino-acids. ex: BLOSUM62 >= 2"""
    # get amino-acids of each vertex in edge
    
    if aa_identical(edge1.v1[3], edge2.v1[3]) and \
       aa_identical(edge1.v2[3], edge2.v2[3]):
        return True
    if aa_identical(edge1.v1[3], edge2.v2[3]) and \
       aa_identical(edge1.v2[3], edge2.v1[3]):
        edge2.reverse()
        return True
    return False


def weight_identical(edge1, edge2, DSS=0.6, DSE=0.75, Dq=35,reverse=False):
    """edge1 and edge2 must have been aligned using residu_identical"""
    if reverse :
        edge2.reverse()
    deltaSS = abs(edge1.w.S1S2 - edge2.w.S1S2)
    #deltaSE = [ edge1.w.S1E2 - edge2.w.S1E2,
    #            edge1.w.S2E1 - edge2.w.S2E1,
    #            edge1.w.S1E2 - edge2.w.S2E1,
    #            edge1.w.S2E1 - edge2.w.S1E2 ]
    deltaSE = np.array([abs(edge1.w.S1E2 - edge2.w.S1E2),
                        abs(edge1.w.S2E1 - edge2.w.S2E1)])
    deltaq = abs(edge1.w.q - edge2.w.q)

    if deltaSS <= DSS and all(deltaSE <= DSE) and deltaq <= Dq:
        return True

    return False


def edges_identical(edge1,edge2) :
    if aa_identical(edge1.v1[3], edge1.v2[3]) and aa_identical(edge2.v1[3], edge2.v2[3]):
        if edges_aa_identical(edge1,edge2) and (weight_identical(edge1,edge2) or weight_identical(edge1, edge2, reverse=True)):
            return True
        else:
            if edges_aa_identical(edge1,edge2) and weight_identical(edge1,edge2):
                return True
    return False


def edges_comp(graph1, graph2):
    common_edges_g1 = {}
    common_edges_g2 = {}
    
    def add_edges(common_edges,edge1,edge2):
        if common_edges.has_key(edge1):
            common_edges[edge1].append(edge2)
        else :
            common_edges[edge1] = [edge2]
            
    for edge1 in graph1:
        for edge2 in graph2:
            if edges_identical(edge1,edge2) : 
                add_edges(common_edges_g1,edge1,edge2)
                add_edges(common_edges_g2,edge2,edge1)
                    
    return (common_edges_g1, common_edges_g2)


def graphlist2dict(graph):
    edges_dict = {}
    for edge in graph:
        if not edges_dict.has_key(edge.v1):
            edges_dict[edge.v1] = set((edge.v2,))
        else:
            edges_dict[edge.v1].add(edge.v2)
        if not edges_dict.has_key(edge.v2):
            edges_dict[edge.v2] = set((edge.v1,))
        else:
            edges_dict[edge.v2].add(edge.v1)
    return edges_dict


def recursive_rm_node(edges_dict, node):
    if len(edges_dict[node]) >= 2:
        return
    if len(edges_dict[node]) == 0:
        edges_dict.pop(node)
        return
    nodep = list(edges_dict[node])[0]
    edges_dict.pop(node)
    edges_dict[nodep].remove(node)
    recursive_rm_node(edges_dict, nodep)


def rm_monogamous(edges_dict):
    nodes = edges_dict.keys()
    for n in nodes:
        if edges_dict.has_key(n) and len(edges_dict[n]) < 2:
            recursive_rm_node(edges_dict,n)


def build_node_edges_dict(common_edges_list): # transfer common edges keys to dict{node:[edges]}
    node_edges_dict = {}
    def add_edge_in_dict(node_edges_dict,node,edge):
        if node_edges_dict.has_key(node) :
            node_edges_dict[node].append(edge)
        else :
            node_edges_dict[node] = [ edge ]
    for edge in common_edges_list :
        add_edge_in_dict(node_edges_dict,edge.v1,edge)
        add_edge_in_dict(node_edges_dict,edge.v2,edge)
    return node_edges_dict


def common_edges_dict_filter(common_edges_dict, node_list): # remove common edges not in a triangle
    new_common_edges_dict = {}
    for e1,e2s in common_edges_dict.iteritems() :
        if e1.v1 in node_list and e1.v2 in node_list :
            new_common_edges_dict[e1] = [ edge for edge in e2s if (edge.v1 in node_list) and (edge.v2 in node_list)  ]
    return new_common_edges_dict
        

def rm_monogamous_from_common_edges(common_edges_dict):
    # filter nodes no clique >= 3
    edges_dict = graphlist2dict(common_edges_dict.keys())
    rm_monogamous(edges_dict)
    # debug_vertex_type(edges_dict)            
    nodes = list(set(edges_dict.keys()))
    return common_edges_dict_filter(common_edges_dict, nodes)


def bfs(cg1,cg2,ce1,ce2,e1,e2,sg1,sg2):
    """
    common graphs 1 {node:[edges]},
    common graphs 2 {node:[edges]},
    common edges dict 1 {e1:[e2s]},
    common edges dict 2 {e2:[e1s]},
    edge 1, edge 2,
    subgraph 1, subgraph 2
    """
    if e1 in sg1 or e2 in sg2 :
        return
    if e2 not in ce1[e1] or e1 not in ce2[e2] :
        return
    if not edges_identical(e1, e2):
        return
    
    sg1.add(e1)
    sg2.add(e2)
    
    for edge1 in cg1[e1.v1] :
        for edge2 in cg2[e2.v1] :
            bfs(cg1,cg2,ce1,ce2,e1,e2,sg1,sg2)
        
    for edge1 in cg1[e1.v2] :
        for edge2 in cg2[e2.v2] :
            bfs(cg1,cg2,ce1,ce2,edge1,edge2,sg1,sg2)
        

def find_subgraphs(common_edges_d1,common_edges_d2,cg1,cg2,sg_min_num = 3) :
    subgraphs1 = []
    subgraphs2 = []
    for e1,e2s in common_edges_d1.iteritems():
        for e2 in e2s :
            print e1
            sg1 = []
            sg2 = []
            print "bfs", e1, e2
            (subg1,subg2) = bfs(cg1,cg2,common_edges_d1,common_edges_d2,e1,e2,sg1,sg2)
            if len(subg1) == len(subg2) and len(subg1) >= sg_min_num :
                subgraphs1.append(sg1)
                subgraphs2.append(sg2)
    return (subgraphs1,subgraphs2)


def build_subgraphs(graph1, graph2):
    print "compare edges"
    (common_edges_g1, common_edges_g2) = edges_comp(graph1, graph2)
    print len(common_edges_g1.keys()),len(common_edges_g2.keys())
    print "rm monogamous from common edges g1"
    common_edges_d1 = rm_monogamous_from_common_edges(common_edges_g1)
    print common_edges_d1
    print "rm monogamous from common edges g2"
    common_edges_d2 = rm_monogamous_from_common_edges(common_edges_g2)
    return ([common_edges_d1.keys()],[common_edges_d2.keys()])
    print "search subgraphs"
    cg1 = build_node_edges_dict(common_edges_g1.keys()) 
    cg2 = build_node_edges_dict(common_edges_g2.keys()) 
    return find_subgraphs(common_edges_d1,common_edges_d2,cg1,cg2)
 
def write_graph(graph, filename, append=False):
    """save graph in text file"""
    mode = 'a' if append else 'w'
    with open(filename, mode) as OUT:
        OUT.write("Node1\tx1\ty1\tz1\tNode2\tx2\ty2\tz2\n")
        for edge in graph:
            OUT.write(edge.tabrepr)
            

def build_graph_sub(pdb_fname, chain = None): 
    structure = load_PDB(pdb_fname)
    accessible_aa = select_aa_acessible(pdb_fname, chain)
    vertices = select_vertices(structure, accessible_aa)
    graph = build_graph(vertices)
    return graph


if __name__ == '__main__':
    print "load graph from pdb"
    graph1 = build_graph_sub('2CPK.pdb', chain='E')
    graph2 = build_graph_sub('3LCK.pdb', chain='A')
    
    print "build subgraphs"
    (sgs1,sgs2) = build_subgraphs(graph1, graph2)

    print "write graph"
    for i, sg in enumerate(sgs1, start=1):
        write_graph(sg, "2CPK_graphs_%02i.tsv" % i)
    for i, sg in enumerate(sgs2, start=1):
        write_graph(sg, "3LCK_graphs_%02i.tsv" % i)

