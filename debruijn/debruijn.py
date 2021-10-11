#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import statistics
import random
import math
from random import randint

import pickle
import networkx as nx
import matplotlib.pyplot as plt

random.seed(9001)

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path

def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', dest='fastq_file', type=isfile,
                        required=True, help="Fastq file")
    parser.add_argument('-k', dest='kmer_size', type=int,
                        default=22, help="K-mer size (default 22)")
    parser.add_argument('-o', dest='output_file', type=str,
                        default=os.curdir + os.sep + "contigs.fasta",
                        help="Output contigs in fasta file")
    parser.add_argument('-f', dest='graphimg_file', type=str,
                        help="Save graph as image (png)")
    return parser.parse_args()

def read_fastq(fastq_file):
    """Read fasta file
    ---
    Read a fastq file to import sequences\n
    Argument:
    - fastq_file (argparse) : fastq file make in arguments of program\n
    Return:
    - Generator of sequnces in fastq
    """
    with open(fastq_file, 'r', encoding ='utf-8') as fastq:
        lines = fastq.readlines()
        for i in range(len(lines)):
            if lines[i].startswith("@"):
                yield lines[i+1][:-1]

def cut_kmer(read, kmer_size):
    """Cut sequence in k-mers
    ---
    Take a sequences and cut it into k-mers a specified lenght\n
    Arguments:
    - read (str): sequence of letters
    - kmer_size (int): size of k-mers wanted (give in argparse)\n
    Return:
    - Generator of k-mers
    """
    nb_kmers = len(read) - kmer_size + 1
    for i in range(nb_kmers):
        kmer = read[i:i + kmer_size]
        yield kmer

def build_kmer_dict(fastq_file, kmer_size):
    """Build dictionnary of k-mers
    ---
    Arguments:
    - fastq_file (argparse): fastq file choosen
    - kmer_size (argparse): lenght of k-mers choosen\n
    Return:
    - kmer_dict (dict): dictionnary of k-mers in the fastq file
    """
    kmers_dict = {}
    for seq in read_fastq(fastq_file):
        for kmer in cut_kmer(seq, kmer_size):
            if kmer not in kmers_dict.keys():
                kmers_dict[kmer]=seq.find(kmer)
            else:
                break
    return kmers_dict

def build_graph(kmer_dict):
    """Produce a directed graph
    ---
    Argument:
    - kmer_dict (dict): dictonnary of k-mers\n
    Return:
    - graph (nx.DiGraph): directed graph for the k-mers
    """
    graph = nx.DiGraph()
    for kmer in kmer_dict.keys():
        graph.add_edge(kmer[:-1],kmer[1:],weight=kmer_dict[kmer])
    return graph

def remove_paths(graph, path_list, delete_entry_node, delete_sink_node):
    """Remove paths in a graph
    ---
    Arguments:
    - graph (nx.Graph): De Bruijn's Graph by Network X
    - path_list (list): list of path between two nodes in graph
    - delete_entry_node (bool): indicate if the entry node need to be deleted
    - delete_sink_node (bool): indicate if the last node need to be deleted\n
    Return:
    - graph (nx.Graph): graph simplified of low informative nodes
    """
    if delete_entry_node is True and delete_sink_node is True:
        graph.remove_nodes_from([node for path in path_list for node in path])
    elif delete_entry_node is True and delete_sink_node is False:
        graph.remove_nodes_from([node for path in path_list for node in path[:-1]])
    elif delete_entry_node is False and delete_sink_node is True:
        graph.remove_nodes_from([node for path in path_list for node in path[1:]])
    else:
        graph.remove_nodes_from([node for path in path_list for node in path[1:-1]])
    return graph

def std(data):
    """Calcul of the standard deviation in a sample
    ---
    Argument:
    - data (list): list of nodes\n
    Return:
    - Standard deviation (int): standard deviation of the data
    """
    mean = sum(data)/len(data)
    ech_mean = sum((l-mean)**2 for l in data) / (len(data)-1)
    return math.sqrt(ech_mean)


def select_best_path(graph, path_list, path_length, weight_avg_list,
                     delete_entry_node=False, delete_sink_node=False):
    """Select the best path in list of path
    ---
    Arguments:
    - graph (nx.DiGraph): De Bruijn's Graph by Network X
    - path_list (list): list of path in the graph
    - path_lenght (list): list of the lenght of paths in path_list
    - weight_avg_list (list): average weight of paths in path_list
    - delete_entry_node (bool): indicate if the entry node need to be deleted
    - delete_sink_node (bool): indicate if the last node need to be deleted\n
    Return:
    - graph (nx.DiGraph): simplified graph
    """
    if std(weight_avg_list) > 0:
        index_keep = weight_avg_list.index(max(weight_avg_list))
    elif std(weight_avg_list) ==0 and  std(path_length) > 0:
        index_keep = path_length.index(max(path_length))
    else:
        index_keep = randint(0,len(path_list))
    del path_list[index_keep]
    remove_paths(graph, path_list, delete_entry_node, delete_sink_node)
    return graph

def path_average_weight(graph, path):
    """Calcul of the mean of weights in a path
    ---
    Arguments:
    - graph (nx.DiGraph): De Bruijn's Graph by Network X
    - path (list): list of nodes in a path\n
    Return:
    - mean weight of the path (int)
    """
    return statistics.mean([d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)])

def solve_bubble(graph, ancestor_node, descendant_node):
    """Resolve the bubble in a graph
    ---
    Select the path with the higher weight or the higher lenght in a bubble
    and remove all others paths\n
    Arguments:
    - graph (nx.DiGraph): De Bruijn's Graph by Network X
    - ancestor_node (node): upstream node in the graph
    - descendent_node (node): downstream node in the graph\n
    Return:
    - graph (DiGraph): De Bruijn Graph without bubbles
    """
    path_list = []
    path_lenght = []
    weight_avg_list = []
    for path in nx.all_simple_paths(graph, ancestor_node,descendant_node):
        path_list.append(path)
        path_lenght.append(len(path))
        weight_avg_list.append(path_average_weight(graph, path))
    select_best_path(graph, path_list, path_lenght, weight_avg_list)
    return graph

def simplify_bubbles(graph):
    """ Simplify bubbles
    ---
    Simplfy bubble in a graph to have the more informative path\n
    Argument:
    - graph (nx.DiGraph): De Bruijn's Graph by Network X\n
    Return:
    - graph (nx.DiGraph): simplified De Bruijn's graph without bubbles
    """
    bubble = False
    for node in graph.nodes:
        if node in graph.nodes():
            list_nodes = list(graph.predecessors(node))
            if len(list_nodes) >1:
                for i, npre in enumerate(list_nodes):
                    for j in range(i+1,len(list_nodes)):
                        node_ancestor = nx.lowest_common_ancestor(graph, npre, list_nodes[j])
                        if node_ancestor is not None:
                            bubble = True
                            break
    if bubble:
        graph = simplify_bubbles(solve_bubble(graph, node_ancestor, node))
    return graph

def solve_entry_tips(graph, starting_nodes):
    """Resolve the entry in a graph
    ---
    Select the entry node of the graph with higher weight or the higher
    lenght of path and remove others.\n
    Arguments:
    - graph (nx.DiGraph): De Bruijn Graph by Network X
    - starting_nodes (list): list of possible entry nodes\n
    Return:
    - graph (nx.DiGraph): De Bruijn's graph with the one entry node
    """
    path_list = []
    path_lenght = []
    weight_avg_path = []
    for nodes in graph.nodes:
        nodes_pred = list(graph.predecessors(nodes))
        if len(nodes_pred) > 1:
            for entries in starting_nodes:
                path_list += list(nx.all_simple_paths(graph, entries, nodes))
    if len(path_list) > 0:
        for path in path_list:
            path_lenght.append(len(path))
            weight_avg_path.append(path_average_weight(graph, path))
        graph = select_best_path(graph, path_list, path_lenght, weight_avg_path,
        delete_entry_node=True)
    return graph

def solve_out_tips(graph, ending_nodes):
    """Resolve the output in a graph
    ---
    Select the output node of the graph with higher weight or
    the higher lenght of path and remove others.\n
    Arguments:
    - graph (nx.DiGraph): De Bruijn Graph by Network X
    - ending_nodes (list): list of possible outputs nodes\n
    Return:
    - graph (nx.DiGraph): De Bruijn's graph with the one output node
    """
    path_list = []
    path_lenght = []
    weight_avg_path = []
    for nodes in graph.nodes:
        successors = list(graph.successors(nodes))
        if len(successors) > 1:
            for outputs in ending_nodes:
                path_list += list(nx.all_simple_paths(graph, nodes, outputs))
    if len(path_list) > 0:
        for path in path_list:
            path_lenght.append(len(path))
            weight_avg_path.append(path_average_weight(graph, path))
        graph = select_best_path(graph, path_list, path_lenght,
        weight_avg_path, delete_sink_node=True)
    return graph

def get_starting_nodes(graph):
    """Get the starting node
    ---
    Select all starting nodes in a graph\n
    Argument:
    - graph (nx.DiGraph): De Bruijn's graph\n
    Return:
    - starting_nodes (list): list of starting nodes in the graph
    """
    starting_nodes = []
    for nodes in graph:
        if len(list(graph.predecessors(nodes))) == 0:
            starting_nodes.append(nodes)
    return starting_nodes

def get_sink_nodes(graph):
    """Get the ending node
    ---
    Select all ending nodes in a graph\n
    Argument:
    - graph (nx.DiGraph): De Bruijn's graph\n
    Return:
    - sink_nodes (list): list of ending nodes in the graph
    """
    sink_nodes = []
    for nodes in graph:
        if len(list(graph.successors(nodes))) == 0:
            sink_nodes.append(nodes)
    return sink_nodes

def concat_path(path_list):
    """Concantenate list of k-mers
    ---
    Argument:
    - path_list (list): list of k-mers\n
    Return:
    - path_str (str): string of the path with k-mers
    """
    path_str=path_list[0]
    for i in path_list[1:]:
        path_str += i[-1]
    return path_str

def get_contigs(graph, starting_nodes, ending_nodes):
    """Build contigs between two nodes
    ---
    Arguments:
    - graph (nx.DiGraph): De Bruijn's graph
    - starting_nodes (list): list of starting nodes in a graph
    - ending_nodes (list): list of ending nodes in a graph\n
    Return:
    - contigs (list): list of sequences (str) possible in the graph
    """
    contigs = []
    for start_node in starting_nodes:
        for end_node in ending_nodes:
            if nx.has_path(graph, start_node,end_node):
                for simple_path in nx.all_simple_paths(graph,start_node,end_node):
                    path = concat_path(simple_path)
                    contigs.append((path,len(path)))
    return contigs

def save_contigs(contigs_list, output_file):
    """Save contigs in a fasta file
    ---
    Arguments:
    - contigs_list (list): list of contigs in a graph
    - output_file (str): name of the fasta file return\n
    Return:
    - fasta (fasta file): write fasta file for contigs
    """
    with open(f"{output_file}.fasta", "w", encoding='utf-8') as fasta:
        # extract number of contig and contig
        for index,contig in enumerate(contigs_list):
            # write a fasta file
            fasta.write(">contig_" + str(index) + " len=" + str(contig[1]) +
            "\n" + fill(contig[0]) + "\n")
    return fasta

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def draw_graph(graph, graphimg_file):
    """Draw the graph
    """
    fig, ax = plt.subplots()
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] > 3]

    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d['weight'] <= 3]
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(graph, pos, edgelist=esmall, width=6, alpha=0.5,
                           edge_color='b', style='dashed')
    plt.savefig(graphimg_file)


def save_graph(graph, graph_file):
    """Save the graph with pickle
    """
    with open(graph_file, "wt") as save:
        pickle.dump(graph, save)


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Build k-mers dictionnary and the graph with it
    print(f"Get sequences and dictionnary of k-mers in ... {args.fastq_file}")
    kmer_dictionnary = build_kmer_dict(args.fastq_file, args.kmer_size)
    print("Build De Bruijn's graph ...")
    graph = build_graph(kmer_dictionnary)
    # Simplify the graph with remove of bubles
    graph = simplify_bubbles(graph)
    # Determine the starting node
    starting_nodes = get_starting_nodes(graph)
    #Determine the ending node
    sink_nodes = get_sink_nodes(graph)
    # Clear the graph of non informative entries
    graph = solve_entry_tips(graph, starting_nodes)
    # Clear the graph of non informative outputs
    graph = solve_out_tips(graph, sink_nodes)
    #Get contigs from the graph and save it in fasta file
    print("Get contigs and build the fasta file ...")
    contigs = get_contigs(graph, starting_nodes, sink_nodes)
    save_contigs(contigs, args.output_file)
    # Fonctions de dessin du graphe
    # A decommenter si vous souhaitez visualiser un petit
    # graphe
    # Plot the graph
    if args.graphimg_file:
        draw_graph(graph, args.graphimg_file)
    # Save the graph in file
    if args.graph_file:
        save_graph(graph, args.graph_file)

if __name__ == '__main__':
    main()
