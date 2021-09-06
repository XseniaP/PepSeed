from tkinter import *
from tkinter import filedialog
from tkinter import scrolledtext
import tkinter as tk
import os
import subprocess
import shlex
from subprocess import check_output
import pathlib
import pandas as pd
import sys
import pathlib

# define a global value as minus infinity
MINUSINF = -1000000000


class Node:
    def __init__(self, name, weight, children, parents):
        self.name = name
        self.weight = weight
        self.children = children
        self.parents = parents

    def add_parent(self, new_parent):
        # self.parents.append(new_parent)
        self.parents.add(new_parent)

    def add_child(self, new_child):
        # self.children.append(new_child)
        self.children.add(new_child)

    def change_weight(self, new_weight):
        self.weight = new_weight


def convert_abc():
    seq_set = []
    arguments = sys.argv
    filepath = ''
    for i in range(len(arguments)):
        if arguments[i] == '-I':
            filepath = arguments[i + 1]
    with open(str(pathlib.Path.cwd()) + "/" + filepath) as fp:
        previous_line = ''
        for line in fp:
            if previous_line.startswith('>'):
                seq_set.append(line.split('\n')[0].split(' ')[0])
            previous_line = line
    fp.close()
    dict_abc = {'R': 'B', 'K': 'B', 'E': 'J', 'D': 'J', 'S': 'O', 'T': 'O', 'L': 'U', 'V': 'U', 'I': 'U', 'Q': 'X',
                'N': 'X', 'W': 'Z', 'F': 'Z', 'A': 'A', 'C': 'C', 'G': 'G', 'H': 'H', 'M': 'M', 'P': 'P', 'Y': 'Y'}
    new_set = []
    for seq in seq_set:
        word = str('')
        for letter in seq:
            word = word + str(dict_abc[letter])
        new_set.append(word)
    return new_set


# extract sssps out of the file sspMatrix.txt created by mapitope
# extract sssps out of the file allPairs.txt created by mapitope
def get_ssps():
    # abc = ['A', 'B', 'X', 'J', 'C', 'H', 'U', 'M', 'Z', 'P', 'O', 'Y', 'G']
    # # pairs = set()
    # pairs = []
    # with open(str(pathlib.Path.cwd()) + "/Results/sspMatrix.txt") as fp:
    #     i = 0
    #     for line in fp:
    #         temp = line.split(' ')
    #         j = 0
    #         for element in temp:
    #             if element == '1':
    #                 # pairs.add(abc[i - 1] + abc[j - 1])
    #                 # pairs.add(abc[j - 1] + abc[i - 1])
    #                 pairs.append(abc[i - 1] + abc[j - 1])
    #                 pairs.append(abc[j - 1] + abc[i - 1])
    #             j += 1
    #         i += 1
    # return pairs
    abc = ['A', 'B', 'X', 'J', 'C', 'H', 'U', 'M', 'Z', 'P', 'O', 'Y', 'G']
    pairs = []
    with open(str(pathlib.Path.cwd()) + "/Results_Mapi/allPairs.txt") as fp:
        for line in fp:
            if line.startswith('Significant: '):
                temp = line.split()
                if len(temp[1]) == 2:
                    pairs.append(temp[1])
                    pairs.append(temp[1][1] + temp[1][0])
    return pairs


# check if the sequence might be in reversed order
def check_reverse(seq, graph):
    rev = False
    count = 0
    rev_count = 0
    for i in range(len(seq) - 1):
        if str(seq[i:i + 2]) in graph:
            if graph[str(seq[i:i + 2])] != 0:
                count += 1
        if str(seq[i + 1: i + 2] + seq[i: i + 1]) in graph:
            if graph[str(seq[i + 1: i + 2] + seq[i: i + 1])] != 0:
                rev_count += 1
    if rev_count > count:
        rev = True
    return rev


# create a seed graph and return the graph and the median weight
def seed_graph_create():
    sequences = convert_abc()
    ssps = sorted(set(get_ssps()))
    graph = dict.fromkeys(ssps, 0)
    for seq in sequences:
        rev = check_reverse(seq, graph)
        prev = ''
        for i in range(len(seq) - 1):
            if not rev:
                temp = str(seq[i:i + 2])
            else:
                # temp = str(seq[i + 1 : i + 2] + seq[i : i + 1])
                # temp = str(seq[-2:])
                # need to check this formula at the end , some issue with i=0
                temp = str(seq[-(i + 1):-i]) + str(seq[-(i + 2):-(i + 1)])
            if temp in ssps:
                if graph[temp] == 0:
                    graph[temp] = Node(name=str(seq[i:i + 2]), weight=1, children=set(), parents=set())
                else:
                    graph[temp].change_weight(graph[temp].weight + 1)
                if prev != '':
                    graph[prev].add_child(temp)
                    graph[temp].add_parent(prev)
                prev = temp
    # need to delete all empty nodes
    # elements = set()
    # for element in graph:
    #     if graph[element] == 0:
    #         elements.add(element)
    # graph.pop(elements)
    count = 0
    sum = 0
    weights = []
    s_set = []
    for element in graph:
        if graph[element] != 0:
            count += 1
            sum += graph[element].weight
            weights.append(graph[element].weight)
    mean = sum / count
    weights.sort()
    median = weights[round(count / 2)]
    for element in graph:
        if graph[element] != 0:
            if graph[element].weight >= mean:
                s_set.append(graph[element].name)
    return graph, median, s_set


# extract graph to csv
def extract_to_csv(graph):
    df = pd.DataFrame(columns={'pair', 'weight', 'parents', 'children'})
    for element in graph:
        if graph[element] != 0:
            df2 = pd.DataFrame([[graph[element].name, graph[element].weight, graph[element].parents, graph[element].children]], columns={'pair', 'weight', 'parents', 'children'})
            df = df.append(df2, ignore_index = True )
    # df = pd.DataFrame(data=graph, index=[0])
    # file_name = str(pathlib.Path.cwd()) + "/graph.xlsx"
    file_name = "graph.csv"
    # df.to_csv(file_name, index = True)


# calculate the weight recursively
def recursive_weight(i, v, S, graph):
    parent_from_set = False
    final_seed =""
    max_v = 0
    # is the condition S==none , s is empty?
    # this is the regular init for each node from the set to be the first node and its weight to be the first weight added
    if (i >= 0) & (S == []):
        return graph[v].weight, v
    # the case when there are still nodes left in the set but current node has no parents to proceed with the search
    if (i >= 0) & (len(S) != 0) & (len(graph[v].parents) == 0):
        # return MINUSINF, v
        return graph[v].weight, v
    # need to change here that also case when the only parent is node itself should use this case
    if (i >= 0) & (len(S) != 0) & (len(graph[v].parents) == 1) & (v in graph[v].parents):
        # return MINUSINF, v
        return graph[v].weight, v
    # more than i nodes not from the SET were used in the path
    if i < 0:
        return MINUSINF, v

    for parent in graph[v].parents:
        S_new = []
        for element in S:
            S_new.append(element)
        if parent in S_new:
            S_new.remove(parent)
            temp, seed = recursive_weight(i, parent, S_new, graph)
        else:
            temp, seed = recursive_weight(i-1, parent, S_new, graph)
        temp += graph[v].weight
        if temp > max_v:
            max_v = temp
            final_seed = seed+v

    # for parent in graph[v].parents:
    #     if parent in S:
    #         parent_from_set = True
    #         break
    # if not parent_from_set:
    #     for parent in graph[v].parents:
    #         S_new = []
    #         for element in S:
    #             S_new.append(element)
    #         # temp, seed = recursive_weight(i - 1, parent, S, graph) + graph[v].weight
    #         temp, seed = recursive_weight(i - 1, parent, S_new, graph)
    #         temp += graph[v].weight
    #         if temp > max_v:
    #             max_v = temp
    # else:
    #     for parent in graph[v].parents:
    #         if parent in S:
    #             S_new = []
    #             # temp  = recursive_weight(i, parent, S.remove(parent), graph) + graph[v].weight
    #             # temp, seed = recursive_weight(i, parent, S.remove(parent), graph) + graph[v].weight
    #             for element in S:
    #                 S_new.append(element)
    #             S_new.remove(parent)
    #             temp, seed = recursive_weight(i, parent, S_new, graph)
    #             temp += graph[v].weight
    #             if temp > max_v:
    #                 max_v = temp
    #                 seed = seed+v
    return max_v, final_seed


# DP algorithm - uses the help function which calculates the weight recursively
def seed_search(graph, s_set):
    # recurrence
    for i in range(len(graph)):
        max = 0
        for v in s_set:
            S = []
            # deep copy of s_set
            for element in s_set:
                S.append(element)
            S.remove(v)
            temp, seed = recursive_weight(i, v, S, graph)
            # temp += graph[v].weight
            if temp > max:
                max = temp
                final_seed = seed
        if max != 0:
            return final_seed, i
        else:
            print("\n no seed found")

