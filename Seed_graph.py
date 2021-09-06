import math
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
# MINUSINF = -1000000000
MINUSINF = -math.inf

class Node:
    def __init__(self, name, weight, children, parents):
        self.name = name
        self.weight = weight
        self.children = children
        self.parents = parents

    def add_parent(self, new_parent):
        self.parents.add(new_parent)

    def add_child(self, new_child):
        self.children.add(new_child)

    def change_weight(self, new_weight):
        self.weight = new_weight

def truncate_carbons(seq, end, start):
    dic = {}
    dic['end'] = '0'
    dic['start'] = '0'
    if (len(seq) == end) and (seq[end-1] == 'C'):
        dic['end'] = '1'
    if (len(seq) == start) and (seq[0] == 'C'):
        dic['start'] = '1'
    if (dic['end'] == '1') and (dic['start'] == '1'):
        seq = seq[1:len(seq)-1]
    elif (dic['end'] == '0') and (dic['start'] == '1'):
        seq = seq[1:len(seq)]
    elif (dic['end'] == '1') and (dic['start'] == '0'):
        seq = seq[0:len(seq)-1]
    return seq

def convert_abc():
    vocab = {}
    seq_set = []
    arguments = sys.argv
    start = 0
    end = 0
    filepath = ''
    for i in range(len(arguments)):
        if arguments[i] == '-I':
            filepath = arguments[i + 1]
    with open(str(pathlib.Path.cwd()) + "/" + filepath) as fp:
        previous_line = ''
        for line in fp:
            if line.startswith('#PAIR_FIXED C* = '):
                start = int(line.split(' = ')[1])
            if line.startswith('#PAIR_FIXED *C = '):
                end = int(line.split(' = ')[1])
            if previous_line.startswith('>'):
                seq_set.append(line.split('\n')[0].split(' ')[0])
            previous_line = line
    fp.close()
    dict_abc = {'R': 'B', 'K': 'B', 'E': 'J', 'D': 'J', 'S': 'O', 'T': 'O', 'L': 'U', 'V': 'U', 'I': 'U', 'Q': 'X',
                'N': 'X', 'W': 'Z', 'F': 'Z', 'A': 'A', 'C': 'C', 'G': 'G', 'H': 'H', 'M': 'M', 'P': 'P', 'Y': 'Y'}
    new_set = []
    if (start == 0) and (end == 0):
        for seq in seq_set:
            word = str('')
            for letter in seq:
                word = word + str(dict_abc[letter])
            new_set.append(word)
    else:
        for seq in seq_set:
            # seq = truncate_carbons(seq, end, start)
            word = str('')
            for letter in seq:
                word = word + str(dict_abc[letter])
            new_set.append(word)
            vocab[word] = seq
    return new_set, vocab


# extract sssps out of the file sspMatrix.txt created by mapitope
# extract sssps out of the file allPairs.txt created by mapitope
def get_ssps():
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

def find_mean(graph):
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
    return mean


# create a seed graph and return the graph and the median weight
def seed_graph_create():
    sequences = convert_abc()[0]
    ssps = sorted(set(get_ssps()))
    graph = dict.fromkeys(ssps, 0)
    for seq in sequences:
        rev = check_reverse(seq, graph)
        prev = ''
        for i in range(len(seq) - 1):
            if not rev:
                temp = str(seq[i:i + 2])
            else:
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
    # delete all empty nodes
    list = []
    for element in graph:
        if graph[element] == 0:
            list.append(element)
    graph_new = {key: graph[key] for key in graph if key not in list}
    s_set = []
    mean = find_mean(graph_new)
    for element in graph_new:
        if graph_new[element] != 0:
            if graph_new[element].weight >= mean:
                s_set.append(graph_new[element].name)
    return graph_new, mean, s_set


# extract graph to csv
def extract_to_csv(graph):
    df = pd.DataFrame(columns={'pair', 'weight', 'parents', 'children'})
    for element in graph:
        if graph[element] != 0:
            df2 = pd.DataFrame([[graph[element].name, graph[element].weight, graph[element].parents, graph[element].children]], columns={'pair', 'weight', 'parents', 'children'})
            df = df.append(df2)
    file_name = "graph.csv"
    df.to_csv(file_name, index = False)


# calculate the weight recursively
def recursive_weight(i, v, S, graph, nodes):
    final_seed = ""
    max_v = MINUSINF
    penalty = 2 * find_mean(graph)

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

    # recurrence
    parents = graph[v].parents - set([v])
    for parent in parents:
    # for parent in graph[v].parents:
        S_new = []
        temp = 0
        # hard copy of the set
        for element in S:
            S_new.append(element)
        nodes_copy = []
        for element in nodes:
            nodes_copy.append(element)
        if parent in S_new:
            S_new.remove(parent)
            nodes_copy.remove(parent)
            temp, seed = recursive_weight(i, parent, S_new, graph, nodes_copy)
        else:
            if parent in nodes_copy:
                nodes_copy.remove(parent)
                temp, seed = recursive_weight(i-1, parent, S_new, graph, nodes_copy)
                temp = temp - penalty
            else:
                seed = ""
        temp += graph[v].weight
        if temp > max_v:
            max_v = temp
            final_seed = seed + v
    if graph[v].weight > max_v:
        max_v = graph[v].weight
        final_seed = v
    return max_v, final_seed


def redefine_seed(final_seed):
    seed = ""
    for i in range(len(final_seed)-1):
        if (i == 2):
            if final_seed[i] != final_seed[i-1]:
                seed = seed + final_seed[i-2: i+2]
            else:
                seed = seed + final_seed[i - 2: i] + final_seed[i + 1]
        elif (i % 2 == 0) and (i > 2):
            if seed[len(seed)-1] == final_seed[i]:
                seed = seed + final_seed[i+1]
            else:
                seed = seed + final_seed[i : i + 2]
    return seed


def convert_to_20(seed, vocab):
    new_seed = ""
    for i in range(len(seed)-1):
        # the seed is not redefined yet (pair is not merged to the next one by last-first letter match)
        # so each consequent pair is reviewed separately
        if i % 2 == 0:
            my_list = list()
            for seq in vocab:
                for j in range(len(seq)-1):
                    if seq[j : j+2] == seed[i : i+2]:
                        my_list.append(vocab[seq][j: j+2])
            my_dict_1 = {}
            my_dict_2 = {}
            for el in my_list:
                if el[0] in my_dict_1:
                    my_dict_1[el[0]] += 1
                else:
                    my_dict_1[el[0]] = 1
                if el[1] in my_dict_2:
                    my_dict_2[el[1]] += 1
                else:
                    my_dict_2[el[1]] = 1
            max1 = 0
            max2 = 0
            for element in my_dict_1:
                if my_dict_1[element] > max1:
                    max1 = my_dict_1[element]
                    letter = element
            new_seed = new_seed + letter

            for element in my_dict_2:
                if my_dict_2[element] > max2:
                    max2 = my_dict_2[element]
                    letter = element
            new_seed = new_seed + letter

            # my_dict_1 = {i[0]: i.count(i[0]) for i in my_list}
            # my_dict_1 = {i: my_list.count(i) for i in my_list}

    return new_seed

# DP algorithm - uses the help function which calculates the weight recursively
def seed_search(graph, s_set, mean):
    # DP target , no more than 4 insersions allowed, each insertion leads to penalty
    abs_max = 0
    k = input("Enter max number of pairs with the weight below mean value on the path: ")
    for i in range(int(k)+1):
        max = 0
        for v in s_set:
            S = []
            temp = 0
            # deep copy of s_set
            for element in s_set:
                S.append(element)
            S.remove(v)
            nodes = set(graph.keys())
            nodes.remove(v)
            temp, seed = recursive_weight(i, v, S, graph, nodes)
            # temp += graph[v].weight
            # temp = temp - mean * i * 1.7
            if temp > max:
                max = temp
                final_seed = seed
        if max > abs_max:
            abs_max = max
            f_final_seed = final_seed
    if (abs_max != 0) and (len(f_final_seed) >= 6):
        vocab = convert_abc()[1]
        f_final_seed = convert_to_20(f_final_seed, vocab)
        seed = redefine_seed(f_final_seed)
        # return final_seed, i
        return seed
    else:
        print("\n no seed found")