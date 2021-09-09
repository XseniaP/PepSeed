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

    def set_name(self, name):
        self.name = name

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
                'N': 'X', 'W': 'Z', 'F': 'Z', 'A': 'A', 'C': 'C', 'G': 'G', 'H': 'H', 'M': 'M', 'P': 'P', 'Y': 'Y', 'q': 'X'}
    new_set = []
    if (start == 0) and (end == 0):
        for seq in seq_set:
            word = str('')
            for letter in seq:
                word = word + str(dict_abc[letter])
            new_set.append(word)
            vocab[word] = seq
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
    mean = (sum / count)
    weights.sort()
    median = weights[round(count * 0.5)]
    return mean


# create a seed graph and return the graph and the median weight
def seed_graph_create():
    reverse_indices = []
    sequences = convert_abc()[0]
    ssps = sorted(set(get_ssps()))
    graph = dict.fromkeys(ssps, 0)
    for seq in sequences:
        rev = check_reverse(seq, graph)
        reverse_indices.append(rev)
        prev = ''
        for i in range(len(seq) - 1):
            if not rev:
                temp = str(seq[i:i + 2])
            else:
                temp = str(seq[-(i + 1):-i]) + str(seq[-(i + 2):-(i + 1)])
            if temp in ssps:
                if graph[temp] == 0:
                    graph[temp] = Node(name=temp, weight=1, children=set(), parents=set())
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
    return graph_new, mean, s_set, reverse_indices


# extract graph to csv
def extract_to_csv(graph):
    df = pd.DataFrame(columns=['pair', 'weight', 'parents', 'children'])
    for element in graph:
        if graph[element] != 0:
            df2 = pd.DataFrame([[graph[element].name, graph[element].weight, graph[element].parents, graph[element].children]], columns=['pair', 'weight', 'parents', 'children'])
            df = df.append(df2)
    file_name = "graph.csv"
    df.to_csv(file_name, index = False)


# calculate the weight recursively
def recursive_weight(i, v, S, graph, nodes):
    final_seed = ""
    max_v = MINUSINF
    penalty = 1.5 * find_mean(graph)

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
    if v == 'q':
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


def convert_to_20(seed, vocab, rev_indices):
    new_seed = ""
    for i in range(len(seed)-1):
        # the seed is not redefined yet (pair is not merged to the next one by last-first letter match)
        # so each consequent pair is reviewed separately
        if i % 2 == 0:
            my_list = list()
            k = 0
            for seq in vocab:
                rev = rev_indices[k]
                for j in range(len(seq)-1):
                    if rev == False:
                        if seq[j: j+2] == seed[i: i+2]:
                            my_list.append(vocab[seq][j: j+2])
                    else:
                        if seq[j + 1] + seq[j] == seed[i: i+2]:
                            my_list.append(vocab[seq][j + 1] + vocab[seq][j])
                k += 1
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

            letter = ""
            for element in my_dict_1:
                if my_dict_1[element] > max1:
                    max1 = my_dict_1[element]
                    letter = element
            new_seed = new_seed + letter

            letter = ""
            for element in my_dict_2:
                if my_dict_2[element] > max2:
                    max2 = my_dict_2[element]
                    letter = element
            new_seed = new_seed + letter

            # my_dict_1 = {i[0]: i.count(i[0]) for i in my_list}
            # my_dict_1 = {i: my_list.count(i) for i in my_list}

    return new_seed


# DP algorithm - uses the help function which calculates the weight recursively
def seed_search(graph, s_set, mean, rev_indices):
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
        print(f_final_seed)
        original_seed = f_final_seed
        f_final_seed = convert_to_20(f_final_seed, vocab, rev_indices)
        print(f_final_seed)
        seed = redefine_seed(f_final_seed)
        print(f_final_seed)
        # return final_seed, i
        return seed, original_seed
    else:
        print("\n no seed found")
        return "", ""


def add_cut_spaces(paths_set, input_alignment_set, output_alignment_set, k, cut_at_start, seq_cut):
    # paths_set_copy = list()
    # input_alignment_set_copy = list()
    # output_alignment_set_copy = list()

    paths_set_k = list()
    input_alignment_set_k = list()
    output_alignment_set_k = list()

    i = 0
    while cut_at_start != 0:
        paths_set_k.append('')
        input_alignment_set_k.append(seq_cut[i])
        output_alignment_set_k.append('')
        cut_at_start = cut_at_start - 1
        i += 1

    for i in range(len(paths_set[k])):
        paths_set_k.append(paths_set[k][i])
        input_alignment_set_k.append(input_alignment_set[k][i])
        output_alignment_set_k.append(output_alignment_set[k][i])

    #
    # for i in range(len(paths_set)):
    #     if i!= k:
    #         paths_set_copy.append(paths_set[i])
    #         input_alignment_set_copy.append((input_alignment_set[i]))
    #         output_alignment_set_copy.append(output_alignment_set[i])
    #     else:
    #         paths_set_copy.append(paths_set_k)
    #         input_alignment_set_copy.append(input_alignment_set_k)
    #         output_alignment_set_copy.append(output_alignment_set_k)
    return paths_set_k, input_alignment_set_k, output_alignment_set_k


# getting the path or set of paths and return the relevant dictionary for each
# which suggests the relevant nodes in the graph
def path_to_graph_dictionary(graph, paths_set, input_alignment_set, output_alignment_set, original_seed, seed):
    dict_abc = {'R': 'B', 'K': 'B', 'E': 'J', 'D': 'J', 'S': 'O', 'T': 'O', 'L': 'U', 'V': 'U', 'I': 'U', 'Q': 'X',
                'N': 'X', 'W': 'Z', 'F': 'Z', 'A': 'A', 'C': 'C', 'G': 'G', 'H': 'H', 'M': 'M', 'P': 'P', 'Y': 'Y', 'q': 'X'}
    my_set = list()
    k = 0
    for path in paths_set:
        cut_at_start = 0
        seq_cut =''
        for i in range(len(seed)):
            if input_alignment_set[k][0] == dict_abc[seed[i]]:
                break
            else:
                seq_cut = seq_cut + dict_abc[seed[i]]
                cut_at_start += 1
        path, input_alignment, output_alignment = add_cut_spaces(paths_set, input_alignment_set, output_alignment_set, k, cut_at_start, seq_cut)
        i = 0
        j = 0
        dic = dict((el, {}) for el in path if el != '')

        while (j < len(input_alignment)) and (i < len(original_seed) - 1):
            # if it's a list position so only 1 letter can be compared
            if j == len(input_alignment) - 1:
                if path[j] == '':
                    i += 2
                    j += 1
                    continue
                if input_alignment[j] == original_seed[i]:
                    dic[path[j]][original_seed[i: i + 2]] = 0
                if input_alignment[j] == original_seed[i + 1]:
                    dic[path[j]][original_seed[i: i + 2]] = 1
                i += 2
                j += 1
            else:
                # if the pair starting in this position is identical to the pair in the original seed
                if input_alignment[j] + input_alignment[j + 1] == original_seed[i: i + 2]:
                    if (path[j] != '') and (path[j + 1] != ''):
                        dic[path[j]][original_seed[i: i + 2]] = 0
                        dic[path[j + 1]][original_seed[i: i + 2]] = 1
                        i += 2
                        j += 1
                        continue
                    elif (path[j] != '') and (path[j + 1] == ''):
                        dic[path[j]][original_seed[i: i + 2]] = 0
                        i += 2
                        j += 1
                        continue
                    elif (path[j] == '') and (path[j + 1] != ''):
                        dic[path[j + 1]][original_seed[i: i + 2]] = 1
                        i += 2
                        j += 1
                        continue
                    elif (path[j] == '') and (path[j + 1] == ''):
                        i += 2
                        j += 1
                        continue
                # if the pairs are not identical
                else:
                    j += 1
        my_set.append(dic)
        k += 1
    print(my_set)
    return my_set
