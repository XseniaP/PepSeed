import math
import pandas as pd
import sys
import pathlib
import networkx as nx
import matplotlib.pyplot as plt

MINUSINF = -math.inf

dict_abc = {'R': 'B', 'K': 'B', 'E': 'J', 'D': 'J', 'S': 'O', 'T': 'O', 'L': 'U', 'V': 'U', 'I': 'U', 'Q': 'X',
                'N': 'X', 'W': 'Z', 'F': 'Z', 'A': 'A', 'C': 'C', 'G': 'G', 'H': 'H', 'M': 'M', 'P': 'P', 'Y': 'Y', 'q': 'X'}

# best alignment of query nodes 1...i (seed letters 1...i) that ends in network vertex "v" (surface graph vertex v)
def qpath_recursive_call(seed, i, v, S, nindel, surface_dict, gap_penalty, match_score):
    final_path = list()
    all_max_paths = dict()
    max_v = MINUSINF
    all_max_paths [max_v] = []
    # all_max_paths[max_v] = set()

    # INIT of the DP
    # case when we found the node matching the first letter of the seed
    if (nindel >= 0) & (i == 0) & (dict_abc[seed[i]] == dict_abc[surface_dict[v].AA_name]):
        # final_path.insert(0, v)
        final_path.append(v)
        all_max_paths[match_score] = [final_path]
        # all_max_paths[match_score] = set()
        # all_max_paths[match_score].add(final_path)
        return match_score, final_path, all_max_paths

    # mismatch of the first letter of the seed with zero penalty
    # if (nindel == 0) & (i == 0) & (dict_abc[seed[i]] != dict_abc[surface_dict[v].AA_name]):
    #     # final_path.insert(0, v)
    #     final_path.append("MISM")
    #     # all_max_paths[MINUSINF] = [final_path]
    #     all_max_paths[0].add(final_path)
    #     return 0, final_path, all_max_paths

    # exceeded # of insertions or deletions allowed
    if (nindel < 0):
        # final_path.insert(0, v)
        final_path.append(v)
        all_max_paths[MINUSINF] = [final_path]
        # all_max_paths[MINUSINF] = set()
        # all_max_paths[MINUSINF].add(final_path)
        return MINUSINF, final_path, all_max_paths

    # if first letter was deletion/mismatch
    if (i < 0) & (nindel >= 0):
        # final_path.append("MISM")
        all_max_paths[0] = [final_path]
        # all_max_paths[0] = set()
        # all_max_paths[0].add(final_path)
        return 0, final_path, all_max_paths

    # match (any letter of the seed other than the first one)
    if (dict_abc[surface_dict[v].AA_name] == dict_abc[seed[i]]) & (i > 0) & (nindel >= 0):
        all_max_paths_match = dict()
        fin_path_match = list()
        max_v_match = MINUSINF
        all_max_paths_match[max_v_match] = []
        # all_max_paths_match[max_v_match] = set()
        for neighbor_node in surface_dict[v].Neighbors_list:
            neighbor = neighbor_node.AA_name + neighbor_node.index
            if neighbor in S:
                S_new = []
                temp_score = 0
                # hard copy of the set
                for element in S:
                    S_new.append(element)
                S_new.remove(neighbor)
                temp_score, temp_path, all_paths = qpath_recursive_call(seed, i - 1, neighbor, S_new, nindel, surface_dict, gap_penalty, match_score)
                temp_score = temp_score + match_score
                temp_path.append(v)

                # all_paths[temp_score] = list(item + [v] for item in all_paths[temp_score - match_score])

                if temp_score > max_v_match:
                    max_v_match = temp_score
                    fin_path_match = temp_path
                    all_max_paths_match[max_v_match] = all_paths[max_v_match - match_score]

                    # all_max_paths_match[max_v_match] = all_paths[temp_score]

                    # all_max_paths_match[max_v_match] = [(list(item)).append(v) for item in set(all_paths[max_v_match - match_score])]
                # elif (temp_score == max_v_match) & (temp_score != MINUSINF):
                #     # all_max_paths_match[max_v_match].append(all_paths[max_v_match - match_score])
                #
                #     all_max_paths_match[max_v_match] = all_max_paths_match[max_v_match] + all_paths[max_v_match - match_score]
                #
                #     # all_max_paths_match[max_v_match].append(all_paths[temp_score])
                #
                #     # all_max_paths_match[max_v_match].append(all_paths[max_v_match - match_score][0])
                #     # all_max_paths_match[max_v_match].append([(list(item)).append(v) for item in set(all_paths[max_v_match - match_score])])
                #     # all_max_paths_match[max_v_match].add(all_paths[max_v_match - match_score])
        if max_v_match > max_v:
            max_v = max_v_match
            final_path = fin_path_match
            # all_max_paths = all_max_paths_match
            all_max_paths[max_v_match] = all_max_paths_match[max_v_match]

        elif (max_v_match == max_v) & (max_v_match != MINUSINF):
            # all_max_paths[max_v].append(all_paths[max_v - gap_penalty][0])
            # all_max_paths[max_v].append(all_paths[max_v_match - gap_penalty])

            # all_max_paths[max_v].append(all_max_paths_match[max_v_match])

            all_max_paths[max_v] = all_max_paths[max_v] + all_max_paths_match[max_v_match]


    # insertion (can't be the last node)
    if (i != len(seed) - 1) & (i >= 0):
        all_max_paths_ins = dict()
        fin_path_ins = list()
        max_v_ins = MINUSINF
        all_max_paths_ins[max_v_ins] = []
        # all_max_paths_ins[max_v_ins] = set()
        for neighbor_node in surface_dict[v].Neighbors_list:
            neighbor = neighbor_node.AA_name + neighbor_node.index
            # if (dict_abc[neighbor_node.AA_name] == dict_abc[seed[i]]) & (nindel >= 0) & (neighbor in S):
            if neighbor in S:
                S_new = []
                temp_score = 0
                # hard copy of the set
                for element in S:
                    S_new.append(element)
                S_new.remove(neighbor)
                temp_score, temp_path, all_paths = qpath_recursive_call(seed, i, neighbor, S_new, nindel-1, surface_dict, gap_penalty, match_score)
                temp_score = temp_score + gap_penalty
                temp_path.append(" INS ")

                # all_paths[temp_score] = list(item + [" INS "] for item in all_paths[temp_score - gap_penalty])

                if temp_score > max_v_ins:
                    max_v_ins = temp_score
                    fin_path_ins = temp_path
                    all_max_paths_ins[max_v_ins] = all_paths[max_v_ins - gap_penalty]

                    # all_max_paths_ins[max_v_ins] = all_paths[temp_score]

                # elif (temp_score == max_v_ins) & (temp_score != MINUSINF):
                #     # all_max_paths_ins[max_v_ins].append(all_paths[max_v_ins - gap_penalty])
                #
                #     all_max_paths_ins[max_v_ins] = all_max_paths_ins[max_v_ins] + all_paths[max_v_ins - gap_penalty]
                #
                #     # all_max_paths_ins[max_v_ins].append(all_paths[temp_score])
                #
                #     # all_max_paths_ins[max_v_ins].append(all_paths[max_v_ins - gap_penalty][0])
                #     # all_max_paths_ins[max_v_ins].add(all_paths[max_v_ins - gap_penalty])
        if max_v_ins > max_v:
            max_v = max_v_ins
            final_path = fin_path_ins
            # all_max_paths = all_max_paths_ins
            all_max_paths[max_v_ins] = all_max_paths_ins[max_v_ins]

        elif (max_v_ins == max_v) & (max_v_ins != MINUSINF):
            # all_max_paths[max_v].append(all_paths[max_v - gap_penalty][0])
            # all_max_paths[max_v].append(all_paths[max_v_ins - gap_penalty])

            # all_max_paths[max_v].append(all_max_paths_ins[max_v_ins])

            all_max_paths[max_v] = all_max_paths[max_v] + all_max_paths_ins[max_v_ins]

    # deletion (skip letter of the seed)
    # if (nindel > 0) & (i > 0):
    if i >= 0:
        S_new = []
        temp_score = 0
        # hard copy of the set
        for element in S:
            S_new.append(element)
        temp_score, temp_path, all_paths = qpath_recursive_call(seed, i-1, v, S_new, nindel - 1, surface_dict, gap_penalty, match_score)
        # if (i != len(seed) - 1):
        temp_score = temp_score + gap_penalty
        temp_path.append(" DEL ")

        # all_paths[temp_score] = list(item + [" DEL "] for item in all_paths[temp_score - gap_penalty])

        # temp_path.insert(0, "")
        if temp_score > max_v:
            max_v = temp_score
            final_path = temp_path
            all_max_paths[max_v] = all_paths[max_v - gap_penalty]

            # all_max_paths[max_v] = all_paths[max_v]

        # elif (temp_score == max_v) & (temp_score != MINUSINF):
        #     # all_max_paths[max_v].append(all_paths[max_v - gap_penalty][0])
        #
        #     # all_max_paths[max_v].append(all_paths[max_v - gap_penalty])
        #
        #     all_max_paths[max_v] = all_max_paths[max_v] + all_paths[max_v - gap_penalty]
        #
        #     # all_max_paths[max_v].append(all_paths[max_v])
        #
        #     # all_max_paths[max_v].add(all_paths[max_v - gap_penalty])

    # # mismatch
    # if (i != 0):
    #     fin_path_mism = list()
    #     max_v_mism = MINUSINF
    #     for neighbor_node in surface_dict[v].Neighbors_list:
    #         neighbor = neighbor_node.AA_name + neighbor_node.index
    #         S_new = []
    #         temp_score = 0
    #         # hard copy of the set
    #         for element in S:
    #             S_new.append(element)
    #         temp_score, temp_path = qpath_recursive_call(seed, i - 1, neighbor, S_new, nindel - 1, surface_dict, penalty)
    #         temp_score = temp_score + penalty
    #         temp_path.append(" MISM ")
    #         # temp_path.insert(0, "")
    #         if temp_score > max_v_mism:
    #             max_v_mism = temp_score
    #             fin_path_mism = temp_path
    #     if max_v_mism > max_v:
    #         max_v = max_v_mism
    #         final_path = fin_path_mism

    return max_v, final_path, all_max_paths


# best alignment of query nodes 1...N that ends in network vertex v
# call it with S = all nodes of the surface graph and i to be the last aa of the query/seed
# nindel is number of deletions/insertions allowed
def qpath_target_call(seed, surface_dict, nindel, gap_penalty, match_score):
    # S = surface_dict.items()
    max_score = 0
    all_paths = dict()
    temp_path = list()
    final_path = list()
    all_paths_scores = dict()
    for node in surface_dict:
        S_new = []
        # for element in S:
        for element in surface_dict:
            S_new.append(element)
        S_new.remove(node)
        if node == 'E213':
            print('found_it')
        temp_score, temp_path, temp_all_max = qpath_recursive_call(seed, len(seed)-1, node, S_new, nindel, surface_dict, gap_penalty, match_score)
        # temp_all_max_unique = set(tuple(i) for i in temp_all_max[temp_score])
        # temp_all_max_unique = set()
        # for el in temp_all_max[temp_score]:
        #     temp_all_max_unique.add(tuple(el))
        # fin = ""
        # for el in temp_path:
        #     fin = fin + " " + el
        # all_paths[fin] = temp_score
        if temp_score in all_paths_scores:
            # all_paths_scores[temp_score] = all_paths_scores[temp_score] + temp_all_max_unique
            # all_paths_scores[temp_score].add(temp_all_max_unique)
            # all_paths_scores[temp_score].append(temp_all_max[temp_score])

            all_paths_scores[temp_score] = all_paths_scores[temp_score] + temp_all_max[temp_score]

        elif (temp_score not in all_paths_scores) & (temp_score != MINUSINF):
            # all_paths_scores[temp_score] = temp_all_max_unique
            all_paths_scores[temp_score] = temp_all_max[temp_score]

        if temp_score > max_score:
            max_score = temp_score
            final_path = temp_path
            all_paths[max_score] = temp_all_max[max_score]
        elif (temp_score == max_score) & (temp_score != MINUSINF):
            all_paths[max_score].append(temp_all_max[max_score])
    if max_score != 0:
        # all_paths_scores.pop(MINUSINF)
        # all_paths.pop(MINUSINF)
        # all_paths_scores_final = dict()
        # for score in all_paths_scores:
        #     all_paths_scores_final[score] = set(all_paths_scores[score])
        # return max_score, final_path, all_paths, all_paths_scores_final
        return max_score, final_path, all_paths, all_paths_scores
    else:
        return {}, 0, "no path found"