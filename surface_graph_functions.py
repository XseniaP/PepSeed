from surface_graph import Surface_Amino_Acid
import pandas as pd

"""""""""
global dicts:
surface_dict - represent the surface, key looks like V11 - the name and index of the AA
path_dict - same format for the AA that part of the path
"""""""""
surface_dict = {}
path_dict = {}


def Initialize_graph(residue_txt):
    """
    This function go over the file with the amino acids and create Surface_Amino_Acid object for each:
    with the index, name(the letter that represent the AA) , empty list for the neighbors , and the
    bool part_of_the_path == False
    :param residue_txt: input tet file
    :return: Initialize the graph
    """
    residue = open(residue_txt, "r")
    for x in residue:
        y = x.split(":")
        #print(y)
        name = y[0]
        #print(name[0])
        #print(name[1:])
        surface_dict[name] = Surface_Amino_Acid(name[0] , name[1:])


def Initialize_Neighbors_list(pairs_distance_txt , distance_param):
    """
    The function go over the pairs_distance_txt and update for each AA in each pair the neighbor list
    if the distance is smaller then distance_param
    :return: Initialize Neighbors list for each AA
    """
    pairs_distance = open(pairs_distance_txt, "r")
    for x in pairs_distance:
        y = x.split()
        if y[2] != 'Distance' and float(y[2]) < distance_param:
            first = y[0]
            second = y[1]
            first = first.split(":")[0]
            second = second.split(":")[0]
            #print(first , second)
            first = surface_dict[first]
            second = surface_dict[second]
            first.add_new_neighbor(second)
            second.add_new_neighbor(first)


"""""""""
important to check if might be another format to this file!!!!
"""""""""

def choose_path_from_pepsurf(significant_path_txt):
    AA_groups_dict = {'R':'B','K':'B','E':'J','D':'J','S':'O','T':'O','L':'U','V':'U','I':'U',
                      'Q':'X','N':'X','W':'Z','F':'Z','A':'A','C':'C','G':'G','H':'H','M':'M',
                      'P':'P','Y':'Y','-':'-'}
    significant_path = open(significant_path_txt, "r")
    i=0
    index = 100                   #arbitrary number
    list_path=[]
    list_alignment_input=[]
    list_alignment_output=[]
    score_list=[]
    for x in significant_path:
        if x.startswith("path:"):
            path = x
            path = path.split(":")
            path=path[1].rstrip()
            path = path.split()
            list_path.append(path)
        if x.startswith("Alignment:"):
            index = i
        if i == index+1 :                ## the input of the alignment
            alignment_input = x
            list_alignment_input.append(list(alignment_input.rstrip().strip(" ")))
        if i == index+2 :                   ## the output of the alignment
            alignment_output = x
            list_alignment_output.append(list(alignment_output.rstrip().strip(" ")))
            index = 0
        i+=1
    i=0
    j=0
    for seq in list_alignment_input:                           ## 20 letters to 13 letters for the input
        for AA in seq:
            list_alignment_input[i][j] = AA_groups_dict[AA]
            j += 1
        i+=1
        j=0
    i=0
    j=0
    for seq in list_alignment_output:                          ## 20 letters to 13 letters for the output
        for AA in seq:
            list_alignment_output[i][j] = AA_groups_dict[AA]
            j += 1
        i+=1
        j=0
    # print(list_alignment_input)
    # print(list_alignment_output)
    # print(list_path)
    for index1 in range(len(list_alignment_input)):          ## calculate the score of each
        #sum=0
        for index2 in range(len(list_alignment_input[index1])):
            if list_alignment_output[index1][index2] == '-':
                list_path[index1].insert(index2,'')
            #if list_alignment_input[index1][index2] == list_alignment_output[index1][index2]:
                #sum+=1
                #print(sum)
            if list_alignment_input[index1][index2] != list_alignment_output[index1][index2]:
                #sum-=1
                #print(sum)
                list_path[index1][index2] = ""
    # print(list_alignment_input)
    # print(list_alignment_output)
    # print(list_path)
    for index in range(len(list_path)):
       list_path[index] = [x[:-1] for x in list_path[index]]
    return (list_path, list_alignment_input, list_alignment_output)
            #print("cc " ,list_path[index1][index2])
        #score_list.append(sum)
    #print(list_alignment_input)
    #print(list_alignment_output)
    #print(list_path)
    #score_list[0] = 4
    # print(score_list)
    # best_indexes = [i for i, x in enumerate(score_list) if x == max(score_list)]
    # print("indexes of the max score: " , best_indexes)
    # new_list_path=[]
    # new_list_input=[]
    # new_list_output=[]
    # for index in best_indexes:
    #     list_path[index] = [x[:-1] for x in list_path[index]]
    #     #list_path[index] =  [x[:-1] for x in list_path[index] if x != ""]
    #     #list_path[index] = []
    #     new_list_input.append(list_alignment_input[index])
    #     new_list_output.append(list_alignment_output[index])
    #     new_list_path.append(list_path[index])
    # print(new_list_path)
    # print(new_list_input)
    # print(new_list_output)
    #return(new_list_path , new_list_input , new_list_output)


def mark_path(path):
    for k in path:
        if k != "":
            #print(k)
            surface_dict[k].insert_to_path()


def update_path():
    for key, value in surface_dict.items():
        if value.check_if_in_path() == True:
          if key not in path_dict:
              path_dict[key] = value
    #print(path_dict)


def extention2(graph_csv ,seed_graph_indexes_dict):
    AA_groups_dict = {'B':['R','K'] , 'J':['E','D'] , 'O':['S','T'] , 'U':['L','V','I'] ,
                      'X': ['Q','N'] , 'Z':['W','F'] , 'A':['A'] , 'C':['C'] , 'G':['G'] ,
                      'H':['H'] , 'M':['M'] , 'P':['P'] , 'Y':['Y']}
    count_dict ={}
    df = pd.read_csv(graph_csv)
    df = df[(df.pair != 'mean') & (df.pair != 'median') & (df.pair.isnull()==False)]
    for index, pairs_dict in seed_graph_indexes_dict.items():
        #print(index)
        for pair, first_or_second in pairs_dict.items():
            all_neighbors_to_look=[]
            neighbors_to_look = df[df["pair"] == pair]
            neighbors_to_look = neighbors_to_look['parents']
            for i in neighbors_to_look:
                neighbors_to_look = i
                for ch in ['{', '}', '"', " ", "'"]:
                    if ch in neighbors_to_look:
                        neighbors_to_look = neighbors_to_look.replace(ch, "")
                neighbors_to_look = neighbors_to_look.split(',')
                if neighbors_to_look != ['set()']:
                    for AA in neighbors_to_look:
                        all_neighbors_to_look.append(AA)
            neighbors_to_look = df[df["pair"] == pair]
            neighbors_to_look = neighbors_to_look['children']
            for i in neighbors_to_look:
                neighbors_to_look = i
                for ch in ['{', '}', '"', " ", "'"]:
                    if ch in neighbors_to_look:
                        neighbors_to_look = neighbors_to_look.replace(ch, "")
                neighbors_to_look = neighbors_to_look.split(',')
                if neighbors_to_look != ['set()']:
                    for AA in neighbors_to_look:
                        all_neighbors_to_look.append(AA)
            #print(all_neighbors_to_look)
            neighbors_to_look_new =[]
            for amino_group in all_neighbors_to_look:
                AA1 = amino_group[0]
                AA2 = amino_group[1]
                for amino_acid1 in AA_groups_dict[AA1]:
                    for amino_acid2 in AA_groups_dict[AA2]:
                        if amino_acid1+amino_acid2 not in neighbors_to_look_new:
                            neighbors_to_look_new.append(amino_acid1+amino_acid2)
                        if amino_acid2+amino_acid1 not in neighbors_to_look_new:
                            neighbors_to_look_new.append(amino_acid2 + amino_acid1)
            #print(neighbors_to_look_new)
            index_neighbors = surface_dict[index].get_neighbors_list()
            #print("index_neighbors :", index_neighbors)
            pair_extention = []
            for pairs in neighbors_to_look_new:
             #   print(pairs)
                first = pairs[0]
                second = pairs[1]
                # if first == 'C':
                #     continue
                for AA in index_neighbors:
                    if first in AA:
                        for AA2 in surface_dict[AA].get_neighbors_list():
                            if second in AA2:
                                #print(AA,AA2)
                                if AA not in pair_extention:
                                    pair_extention.append(AA)
                                if AA2 not in pair_extention:
                                    pair_extention.append(AA2)
            #print("pair_extention ", pair_extention)

            for residue in pair_extention:
                if residue in count_dict.keys():
                    count_dict[residue] += 1
                if residue not in count_dict.keys():
                    count_dict[residue] = 1
    #print("count_dict: ",count_dict)
    return (count_dict)






def extention(graph_csv ,seed_graph_indexes_dict):
    AA_groups_dict = {'B':['R','K'] , 'J':['E','D'] , 'O':['S','T'] , 'U':['L','V','I'] ,
                      'X': ['Q','N'] , 'Z':['W','F'] , 'A':['A'] , 'C':['C'] , 'G':['G'] ,
                      'H':['H'] , 'M':['M'] , 'P':['P'] , 'Y':['Y']}
    count_dict ={}
    df = pd.read_csv(graph_csv)
    #del df['Unnamed: 0']
    df = df[(df.pair != 'mean') & (df.pair != 'median') & (df.pair.isnull()==False)]
    for index, pairs_dict in seed_graph_indexes_dict.items():
        #print(index)
        for pair, first_or_second in pairs_dict.items():

            if first_or_second == 0:
                #print(index  ," parents")
                neighbors_to_look_new =[]
                neighbors_to_look = df[df["pair"] == pair]
                neighbors_to_look = neighbors_to_look['parents']
                for i in neighbors_to_look:
                    neighbors_to_look = i
                    for ch in ['{', '}', '"', " ", "'"]:
                        if ch in neighbors_to_look:
                            neighbors_to_look = neighbors_to_look.replace(ch, "")
                    neighbors_to_look = neighbors_to_look.split(',')
                    if neighbors_to_look == ['set()']:
                        neighbors_to_look = []
                    neighbors_to_look_new=[]
                   # print("before convert to 20 AA: ",neighbors_to_look)
                    for amino_group in neighbors_to_look:
                        first = amino_group[0]
                        second = amino_group[1]
                        for amino_acid1 in AA_groups_dict[first]:
                            for amino_acid2 in AA_groups_dict[second]:
                                neighbors_to_look_new.append(amino_acid1+amino_acid2)
                  #  print("after convert to 20 AA: ",neighbors_to_look_new)

                index_neighbors = surface_dict[index].get_neighbors_list()
                #print("index_neighbors :", index_neighbors)
                pair_extention = []
                for pairs in neighbors_to_look_new:
                 #   print(pairs)
                    first = pairs[0]
                    second = pairs[1]
                    # if first == 'C':
                    #     continue
                    for AA in index_neighbors:
                        if second in AA:
                            for AA2 in surface_dict[AA].get_neighbors_list():
                                if first in AA2:
                                    # print(AA,AA2)
                                    if AA not in pair_extention:
                                        pair_extention.append(AA)
                                    if AA2 not in pair_extention:
                                        pair_extention.append(AA2)
                #print("pair_extention ", pair_extention)

            if first_or_second == 1:
                #print(index, " children")
                neighbors_to_look = df[df["pair"] == pair]
                neighbors_to_look = neighbors_to_look['children']
                for i in neighbors_to_look:
                    neighbors_to_look = i
                    for ch in ['{', '}', '"'," " ,"'"]:
                        if ch in neighbors_to_look:
                            neighbors_to_look = neighbors_to_look.replace(ch,"")
                    neighbors_to_look = neighbors_to_look.split(',')
                    #print(neighbors_to_look)
                    if neighbors_to_look == ['set()']:
                        neighbors_to_look = []

                    neighbors_to_look_new=[]
                   # print("before convert to 20 AA: ",neighbors_to_look)
                    for amino_group in neighbors_to_look:
                        first = amino_group[0]
                        second = amino_group[1]
                        for amino_acid1 in AA_groups_dict[first]:
                            for amino_acid2 in AA_groups_dict[second]:
                                neighbors_to_look_new.append(amino_acid1+amino_acid2)
                  #  print("after convert to 20 AA: " ,neighbors_to_look_new)


                index_neighbors = surface_dict[index].get_neighbors_list()
                #print("index_neighbors :",index_neighbors)
                pair_extention = []
                for pairs in neighbors_to_look_new:
                 #   print(pairs)
                    first = pairs[0]
                    second = pairs[1]
                    # if second == 'C':
                    #     continue
                    for AA in index_neighbors:
                        if first in AA:
                            for AA2 in surface_dict[AA].get_neighbors_list():
                                if second in AA2:
                                    #print(AA,AA2)
                                    if AA not in pair_extention:
                                        pair_extention.append(AA)
                                    if AA2 not in pair_extention:
                                        pair_extention.append(AA2)
                #print("pair_extention " ,pair_extention)

        #### to see where it is better to count!!!!!!!!!!

        #for residue in pair_extention:
         #   if residue in count_dict.keys():
          #      count_dict[residue]+=1
           # if residue not in count_dict.keys():
            #    count_dict[residue]=1
        #print(count_dict)

            for residue in pair_extention:
                if residue in count_dict.keys():
                    count_dict[residue]+=1
                if residue not in count_dict.keys():
                    count_dict[residue]=1
    #print(count_dict)
    return(count_dict)


def get_paths_and_return_best_epitope(paths_list ,seed_graph_indexes_dict ,graph_csv , extention_param):
    list_scores=[]
    list_dict=[]
    epitope_list=[]
    for i in range(len(seed_graph_indexes_dict)):
        #print("this is the ", i + 1, "path")
        epitope = []
        score = 0
        epitope_dict = extention2(graph_csv, seed_graph_indexes_dict[i])
        for AA in paths_list[i]:
            if AA != '':
                epitope.append(AA)
        for key,value in epitope_dict.items():
            if epitope_dict[key] > extention_param:
                #print("key =",key)
                score+=1
                #print("score =",score)
                if key not in epitope:
                    epitope.append(key)
            if epitope_dict[key] == extention_param:
                score+=0.5
                #print("key =",key)
                #print("score =",score)
                if key not in epitope:
                    epitope.append(key)
        #print("path score :" ,score)
        list_scores.append(score)
        list_dict.append(epitope_dict)
        epitope_list.append(epitope)
        #print("epitope: " ,epitope)
        #print("  ")
    #print(list_scores)
    s = [i[0] for i in sorted(enumerate(list_scores), key=lambda x: x[1])]
    #print(s)
    print("\nthe results sorted from the best: ")
    for index in reversed(s):
        print("the path is: ", paths_list[index])
        print("path score :" ,list_scores[index])
        print("dict is: " , list_dict[index])
        print("epitope: " ,epitope_list[index])
        print("  ")

