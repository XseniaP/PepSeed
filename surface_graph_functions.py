from surface_graph import Surface_Amino_Acid

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
    for index1 in range(len(list_alignment_input)):          ## calculate the score of each
        sum=0
        for index2 in range(len(list_alignment_input[index1])):
            if list_alignment_input[index1][index2] == list_alignment_output[index1][index2]:
                sum+=1
            else:
                sum-=1
                list_path[index1][index2] = ""
        score_list.append(sum)
    print(list_alignment_input)
    print(list_alignment_output)
    print(list_path)
    #score_list[0] = 4
    print(score_list)
    best_indexes = [i for i, x in enumerate(score_list) if x == max(score_list)]
    print("indexes of the max score: " , best_indexes)
    new_list_path=[]
    new_list_input=[]
    new_list_output=[]
    for index in best_indexes:
        list_path[index] = [x[:-1] for x in list_path[index]]
        #list_path[index] =  [x[:-1] for x in list_path[index] if x != ""]
        #list_path[index] = []
        new_list_input.append(list_alignment_input[index])
        new_list_output.append(list_alignment_output[index])
        new_list_path.append(list_path[index])
    print(new_list_path)
    print(new_list_input)
    print(new_list_output)
    return(new_list_path , new_list_input , new_list_output)


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





