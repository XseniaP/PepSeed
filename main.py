# This is a Main module
import Mapitope_run
import PepSurf_run
import Seed_graph
import Cluster_extend
import surface_graph_functions
import pathlib
from surface_graph_functions import Initialize_Neighbors_list
from surface_graph_functions import choose_path_from_pepsurf
from surface_graph_functions import Initialize_graph
from surface_graph_functions import get_paths_and_return_best_epitope


# -U TEST4 -P 2ghw.pdb -S 2ghw.txt -I 17b.txt -C A -D 7.0 -V 5.0 -R rasmol.txt
# -U TEST6 -P 2ghw.pdb -S 2ghw.txt -I 17b.txt -C A -D 9.0 -V 3.0 -R rasmol.txt -F 3
#  -U TEST4 -P 2ghw.pdb -S 2ghw.txt -I 17b.txt -C A -D 7.0 -V 5.0 -R rasmol.txt
# -U Results_Mapi -P 2ghw.pdb -S 2ghw.txt -I 17b.txt -C A -D 9.0 -V 3.0 -R rasmol.txt -F 3

# -U Results_Mapi -P 2ghw.pdb -S 2ghw.txt -I 17b.txt -C A -D 9.0 -V 2 -R rasmol.txt -F 0
# -U Results_Mapi -P 1e6j_P.pdb -S 1e6j_P.txt -I 13b5.txt -C P -D 17.0 -V 2.5 -R rasmol.txt -F 0

# -U Results_Mapi -P 1pc3_model1_A.pdb -S 1pc3_model1_A.txt -I fasta_all_motifs_p4_36.txt -C A -D 6.0 -V 11.4 -R rasmol.txt -F 0

# -U Results_Mapi -P 6ch7_Q.pdb -S 6ch7_Q.txt -I fatsa_BG18_top10motifs.txt -C Q -D 8.0 -V 3.64 -R rasmol.txt -F 0

def main_func():
#   Running mapi with the given arguments
    Mapitope_run.mapi_run()

#   Ksenia's: create a seed graph
    graph, mean, s_set, rev_indices = Seed_graph.seed_graph_create()

#   Extract graph to scv
    Seed_graph.extract_to_csv(graph)

#   Ksenia's: find a seed
    seed, original_seed = Seed_graph.seed_search(graph, s_set, mean, rev_indices)
    print(seed)

#   Ksenia's: running Pepsurf with the seed found and given arguments to find the preliminary cluster
    # seed = 'LRTQRNRP'
    # seed = 'TNRLRN'
    print(seed)
    PepSurf_run.perpsurf_run(seed)


    #   Sapir's: functions calls to create a surface graph
    residue_txt = str(pathlib.Path.cwd()) + "/Results_Mapi/surface.txt"
    # residue_txt = "surface.txt"
    Initialize_graph(residue_txt)


    pairs_distance_txt = str(pathlib.Path.cwd()) + "/Results_Mapi/pairsDistance.txt"
    # pairs_distance_txt = "pairsDistance.txt"
    D_param = 8
    Initialize_Neighbors_list(pairs_distance_txt, D_param)


#   Sapir's: get the highest score paths
    significant_path_txt = str(pathlib.Path.cwd()) + "/RESULTS/sigPathsAln/0_significantPaths.txt"
    # significant_path_txt = '0_significantPaths.txt'
    paths_set, input_alignment_set, output_alignment_set = choose_path_from_pepsurf(significant_path_txt)
    print(paths_set)
    print(input_alignment_set)
    print(output_alignment_set)

    # paths_set = [['N183', 'W184', 'T188', '', 'A194', '', '', 'C198'], ['N74', '', '', 'A77', 'A78', 'D81', 'E79'], ['N120', '', 'T48', '', 'A47', 'D51', ''], ['Q219', '', 'T216', 'A217', '', 'D197', '', 'C198'], ['N21', '', 'T19', 'A22', '', 'E29', 'E28'], ['N193', '', 'T171', 'A174', '', 'D166', 'D163'], ['', 'T210', 'A208', 'A209', 'E213', 'E212'], ['', '', 'T186', '', 'A177', '', 'E180']]
    # input_alignment_set = [['X', 'Z', 'O', 'A', 'A', 'J', 'J', 'C'], ['X', 'Z', 'O', 'A', 'A', 'J', 'J'], ['X', 'Z', 'O', 'A', 'A', 'J', 'J'], ['X', 'Z', 'O', 'A', 'A', 'J', 'J', 'C'], ['X', 'Z', 'O', 'A', 'A', 'J', 'J'], ['X', 'Z', 'O', 'A', 'A', 'J', 'J'], ['Z', 'O', 'A', 'A', 'J', 'J'], ['X', 'Z', 'O', 'A', 'A', 'J', 'J']]
    # output_alignment_set = [['X', 'Z', 'O', 'U', 'A', 'X', 'X', 'C'], ['X', 'U', '-', 'A', 'A', 'J', 'J'], ['X', 'M', 'O', 'G', 'A', 'J', 'X'], ['X', 'M', 'O', 'A', 'U', 'J', 'B', 'C'], ['X', 'U', 'O', 'A', 'U', 'J', 'J'], ['X', 'U', 'O', 'A', 'B', 'J', 'J'], ['M', 'O', 'A', 'A', 'J', 'J'], ['J', 'Y', 'O', 'B', 'A', 'X', 'J']]

#   Ksenia's: return set of dictionaries, 1 for each path
    seed_graph_indexes_dict = Seed_graph.path_to_graph_dictionary(graph, paths_set, input_alignment_set, output_alignment_set, original_seed, seed)

#   Sapir's: locating the cluster found by pepsurf on the surface graph
#   Running cluster extension
    graph_csv = str(pathlib.Path.cwd()) + "/graph.csv"
    extention_param = 4
    get_paths_and_return_best_epitope(paths_set, seed_graph_indexes_dict, graph_csv, extention_param)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main_func()


