# This is a Main module
import Mapitope_run
import PepSurf_run
import Seed_graph
import pathlib
from surface_graph_functions import Initialize_Neighbors_list
from surface_graph_functions import choose_path_from_pepsurf
from surface_graph_functions import Initialize_graph
from surface_graph_functions import get_paths_and_return_best_epitope

#  3 examples parameters (configurations) to run the
# -U Results_Mapi -P 1e6j_P.pdb -S 1e6j_P.txt -I 13b5.txt -C P -D 6 -V 2.5 -R rasmol.txt -F 0
# -U Results_Mapi -P 6ch7_Q.pdb -S 6ch7_Q.txt -I fatsa_BG18_top10motifs.txt -C Q -D 5.0 -V 3.64 -R rasmol.txt -F 0
# -U Results_Mapi -P 2ghw.pdb -S 2ghw.txt -I 17bcopy.txt -C A -D 6.0 -V 2 -R rasmol.txt -F 0

# Default user input (stdin)
# Enter max number of pairs with the weight below mean value on the path: 0
# Enter distance which defines the neighbors: 4
# Enter the extension parameter: 1

def main_func():
#   Rnning Mapitope in order to get the file with the SSPs
    Mapitope_run.mapi_run()

#   Moving to the ABC of 13 letters and Creating a seed graph
    graph, mean, s_set, rev_indices = Seed_graph.seed_graph_create()
#   Displaying an image of the seed graph
    Seed_graph.print_seed_graph(graph, mean)

#   Extract graph (dictionary) to csv file
    Seed_graph.extract_to_csv(graph)

#   Find a seed using a DP algorithm
    seed, original_seed, all_seeds = Seed_graph.seed_search(graph, s_set, mean, rev_indices)

#   Running QPath with the seed as an input to find the preliminary path on the antigen surface
    PepSurf_run.perpsurf_run(seed)

#   Create a surface graph
    residue_txt = str(pathlib.Path.cwd()) + "/Results_Mapi/surface.txt"
    Initialize_graph(residue_txt)

    pairs_distance_txt = str(pathlib.Path.cwd()) + "/Results_Mapi/pairsDistance.txt"


    D_param = int(input("Enter distance which defines the neighbors: "))
    # D_param = 4
    Initialize_Neighbors_list(pairs_distance_txt, D_param)

#   Get the highest score paths
    significant_path_txt = str(pathlib.Path.cwd()) + "/RESULTS/sigPathsAln/0_significantPaths.txt"
    paths_set, input_alignment_set, output_alignment_set = choose_path_from_pepsurf(significant_path_txt)
    print(paths_set)

#   Return set of dictionaries, 1 for each path
    seed_graph_indexes_dict = Seed_graph.path_to_graph_dictionary(graph, paths_set, input_alignment_set,
                                                                  output_alignment_set, original_seed, seed)

#   Locating the cluster found by PepSurf on the surface graph and Running cluster extension
    graph_csv = str(pathlib.Path.cwd()) + "/graph.csv"
    extension_param = int(input("Enter the extension parameter: "))
    # extention_param = 1
    get_paths_and_return_best_epitope(paths_set, seed_graph_indexes_dict, graph_csv, extension_param)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main_func()


