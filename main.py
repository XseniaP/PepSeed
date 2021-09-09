# This is a Main module
import Mapitope_run
import PepSurf_run
import Seed_graph
import Cluster_extend
import surface_graph
import surface_graph_functions

def main_func():
#   Running mapi with the given arguments
    Mapitope_run.mapi_run()

#   Ksenia's: create a seed graph
    graph, median, s_set = Seed_graph.seed_graph_create()

#   Ksenia's: find a seed
    seed = Seed_graph.seed_search(graph, s_set)

#   Ksenia's: running Pepsurf with the seed found and given arguments to find the preliminary cluster
    PepSurf_run.perpsurf_run(seed)

    #   Sapir's: functions calls to create a surface graph
    residue_txt = "surface.txt"
    Initialize_graph(residue_txt)

    pairs_distance_txt = "pairsDistance.txt"
    D_param = 10
    Initialize_Neighbors_list(pairs_distance_txt , D_param)

    significant_path_txt = '0_significantPaths.txt'

    pepsurf_result = choose_path_from_pepsurf(significant_path_txt)

    ## ksenia, this is the input for the dict creation:
    ## 0 - path
    ## 1 - input
    ## 2 - output
    pepsurf_result[0]
    pepsurf_result[1]
    pepsurf_result[2]

    graph_csv = "our_graph.csv"

    seed_graph_indexes_dict = "here put the function that create the dict"

    #   Sapir's: locating the cluster found by pepsurf on the surface graph

    #   Running cluster extension
    extention_param = 3
    get_paths_and_return_best_epitope(pepsurf_result[0] ,seed_graph_indexes_dict ,graph_csv ,extention_param)



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main_func()


