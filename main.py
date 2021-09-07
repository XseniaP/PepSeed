# This is a Main module
import Mapitope_run
import PepSurf_run
import Seed_graph
import Cluster_extend
import surface_graph_functions
import pathlib

def main_func():
#   Running mapi with the given arguments
#     Mapitope_run.mapi_run()

#   Ksenia's: create a seed graph
    graph, mean, s_set = Seed_graph.seed_graph_create()

#   Extract graph to scv
#     Seed_graph.extract_to_csv(graph)

#   Ksenia's: find a seed
    seed, original_seed = Seed_graph.seed_search(graph, s_set, mean)
    print(seed)

#   Ksenia's: running Pepsurf with the seed found and given arguments to find the preliminary cluster
#     seed = "QFTAAE"
#     seed = "QMQFTAAEMGQYT"
#     seed = "QFTAAE"
#     seed = "QFTAAEEC"
    print(seed)
    PepSurf_run.perpsurf_run(seed)

#   Sapir's: get the highest score paths
    filepath = str(pathlib.Path.cwd()) + "/RESULTS/sigPathsAln/0_significantPaths.txt"
    paths_set, input_alignment_set, output_alignment_set = surface_graph_functions.choose_path_from_pepsurf(filepath)
    print(paths_set)
    print(input_alignment_set)
    print(output_alignment_set)

#   Ksenia's: return set of dictionaries, 1 for each path
    set_of_dicts = Seed_graph.path_to_graph_dictionary(graph, paths_set, input_alignment_set, output_alignment_set, original_seed)

#   Sapir's: functions calls to create a surface graph

#   Sapir's: locating the cluster found by pepsurf on the surface graph

#   Running cluster extension


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main_func()


