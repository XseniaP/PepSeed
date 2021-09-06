# This is a Main module
import Mapitope_run
import PepSurf_run
import Seed_graph
import Cluster_extend

def main_func():
#   Running mapi with the given arguments
    Mapitope_run.mapi_run()

#   Ksenia's: create a seed graph
    graph, mean, s_set = Seed_graph.seed_graph_create()
#   Extract graph to scv
#     Seed_graph.extract_to_csv(graph)
#   Ksenia's: find a seed
    seed = Seed_graph.seed_search(graph, s_set, mean)
    print(seed)
#   Ksenia's: running Pepsurf with the seed found and given arguments to find the preliminary cluster
#     seed = "QFTAAE"
#     seed = "QMQFTAAEMGQYT"
#     seed = "QFTAAE"
    # seed = "NFTAAD"
    print(seed)
    PepSurf_run.perpsurf_run(seed)

#   Sapir's: functions calls to create a surface graph

#   Sapir's: locating the cluster found by pepsurf on the surface graph

#   Running cluster extension


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main_func()


