# This is a Main module
import Mapitope_run
import PepSurf_run
import Seed_graph
import Cluster_extend
import surface_graph_functions
import pathlib

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
#     seed = "QFTAAE"
#     seed = "QMQFTAAEMGQYT"
#     seed = "QFTAAE"
#     seed = "QFTAAEEC"
#     seed = 'NSPRVIR'
#     seed = 'MKLPRTNQ'
#     seed = 'RLVRPTQ'
#     seed = 'FYINAL'
#     seed = 'FYIPAEI'
#     seed = 'LRLRNRPTQ'
    # seed = 'LRNRPTQ'
    # seed = 'LRNTQRP'
    # seed = 'LRTQRNRP'
    # seed = 'TNRLRN'
    print(seed)
    PepSurf_run.perpsurf_run(seed)

#   Sapir's: get the highest score paths
    filepath = str(pathlib.Path.cwd()) + "/RESULTS/sigPathsAln/0_significantPaths.txt"
    paths_set, input_alignment_set, output_alignment_set = surface_graph_functions.choose_path_from_pepsurf(filepath)
    print(paths_set)
    print(input_alignment_set)
    print(output_alignment_set)

    # paths_set = [['N183', 'W184', 'T188', '', 'A194', '', '', 'C198'], ['N74', '', '', 'A77', 'A78', 'D81', 'E79'], ['N120', '', 'T48', '', 'A47', 'D51', ''], ['Q219', '', 'T216', 'A217', '', 'D197', '', 'C198'], ['N21', '', 'T19', 'A22', '', 'E29', 'E28'], ['N193', '', 'T171', 'A174', '', 'D166', 'D163'], ['', 'T210', 'A208', 'A209', 'E213', 'E212'], ['', '', 'T186', '', 'A177', '', 'E180']]
    # input_alignment_set = [['X', 'Z', 'O', 'A', 'A', 'J', 'J', 'C'], ['X', 'Z', 'O', 'A', 'A', 'J', 'J'], ['X', 'Z', 'O', 'A', 'A', 'J', 'J'], ['X', 'Z', 'O', 'A', 'A', 'J', 'J', 'C'], ['X', 'Z', 'O', 'A', 'A', 'J', 'J'], ['X', 'Z', 'O', 'A', 'A', 'J', 'J'], ['Z', 'O', 'A', 'A', 'J', 'J'], ['X', 'Z', 'O', 'A', 'A', 'J', 'J']]
    # output_alignment_set = [['X', 'Z', 'O', 'U', 'A', 'X', 'X', 'C'], ['X', 'U', '-', 'A', 'A', 'J', 'J'], ['X', 'M', 'O', 'G', 'A', 'J', 'X'], ['X', 'M', 'O', 'A', 'U', 'J', 'B', 'C'], ['X', 'U', 'O', 'A', 'U', 'J', 'J'], ['X', 'U', 'O', 'A', 'B', 'J', 'J'], ['M', 'O', 'A', 'A', 'J', 'J'], ['J', 'Y', 'O', 'B', 'A', 'X', 'J']]

    # paths_set =  [['Q396', '', '', 'R395', 'V394', '', 'R495'], ['', 'S461', 'P462', '', '', 'L443', 'R444'], ['', '', 'P450', 'R449', '', '', 'R342']]
    # input_alignment_set = [['X', 'O', 'P', 'B', 'U', 'U', 'B'], ['X', 'O', 'P', 'B', 'U', 'U', 'B'], ['X', 'O', 'P', 'B', 'U', 'U', 'B']]
    # output_alignment_set = [['X', 'A', '-', 'B', 'U', 'Y', 'B'], ['B', 'O', 'P', '-', 'Z', 'U', 'B'], ['J', '-', 'P', 'B', '-', 'Z', 'B']]

    # paths_set = [['R444', 'L443', '', '', 'P462', 'S461', '', '', 'C467']]
    # input_alignment_set = [['B', 'U', 'U', 'B', 'P', 'O', 'X', 'B', 'C']]
    # output_alignment_set = [['B', 'U', 'Z', '-', 'P', 'O', 'B', '-', 'C']]

    # paths_set = [['G110', '', 'W112', 'F113', 'H114', '', 'R98', 'D118', ''], ['W120', 'G121', '', '', '', '', '', 'R39', 'E47', 'A45'], ['', 'W112', '', 'H33', 'D32', 'R29', 'D27', ''], ['W120', 'G121', 'Q122']]
    # input_alignment_set = [['G', 'X', 'Z', 'Z', 'H', 'J', 'B', 'J', 'A'], ['Z', 'G', 'X', 'Z', 'Z', 'H', 'J', 'B', 'J', 'A'], ['X', 'Z', 'Z', 'H', 'J', 'B', 'J', 'A'], ['Z', 'G', 'X']]
    # output_alignment_set =[['G', 'J', 'Z', 'Z', 'H', 'X', 'B', 'J', 'G'], ['Z', 'G', 'J', '-', 'Y', '-', 'X', 'B', 'J', 'A'], ['B', 'Z', 'U', 'H', 'J', 'B', 'J', 'O'], ['Z', 'G', 'X']]

#   Ksenia's: return set of dictionaries, 1 for each path
    set_of_dicts = Seed_graph.path_to_graph_dictionary(graph, paths_set, input_alignment_set, output_alignment_set, original_seed, seed)

#   Sapir's: functions calls to create a surface graph

#   Sapir's: locating the cluster found by pepsurf on the surface graph

#   Running cluster extension


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main_func()


