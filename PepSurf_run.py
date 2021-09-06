import sys
import subprocess
import pathlib


def create_seed_txt(seed):
    f = open("seed.txt", "w+")
    f.write("#LIBRARY_TYPE = NNK\r\n#STOP_CODON_MODIFICATION TAG = GLU\r\n#PEPTIDES_START\r\n>\r\n"+seed+"\r\n#PEPTIDES_END")
    f.close()


# arguments
# -I seed.txt -P 1e6j.pdb -C P -S 1e6j_P.txt
def perpsurf_run(seed):
    my_list = sys.argv
    exec_name = str(pathlib.Path.cwd()) + "/Pepsurf"
    arguments = []
    arguments.append(exec_name)
    arguments.append('-I')
    create_seed_txt(seed)
    arguments.append('seed.txt')
    for i in range(len(my_list)):
        if (my_list[i] == '-P') or (my_list[i] == '-C') or (my_list[i] == '-S'):
            arguments.append(my_list[i])
            arguments.append(my_list[i+1])
    print(arguments[0]+'\n')
    res = subprocess.run(arguments, capture_output=True, text=True)