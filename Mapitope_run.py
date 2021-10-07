import sys
import subprocess
import pathlib
import os


# Mapitope executable (under the name: mapi_exec) should be located in the same folder with
# all the input files and main.py code
# input files:
# mapi_exec arguments:
# -U Results -P 1e6j_P.pdb -S 1e6j_P.txt -I 13b5.txt -C P -D 17.0 -V 2.5 -R rasmol.txt -F 0
def mapi_run():
    arguments = sys.argv
    exec_name = str(pathlib.Path.cwd()) + "/mapi_exec"
    arguments[0] = exec_name
    # arguments[0] = "Results_Mapi"
    print(arguments[0]+'\n')
    # print(arguments[4] + '\n')
    # arguments2 = [None] * 5
    # arguments2[0] = str(pathlib.Path.cwd()) + "/Surface_Racer_5.0"
    # arguments2[1] = arguments[6].split(".")[0] + ".pdb"
    # os.startfile(str(pathlib.Path.cwd()) + "/Surface_Racer_5.0")
    # subprocess.call(str(pathlib.Path.cwd()) + "/Surface_Racer_5.0")
    # name = str(pathlib.Path.cwd()) + "/surfrace"
    # arguments2[0] = name
    # arguments2[1] = 2
    # arguments2[2] = arguments[4]
    # arguments2[3] = 1.4
    # arguments2[4] = 1
    # arguments2 = [name, 2, arguments[4], 1.4, 1]
    # print (arguments2)
    # subprocess.run(arguments2)
    # subprocess.run([str(pathlib.Path.cwd()) + "/Surface_Racer_5.0"])
    res = subprocess.run(arguments, capture_output=True, text=True)

