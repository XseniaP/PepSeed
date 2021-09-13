import sys
import subprocess
import pathlib


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
    res = subprocess.run(arguments, capture_output=True, text=True)

