import sys
import subprocess
import pathlib

# Mapitope executable (mapi_exec) should be located in the same folder with all the input files and main.py code
# input files:
# mapi_exec arguments:
def mapi_run():
    arguments = sys.argv
    exec_name = str(pathlib.Path.cwd()) + "/mapi_exec"
    arguments[0] = exec_name
    print(arguments[0]+'/n')
    res = subprocess.run(arguments, capture_output=True, text=True)
    print("stdout:", res.stdout)
    print("stderr:", res.stderr)

