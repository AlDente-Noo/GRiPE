import os
import sys
import subprocess

from utils import *

def main(filename):
    print(filename)
    f = open(filename)
    for l in f:
        print(l.strip())
        proc = subprocess.run(['java', '-classpath', JAVA_CLASS_PATH, JAVA_CLASS, l.strip()], capture_output=True)
        output = proc.stdout.decode('utf-8')
        print(output)
        output_lines = output.splitlines()
        if not output_lines[0].startswith('output'):
            print(output_lines)
            print('Error during GRiP simulation.')
            sys.exit()
        id = output_lines[0].split(' ')[-1]
        prefix = output_lines[1].split(' ')[-1]
        dir = output_lines[2].split(' ')[-1]
        subprocess.run(['python', DATAPROC_DIR + SITE_OCCUPANCY_SCRIPT, '-i', str(id), '-p', str(prefix),
                        '-d', str(dir)])
    f.close()

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Error: no filename is given.")
        #sys.exit()
    if len(sys.argv) > 2:
        print("Warning: too many arguments, using only the first one.")
    main(sys.argv[1])
    #stream = os.popen('echo Returned output')
    #output = stream.read()
    #print(type(output.splitlines()))
    #stream = subprocess.run(['echo', 'Returned output'], capture_output=True)
    #print(stream.stdout, type(stream.stdout))
