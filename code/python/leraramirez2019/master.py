#!/usr/bin/env python
"""
"""
import sys, os
from leraramirez2019 import scan
import argparse

def report_fromfile(report_folder,in_file,targets):

    report = os.path.join(report_folder, "report")

    with open(in_file) as ins:
        for line in ins:
            if "#" not in line:
                print "-----> Reporting",line.strip()
                scan.main([report + ' ' + line] + targets, njobs=nb_jobs)

def get_fromdir(parent,summary_dir,targets):

    for python_file in os.listdir(parent):

        if python_file.startswith("get_") and python_file.endswith(".py"):

            get_what = python_file[4:-3]+".txt"
            print "-----> Getting", python_file[4:-3]

            out_file = os.path.join(summary_dir,get_what)
            path_script = os.path.join(parent,python_file)

            if os.path.isfile(out_file):
                os.remove(out_file)

            # Call scan with the same python interpreter as the one from the file
            scan.main([sys.executable + ' ' + path_script + ' > ' + get_what] + targets,njobs=nb_jobs)

            if get_what=="parameters.txt":

                scan.main(['/bin/cat ' + get_what + ' >>'+out_file] + targets[0:1])
                scan.main(['/usr/bin/tail -n1 ' + get_what + ' >> '+out_file] + targets[1:])

            else:
                scan.main(['/bin/cat ' + get_what + ' >>' + out_file] + targets)

def main(args):

    paths = []

    for arg in args:

        if os.path.isdir(arg):
            # Absolute paths
            paths.append(os.path.abspath(arg))
        else:
            sys.stderr.write("  Error: unexpected argument `%s'\n" % arg)
            sys.exit()

    if not paths:
        sys.stderr.write("  Error: you must specify directories\n")
        sys.exit()

    for p in paths:

        print "-----> Called master.py for",p

        scan_folder = os.path.join(p,'scan')

        if not os.path.isdir(scan_folder):
            print p, "skipped"
            continue

        targets = [os.path.join(scan_folder,i) for i in os.listdir(scan_folder) if "DS_Store" not in i]
        targets.sort()

        summary_dir = os.path.join(p,'..', '..', 'data_summary')

        run_dir_name=os.path.basename(os.path.normpath(p))

        if not os.path.isdir(summary_dir):
            os.mkdir(summary_dir)

        # Create the folder ../data_summary/p
        summary_dir = os.path.join(summary_dir, run_dir_name)

        if not os.path.isdir(summary_dir):
            os.mkdir(summary_dir)

        report_fromfile(bin_folder,os.path.join(p,"..","report_instructions.txt"),targets)

        get_fromdir(os.path.join(p,".."),summary_dir,targets)

# ------------------------------------------------------------------------

if __name__ == "__main__":

    bin_folder = "/Users/Manu/Documents/Projects/TheoryPaper/rigid_crosslinkers/code/cytosim/bin"
    nb_jobs = 20

    parser = argparse.ArgumentParser(description=__doc__)

    parser.add_argument('--nb_jobs',type=int)
    parser.add_argument('--bin_folder',type=str)

    args,unknown = parser.parse_known_args()

    if args.nb_jobs is not None:
        nb_jobs = args.nb_jobs
    if args.bin_folder is not None:
        bin_folder = args.bin_folder


    if len(unknown) == 0 or unknown[0] == 'help':
        parser.print_help()
    else:
        main(unknown)