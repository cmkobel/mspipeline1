#!/usr/bin/env python

"""
This script makes snakemake able to read out the last n_lines of a failed job stderr-log.
"""

import re
import glob
import sys


# Documentation from snakemake:
# https://snakemake.readthedocs.io/en/stable/executing/cli.html
#--log-handler-script
# Provide a custom script containing a function ‘def log_handler(msg):’. Snakemake will call this function for every logging output (given as a dictionary msg)allowing to e.g. send notifications in the form of e.g. slack messages or emails.


n_lines = 20


def log_handler(msg):

    if msg["level"] == "error":
        print("=== ", msg) # debug

        # Many different messages might come out of on level == error. By using the try, we make sure to only use results from successfull parsings.
        try: 

            # Extract values from output string.
            m = re.search('Error executing rule (?P<name>.+) on cluster \(jobid: (?P<jobid>\d+), external: (?P<external>\d+), jobscript: (?P<jobscript>.+)\)\. For error details see the cluster log and the log files of the involved rule\(s\)\.', msg['msg'])

            name = m.groupdict()['name']
            jobid = m.groupdict()['jobid']
            external = m.groupdict()['external'] # "external" is the external jobid
            #print("=== ", msg) # debug
            #print("=== parsed", name, jobid, external, "===")
        
            # print the last few lines of the associated log file
            stderr_file = glob.glob("logs/" + external + "-" + jobid + "-" + name + ".err.log")[0]

            #print("Printing the last " + str(n_lines) + "lines from " + stderr_file + " below:")
            print("\033[91mDebug: Printing the last", str(n_lines), "lines from: ", stderr_file)

            # stderr_file = "logs/7641656-2-fail_tester.err.log"
            # TODO: Print only the last ~10 lines by first counting the number of lines.
            with open(stderr_file, 'r') as stderr_open:
                full = stderr_open.readlines()
                len_full = len(full)
                for i, line in enumerate(full[-n_lines:]):
                    #print("    stderr", i+len_full-n_lines+1, line, end ='')
                    print(f"    stderr {i+len_full-n_lines+1} |", line, end = '', file = sys.stderr)

        except:
            pass

    
    # I can't come up with anything relevant to add from these messages. 
    # if msg["level"] == "progress":
    #     print("%%%", msg["done"], msg["total"])
    #     print(msg)