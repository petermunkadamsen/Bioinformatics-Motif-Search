#!/usr/bin/env python3
import sys
import re

#function at exits program with a usage message
def usage(mesg=None):
    if mesg is not None:
        print("Error: ",  mesg)
    # print("Usage: motif_search.py <sequencefile> <signaldescriptionfile> <deviation_number>")
    sys.exit(1)

# exit the program if the correct arguments are not given

#Julies input: python3 /Users/juliehaugaard/Desktop/unix_og_python/projekt/scr/motif_search_copy.py /Users/juliehaugaard/Desktop/unix_og_python/projekt/motif.fsa /Users/juliehaugaard/Desktop/unix_og_python/projekt/signal_description.tsv 16



def read_file(filename):
    """takes a filename as a parameter and returns a dict {header: sequence}"""
    # make the dict to store sequences
    sequences = dict()

    # define default values
    fasta = False
    header  = None
    seq = ""

    # open the file
    try:
        with open(filename, "r") as infile:
            #read lines
            lines = [line.strip() for line in infile if line.strip()] # removes if line is empty

            #see if there is a ">" at any point 
            for line in lines:
                if line.startswith(">"):
                    fasta = True
                    
            #if it is a fasta file
            if fasta: 
                for line in lines: 
                    if line.startswith(">"):
                        if header and seq: 
                            sequences[header] = seq
                        header = line
                        seq = ""
                    else: 
                        seq += "".join(line.split()).upper()
                if header:
                    sequences[header] = seq # remember the last header and sequence
            else: # if the file is not a fasta file 
                lineno = 0
                for line in lines:
                    lineno += 1
                    header = f">line{lineno}" # header becomes: line1, line2, line3 etc.
                    sequences[header] = line.upper()
            
            if not sequences:
                usage(f"File '{filename}' is empty or contains no valid sequences.")

    except IOError as err:
        usage(f"Can't read file '{filename}', reason: {err}")
    except Exception as err:
        usage(f"Unexpected error: {err}")
    return sequences 



def read_signal_file(signal_description_file):
    """takes a filename as a parameter and returns lists of dicts"""
    # make the list to store the dicts
    signal_description = list() 

    # open the file
    try: 
        with open(signal_description_file, "r") as infile:
            for line in infile:
                line = line.strip()
                if not line or line.startswith("#"): # ignored by program
                    continue
                if line.startswith("*"):
                    gap = re.search(r'\*\s+(\d+)\-(\d+)', line)
                    if gap: 
                        min_gap = gap.group(1)
                        max_gap = gap.group(2)
                        signal_description.append({
                            "type": "GAP",
                            "range": (min_gap, max_gap)
                        })
                    else: 
                        usage(f"Invalid GAP line format: {line}")
                else:
                    column = line.split("\t")
                    if len(column) != 2:
                        usage(f"Invalid MATCH line format: {line}")
                    signal = list(column[0]) # multiple allowed letters like 'AT' become ['A', 'T']
                    penalty = int(column[1])
                    signal_description.append({
                        "type": "MATCH",
                        "allowed": signal,
                        "penalty": penalty
                    })
            if not signal_description:
                usage(f"File '{filename}' is empty or contains no valid signals.")

    except IOError as err:
        usage(f"Can't read file '{signal_description_file}', reason: {err}")
    except Exception as err:
        usage(f"Unexpected error: {err}")
    return signal_description



def motif_search(sequences, signal_description, deviation):
    results = {}
    for header, seq in sequences.items():
        seq = seq.upper()
        matches = []
        seq_len = len(seq)

        for start in range(len(seq)):
            seq_pos = start
            signal_description_index = 0
            penalty = 0

            while signal_description_index < len(signal_description):
                signal = signal_description[signal_description_index]

                if signal["type"] == "MATCH":
                    if seq_pos >= seq_len:
                        break
                    if seq[seq_pos] not in signal["allowed"]:
                        penalty += signal["penalty"]
                    seq_pos += 1
                    signal_description_index += 1

                elif signal["type"] == "GAP":
                    min_gap = int(signal["range"][0])
                    max_gap = int(signal["range"][1])
                    if seq_pos + min_gap > seq_len:
                        break
                    seq_pos += min_gap
                    signal_description_index += 1

                if penalty > deviation:
                    break

            if signal_description_index == len(signal_description) and penalty <= deviation:
                match = seq[start:seq_pos]
                matches.append((start, seq_pos, match, penalty))

        results[header] = matches 

    return results


# main program - flyt dette ind i en if __name__ == "__main__" blok

if __name__ == "__main__":
    # exit the program if the correct arguments are not given
# exit the program if the correct arguments are not given
    if len(sys.argv) != 4:
        usage("Usage: motif_search.py <sequencefile> <signaldescriptionfile> <deviation_number>")
    if len(sys.argv == 4):
        filename = sys.argv[1]
        signal_description_file = sys.argv[2]
        deviation_number = int(sys.argv[3])
    else:
        filename = input("Give me a file containing one or more sequences: ")
        signal_description_file = input("Give me a signal description file: ")
        deviation_number = int(input("What should the deviation number be? "))

    sequences = read_file(filename)
    print(sequences)

    signals = read_signal_file(signal_description_file)
    print(signals)

    results = motif_search(sequences, signals, deviation_number)

    for header, hits in results.items():
        if hits:  # kun hvis der er matches i listen
            print(f"\nMatches in {header}:")
            for start, end, subseq, pen in hits:
                print(f"  {start}-{end}: {subseq} (penalty: {pen})")
