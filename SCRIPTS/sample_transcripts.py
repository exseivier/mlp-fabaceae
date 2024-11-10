import os
import sys
import pandas as pd
import numpy as np
import subprocess as sp
import json as js
import random as rnd

# Helper functions.

def save_matrix(matrix, fileout):
    """
        Input:  Matrix of data.
        Output: Save to file. (0)
    """
    matrix.to_csv(fileout, index = False, header = False)

def create_dictionary(kmer_size):
    """
        Creates all posible kamers of size kmer_size.
        Input:  kmer_size.
    """
    kmers = []
    opt1 = ['A', 'C', 'G', 'T']
    opt2 = ['A', 'C', 'G', 'T']
    for _ in range(kmer_size-1):
        for X in (opt1):
            for Y in (opt2):
                kmers.append(Y+X)

        opt2 = kmers
        kmers = []

    d = {}
    for k in opt2:
        d[k] = 0

    return d

def reset_dic(dic):
    """

    """
    for key in dic.keys():
        dic[key] = 0
    return dic

def calculate_gc_mfe_heat(seqs):
    """
        Calculate GC%, Minimum Free Energy and Specifric Heat
    """
    thermos = {}
    colnames = ["%GC", "Length", "MFE", "MFE_by_nt", "t55", "t56", "t57", "t58", "t59", "t60", "t61", "t62", "t63", "t64", "t65"]
    for h, s in seqs.items():
        print(f"Processing {h}...")
        thermos[h] = [(s.count("G") + s.count("C")) / len(s)]
        thermos[h].append(len(s))
        p = sp.Popen(["RNAfold", "--noPS"], stdin = sp.PIPE, stdout = sp.PIPE, stderr = sp.PIPE)
        out, err = p.communicate(input = s.encode())
        #print(out.decode().split(" ")[-1].strip().strip("(").strip(")"))
        out = float(out.decode().split(" ")[-1].strip().strip("(").strip(")"))
        thermos[h].append(out)
        thermos[h].append(out/len(s))
        p.kill()

        p = sp.Popen(["RNAheat", "--Tmin=55", "--Tmax=65"], stdin = sp.PIPE, stdout = sp.PIPE, stderr = sp.PIPE)
        out, err = p.communicate(input = s.encode())
        for item in out.decode().split("\n")[:-1]:
            thermos[h].append(item.split("\t")[1])

    thermos = pd.DataFrame(thermos)
    thermos = thermos.transpose()
    thermos.columns = colnames

    return thermos
        

def calculate_kmers_freq(seqs, k_dic, k_size, true = "True"):
    """
        Input:  Sequences with IDs (dictionary).
        Output: Matrix of scaled and normalised frequencies.
                scaling = kmer(i) / total_kmers.
                nomr_data = min-max_Normalization.
    """
    res_dict = {}
    for head, seq in seqs.items():
        k_dic = reset_dic(k_dic)
        res_dict[head] = []
        for i in range(len(seq) -k_size +1):
            seq_slice = seq[i:i + k_size]
            if seq_slice.count("N") > 0\
                    or seq_slice.count("Y") > 0\
                    or seq_slice.count("R") > 0\
                    or seq_slice.count("M") > 0\
                    or seq_slice.count("K") > 0\
                    or seq_slice.count("S") > 0\
                    or seq_slice.count("W") > 0:
                continue
            k_dic[seq_slice] += 1


#        xmax = pd.Series(k_dic).max()
#        xmin = pd.Series(k_dic).min()
        for key, value in k_dic.items():
#            value = (value - xmin) / (xmax - xmin)
            res_dict[head].append(value)

        if true in ("True"):
            res_dict[head].append(1)
        else:
            res_dict[head].append(0)
        
    df = pd.DataFrame(res_dict)
    df = df.transpose()
    
    return df

def load_fasta(filename):
    print(f"Loading {filename}")
    if filename == "":
        print("Filename is empty")
        print("An empty string is returned")
        return ""
        
    with open(filename, "r") as fh:
        data = fh.read().split(">")
        fh.close()
        data = data[1:]
    
    return data

def sample_fasta(filename, fraction):
    seqs = load_fasta(filename)
    len_seqs = len(seqs)
    num_seqs = len_seqs * fraction
    sample_seqs = rnd.sample(seqs, k = int(round(num_seqs)))
    return sample_seqs

def subseq_sampler(seqs, size_range, num_seqs):
    """
        Takes a random sample of genomic secuences
        of determined size and number.
        Requires
        seqs:            sequences (dictionary)
        siaze_range:     Range of sizes (str) example: "100-200"
        num_seqs:        Number of sequences (int)
    """
    index = np.random.choice(range(len(seqs)), size = num_seqs, replace = True)
    index.sort()
    size_range = size_range.split("-")
    size_range = [int(x) for x in size_range]
    sample = {}
    for i, v in enumerate(index):
        seq = seqs[v].split("\n")
        head = seq[0]
        seq = ''.join(seq[1:])
        lenseq = len(seq)
        rnd_size = np.random.randint(size_range[0], size_range[1]+1)
        
        if rnd_size >= lenseq:
            rnd_position = 0
            rnd_size = lenseq
        else:
            rnd_position = np.random.randint(0, lenseq -rnd_size +1)

        sample[f'RndSeq_{i}'] = seq[rnd_position : rnd_position +rnd_size]
        if len(sample[f'RndSeq_{i}']) <= 0:
            print("Warning! Subsequence is equal 0!")
            print(f'Seq length: {len(seq)}')
            print(f'Random size: {rnd_size}')
            print(f'Empty sequence: {head} - {seq[rnd_position : rnd_position +rnd_size]}')
            print(f'coordinates: {rnd_position} - {rnd_position +rnd_size}')
        sample[f'RndSeq_{i}'].upper()
    
    return sample

def load_json(jsonfile):
    with open(jsonfile, "rb") as jfh:
        data = js.load(jfh)
        jfh.close()
    return data


# Classes

class GET_FILES():
    file = None
    data = None
    counter = 0
    len_data = 0
    def __init__(self, path, fmt = 'json'):
        self.path = path
        if fmt == 'json':
            self.data = load_json(f'{self.path}/dataset_catalog.json')
            self.len_data = len(self.data['assemblies']) -1
        elif fmt == 'tsv':
            print("Under construction!")
        elif fmt == 'csv':
            print("Under construction!")
        else:
            raise Exception("ErrorIO: Unknown format\nSupported formats: [json | tsv | csv]")
        
    def __str__(self):
        return f"Files uploaded from {self.path}"
    
    def next(self):
        if self.counter == self.len_data:
            print("It is the end of data")
            print("next will return empty string")
            print("If you want to go start, use rewind method")
            return ""
        self.counter += 1
        return f"{self.path}/{self.data['assemblies'][self.counter]['files'][0]['filePath']}"
    
    def rewind(self):
        print("Returning to the start position")
        self.counter = 0
    
    def total_files(self):
        return self.len_data


def main():

    if sys.argc != 3:
        print("Usage: script <transcripts_path> <mirna_seqs>")
    transcripts_path = sys.argv[1]
    mirna_seqs = sys.argv[2]
    files = GET_FILES(transcripts_path, fmt = 'json')
    print(f'Total files: {files.total_files()}')
    sample = []
    for i in range(files.total_files()):
        print(f'Processing {i}')
        sample.extend(sample_fasta(files.next(), 0.1))

    spt_sample = subseq_sampler(sample, "100-300", 4500)
    print(len(spt_sample))

    miRNA_list = load_fasta(mirna_seqs)
    miRNA_dict = {}
    for item in miRNA_list:
        item = item.split("\n")
        head = item[0]
        seq = ''.join(item[1:]).upper()
        miRNA_dict[head] = seq
        if 0 in [len(x) for x in spt_sample.values()]:
            print("There is an empty sequence")

    print(len(spt_sample))
    print(len(miRNA_dict))

    k_dic = create_dictionary(5)
    mat_0 = calculate_kmers_freq(spt_sample, k_dic, 5, "False")
    mat_1 = calculate_kmers_freq(miRNA_dict, k_dic, 5, "True")
    mat = pd.concat([mat_0, mat_1], axis = 0)
    colnames = list(k_dic.keys())
    colnames.append("Train")
    mat.columns = colnames
    mat.reset_index(inplace = True, drop = True)
    print(mat)
    
    gc_mfe_sha_0 = calculate_gc_mfe_heat(spt_sample)
    gc_mfe_sha_1 = calculate_gc_mfe_heat(miRNA_dict)
    gc_mfe_sha = pd.concat([gc_mfe_sha_0, gc_mfe_sha_1], axis = 0)
    gc_mfe_sha.reset_index(inplace = True, drop = True)
    print(gc_mfe_sha)
    
    dataset = pd.concat([mat, gc_mfe_sha], axis = 1)
    print(dataset)

if __name__ == "__main__":
    main()