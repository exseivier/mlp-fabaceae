#!/usr/bin/env python3
import numpy as np
import pandas as pd
import os
import random as rnd
import subprocess as sp
import json as js
import gc
from tqdm import tqdm


# Helper functions.

def save_matrix(matrix, fileout):
    """
        Input:  Matrix of data.
        Output: Save to file. (0)
    """
    matrix.to_csv(fileout, index = False, header = False)

def create_dictionary(kmer_size, alpha):
    """
        Creates all posible kamers of size kmer_size.
        Input:  kmer_size.
    """
    kmers = []
    opt1 = alpha.copy()
    opt2 = alpha.copy()
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

def calculate_gc_mfe_heat(seqs, k_dic_structure):
    """
        Calculate GC%, Minimum Free Energy and Specifric Heat
    """
    thermos = {}
    colnames = ["%GC", "MFE_by_nt", "MFE", "len"]
    colnames.extend(list(k_dic_structure.keys()))
    for h, s in tqdm(seqs.items()):     
        thermos[h] = [(s.count("G") + s.count("C")) / len(s)]
        p = sp.Popen(["RNAfold", "--noPS"], stdin = sp.PIPE, stdout = sp.PIPE, stderr = sp.PIPE)
        out, err = p.communicate(input = s.encode())
        mfe = float(out.decode().split(" ")[-1].strip().strip("(").strip(")"))
        thermos[h].append(mfe/len(s))
        thermos[h].append(mfe)
        thermos[h].append(len(s))
        struct = calculate_kmers_freq({h : out.decode().split("\n")[1].split(" ")[0]}, k_dic_structure, 5, true = "False")
        struct.pop(len(struct.iloc[0, :]) -1)
        thermos[h].extend(list(struct.iloc[0, :]))
        p.kill()
        
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
    skips = 0
    for head, seq in seqs.items():
        k_dic = reset_dic(k_dic)
        for i in range(len(seq) -k_size +1):
            seq_slice = seq[i:i + k_size].upper()
            try:
                k_dic[seq_slice] += 1
            except Exception as e:
                skips += 1
                print(f'KeyError: {seq_slice} - skips: {skips}')
        
        total_kmers = sum(k_dic.values())
        res_dict[head] = []
        if total_kmers >= len(seq)/2 and total_kmers > 0:
            for key, value in k_dic.items():
                res_dict[head].append(value / total_kmers)
        elif total_kmers < len(seq)/2 or total_kmers <= 0:
            for key, value in k_dic.items():
                res_dict[head].append(np.nan)

        if true in ("True"):
            res_dict[head].append(1)
        else:
            res_dict[head].append(0)
        
    df = pd.DataFrame(res_dict)
    df = df.transpose()
    
    return df

def load_fasta(filename):
#    print(f"Loading {filename}")
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
    files = None
    data = None
    counter = 0
    len_data = 0
    def __init__(self, path, fmt = 'json'):
        self.path = path
        if fmt == 'json':
            self.files = []
            self.data = load_json(f'{self.path}/dataset_catalog.json')
            for items in self.data['assemblies']:
                for f in items['files']:
                    self.files.append(f['filePath'])
            self.files = self.files[2:]
            self.len_data = len(self.files)

        elif fmt == 'txt':
            with open(f'{self.path}/dataset_catalog.txt') as fh:
                self.files = [x for x in fh.read().split("\n") if len(x) > 0]
                self.len_data = len(self.files)
                fh.close()
        else:
            raise Exception("ErrorIO: Unknown format\nSupported formats: [json | txt]")
        
    def __str__(self):
        return f"Files uploaded from {self.path}"
    
    def next(self):
        if self.counter == self.len_data:
            print("It is the end of data")
            print("next will return empty string")
            print("If you want to go start, use rewind method")
            return ""
        self.counter += 1
        return f'{self.path}/{self.files[self.counter -1]}'
        
    
    def rewind(self):
        print("Returning to the start position")
        self.counter = 0
    
    def total_files(self):
        return self.len_data

def cleaning_sequences(seqs):
    unwanted_letters = ["N","Y","R","M","K","S","W","H","D","B","V"]
    unwanted_letters.extend(list("".join(unwanted_letters).lower()))

    for head, seq in list(seqs.items()):
        counting = 0
        for letter in unwanted_letters:
            counting += seq.count(letter)

        if counting >= round(len(seq)/100):
            seqs.pop(head)

def create_fasta_subsamples(path):

    """

    """
    i = 0
    files = os.listdir(path)
    for f in files:
        i += 1
        with open(f'{path}/{f}') as fh:
            seqs = fh.read()

        seqs = seqs.split(">")
        
        for j in range(10):
            with open(f"../fasta/subsamples/dataset-{i}-{j}.fa", "w") as fh:
                fh.write(">"+">".join(rnd.sample(seqs, k = int(len(seqs) / 10))))



def create_fasta_dataset(rna_fasta_path, out_file, true_mirna_fasta, fmto = "txt"):

    """

    """

    files = GET_FILES(rna_fasta_path, fmt = fmto)
    print(f'Loading: {files.total_files()} files...')
    sample = []
    for i in tqdm(range(files.total_files())):
        sample.extend(sample_fasta(files.next(), 0.1))

    spt_sample = subseq_sampler(sample, "80-250", 4500)
    print(len(spt_sample), type(spt_sample))

    miRNA_list = load_fasta(true_mirna_fasta)
    miRNA_dict = {}
    for item in miRNA_list:
        item = item.split("\n")
        head = item[0]
        seq = ''.join(item[1:]).upper()
        if len(seq) <= 0:
            continue
        miRNA_dict[head] = seq
    if 0 in [len(x) for x in spt_sample.values()]:
        print("There is an empty sequence")

    # Cleaning the sequences
    print(f'Before cleaning. Got {len(spt_sample)} random sequences, {type(spt_sample)}')
    print(f'Before cleaning. Got {len(miRNA_dict)} true miRNA sequences, {type(miRNA_dict)}')

   
    cleaning_sequences(miRNA_dict)
    cleaning_sequences(spt_sample)

    print(f'After cleaning. Got {len(spt_sample)} random sequences, {type(spt_sample)}')
    print(f'After cleaning. Got {len(miRNA_dict)} true miRNA sequences, {type(miRNA_dict)}')

    seqs_list = []
    with open(out_file, "w") as fh:
        for head, seq in miRNA_dict.items():
            seqs_list.append(f'>v_{head}\n{seq}\n')
        for head, seq in spt_sample.items():
            seqs_list.append(f'>f_{head}\n{seq}\n')

        fh.write("".join(seqs_list))
    
    fh.close()

def create_dataset(rna_fasta_path, out_file, true_mirna_fasta, fmto = "txt"):

    files = GET_FILES(rna_fasta_path, fmt = fmto)
    print(f'Loading: {files.total_files()} files...')
    sample = []
    for i in tqdm(range(files.total_files())):
        sample.extend(sample_fasta(files.next(), 0.1))

    spt_sample = subseq_sampler(sample, "60-200", 4500)
    print(len(spt_sample), type(spt_sample))

    miRNA_list = load_fasta(true_mirna_fasta)
    miRNA_dict = {}
    for item in miRNA_list:
        item = item.split("\n")
        head = item[0]
        seq = ''.join(item[1:]).upper()
        if len(seq) <= 0:
            continue
        miRNA_dict[head] = seq
    if 0 in [len(x) for x in spt_sample.values()]:
        print("There is an empty sequence")

    # Cleaning the sequences
    print(f'Before cleaning. Got {len(spt_sample)} random sequences, {type(spt_sample)}')
    print(f'Before cleaning. Got {len(miRNA_dict)} true miRNA sequences, {type(miRNA_dict)}')

   
    cleaning_sequences(miRNA_dict)
    cleaning_sequences(spt_sample)

    print(f'After cleaning. Got {len(spt_sample)} random sequences, {type(spt_sample)}')
    print(f'After cleaning. Got {len(miRNA_dict)} true miRNA sequences, {type(miRNA_dict)}')

    k_dic_structure = create_dictionary(5, [".", ")", "("])
    k_dic_sequence = create_dictionary(5, ["A", "C", "G", "T"])
    print(list(k_dic_structure.items())[0:10])
    print(list(k_dic_sequence.items())[0:10])

    print("Calculating kmers in random sample...")
    mat_0 = calculate_kmers_freq(spt_sample, k_dic_sequence, 5, "False")
    print("Calculating kmers in true miRNA sequences...")
    mat_1 = calculate_kmers_freq(miRNA_dict, k_dic_sequence, 5, "True")
    mat = pd.concat([mat_0, mat_1], axis = 0)
    colnames = list(k_dic_sequence.keys())
    colnames.append("Train")
    mat.columns = colnames
    mat.reset_index(inplace = True, drop = True)
    
    print("Getting 2D structure information from random sample...")
    gc_mfe_sha_0 = calculate_gc_mfe_heat(spt_sample, k_dic_structure)
    print("Getting 2D structure information from true miRNA sequences...")
    gc_mfe_sha_1 = calculate_gc_mfe_heat(miRNA_dict, k_dic_structure)
    gc_mfe_sha = pd.concat([gc_mfe_sha_0, gc_mfe_sha_1], axis = 0)
    gc_mfe_sha.reset_index(inplace = True, drop = True)
    print(gc_mfe_sha)
    
    dataset = pd.concat([mat, gc_mfe_sha], axis = 1)
    print(dataset)
    
    print("Saving dataset...")
    dataset = dataset.dropna()
    dataset.to_csv(out_file, header = True, index = False)

def testing_models(models_path, datasets_path):

    """
        Test every model found in models_path with all testing
        datasets found in dataset_path
        Requires:
            models_path     [STR] Full path where keras models are.
            datasets_path   [STR] Full path where datasets are.
    """

    import tensorflow as tf
    models_files = os.listdir(models_path)
    datasets_files = os.listdir(datasets_path)
    res = {"dataset":[], "model":[], "good":[], "bad":[], "total":[], "accuracy":[]}
    print(f'Processing {len(models_files)} models and {len(datasets_files)} datasets')

    for model in tqdm(models_files):
        model = f'{models_path}/{model}'
        for ds_file in datasets_files:
            if model[-5:] == "keras" and ds_file.split(".")[-1] == "csv":
                ds_file = f'{datasets_path}/{ds_file}'
                
                lmodel = tf.keras.models.load_model(model)
                
                test_x = pd.read_csv(ds_file)
                test_label = test_x.pop("Train")
                test_x = tf.keras.utils.normalize(test_x)
                result = lmodel.predict(test_x, verbose = False)
                
                good = 0
                bad = 0
                for i in range(len(result)):
                    if (result[i].round() == test_label.iloc[i]).any():
                        good += 1
                    else:
                        bad += 1
                
                res["dataset"].append(ds_file)
                res["model"].append(model)
                res["good"].append(good)
                res["bad"].append(bad)
                res["total"].append(good + bad)
                res["accuracy"].append(round((good / (good + bad)) * 100, 4))

    return res


def predict(exon_fasta, models_path, out_dir, cutoff = 0.85):
    """
        Predicts miRNAs from sequences assembled from deep sequencing of RNA
        
        Requires:
            exon_fasta      Full path and file name of fasta file which contains
                            the input sequences.
            model           Full path and file name of keras file of the model.
            out_dir         Output directory name.
    """
    import tensorflow as tf
    transcripts = load_fasta(exon_fasta)
    models = os.listdir(models_path)
    
    if len(models) <= 0:
        raise Exception("Error: list of models is empty")

    lmodel_list = []
    for model in models:
        lmodel_list.append(tf.keras.models.load_model(f'{models_path}/{model}'))
    k_dic_seq = create_dictionary(5, ["A", "C", "G", "T"])
    k_dic_str = create_dictionary(5, [".", "(", ")"])
    
    spt_seqs = {}
    counter = 0
    positive = 0
    negative = 0
    total = len(transcripts)
    for t in transcripts:
        counter += 1        
        t = t.split("\n")
        t_head = t[0]
        t_seq = ''.join(t[1:])
#        if len(t_seq) <= 50 or len(t_seq) > 250:
#            continue
        
        spt_seqs[t_head] = t_seq

    with open(f"{out_dir}/predicted.txt", "w") as fp, open(f"{out_dir}/discarded.txt", "w") as fd:
        print("Predicting ...")
        total = 0
        positive = 0
        negative = 0
        num_of_trans = len(transcripts)
        counter = 0
        mat = calculate_kmers_freq(spt_seqs, k_dic_seq, 5, "False")
        colnames = list(k_dic_seq.keys())
        colnames.append("Train")
        try:
            mat.columns = colnames
        except Exception as e:
            print(mat)
            print(colnames)
            print(mat.columns)
            print(spt_seqs)
            print(t_seq)
        gc_mfe_sha = calculate_gc_mfe_heat(seqs=spt_seqs, k_dic_structure=k_dic_str)
        mat = pd.concat([mat, gc_mfe_sha], axis = 1)
        mat.pop("Train")
        mat = tf.keras.utils.normalize(mat.astype(float))
        results = [lmodel.predict(mat, verbose = False) for lmodel in lmodel_list]
        results = (sum(results)/len(results))
        results = [1 if x >= cutoff else 0 for x in results]
        print(results[0:10])
        
        heads = list(spt_seqs.keys())
        for i in range(len(heads)):
            if results[i] > 0:
                fp.write(f"{heads[i]} positive miRNA.\n")
                positive += 1
                print(f"{heads[i]} Positive (+) - {counter} / {num_of_trans} | +{positive} / {num_of_trans} | -{negative} / {num_of_trans}") 
            else:
                negative += 1
                fd.write(f"{heads[i]} was discarded as miRNA.\n")
                print(f"{heads[i]} Negative (-) - {counter} / {num_of_trans} | +{positive} / {num_of_trans} | -{negative} / {num_of_trans}")
            
        fp.write(f'Total transcripts: {total}')
        fp.write(f'Putative miRNA encoding gene: {positive}')
        fp.write(f'No miRNA encoding gene: {total - positive}')
        fd.close()
        fp.close()


def predict_and_stats(exon_fasta, models_path, model_name, cutoff = 0.85):
    """
        Predicts miRNAs from sequences assembled from deep sequencing of RNA
        and calculates Specificity, Sensitivity, PPV and NPV.
        
        Requires:
            exon_fasta      Full path and file name of fasta file which contains
                            the input sequences.
            model           Full path and file name of keras file of the model.
    """
    import tensorflow as tf
    transcripts = load_fasta(exon_fasta)
    models = os.listdir(models_path)
    
    if len(models) <= 0:
        raise Exception("Error: list of models is empty")

    lmodel_list = []
    for model in models:
        lmodel_list.append(tf.keras.models.load_model(f'{models_path}/{model}'))
    k_dic_seq = create_dictionary(5, ["A", "C", "G", "T"])
    k_dic_str = create_dictionary(5, [".", "(", ")"])
    
    spt_seqs = {}
    counter = 0
    positive = 0
    negative = 0
    total = len(transcripts)
    for t in transcripts:
        counter += 1        
        t = t.split("\n")
        t_head = t[0]
        t_seq = ''.join(t[1:])
        spt_seqs[t_head] = t_seq

    print("Predicting ...")
    total = 0
    true_positive = 0
    true_negative = 0
    false_positive = 0
    false_negative = 0
    num_of_trans = len(transcripts)
    counter = 0
    mat = calculate_kmers_freq(spt_seqs, k_dic_seq, 5, "False")
    colnames = list(k_dic_seq.keys())
    colnames.append("Train")
    try:
        mat.columns = colnames
    except Exception as e:
        print(mat)
        print(colnames)
        print(mat.columns)
        print(spt_seqs)
        print(t_seq)
    gc_mfe_sha = calculate_gc_mfe_heat(seqs=spt_seqs, k_dic_structure=k_dic_str)
    mat = pd.concat([mat, gc_mfe_sha], axis = 1)
    mat.pop("Train")
    mat = tf.keras.utils.normalize(mat.astype(float))
    results = [lmodel.predict(mat, verbose = False) for lmodel in lmodel_list]
    results = (sum(results)/len(results))
    results = [1 if x >= cutoff else 0 for x in results]
        
    heads = list(spt_seqs.keys())
    for i in range(len(heads)):
        if results[i] > 0:
            if heads[i].split("_")[0] == "v":
                true_positive += 1
            elif heads[i].split("_")[0] == "f":
                false_positive += 1
        else:
            if heads[i].split("_")[0] == "v":
                false_negative += 1
            elif heads[i].split("_")[0] == "f":
                true_negative += 1

    #   Returns Sensitivity, Specificity, PPV and NPV in that order
    return true_positive / (true_positive + false_negative), \
            true_negative / (true_negative + false_positive), \
            true_positive / (true_positive + false_positive), \
            true_negative / (false_negative + true_negative)

