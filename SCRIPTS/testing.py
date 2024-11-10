#!/usr/bin/env python3

print("Loading os")
import os
print("Loading tensorflow")
import tensorflow as tf
print("Loading pandas")
import pandas as pd
print("Loading Te Quiero DeMasiado")
from tqdm import tqdm
print("Loading subprocess")
import subprocess as sp

def testing_models(models_path, datasets_path):
    models_files = os.listdir(models_path)
    datasets_files = os.listdir(datasets_path)
    res = {"dataset":[], "model":[], "good":[], "bad":[], "total":[], "accuracy":[]}
    print(f'Processing {len(models_files)} models and {len(datasets_files)} datasets')
    for model in tqdm(models_files):
        model = f'{models_path}/{model}'
        for ds_file in datasets_files:
            if model[-5:] == "keras" and ds_file.split(".")[-1] == "csv":
                ds_file = f'{datasets_path}/{ds_file}'
                #print(f"Processing {model}")
                lmodel = tf.keras.models.load_model(model)
                #print(f'Processing {ds_file}')
                test_x = pd.read_csv(ds_file)
                test_label = test_x.pop("Train")
                test_x = tf.keras.utils.normalize(test_x)
                result = lmodel.predict(test_x, verbose = False)
                #print(len(result))
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

                #print(f'\tTotal predictions: {good + bad}')
                #print(f'\tGood predictions: {good}')
                #print(f'\tBad preditions: {bad}')
                #print(f'\t\tAccuracy {good / (good+bad)}')
                #print("-----------------------------------------")

    return res

#######################################################
##
##    Helper functions
##

def split_seqs(seqs, size, step):
    """
        Splits DNA sequence.
        seqs:    Dictionary with sequences
        size:    Window size
        step:    Window step
    """
    spt_dna = {}
    for head, seq in seqs.items():
        lenseq = len(seq)
        if size >= lenseq:
            stop = lenseq
        else:
            stop = lenseq -size +1
        for i in range(0, stop, step):
            if len(seq[i: i +size]) < size -1:
                continue
            spt_dna[f'{head}_frag_{i}'] = seq[i : i +size]

    return spt_dna

def reset_dic(dic):
    """

    """
    for key in dic.keys():
        dic[key] = 0
    return dic

def select_the_biggest_sequence(seqs):
    """
        Will select the biggest transcript-variant sequence.
        Require:
            seqs:    A vector containing items of header and sequence
    """
    filtered = {}
    for i in set([item.split("\n")[0].split("|")[0] for item in seqs]):
        filtered[i] = ""

    for item in seqs:
        item = item.split("\n")
        head = item[0].split("|")[0]
        seq = ''.join(item[1:])
        filtered[head] = seq if len(seq) > len(filtered[head]) else filtered[head]
    out_vector = []
    for head, seq in filtered.items():
        out_vector.append(f'{head}\n{seq}\n')

    return out_vector

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
    for h, s in seqs.items():
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
    for head, seq in seqs.items():
        k_dic = reset_dic(k_dic)
        for i in range(len(seq) -k_size +1):
            seq_slice = seq[i:i + k_size].upper()
            if seq_slice.count("N") > 0\
                    or seq_slice.count("Y") > 0\
                    or seq_slice.count("R") > 0\
                    or seq_slice.count("M") > 0\
                    or seq_slice.count("K") > 0\
                    or seq_slice.count("S") > 0\
                    or seq_slice.count("W") > 0:
                continue
            try:
                k_dic[seq_slice] += 1
            except Exception as e:
                print(seq_slice)
                print(seq)
                raise Exception("Check this out!!!")

        total_kmers = sum(k_dic.values())
        if total_kmers < len(seq)/2:
            continue

        res_dict[head] = []
        for key, value in k_dic.items():
            try:
                res_dict[head].append(value / total_kmers)
            except Exception as e:
                print(seq)
                print(total_kmers)
                raise Exception(e)

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


def predict(exon_fasta, model, out_dir):

    transcripts = load_fasta(exon_fasta)
    lmodel = tf.keras.models.load_model(model)
    k_dic_seq = create_dictionary(5, ["A", "C", "G", "T"])
    k_dic_str = create_dictionary(5, [".", "(", ")"])

    print(f'Before filtering {len(transcripts)}')
    transcripts = select_the_biggest_sequence(transcripts)
    print(f'After filtering {len(transcripts)}')
    print(transcripts[0:5])

    with open(f"{out_dir}/predicted.txt", "a") as fp, open(f"{out_dir}/discarded.txt", "a") as fd:
        total = 0
        positive = 0
        negative = 0
        num_of_trans = len(transcripts)
        counter = 0
        for t in transcripts:
            total += 1
            counter += 1        
            t = t.split("\n")
            t_head = t[0]
            t_seq = ''.join(t[1:])
            #spt_seqs = split_seqs({t_head : t_seq}, 150, 100)
            spt_seqs = {t_head : t_seq}
            if len(t_seq) <= 50 or len(t_seq) > 200:
                continue
            #fh.write(f"Testing {t_head}...\n")
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
            #print(mat.astype(float))
            mat = tf.keras.utils.normalize(mat.astype(float))
            results = lmodel.predict(mat, verbose = False)
            if (sum(results.round()) > 0).any():
                #print(results.round())
                fp.write(f"{t_head} is a putative miRNA encoding gene!\n")
                positive += 1
                print(f"{t_head} Positive (+) - {counter} / {num_of_trans} | +{positive} / {num_of_trans} | -{negative} / {num_of_trans}") 
            else:
                negative += 1
                fd.write(f"{t_head} was discarded as miRNA encoding gene!\n")
                print(f"{t_head} Negative (-) - {counter} / {num_of_trans} | +{positive} / {num_of_trans} | -{negative} / {num_of_trans}")
            
        fp.write(f'Total transcripts: {total}')
        fp.write(f'Putative miRNA encoding gene: {positive}')
        fp.write(f'No miRNA encoding gene: {total - positive}')
        fd.close()
        fp.close()
