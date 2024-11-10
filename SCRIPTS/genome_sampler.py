def load_fasta(filename):
    """
        Input: File name (str).
        Output: File handler (*FH).
    """
    seqs = {}
    with open(filename, "r") as fh:
        for line in fh:
            line = line.strip()
            if len(line) <= 0:
                continue
            if line[0] == ">":
                header = line[1:]
                seqs[header] = []
            else:
                seqs[header].append(line.upper())

        fh.close()
    for key in seqs.keys():
        seqs[key] = "".join(seqs[key])

    return seqs

def genome_sampler(seqs, size_range, num_seqs):
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
    keys = list(seqs.keys())
    sample = {}
    for i, v in enumerate(index):
        lenseq = len(seqs[keys[v]])
        rnd_size = np.random.randint(size_range[0], size_range[1]+1)
        rnd_position = np.random.randint(0, lenseq -rnd_size +1)
        sample[f'RndSeq_{i}'] = seqs[keys[v]][rnd_position : rnd_position +rnd_size]

    return sample