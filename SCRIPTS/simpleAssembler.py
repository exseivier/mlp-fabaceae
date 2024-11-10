#!/usr/bin/env python3

import pandas as pd
import pysam as ps
from tqdm import tqdm


class BAMHOLDER():
    
    bamfile = None
    chrs_number = None
    headers = None
    chr_count = None

    samiter = None
    samiter_len = None
    samiter_count = None

    def __init__(self, filename):
        self.bamfile = ps.AlignmentFile(filename, "rb")
        self.headers = self.bamfile.header.as_dict()["SQ"]
        self.chrs_number = len(self.headers)
        self.chr_count = 0
        self.samiter_count = 0

    def __str__(self):
        if self.bamfile != None:
            return self.bamfile
        else:
            return "bamfile is empty"

    def get_headers(self):
        return self.headers

    def get_chrs_number(self):
        return self.chrs_number

    def get_bamfile(self):
        return self.bamfile

    def get_aligns_by_chr(self):
        if self.chr_count == self.chrs_number:
            print(f"It is the end of the BAM file\nUse rewind function to initialize a new iterator")
            return "-101"
        else:
            self.samiter = self.bamfile.fetch(self.headers[self.chr_count]["SN"], \
                    0, \
                    self.headers[self.chr_count]["LN"], \
                    multiple_iterators = True)
            self.samiter_len = self.bamfile.count(self.headers[self.chr_count]["SN"])
            self.chr_count += 1
            self.samiter_count = 0
    
    def get_aligns_by_segment(self, chr_num = 0, start = 0, end = 20000000):
        if self.chr_count == self.chrs_number:
            print(f"It is the end of the BAM file\nUse rewind function to initialize a new iterator")
            return "-101"
        else:
            self.samiter = self.bamfile.fetch(self.headers[self.chr_count]["SN"], \
                    start, \
                    end, \
                    multiple_iterators = True)
            self.samiter_len = len(list(self.bamfile.fetch(self.headers[self.chr_count]["SN"], \
                    start, \
                    end, \
                    multiple_iterators = True)))
            self.chr_count = 0
            self.samiter_count = 0
 

    def get_align(self):
        if self.samiter_count == self.samiter_len:
            print("It is the end of chromosome\nUse get_aligns_by_chr function to continue to the next chromosome")
            return "-101"
        else:
            self.samiter_count += 1
            return next(self.samiter)


class SIMPLEASSEMBLY():
    dgtf = None # Dictionary with information used to create a gtf file.
    bamh = None # This variable will be the BAM-file holder class.

    def __init__(self, bamfilename):
        self.bamh = BAMHOLDER(bamfilename)
        self.dgtf = {}
        self.dgtf_filt = {}

    def __str__(self):
        return "Hello world!"

    def assembly(self, cutoff):
        assemblies_count = 0
        mapped_reads_count = 0
        for i in tqdm(range(self.bamh.chrs_number)):
            assemblies_count += 1
            self.bamh.get_aligns_by_chr()
            assembly_name = f'ASSEMBLY_{assemblies_count}'           
            self.dgtf[self.bamh.headers[i]["SN"]] = {assembly_name : {"start" : 1, "end" : 0, "mapped" : 0}}
            for j in range(self.bamh.samiter_len):
                ii = self.bamh.get_align()
                if len(self.dgtf) > 0:
                    if ii.reference_start +1 - self.dgtf[ii.reference_name][assembly_name]["end"] <= cutoff:
                        mapped_reads_count += 1
                        self.dgtf[ii.reference_name][assembly_name]["end"] = ii.reference_end
                        self.dgtf[ii.reference_name][assembly_name]["mapped"] = mapped_reads_count
                        
                    else:
                        assemblies_count += 1
                        mapped_reads_count = 1
                        assembly_name = f"ASSEMBLY_{assemblies_count}"
                        self.dgtf[ii.reference_name].update({assembly_name : \
                                {"start" : ii.reference_start +1, \
                                "end" : ii.reference_end, \
                                "mapped" : mapped_reads_count}})
                else:
                    mapped_reads_count += 1
                    self.dgtf[ii.reference_name] = {assembly_name : \
                            {"start" : ii.reference_start +1, "end" : ii.reference_end, "mapped" : mapped_reads_count}}

    def join_assembly(self, select_range = "20-30", cov = 200):

        srange = [int(x) for x in select_range.split("-")]
        for chromosome, values in tqdm(self.dgtf.items()):
            if len(values) == 0:
                continue

            data = pd.DataFrame(values)
            data = data.transpose()

            keep_coverage = data["mapped"] >= cov
            keep_size = (data["end"] - data["start"]  +1 >= srange[0]) & (data["end"] - data["start"] +1 <= srange[1])
            keep = keep_coverage & keep_size

            data = data[keep]
            data_copy = data.copy()
            data_copy.index = [f'{x}_b' for x in data_copy.index]
            data["start"] -= 110
            data.loc[(data["start"] <= 0), "start"] = 1                       ####    BUG     ####  Negative start values
            data_copy["end"] += 110
            data = pd.concat([data, data_copy], axis=0)
            data = data.sort_values("start")

            data = data.transpose()
            data = data.to_dict()
            
            self.dgtf[chromosome] = data

    def filter(self, select_range = "60-130", cov = 200):

        """

        """
        srange = [int(x) for x in select_range.split("-")]

        for chrs, vals in tqdm(self.dgtf.items()):
            if len(vals) == 0:
                continue
        
            data = pd.DataFrame(vals)
            data = data.transpose()

            keep_cov = data["mapped"] >= cov
            keep_size = (data["end"] - data["start"]  +1 >= srange[0]) \
                    & (data["end"] - data["start"] +1 <= srange[1])
            keep = keep_cov & keep_size

            data = data[keep]
            data.sort_values(inplace = True, by = "start")
            data = data.transpose()
            data = data.to_dict()

            self.dgtf[chrs] = data

    def join_filt(self, srange = "20-60", cov = 50):

        """
            Use this function after assembly
        """
        srange = [int(x) for x in srange.split("-")]
        for key in tqdm(self.dgtf.keys()):
            
            # Getting the transpose dataframe of dictionary dgtf
            data = pd.DataFrame(self.dgtf[key]).transpose()

            # Filtering by coverage
            data = data[data["mapped"] >= cov]

            # Filtering by sequence size
            keep_size = (data["end"] - data["start"] +1 >= srange[0]) & \
                    (data["end"] - data["start"] +1 <= srange[1])
            data_t = data[keep_size]
            data_f = data[~keep_size]

            # If sequence size is higher or equal than 20 and
            # lower or equal than 60. Create a copy of the coordinates
            # and one copy subtract 110 to start and the other copy
            # sum 110 to end.
            data_t_copy = data_t.copy()
            data_t_copy.index = [f'{x}_b' for x in data_t_copy.index]

            # I going to try the approach I used in join_assembly function.
            data_t.loc[:,"start"] = data_t["start"] - 110
            #data_t["start"] -= 110
            data_t.loc[data_t["start"] <= 0, "start"] = 1
            data_t_copy.loc[:,"end"] = data_t_copy["end"] + 110
            #data_t_copy["end"] += 110
            data_t = pd.concat([data_t, data_t_copy], axis = 0)
            data = pd.concat([data_t, data_f], axis = 0)

            # Filter sequences by size
            data = data[(data["end"] - data["start"] +1 >= 20) & (data["end"] - data["start"] + 1 <= 500)]

            self.dgtf_filt[key] = data.transpose().to_dict()


    def to_gtf(self, gtffile, nmap_reads = 200):
        with open(gtffile, "w") as fh:
            alist = []
            for chromosome, assemblies in self.dgtf.items():
                for assembly, coords in assemblies.items():
                    if coords["mapped"] >= nmap_reads: 
                        alist.append(f'{chromosome}\tsimpleAssembly\texon\t{coords["start"]}\t{coords["end"]}\t100.0\t+\t0\tgene_id "{assembly}"; transcript_id "t{assembly}"; count "{coords["mapped"]}"')
            alist = "\n".join(alist)
            fh.write(alist)

        fh.close()


