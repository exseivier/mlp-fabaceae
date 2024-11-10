from sys import argv
import pandas as pd

def blast2gtf():
    """

    """
    if len(argv) != 3:
        print("ArgsError: One or some arguments are missing!")
        return False
    input_file = argv[1]
    output_file = argv[2]
    df = pd.read_csv(input_file, sep = "\t", header = None)
    df = df[(df.iloc[:,6] == 100) & (df.iloc[:, 7] == 100)]
    df.reset_index(inplace = True, drop = True)
    strand = ["+" if x == "plus" else "-" for x in df[3]]
    df[3] = strand
    df[4] = 0
    array_out = []
    counter = 0

    with open(output_file, "w") as fout:

        for i in df.index:
            chrname = df.iloc[i, 0]
            prog = "blast2gtf"
            feature = "miRNA_stemloop"
            start = df.iloc[i, 1]
            end = df.iloc[i, 2]
            score = df.iloc[i, 6]
            strand = df.iloc[i, 3]
            frame = df.iloc[i, 4]
            desc = f'gene_id "{df.iloc[i, 5]}_{counter}"; transcript_id "t{df.iloc[i, 5]}_{counter}"'
            array_out.append(f'{chrname}\t{prog}\t{feature}\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{desc}')
        
            counter += 1

        array_out = "\n".join(array_out)
        array_out = f'{array_out}\n'
        fout.write(array_out)

    return True

if __name__ == "__main__":
    if blast2gtf():
        print("Success!")
    else:
        print("Failure!")

