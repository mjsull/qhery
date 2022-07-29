import sys
import glob
import os

qld_ids = {}
for table in glob.glob("/data/Sars-Cov-2/FSS/sequence*/*.csv"):
    with open(table) as f:
        f.readline()
        for line in f:
            labno, qldid = line.split(',')[:2]
            qld_ids[labno.replace('-', '')] = qldid

b_comments = {
    "Sotrovimab":[
        "S:P337L",
        "S:P337R",
        "S:P337H",
        "S:P337T",
        "S:E340K",
        "S:E340A",
        "S:E340V",
        "S:E340G",
        "S:3S71L",
        "S:S371F"
    ],
    "Evusheld":[
        "S:E484A",
        "S:F486V",
        "S:Q498R",
        "S:E990A"
    ],
    "Remdesivir":[
        "RdRP:E802A",
        "RdRP:E802D"
    ]
}



gene_translate = {
    "S": "Spike",
    "RdRP": "RdRp",
    "3CL":"3CL"
}

rows = [["Spike", "Sotrovimab"],["", "Evusheld"],["RdRp", "Remdesivir"], ["3CL", "Paxlovid"]]





for table in glob.glob(os.path.join(sys.argv[1], "*.final.tsv")):
    mut_dict = {}
    drug_comments = {
        "Sotrovimab": "A",
        "Evusheld": "A",
        "Remdesivir": "A",
        "Paxlovid": "A"
    }

    with open(table) as f:
        headers = f.readline().rstrip().split("\t")
        resistance_headers = headers[9:]
        muts = []
        mut_dict = {}
        for line in f:
            mutation, alt_names, in_sample, in_variant, covered, resistance_mutation, in_codon, codons, codon_depth = line.split("\t")[:9]
            crosses = set()
            if covered == "True" and (in_sample == "True" or in_codon == "True") and in_variant == "False":
                for i in b_comments:
                    if mutation in b_comments[i]:
                        drug_comments[i] = "B"
                gene, mut = mutation.split(":")
                resistance_info = line.split("\t")[9:]
                mutant_base = ""
                for j in mutation:
                    if j.isdigit():
                        mutant_base = ""
                    else:
                        mutant_base += j
                codon_freq = 0
                for j in codons.split("("):
                    if j.split(")")[0] == mutant_base:
                        codon_freq = float(j.split(":")[1].split(';')[0][:-1])
                        break
                if codon_freq < 90:
                    mut += "*"
                for num, k in enumerate(resistance_info):
                    if resistance_headers[num] == "Sotrovimab_in_epitope" and k == "True":
                        crosses.add("†")
                    elif resistance_headers[num] == "Sotrovimab_average_fold_reduction" and float(k) >= 2:
                        crosses.add("††")
                    elif resistance_headers[num] in ["Cilgavimab_in_epitope", "Tixagevimab_in_epitope"] and k == "True":
                        crosses.add("†††")
                    elif resistance_headers[num] in ["Cilgavimab_average_fold_reduction", "Tixagevimab_average_fold_reduction"] and float(k) >= 2:
                        crosses.add("††††")
                crosses = list(crosses)
                crosses.sort()
                mut += " ".join(crosses)
                if not gene_translate[gene] in mut_dict:
                    mut_dict[gene_translate[gene]] = []
                mut_dict[gene_translate[gene]].append(mut)

            if covered == "False":
                for i in b_comments:
                    if mutation in b_comments[i]:
                        if drug_comments[i] == "A":
                            drug_comments[i] = "D"
        sample_name = qld_ids[os.path.basename(table).split('.')[0]]
        for i in rows:
            if i[0] in mut_dict:
                muts = mut_dict[i[0]]
            else:
                muts = ""
            sys.stdout.write("\t".join([sample_name, i[0], ', '.join(muts), i[1], drug_comments[i[1]]]) + "\n")
            sample_name = ''









