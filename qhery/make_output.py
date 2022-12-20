import subprocess
import os
import pysam
from collections import defaultdict


def get_mutant_name_ref(mut):
    if "_" in mut:
        return mut.split("_")[0][:-1] + "ins"
    else:
        return mut

def get_fasta_depth(fasta_coverage, position):
    for i in range(position, position + 3):
        if fasta_coverage[i] == "n":
            return("not_in_consensus")
    return("in_consensus")


def make_final_tables(
    sample_muts,
    resistance_muts,
    variant_muts,
    epitopes,
    pipeline_folder,
    sample_name,
    bam_file,
    aa_to_nuc_dict,
    fasta_coverage,
    min_fold_reduction_to_report=2,
):
    with open(os.path.join(pipeline_folder, sample_name + ".full.tsv"), "w") as full_table, open(
        os.path.join(pipeline_folder, sample_name + ".final.tsv"), "w"
    ) as final_table:
        headers = [
            "Mutation",
            "alt_names",
            "in_sample",
            "in_variant",
            "covered",
            "resistance_mutation",
            "codon_present",
            "codons",
            "codon_depth",
        ]
        min_depth = 20
        drug_list = list(resistance_muts)
        drug_list.sort()
        for i in drug_list:
            headers.append("{}_average_fold_reduction".format(i))
            headers.append("{}_fold_reductions".format(i))
            headers.append("{}_in_epitope".format(i))
        full_table.write("\t".join(headers) + "\n")
        final_table.write("\t".join(headers) + "\n")
        all_mutations = set()
        alt_mutation_names = {}
        for i in sample_muts:
            gene, mutation = i.split(":")
            if "_" in mutation:
                name = gene + ":" + mutation.split("_")[0][:-1] + "ins"
                all_mutations.add(name)
                if name in alt_mutation_names:
                    alt_mutation_names[name].append(i)
                else:
                    alt_mutation_names[name] = [i]
            else:
                all_mutations.add(i)
        resistance_genes = set()
        for i in resistance_muts:
            for j in resistance_muts[i]:
                gene, mutation = j.split(":")
                resistance_genes.add(gene)
                all_mutations.add(j)
        all_mutations = list(all_mutations)
        all_mutations.sort(
            key=lambda x: (
                x.split(":")[0],
                int("0" + "".join([n for n in x.split("-")[0] if n.isdigit()])),
            )
        )
        for i in all_mutations:
            columns = [i]
            if i in alt_mutation_names:
                columns.append(",".join(alt_mutation_names[i]))
            else:
                columns.append("-")
            if i in alt_mutation_names or i in sample_muts:
                columns.append("True")
            else:
                columns.append("False")
            if i in variant_muts:
                columns.append("True")
            else:
                columns.append("False")
            gene, mutation = i.split(":")
            if "-" in mutation:
                start = int("".join([x for x in mutation.split("-")[0] if x.isdigit()]))
                stop = int("".join([x for x in mutation.split("-")[1] if x.isdigit()]))
            else:
                start = int("".join([x for x in mutation if x.isdigit()]))
                stop = start
            if not bam_file is None:
                codon_usage, codon_depth = get_codons(bam_file, aa_to_nuc_dict[gene][start])
            elif not fasta_coverage is None:
                codon_usage = ""
                codon_depth = get_fasta_depth(fasta_coverage, aa_to_nuc_dict[gene][start])
            else:
                codon_usage, codon_depth = "", "NA"
            if codon_depth == "NA":
                coverage_status = "NA"
            elif codon_depth == "not_in_consensus":
                coverage_status = "False"
            elif codon_depth == "in_consensus" or codon_depth >= min_depth:
                coverage_status = "True"
            else:
                coverage_status = "False"
            columns.append(coverage_status)
            mutant_base = ""
            for j in i:
                if j.isdigit():
                    mutant_base = ""
                else:
                    mutant_base += j
            codon_freq = 0
            for j in codon_usage.split("("):
                if j.split(")")[0] == mutant_base:
                    codon_freq = float(j.split(":")[1].split(";")[0][:-1])
                    break
            rx_columns = []
            above_min_fold_reduction = "False"
            for j in drug_list:
                if i in resistance_muts[j]:
                    tot_res = 0
                    col = ""
                    for k in resistance_muts[j][i]:
                        tot_res += k[1]
                        col += "," + k[0] + str(k[1])
                    col = col[1:]
                    rx_columns += [
                        "{:.2f}".format(tot_res / len(resistance_muts[j][i])),
                        col,
                    ]
                    if tot_res / len(resistance_muts[j][i]) >= min_fold_reduction_to_report:
                        above_min_fold_reduction = "True"
                else:
                    rx_columns += ["0", "-"]
                in_epitope = "False"
                if gene == "S":
                    for num in range(start, stop + 1):
                        if num in epitopes[j]:
                            in_epitope = "True"
                            break
                rx_columns.append(in_epitope)
            columns.append(above_min_fold_reduction)
            if codon_depth == "NA" or codon_depth == "in_consensus" or codon_depth == "not_in_consensus":
                codon_present = "NA"
            elif codon_freq > 5 and codon_depth >= min_depth:
                codon_present = "True"
            else:
                codon_present = "False"
            if codon_usage == "":
                codon_usage = "NA"
            columns += [codon_present, codon_usage, str(codon_depth)]
            columns += rx_columns
            outline = "\t".join(columns) + "\n"
            full_table.write(outline)
            if columns[2] == "True" and columns[3] == "False" and columns[0].split(":")[0] in resistance_genes:
                final_table.write(outline)
            elif columns[4] == "False" and columns[5] == "True":
                final_table.write(outline)


def get_nuc_aa_translations():
    orf_coords = {
        "nsp1": (266, 805),
        "nsp2": (806, 2719),
        "PLpro": (2720, 8554),
        "nsp4": (8555, 10054),
        "_3CLpro": (10055, 10972),  # 3C
        "nsp6": (10973, 11842),  #
        "nsp7": (11843, 12091),
        "nsp8": (12092, 12685),
        "nsp9": (12686, 13024),
        "nsp10": (13025, 13441),
        "RdRP": (13442, 13483),  # rdrp
        "nsp13": (16237, 18039),
        "nsp14": (18040, 19620),
        "nsp15": (19621, 20658),
        "nsp16": (20659, 21552),
    }
    nuc_to_aa_dict = defaultdict(lambda: set())
    aa_to_nuc_dict = {}
    for i in orf_coords:
        aa_to_nuc_dict[i] = {}
        for num in range(orf_coords[i][0], orf_coords[i][1] + 1):
            aa_pos = (num - orf_coords[i][0]) // 3 + 1
            if not aa_pos in aa_to_nuc_dict[i] or num - 1 < aa_to_nuc_dict[i][aa_pos]:
                aa_to_nuc_dict[i][aa_pos] = num - 1
            nuc_to_aa_dict[num].add((i, aa_pos))
    for num in range(13483, 16237):
        aa_pos = (num - 13483) // 3 + 15
        if not aa_pos in aa_to_nuc_dict["RdRP"] or num - 1 < aa_to_nuc_dict["RdRP"][aa_pos]:
            aa_to_nuc_dict["RdRP"][aa_pos] = num - 1
        nuc_to_aa_dict[num].add(("RdRP", aa_pos))
    with open("{}/Sars_cov_2.ASM985889v3.101.gff3".format(os.path.join(os.path.dirname(__file__), "data"))) as f:
        for line in f:
            if line.startswith("#"):
                continue
            (
                ref,
                method,
                feature,
                start,
                stop,
                frame,
                strand,
                qual,
                details,
            ) = line.rstrip().split("\t")
            start, stop = int(start), int(stop)
            if feature == "gene":
                for i in details.split(";"):
                    if i.startswith("Name="):
                        gene = i.split("=")[1]
                if not gene in aa_to_nuc_dict:
                    aa_to_nuc_dict[gene] = {}
                for num in range(start, stop + 1):
                    aa_pos = (num - start) // 3 + 1
                    if not aa_pos in aa_to_nuc_dict[gene] or num - 1 < aa_to_nuc_dict[gene][aa_pos]:
                        aa_to_nuc_dict[gene][aa_pos] = num - 1
                    nuc_to_aa_dict[num].add((gene, aa_pos))
    return (nuc_to_aa_dict, aa_to_nuc_dict)


def make_epitope_graphs(bam, epitopes, pipeline_folder, sample_name):
    SPIKE_PROTEIN_START = 21563
    for i in epitopes:
        positions = []
        for pos in epitopes[i]:
            positions.append(str(SPIKE_PROTEIN_START + pos * 3 - 3))
            positions.append(str(SPIKE_PROTEIN_START + pos * 3 - 2))
            positions.append(str(SPIKE_PROTEIN_START + pos * 3 - 1))
        if positions != []:
            cmd = "bammix -b {} -p {} -o {} ".format(
                bam,
                " ".join(positions),
                os.path.join(pipeline_folder, sample_name + "." + i),
            )
            subprocess.Popen(cmd, shell=True).wait()


def get_aa(codon):
    codon2aa = {
        "AAA": "K",
        "AAC": "N",
        "AAG": "K",
        "AAT": "N",
        "ACA": "T",
        "ACC": "T",
        "ACG": "T",
        "ACT": "T",
        "AGA": "R",
        "AGC": "S",
        "AGG": "R",
        "AGT": "S",
        "ATA": "I",
        "ATC": "I",
        "ATG": "M",
        "ATT": "I",
        "CAA": "Q",
        "CAC": "H",
        "CAG": "Q",
        "CAT": "H",
        "CCA": "P",
        "CCC": "P",
        "CCG": "P",
        "CCT": "P",
        "CGA": "R",
        "CGC": "R",
        "CGG": "R",
        "CGT": "R",
        "CTA": "L",
        "CTC": "L",
        "CTG": "L",
        "CTT": "L",
        "GAA": "E",
        "GAC": "D",
        "GAG": "E",
        "GAT": "D",
        "GCA": "A",
        "GCC": "A",
        "GCG": "A",
        "GCT": "A",
        "GGA": "G",
        "GGC": "G",
        "GGG": "G",
        "GGT": "G",
        "GTA": "V",
        "GTC": "V",
        "GTG": "V",
        "GTT": "V",
        "TAA": "_",
        "TAC": "Y",
        "TAG": "_",
        "TAT": "T",
        "TCA": "S",
        "TCC": "S",
        "TCG": "S",
        "TCT": "S",
        "TGA": "_",
        "TGC": "C",
        "TGG": "W",
        "TGT": "C",
        "TTA": "L",
        "TTC": "F",
        "TTG": "L",
        "TTT": "F",
    }
    if codon.upper() in codon2aa:
        return codon2aa[codon.upper()]
    else:
        return "X"


def get_codons(bam_file, start_position):
    samfile = pysam.AlignmentFile(bam_file, "rb")
    codonfreq = defaultdict(lambda: 0)
    for pileupcolumn in samfile.pileup("MN908947.3", start_position):
        if pileupcolumn.pos != start_position:
            continue
        depth = 0
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                try:
                    codon = (
                        pileupread.alignment.query_sequence[pileupread.query_position]
                        + pileupread.alignment.query_sequence[pileupread.query_position + 1]
                        + pileupread.alignment.query_sequence[pileupread.query_position + 2]
                    )
                    codonfreq[codon] += 1
                    depth += 1
                except IndexError:
                    pass
            elif pileupread.is_del:
                codonfreq["del"] += 1
                depth += 1
        outstring = ""
        codonfreqlist = list(codonfreq)
        codonfreqlist.sort(key=lambda x: codonfreq[x], reverse=True)
        for i in codonfreq:
            if i == "del":
                outstring += "{}:{:.1%};".format(i, codonfreq[i] / depth)
            else:
                outstring += "{}({}):{:.1%};".format(i, get_aa(i), codonfreq[i] / depth)
            depth += codonfreq[i]
        return [outstring, depth]