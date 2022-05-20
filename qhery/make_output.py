import subprocess
import os
from collections import defaultdict


def get_mutant_name_ref(mut):
    if "_" in mut:
        return(mut.split("_")[0][:-1] + "ins")
    else:
        return(mut)



def make_final_tables(sample_muts, resistance_muts, variant_muts, epitopes, uncovered, pipeline_folder, sample_name, min_fold_reduction_to_report=2):
    with open(os.path.join(pipeline_folder, sample_name + '.full.tsv'), 'w') as full_table, \
        open(os.path.join(pipeline_folder, sample_name + '.final.tsv'), 'w') as final_table:
        headers = ["Mutation", "alt_names", "in_sample", "in_variant", "covered", "resistance_mutation"]
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
            key=lambda x: (x.split(':')[0], int('0' + ''.join([n for n in x.split('-')[0] if n.isdigit()]))))
        for i in all_mutations:
            columns = [i]
            if i in alt_mutation_names:
                columns.append(','.join(alt_mutation_names[i]))
            else:
                columns.append('-')
            if i in alt_mutation_names or i in sample_muts:
                columns.append("True")
            else:
                columns.append("False")
            if i in variant_muts:
                columns.append("True")
            else:
                columns.append("False")
            gene, mutation = i.split(':')
            coverage_status = "True"
            if "-" in mutation:
                start = int(''.join([x for x in mutation.split('-')[0] if x.isdigit()]))
                stop = int(''.join([x for x in mutation.split('-')[1] if x.isdigit()]))
            else:
                start = int(''.join([x for x in mutation if x.isdigit()]))
                stop = start
            for pos in range(start, stop+1):
                if uncovered is None:
                    coverage_status = "Unknown"
                elif (gene, pos) in uncovered:
                    coverage_status = "False"
                    break
            columns.append(coverage_status)
            rx_columns = []
            above_min_fold_reduction = "False"
            for j in drug_list:
                if i in resistance_muts[j]:
                    tot_res = 0
                    col = ''
                    for k in resistance_muts[j][i]:
                        tot_res += k[1]
                        col += ',' + k[0] + str(k[1])
                    col = col[1:]
                    rx_columns += ["{:.2f}".format(tot_res/len(resistance_muts[j][i])), col]
                    if tot_res/len(resistance_muts[j][i]) >= min_fold_reduction_to_report:
                        above_min_fold_reduction = "True"
                else:
                    rx_columns += ["0", "-"]
                in_epitope = "False"
                if gene == "S":
                    for num in range(start, stop+1):
                        if num in epitopes[j]:
                            in_epitope = "True"
                            break
                rx_columns.append(in_epitope)
            columns.append(above_min_fold_reduction)
            columns += rx_columns
            outline = "\t".join(columns) + "\n"
            full_table.write(outline)
            if columns[2] == "True" and columns[3] == "False" and columns[0].split(":")[0] in resistance_genes:
                final_table.write(outline)
            elif columns[4] == "False" and columns[5] == "True":
                final_table.write(outline)


def get_coverage_genes(bam_file, pipeline_folder, sample_name, min_cov=20):
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
        "nsp16": (20659, 21552)
    }
    coord_dict = defaultdict(lambda: set())
    for i in orf_coords:
        for num in range(orf_coords[i][0], orf_coords[i][1]+1):
            coord_dict[num].add((i, (num - orf_coords[i][0]) //3+1))
    for num in range(13483, 16237):
        coord_dict[num].add(("RdRP",  (num - 13483) //3+15))
    with open("{}/Sars_cov_2.ASM985889v3.101.gff3".format(os.path.join(os.path.dirname(__file__), 'data'))) as f:
        for line in f:
            if line.startswith("#"):
                continue
            ref, method, feature, start, stop, frame, strand, qual, details = line.rstrip().split("\t")
            start, stop = int(start), int(stop)
            if feature == "gene":
                for i in details.split(';'):
                    if i.startswith("Name="):
                        gene = i.split('=')[1]
                for num in range(start, stop + 1):
                    coord_dict[num].add((gene, (num - start) // 3 + 1))
    cov_file = os.path.join(pipeline_folder, sample_name + '.cov')
    subprocess.Popen("samtools depth -aa {} > {}".format(bam_file, cov_file), shell=True).wait()
    uncovered_aa = set()
    with open(cov_file) as f:
        for line in f:
            ref, pos, depth = line.split()
            if int(depth) < min_cov:
                uncovered_aa = uncovered_aa.union(coord_dict[int(pos)])
    return(uncovered_aa)



def make_alignment_files(fasta, pipeline_folder, sample_name):
    dirname = os.path.dirname(__file__)
    data_dir = os.path.join(dirname, 'data')
    subprocess.Popen("blastx -query {} -subject {} -out {}".format(
        fasta, os.path.join(data_dir, "proteins.faa"), os.path.join(pipeline_folder, sample_name + '.blastx_alignment.txt')),
        shell=True).wait()

def make_epitope_graphs(bam, epitopes, pipeline_folder, sample_name):
    SPIKE_PROTEIN_START = 21563
    for i in epitopes:
        positions = []
        for pos in epitopes[i]:
            positions.append(str(SPIKE_PROTEIN_START+pos*3-3))
            positions.append(str(SPIKE_PROTEIN_START+pos*3-2))
            positions.append(str(SPIKE_PROTEIN_START+pos*3-1))
        if positions != []:
            cmd = "bammix -b {} -p {} -o {} ".format(bam, " ".join(positions), os.path.join(pipeline_folder, sample_name + '.' + i))
            subprocess.Popen(cmd, shell=True).wait()