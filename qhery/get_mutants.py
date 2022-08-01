import sys
import os
import subprocess
import codecs

class mutantFinder:
    def __init__(self, vcf_file, working_dir, sample_name):
        self.vcf_file = vcf_file
        self.sample_name = sample_name
        self.working_dir = working_dir
        dirname = os.path.dirname(__file__)
        self.data_dir = os.path.join(dirname, 'data')
        self.bcftools_vcf = os.path.join(working_dir, "{}.bcftools.vcf".format(self.sample_name))
        self.csq_file = os.path.join(self.working_dir, "{}.csq.tsv".format(self.sample_name))
        self.populate_mature_proteins()

    def populate_mature_proteins(self):
        self.mature_protein_coords = {}
        self.mature_proteins_lookup = {}
        with open(os.path.join(self.data_dir, "nsp_coords.tsv")) as f:
            for line in f:
                parent, name, start, stop = line.split()
                start, stop = int(start), int(stop)
                if not parent in self.mature_protein_coords:
                    self.mature_protein_coords[parent] = {}
                    self.mature_proteins_lookup[parent] = {}
                self.mature_protein_coords[parent][name] = (start, stop)
        offset = {}
        with open(os.path.join(self.data_dir, "Sars_cov_2.ASM985889v3.101.gff3")) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                else:
                    ref, method, feature, start, stop, frame, strand, qual, extra = line.rstrip().split("\t")
                    if feature == "mRNA":
                        for i in extra.split(';'):
                            if i.startswith("transcript_id="):
                                offset[i.split('=')[1]] = int(start)
        for i in self.mature_protein_coords:
            for j in self.mature_protein_coords[i]:
                start, stop = self.mature_protein_coords[i][j]
                start = (start-offset[i]) //3 + 1
                stop = (stop-offset[i]) //3 + 1
                for num1, num2, in enumerate(range(start, stop+1)):
                    self.mature_proteins_lookup[i][num2] = [j, num1+1]





    def convert_vcf(self, vcf_in=None, vcf_out=None):
        if vcf_in is None:
            vcf_in = self.vcf_file
        if vcf_out is None:
            vcf_out = self.bcftools_vcf
        with open(vcf_in) as f, open(vcf_out, 'w') as o:
            for line in f:
                if line.startswith("##source"):
                    o.write(line)
                    o.write('##contig=<ID=MN908947.3,length=29903,assembly=MN908947.3,species="SARS-CoV-2",taxonomy=x>\n')
                elif line.startswith("##"):
                    o.write(line)
                else:
                    o.write("\t".join(line.rstrip().split("\t")[:8]) + "\n")

    def run_bcf_csq(self, vcf_in=None, csq_out=None):
        if vcf_in is None:
            vcf_in = self.bcftools_vcf
        if csq_out is None:
            csq_out = self.csq_file
        subprocess.Popen("bcftools csq -f {}/nCoV-2019.reference.fasta -g {}/Sars_cov_2.ASM985889v3.101.gff3  {} -O t -o {}".format(
            self.data_dir, self.data_dir, vcf_in, csq_out), shell=True).wait()




    def parse_csq(self, csq_file=None):
        mut_list = []
        if csq_file is None:
            csq_file = self.csq_file
        with codecs.open(csq_file, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if not line.startswith("CSQ\t"):
                    continue
                csq, sample, haplo, chrom, pos, consequence = line.rstrip().split("\t")
                mut_type, gene, acc, gene_type, strand, prot_mut, nuc_mut = consequence.split("|")
                if mut_type.startswith("*"):
                    mut_type = mut_type[1:]
                if ">" in prot_mut:
                    ref, alt = prot_mut.split(">")
                else:
                    ref = prot_mut
                    alt = prot_mut
                ref_pos, ref_aa = "", ""
                for i in ref:
                    if i.isdigit():
                        ref_pos += i
                    else:
                        ref_aa += i
                ref_pos = int(ref_pos)
                alt_aa = ""
                for i in alt:
                    if not i.isdigit():
                        alt_aa += i
                if gene == "ORF1ab":
                    gene, ref_pos = self.mature_proteins_lookup[acc][ref_pos]
                if mut_type == "synonymous":
                    pass
                elif mut_type == "missense" and len(ref_aa) == 1 and len(alt_aa) == 1:
                    mut = gene + ":" + ref_aa + str(ref_pos) + alt_aa
                    mut_list.append(mut)
                elif mut_type == "missense" and len(ref_aa) == len(alt_aa):
                    for i in range(len(ref_aa)):
                        mut = gene + ':' + ref_aa[i] + str(ref_pos+i) + alt_aa[i]
                        mut_list.append(mut)
                elif mut_type == "inframe_deletion":
                    for i in range(len(alt_aa)):
                        if alt_aa[i] != ref_aa[i]:
                            mut = gene + ':' + ref_aa[i] + str(ref_pos+i) + alt_aa[i]
                            mut_list.append(mut)
                    del_start = len(alt_aa)
                    if len(ref_aa) == del_start + 1:
                        mut = gene + ':' + ref_aa[del_start:] + str(ref_pos + del_start) + "∆"
                    else:
                        mut = gene + ':' + ref_aa[del_start:] + str(ref_pos + del_start) + '-' + str(ref_pos + len(ref_aa)-1) + "∆"
                    mut_list.append(mut)
                elif mut_type in "inframe_insertion":
                    mut = gene + ':' + ref_aa + str(ref_pos) + alt_aa[0] + '_' + alt_aa[1:]
                    mut_list.append(mut)
                elif mut_type == "frameshift":
                    pass
                else:
                    sys.stderr.write("{}\nI haven't seen this mutation before, what do I do?!?!?!?!?!!?!!!  ...ignoring....\n".format(mut_type))
        mut_list = set(mut_list)
        mut_list = list(mut_list)
        return(mut_list)

    def run_lofreq(self, bam_file, min_cov=20):
        self.indel_qual_bam = os.path.join(self.working_dir, self.sample_name + ".indelbq.bam")
        self.lofreq_vcf = os.path.join(self.working_dir, self.sample_name + ".lofreq.vcf")
        subprocess.Popen("lofreq indelqual -u 5 {} > {}".format(bam_file, self.indel_qual_bam), shell=True).wait()
        subprocess.Popen("lofreq call --call-indels -C 20 -f {}/nCoV-2019.reference.fasta {} > {}".format(self.data_dir, self.indel_qual_bam, self.lofreq_vcf), shell=True).wait()



    def recover_low_freq(self, bam_file):
        self.run_lofreq(bam_file)
        converted_vcf = os.path.join(self.working_dir, self.sample_name + ".lofreq.csq.vcf")
        lofreq_csq = os.path.join(self.working_dir, self.sample_name + ".lofreq.csq.tsv")
        self.convert_vcf(vcf_in=self.lofreq_vcf, vcf_out=converted_vcf)
        self.run_bcf_csq(vcf_in=converted_vcf, csq_out=lofreq_csq)
        return(self.parse_csq(csq_file=lofreq_csq))


