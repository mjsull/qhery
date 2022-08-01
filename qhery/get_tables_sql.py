import sqlite3
import sys, os
import subprocess

class covid_drdb:
    def __init__(self, drug_list, database_folder):
        self.drug_list = drug_list
        self.database_folder = database_folder
        self.mutant_synonyms = {

        }


    def download_latest(self):
        url = subprocess.check_output("curl -s https://api.github.com/repos/hivdb/covid-drdb-payload/releases/latest | "
                         "grep -P \"browser_download_url.*covid-drdb-\d{8}.db\" | "
                         "cut -d : -f 2,3 | tr -d \\\"", shell=True)
        url = url.decode()[1:-1]
        sys.stdout.write("Latest covid-drdb is {}.\n".format(url))
        if os.path.exists(os.path.join(self.database_folder, url.split('/')[-1])):
            sys.stdout.write("Latest version already downloaded.\n")
        else:
            subprocess.Popen("wget -P {} {}".format(self.database_folder, url), shell=True).wait()
        self.database = os.path.join(self.database_folder, url.split('/')[-1])

    def get_database(self):
        db_list = []
        for i in os.listdir(self.database_folder):
            if i.startswith("covid-drdb-") and i.endswith(".db"):
                db_list.append(i)
        db_list.sort()
        if db_list == []:
            sys.stderr.write("No databases found in database directory.")
            sys.exit(0)
        self.database = os.path.join(self.database_folder, db_list[-1])

    def connect(self):
        self.con = sqlite3.connect(self.database)

    def get_ref(self):
        self.ref_aa = {}
        for row in self.con.execute("SELECT gene, position, amino_acid FROM ref_amino_acid"):
            gene, position, aa = row
            if not gene in self.ref_aa:
                self.ref_aa[gene] = {}
            self.ref_aa[gene][position] = aa

    def get_epitopes(self):
        self.epitopes = {}
        for i in self.drug_list:
            self.epitopes[i] = []
            for row in self.con.execute('SELECT position FROM antibody_epitopes WHERE ab_name = "{}"'.format(i)):
                self.epitopes[i].append(row[0])


    def get_fold_resistance(self):
        self.resistances = {}
        for i in self.drug_list:
            self.resistances[i] = {}
            for j in self.drug_synonyms[i]:
                for row in self.con.execute('SELECT iso_name, fold_cmp, fold FROM rx_fold WHERE rx_name = "{}"'.format(j)):
                    if row[0] in self.iso_dict:
                        if ":" in self.iso_dict[row[0]]:
                            mutation_name = self.iso_dict[row[0]]
                        else:
                            mutation_name = "S:" + self.iso_dict[row[0]]
                        if mutation_name in self.resistances[i]:
                            self.resistances[i][mutation_name].append([row[1], row[2]])
                        else:
                            self.resistances[i][mutation_name] = [[row[1], row[2]]]


    def add_local_resitances(self):
        dirname = os.path.dirname(__file__)
        data_dir = os.path.join(dirname, 'data')
        with open(os.path.join(data_dir, "resistance_table.tsv")) as f:
            for line in f:
                rx, gene, refaa, pos, mutaa, symbol, fold_change = line.rstrip().split("\t")
                mut_name = gene + ":" + refaa + pos + mutaa
                if not rx in self.resistances:
                    continue
                if not mut_name in self.resistances[rx]:
                    self.resistances[rx][mut_name] = [[symbol, float(fold_change)]]



    def list_resistances(self):
        out_list = []
        for i in self.resistances:
            for j in self.resistances[i]:
                max_resist = 0
                fold_change = ""
                for k in self.resistances[i][j]:
                    if j[1] > max_resist:
                        max_resist = j[1]
                    fold_change += ',' + j[0] + str(j[1])
                out_list.append([i, j, max_resist, fold_change[1:]])
        return(out_list)


    def get_variant_mutations(self, variant):
        self.get_ref()
        self.variant_mutations = {}
        consensus_available = None
        for row in self.con.execute('SELECT var_name FROM variant_synonyms WHERE synonym = "{}"'.format(variant)):
            variant = row[0]
        for row in self.con.execute('SELECT consensus_availability FROM variants WHERE var_name = "{}"'.format(variant)):
            consensus_available = row[0]
        if consensus_available is None:
            sys.stderr.write("Variant not in database.\n")
        elif consensus_available == 0:
            sys.stderr.write("Consensus for variant not available.\n")
        mut_list = []
        for row in self.con.execute('SELECT var_name, gene, position, amino_acid FROM variant_consensus WHERE var_name = "{}"'.format(variant)):
            var, gene, pos, aa = row
            if var == variant:
                mut_list.append([gene, pos, self.ref_aa[gene][pos], aa])
        mut_list.sort()
        lastdel = None
        out_list = []
        for i in mut_list:
            if i[3] == "del" and not lastdel is None and lastdel[0] == i[0] and lastdel[1] + lastdel[4] == i[1]-1:
                lastdel[4] += 1
            else:
                if not lastdel is None:
                    refseq = ""
                    for j in range(lastdel[1], lastdel[1]+lastdel[4]+1):
                        refseq += self.ref_aa[lastdel[0]][j]
                    if lastdel[4] == 0:
                        thepos = lastdel[1]
                    else:
                        thepos = "{}-{}".format(lastdel[1], lastdel[1] + lastdel[4])
                    out_list.append("{}:{}{}∆".format(lastdel[0], refseq, thepos))
                    lastdel = None
                if i[3] == "del":
                    lastdel = i + [0]
                else:
                    out_list.append("{}:{}{}{}".format(i[0], i[2], i[1], i[3]))
        return(out_list)






    def get_single_mutations(self):
        self.iso_dict = {}
        for row in self.con.execute('SELECT iso_name, single_mut_name FROM isolate_mutations_single_s_mut_view'):
            iso_name, single_mut_name = row
            self.iso_dict[iso_name] = single_mut_name

    def get_synonyms(self):
        self.drug_synonyms = {}
        for i in self.drug_list:
            self.drug_synonyms[i] = [i]
            for row in self.con.execute('SELECT abbreviation_name FROM antibodies WHERE ab_name = "{}"'.format(i)):
                self.drug_synonyms[i].append(row[0])
            for row in self.con.execute('SELECT synonym FROM antibody_synonyms WHERE ab_name = "{}"'.format(i)):
                self.drug_synonyms[i].append(row[0])

    def list_rx(self):
        rx_list = set()
        for row in self.con.execute('SELECT ab_name FROM antibodies'):
            rx_list.add(row[0])
        dirname = os.path.dirname(__file__)
        data_dir = os.path.join(dirname, 'data')
        with open(os.path.join(data_dir, "resistance_table.tsv")) as f:
            for line in f:
                rx, gene, refaa, pos, mutaa, symbol, fold_change = line.split("\t")
                rx_list.add(rx)
        rx_list = list(rx_list)
        rx_list.sort()
        for i in rx_list:
            sys.stdout.write(i + "\n")

