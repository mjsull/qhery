import get_mutants
import get_tables_sql
import make_output
import argparse
import sys
import os








parser = argparse.ArgumentParser(prog='covresid')
subparsers = parser.add_subparsers(dest="subparser_name")

list_parser = subparsers.add_parser("list_rx", help="Print list of treatements with recorded mutations and exit.")
list_parser.add_argument("--database_dir", "-d", help="Directory with latest Stanford resistance database.")

run_parser = subparsers.add_parser("run", help="Run CoViD resistance identifier.")
run_parser.add_argument("--sample_name", "-n", help="Sample name.")
run_parser.add_argument("--vcf", "-v", help="vcf file")
run_parser.add_argument("--bam", "-b", help="bam file")
run_parser.add_argument("--database_dir", "-d", help="Directory with latest Stanford resistance database.")
run_parser.add_argument("--pipeline_dir", "-p", help="Pipeline to run program in.")
run_parser.add_argument("--lineage", "-l", help="Lineage report of variants.")
run_parser.add_argument("--rx_list", "-rx", nargs="+", help="List of drugs to analyze.")
run_parser.add_argument("--fasta", "-f", help="Consensus fasta.")

mut_parser = subparsers.add_parser("mutations", help="Run CoViD resistance identifier.")
mut_parser.add_argument("--vcf", "-v", help="List of VCF files")
mut_parser.add_argument("--database_dir", "-d", help="Directory with latest Stanford resistance database.")
mut_parser.add_argument("--pipeline_dir", "-p", help="Pipeline to run program in.")
mut_parser.add_argument("--lineage", "-l", help="Lineage report of variants.")
mut_parser.add_argument("--keep_lineage", "-k", help="Report lineage mutations as well.")
mut_parser.add_argument("--sample_name", "-n", help="Sample name.")


args = parser.parse_args()

if args.subparser_name == "list_rx":
    gt = get_tables_sql.covid_drdb([], args.database_dir)
    gt.get_database()
    gt.connect()
    sys.stdout.write("The following MABs have escape information.")
    gt.list_rx()

elif args.subparser_name == "mutations":
    gt = get_tables_sql.covid_drdb([], args.database_dir)
    gt.get_database()
    gt.connect()
    if not args.lineage is None:
        mut_list_var = gt.get_variant_mutations(args.lineage)
    else:
        mut_list_var = []
    mf = get_mutants.mutantFinder(args.vcf, args.pipeline_dir, args.sample_name)
    mf.convert_vcf()
    mf.run_bcf_csq()
    mut_list_sample = mf.parse_csq()
    all_muts = list(set(mut_list_var + mut_list_sample))
    all_muts.sort(key=lambda x: (x.split(':')[0].split('-')[0],  int(''.join([n for n in x if n.isdigit()]))))
    with open(os.path.join(args.pipeline_dir, "{}.mutations.txt".format(args.sample_name)), "w") as o:
        for i in all_muts:
            if not i in mut_list_var or not i in mut_list_sample:
                o.write("{}\n".format(i))

#    print(set(mut_list_var) - set(mut_list_sample))
#    print(set(mut_list_sample) - set(mut_list_var))
    

elif args.subparser_name == "run":
    if os.path.exists(args.pipeline_dir) and os.path.isdir(args.pipeline_dir):
        pass
    elif os.path.exists(args.pipeline_dir):
        sys.stderr.write("Pipeline dir exists and is not a directory.\n")
        sys.exit(1)
    else:
        os.makedirs(args.pipeline_dir)
    gt = get_tables_sql.covid_drdb(args.rx_list, args.database_dir)
    gt.download_latest()
    gt.connect()
    gt.get_epitopes()
    make_output.make_epitope_graphs(args.bam, gt.epitopes, args.pipeline_dir, args.sample_name)
    gt.get_synonyms()
    gt.get_single_mutations()
    gt.get_fold_resistance()
    gt.add_local_resitances()
    mut_list_var = gt.get_variant_mutations(args.lineage)
    mf = get_mutants.mutantFinder(args.vcf, args.pipeline_dir, args.sample_name)
    mf.convert_vcf()
    mf.run_bcf_csq()
    mut_list_sample = mf.parse_csq()
    if not args.bam is None:
        uncovered = make_output.get_coverage_genes(args.bam, args.pipeline_dir, args.sample_name)
        mut_list_lofreq = mf.recover_low_freq(args.bam)
    mut_list_sample = list(set(mut_list_lofreq + mut_list_sample))
    make_output.make_final_tables(mut_list_sample, gt.resistances, mut_list_var, gt.epitopes, uncovered, args.pipeline_dir, args.sample_name)
    make_output.make_alignment_files(args.fasta, args.pipeline_dir, args.sample_name)
    make_output.make_epitope_graphs(args.bam, gt.epitopes, args.pipeline_dir, args.sample_name)





