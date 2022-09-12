#!/usr/bin/env python
import argparse
import sys
import os

try:
    import qhery.get_mutants as get_mutants
    import qhery.get_tables_sql as get_tables_sql
    import qhery.make_output as make_output
except ModuleNotFoundError:
    import get_mutants
    import get_tables_sql
    import make_output



def main(args=None):
    parser = argparse.ArgumentParser(prog="qhery")
    subparsers = parser.add_subparsers(dest="subparser_name")

    list_parser = subparsers.add_parser(
        "list_rx",
        help="Print list of treatements with recorded mutations and exit.",
    )
    list_parser.add_argument(
        "--database_dir",
        "-d",
        help="Directory with latest Stanford resistance database.",
    )

    run_parser = subparsers.add_parser("run", help="Run CoViD resistance identifier.")
    run_parser.add_argument("--sample_name", "-n", required=True, help="Sample name.")
    run_parser.add_argument("--vcf", "-v", help="vcf file")
    run_parser.add_argument("--bam", "-b", help="bam file")
    run_parser.add_argument(
        "--database_dir",
        "-d",
        required=True,
        help="Directory with latest Stanford resistance database.",
    )
    run_parser.add_argument(
        "--download",
        help="Download the latest database.",
        action="store_true",
        default=False,
    )
    run_parser.add_argument("--pipeline_dir", "-p", required=True, help="Pipeline to run program in.")
    run_parser.add_argument("--lineage", "-l", help="Lineage report of variants.")
    run_parser.add_argument("--rx_list", "-rx", nargs="+", help="List of drugs to analyze.")
    run_parser.add_argument("--fasta", "-f", help="Consensus fasta.")
    run_parser.add_argument("--download_nextclade_data", "-dn", action="store_true", help="download nextclade data.")
    run_parser.add_argument("--nextclade_data", "-nd", help="directory of sars-cov-2 nextclade data, "
                            "can be downloaded with "
                            "\"nextclade dataset get --name 'sars-cov-2' --output-dir 'data/sars-cov-2'\"")

    mut_parser = subparsers.add_parser("mutations", help="List mutations without resistance information.")
    mut_parser.add_argument("--vcf", "-v", help="List of VCF files")
    mut_parser.add_argument(
        "--database_dir",
        "-d",
        help="Directory with latest Stanford resistance database.",
    )
    mut_parser.add_argument("--pipeline_dir", "-p", help="Pipeline to run program in.")
    mut_parser.add_argument("--lineage", "-l", help="Lineage report of variants.")
    mut_parser.add_argument("--sample_name", "-n", help="Sample name.")
    mut_parser.add_argument("--bam", "-b", help="bam file")

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
        mf = get_mutants.mutantFinder(args.pipeline_dir, args.sample_name)
        mf.convert_vcf(args.vcf)
        mf.run_bcf_csq()
        mut_list_sample = mf.parse_csq()
        if not args.bam is None:
            if not os.path.exists(args.bam + '.bai'):
                sys.stderr.write("Please sort and index your bam file before running ")
                sys.exit(1)
            mut_list_lofreq = mf.recover_low_freq(args.bam)
        else:
            mut_list_lofreq = []
        all_muts = list(set(mut_list_sample + mut_list_lofreq))
        all_muts.sort(
            key=lambda x: (
                x.split(":")[0].split("-")[0],
                int("".join([n for n in x if n.isdigit()])),
            )
        )
        with open(
            os.path.join(args.pipeline_dir, "{}.mutations.txt".format(args.sample_name)),
            "w",
        ) as o:
            for i in all_muts:
                if "_" in i:
                    mut = i.split("_")[0][:-1] + "ins"
                else:
                    mut = i
                if not mut in mut_list_var:
                    o.write("{}\n".format(i))

    elif args.subparser_name == "run":
        if os.path.exists(args.pipeline_dir) and os.path.isdir(args.pipeline_dir):
            pass
        elif os.path.exists(args.pipeline_dir):
            sys.stderr.write("Pipeline dir exists and is not a directory.\n")
            sys.exit(1)
        else:
            os.makedirs(args.pipeline_dir)
        if args.rx_list is None:
            args.rx_list = [
                "Sotrovimab",
                "Paxlovid",
                "Remdesivir",
                "Tixagevimab",
                "Cilgavimab",
                "Evusheld",
            ]
        gt = get_tables_sql.covid_drdb(args.rx_list, args.database_dir)
        if args.download:
            gt.download_latest()
        else:
            gt.get_database()
        gt.connect()
        gt.get_epitopes()
        # if not args.bam is None:
        #     make_output.make_epitope_graphs(args.bam, gt.epitopes, args.pipeline_dir, args.sample_name)
        gt.get_synonyms()
        gt.get_single_mutations()
        gt.get_fold_resistance()
        gt.add_local_resitances()
        mut_list_var = gt.get_variant_mutations(args.lineage)
        mf = get_mutants.mutantFinder(args.pipeline_dir, args.sample_name)
        fasta_coverage = None
        if args.vcf is None and args.fasta is None:
            sys.stderr.write("Please provide either a FASTA or a VCF file.")
        if not args.fasta is None:
            if not args.nextclade_data is None:
                nextclade_json = mf.run_nextclade(args.fasta, args.pipeline_dir, args.sample_name, args.nextclade_data)
            elif args.download_nextclade_data:
                nextclade_json = mf.run_nextclade(args.fasta, args.pipeline_dir, args.sample_name)
            else:
                sys.stderr.write("To use a fasta file please provide qhery with a nextclade "
                                 "data directory or set --download_nextclade_data\n")
                sys.exit(1)
            mut_list_sample, fasta_coverage, variant = mf.parse_nextclade(nextclade_json)
            if args.lineage is None:
                mut_list_var = gt.get_variant_mutations(variant)
        if not args.vcf is None:
            if not args.fasta is None:
                sys.stdout.write("VCF provided, will not use consensus sequence to determine mutations.\n")
            mf.convert_vcf(args.vcf)
            mf.run_bcf_csq()
            mut_list_sample = mf.parse_csq()
        mut_list_sample.sort(
            key=lambda x: (
                x.split(":")[0],
                int("0" + "".join([n for n in x.split("-")[0] if n.isdigit()])),
            )
        )
        print(mut_list_sample)
        nuc_to_aa_dict, aa_to_nuc_dict = make_output.get_nuc_aa_translations()
        if not args.bam is None:
            mut_list_lofreq = mf.recover_low_freq(args.bam)
        else:
            mut_list_lofreq = []
        mut_list_sample = list(set(mut_list_lofreq + mut_list_sample))
        make_output.make_final_tables(
            mut_list_sample,
            gt.resistances,
            mut_list_var,
            gt.epitopes,
            args.pipeline_dir,
            args.sample_name,
            args.bam,
            aa_to_nuc_dict,
            fasta_coverage
        )
        if not args.fasta is None:
            make_output.make_alignment_files(args.fasta, args.pipeline_dir, args.sample_name)


main("command_line")