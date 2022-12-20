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
    list_parser.add_argument(
        "--download",
        help="Download the latest database.",
        action="store_true"
    )
    list_parser.add_argument(
        "--details",
        help="List rx type and synonyms.",
        action="store_true"
    )
    list_parser.add_argument(
        "--include_unapproved",
        "-a",
        help="List rx that have not been give EUA.",
        action="store_false"
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
    run_parser.add_argument("--rx_list", "-rx", nargs="+", help="List of drugs to analyze, if empty will run on all approved drugs.")
    run_parser.add_argument("--fasta", "-f", help="Consensus fasta.")
    run_parser.add_argument("--download_nextclade_data", "-dn", action="store_true", help="download nextclade data.")
    run_parser.add_argument("--nextclade_data", "-nd", help="directory of sars-cov-2 nextclade data, "
                            "can be downloaded with "
                            "\"nextclade dataset get --name 'sars-cov-2' --output-dir 'data/sars-cov-2'\"")


    args = parser.parse_args()
    gt = get_tables_sql.covid_drdb(args.database_dir)
    if args.download:
        gt.download_latest()
    else:
        gt.get_database()
    gt.connect()
    if args.subparser_name == "list_rx":
        gt.list_rx(args.include_unapproved, args.details)
    elif args.subparser_name == "run":
        if os.path.exists(args.pipeline_dir) and os.path.isdir(args.pipeline_dir):
            pass
        elif os.path.exists(args.pipeline_dir):
            sys.stderr.write("Pipeline dir exists and is not a directory.\n")
            sys.exit(1)
        else:
            os.makedirs(args.pipeline_dir)
        print(args)
        if args.rx_list is None:
            rx_list = gt.get_rx(True)
        else:
            temp_rx_list = gt.get_rx(False)
            rx_list = []
            for i in temp_rx_list:
                add_it = False
                for j in args.rx_list:
                    if i.name.lower() == j.lower() or j.lower() in i.synonyms_lower():
                        add_it = True
                        break
                if add_it:
                    rx_list.append(i)
        gt.drug_list = rx_list
        gt.connect()
        gt.get_ref()
        gt.get_epitopes()
        gt.get_single_mutations()
        print('dong')
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
            gt.epitope_dict(),
            args.pipeline_dir,
            args.sample_name,
            args.bam,
            aa_to_nuc_dict,
            fasta_coverage
        )



main("command_line")