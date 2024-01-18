---
title: 'qhery: uncovering resistance mutations in SARS-CoV-2'
tags:
  - Python
  - Bioinformatics
  - Viral Genomics
  - Drug Resistance
  - sequence analysis
authors:
  - name: Mitchell John Sullivan
    orcid: 0000-0003-0872-7098
    affiliation: 1
  - name: Alyssa Pyke
    affiliation: 1
  - name: Wytamma Wirth
    affiliation: 2
  - name: Amy Jennison
    affiliation: 1
  - name: Sanmarie Schlebusch
    affiliation: 1
affiliations:
  - name: Q-PHIRE genomics, Forensic and Scientific Services, Queensland Health
    index: 1
  - name: Peter Doherty Institute for Infection and Immunity, University of Melbourne, Melbourne, VIC, Australia
    index: 2
date: 17 January 2024
bibliography: paper.bib

---

# Summary

Antivirals and neutralising monoclonal antibodies (mAbs) are now extensively
used for the prevention and treatment of COVID-19. Unfortunately, the virus
SARS-CoV-2 is rapidly changing and able to escape treatment[@cox_sars-cov-2_2023]. SARS-CoV-2 has been
shown to adapt to escape binding of Sotrovimab in patients[@rockett_resistance_2022]. Mutations that reduce
the effectiveness of antiviral treatments have also been observed[@iketani_functional_2022;iketani_multiple_2023].

Two tools exist for the analysis of resistance mutations in SARS-CoV-2.
The Stanford Coronavirus Antiviral & Resistance database (CoV-RDB)[@tzou_coronavirus_2022] has comprehensively
curated published data on the neutralising susceptibility of SARS-CoV-2 variants and
spike mutations to mAbs. It also includes mutations that grant resistance to commonly used antivirals.
Included on the website is a sequence analysis program that when provided with a consensus sequence or raw reads can
annotate mutations. These mutations are then checked against the Cov-RDB to determine whether mutations may cause a
decrease in mAbs susceptibility. Sabres is a command-line tool that compares variants found in a VCF file to an
internally
maintained database of resistance mutations [@fong_sabres_2023].

Qhery is a Python application for the detection of mutations in SARS-CoV-2 that
may confer resistance to treatment. Qhery is freely available under a GPL license
for macOS and GNU/Linux from https://github.com/mjsull/qhery/

# Statement of Need

Coronavirus Resistance Database (CoV-RDB) is an excellent tool for the analysis of resistance mutations. However, it's
analysis pipeline is currently only available via web portal. This limits the number of sequences that can be processed
efficiently, and it's ability to be incorporated into an automated monitoring platform such as
Austrakka[@hoang_austrakka_2022]. There may also privacy concerns with uploading private sequencing information to a
third-party website. SABRES is also a useful tool,
however it's internal database is updated infrequently. It is also unable to detect low frequency variants natively. Qhery
removes
these limitations by allowing the user to compare SARS-CoV-2 genomic information to the frequently updated CoV-RDB
locally.

# Implementation

Qhery can take a combination of variant call format (vcf) file, a FASTA file of the consensus sequence or a BAM file of raw reads aligned to the SARS-CoV-2 reference (MN908947.3) and converts them into amino acid changes. If a BAM file is
provided QHERY identifies low frequency variants using lofreq[@wilm_lofreq_2012] which is then converted to Amino acid changes using
bcftools csq[@danecek_bcftoolscsq_2017]. Provided VCF files are also converted to amino acid changes using bcftools csq. Finally, nextclade
is used
to determine amino acid changes from a FASTA file. Nextclade can also be used to determine the lineage of the query if
not provided by the user. For a list of provided treatments, mutations that increase resistance to treatment are
retrieved
from the CoV-RDB database. This database can be automatically downloaded or provided to the software. Identified
resistance
mutations are then checked against the list of observed mutations for the query and reported in table. In addition,
codon
frequencies in the BAM file are examined using pysam (https://github.com/pysam-developers/pysam) which is a wrapper around HTSlib [@bonfield_htslib_2021] and 
samtools [@danecek_twelve_2021] and if a resistance mutation is present above a user defined
threshold that mutation is also reported. This information is compiled into a table that details whether a mutation
is present in the sample, the codon frequency at that site, and whether that mutation is a lineage defining mutation.
The fold change in resistance for each of the treatments specified and whether that mutation falls in the epitope or
binding pocket of the treatment is also included.



![Flowchart of qhery](https://github.com/mjsull/qhery/blob/main/paper/flowchart.svg?raw=true)

## Figure 1: Flowchart of Qhery

# Availability

Qhery is a Python script available under a GPL license. It runs on macOS, GNU/Linux
and Microsoft Windows operating systems. Qhery can be used to assess SARS-CoV-2 for
mutations that confer resistance to a wide variety of treatments.

# Funding

This work was funded by the Queensland Government.

# Acknowledgements

We would like to thank

# References