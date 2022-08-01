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
  - name: Amy Jennison
    affiliation: 1
  - name: Sanmarie Schlebusch
    affiliation: 1
affiliations:
  - name: Q-PHIRE genomics, Forensic and Scientific Services, Queensland Health
    index: 1
date: 06 June 2022
bibliography: paper.bib

---

# Summary
Qhery is a Python application for the detection of mutations in SARS-CoV-2 that
may confer resistance to treatment. Qhery is freely available under a GPL license
for macOS, GNU/Linux and Windows from https://github.com/mjsull/qhery/


# Introduction 

Antivirals and neutralising monocolanl antibodies (mAbs) are now extensively
used for the prevention and treatment of COVID-19. Unfortunately, the virus
SARS-CoV-2 is able to rapidly adapt to escape treatment[@RN1]. SARS-CoV-2 has been
shown to adapt to escape binding of Sotrovimab in patients. It has also been shown
to adapt to antiviral treatments.

The Stanford Coronavirus Antiviral & Resistance database (CoV-RDB)[@RN2] has comprehensively 
curated published data on the neutralising susceptibility of SARS-CoV-2 variants and
spike mutations to mAbs. It also includes a sequence analysis program that when provided
with a consensus sequence can annotate mutations. These mutations are then checked against
the Cov-RDB to determine whether mutations may cause a decrease in Mabs susceptibility.

# Statement of Need
Coronavirus Resistance Database (CoV-RDB) is an excellent tool for the analysis of

limits processing large amounts of data. Privacy concerns. Limited to monoclonal antibodies.



# Implementation

Qhery is a pythong script available under a GPL license. It runs on macOS, GNU/Linux
and Microsoft Windows operating systems. Qhery can be used to assess SARS-CoV-2 for
mutations that confer resistance to a wide variety of treatments. 

Qhery takes a variant 


![Flowchart of qhery](https://github.com/mjsull/qhery/blob/main/paper/flowchart.svg?raw=true)

Legend goes here

More text



# Funding
This work was funded by the Queensland Government.

# References