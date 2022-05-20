# qhery

## 
![release](https://img.shields.io/github/v/release/mjsull/qhery) ![license](https://img.shields.io/badge/license-GPLv3-green)



### Installation 

#### git

While in development qhery can only be installed by downloading from git

```
git clone https://github.com/mjsull/qhery.git
```



#### requirements:

- **Python >= 3.9.12**
- **bcftools >= 1.10.2**
- **curl >= 7.83.1**
- **wget >= 1.20.3**

##### Optional requirements:
- **ncbi-blast+ >= 2.9.0+** - This will generate a BLASTx alignment of the genome for visualization
- **lofreq >= 2.1.5** - if provided with a BAM file **qhery** will look for minor alleles in the alignment with lofreq
- **samtools >= 1.7** - samtools is used to determine the depth of sequence along the genome, and which resistance mutations
cannot be reported on due to lack of coverage.

### Example usage: 

`qhery run --database_dir database_dir --vcf sample.vcf --pipeline_dir output_dir --lineage Omicron/BA.1 --sample_name mysample --rx_list Sotrovimab`

Determines the amino acid changes caused by the mutations listed in sample.vcf and then compares them to a list of mutations that cause a reduction in Sotrovimab binding.

`qhery run --database_dir database_dir --vcf sample.vcf --pipeline_dir output_dir --lineage Omicron/BA.1 --sample_name mysample --rx_list Sotrovimab Remdesivir --fasta sample.consensus.fasta --bam sample.primertrimmed.rg.sorted.bam` 

Determines the amino acid changes caused by the mutations listed in sample.vcf. Additionally will use lofreq to find minor alleles in the BAM file. Finally they are  compared to a list of mutations that cause a reduction in Sotrovimab binding or reduction in Remdesivir efficiency.
 

`qhery list_rx`

List treatments for which resistance information exists.


### Example output:

| Mutation    | alt_names   | in_sample | in_variant | covered  | resistance_mutation | Remdesivir_average_fold_reduction | Remdesivir_fold_reductions | Remdesivir_in_epitope | Sotrovimab_average_fold_reduction | Sotrovimab_fold_reductions | Sotrovimab_in_epitope |
|-------------|-------------|-----------|------------|----------|---------------------|-----------------------------------|----------------------------|-----------------------|-----------------------------------|----------------------------|-----------------------|
| E:T9I       | -           | True      | True       | True     | False               | 0                                 | -                          | False                 | 0                                 | -                          | False                 |
| M:D3G       | -           | True      | True       | True     | False               | 0                                 | -                          | False                 | 0                                 | -                          | False                 |
| N:ERS31-33∆ | -           | True      | True       | False    | False               | 0                                 | -                          | False                 | 0                                 | -                          | False                 |
| ORF3a:L52F  | -           | True      | False      | True     | False               | 0                                 | -                          | False                 | 0                                 | -                          | False                 |
| RdRP:802D   | -           | False     | False      | True     | True                | 2.54                              | =2.54                      | False                 | 0                                 | -                          | False                 |
|S:R214ins    | S:R214R_EPE | True      | True       | True     | True                | 0                                 | -                          | False                 | 3.00                              | =3.0                       | False
| S:P337T     | -           | True      | False      | True     | True                | 0                                 | -                          | False                 | 8.00                              | =5.4,=10.6                 | True                  |

### Columns
| column | header                     | description                                                                                                    |
|--------|----------------------------|----------------------------------------------------------------------------------------------------------------|
| 1      | mutation                   | The mutation name. Gene name comes before the colon, then reference amino acid, position and sample amino acid |
| 2      | alt_names                  | Discrepency between database mutation name and csq mutation name                                               |
| 3      | in_sample                  | Is mutation in the query                                                                                       |
| 4      | in_variant                 | Is mutation a lineage defining mutation                                                                        |
| 5      | covered                    | Is the mutation covered by 20 or more reads                                                                    |
| 6      | resistance_mutation        | Is there evidence the mutation may confer some resistance to one of the treatments listed                      |
| 7      | rx1_average_fold_reduction | Average fold reduction of listed fold reductions                                                               |
| 8      | rx1_fold_reductions        | Fold reductions listed in the database                                                                         |
| 9      | rx1_in_epitope             | Is the mutation in the epitope of this treatment (MABs only)                                                   |
| 10     | rx2_average_fold_reduction | The previous 3 columns repeate for each treatment provided by the user                                         |
| 11     | rx2_fold_reductions        | ...                                                                                                            |
| 12     | rx2_in_epitope             | ...                                                                                                            |

## TODO
1. Need a phasing step between lofreq and bcftools csq (or to switch to a vcf caller that does phasing)
2. Add allele frequency information to output table


## Arguments:


``-h``, ``--help``

show this help message and exit

### subcommands

### list_rx
___

List all drugs for which resistance information is available.

Takes no arguments



<hr style="border:1px solid gray">

### run
___

Determines mutations in samples and then checks against resistance data.

#### arguments


``-n``, ``--sample_name`` <sample_name>

Sample name, output files will be prefixed with this.

``-v``, ``--vcf`` <sample.vcf>

vcf file, variants called against the Wuhan-Hu-1 reference (MN908947.3)

``-b``, ``--bam`` <sample.sorted.bam>

Sorted bam file. File of read alignments for the sample mapped against the Wuhan-Hu-1 reference (MN908947.3)

``-d``, ``--database_dir`` <path/to/database_dir>

Directory with the latest version of the Stanford [resistance database](https://github.com/hivdb/covid-drdb-payload/releases/latest).
If the latest version is not in this folder it will be downloaded to this location.

``-p``, ``--pipeline_dir`` <path/to/pipeline_dir>

All script output and intermediated files will be put here. Script will create a directory if none exists.

``-l``, ``--lineage`` <BA.1>

Lineage of the query (BA.1/BA.2/BA.3/Delta etc.)

``-rx``, ``--rx_list`` \<Sotrovimab Remdesevir>

List of treatments to interrogate. 

``--fasta``, ``--fasta`` <sample.fasta>

Fasta file of the consensus sequence of the sample, only used to generate a BLASTx alignment for double checking mutations.

<hr style="border:1px solid gray">

#### mutations
___

Only list mutations and not resistance information.

``-n``, ``--sample_name`` <sample_name>

Sample name, output files will be prefixed with this.

``-v``, ``--vcf`` <sample.vcf>

vcf file, variants called against the Wuhan-Hu-1 reference (MN908947.3)

``-b``, ``--bam`` <sample.sorted.bam>

Sorted bam file. File of read alignments for the sample mapped against the Wuhan-Hu-1 reference (MN908947.3)

``-d``, ``--database_dir`` <path/to/database_dir>

Directory with the latest version of the Stanford [resistance database](https://github.com/hivdb/covid-drdb-payload/releases/latest).
If the latest version is not in this folder it will be downloaded to this location.

``-p``, ``--pipeline_dir`` <path/to/pipeline_dir>

All script output and intermediated files will be put here. Script will create a directory if none exists.

``-l``, ``--lineage`` <BA.1>

Lineage of the query (BA.1/BA.2/BA.3/Delta etc.)


``-k``, ``--keep_lineage``

report lineage defining mutations as well


<hr style="border:1px solid gray">


## Method

A flowchart of how ``qhery run`` works

![flowchart](https://github.com/mjsull/qhery/blob/master/paper/flowchart.svg?raw=true)
