# PathotypeR
DEC pathotype assignment of _E. coli_ genomes in R.

Diagnostic microbiology has developed several schemes to subtype _Escherichia coli_. These are useful for understanding the epidemiology and pathogenesis of particular strains. Diarrheagenic _E. coli_ (DEC) pathotypes group _E. coli_ strains which possess similar virulence factors (VF) and cause diseases with similar pathology. Whole-genome sequencing can predict the presence of VF genes in an _E. coli_ isolate with high accuracy, permitting the assignment of DEC pathotypes without additional molecular, biochemical, or phenotypic assays. 

`PathotypeR` assigns _E. coli_ genomes a DEC pathotype based on the presence/absence of specific VF genes. Inputs the output of [`AMRFinderPlus`](https://github.com/ncbi/amr).

`PathotypeR` includes two functions:
1. `pathotypeR()`: Function that quantifies each samples' VFs and assigns a DEC pathotype. Can also return total VF count per sample, VF presence/absence, the prevalence of each pathotype, and the prevalence of each VF.
2. `amrfinder_process()`: Merges `AMRFinderPlus` output files into a single dataframe. (Called by `pathotypeR()` but can be used on own.)

Genomes are assigned a DEC pathotype based on the presence/absence of specific VF genes. Namely:
- Shiga toxin-producing _E. coli_ (STEC): _stx1_ and/or _stx2_ (without _eae_)
- Enteropathogenic _E. coli_ (EPEC): _eae_ and/or _bfpA_ (without _stx1_ and/or _stx2_)
- Enterohaemorrhagic _E. coli_ (EHEC): _stx1_ and/or _stx2_, and _eae_
- Enteroinvasive _E. coli_ (EIEC): _ipaH_
- Enterotoxigenic _E. coli_ (ETEC): _ltcA_ and/or _sta1_
- Enteroaggregative _E. coli_ (EAEC): _aatA_ and/or _aaiC_ and/or _aggR_
- Diffusely adherent _E. coli_ (DAEC): _afaC_ and/or _afaE_
- none: does not encode any of the above VF genes

Hybrid strains contain genes associated with multiple DEC pathotypes and will be reported as all pathotypes detected (e.g., STEC-EAEC, EPEC-ETEC). 

NOTE: PathotypeR does NOT assign pathotypes based on collection site or association with disease, such as: extraintesinal pathogenic _E. coli_ (ExPEC), uropathogenic _E. coli_ (UPEC), neonatal meningitis-associated _E. coli_ (NMEC), and sepsis-associated _E. coli_ (SEPEC).

# Content

[Installation](#install)

[Pre-process](#pre-processing)

[Functions](#functions)

[References](#references)


## Install

Install directly from GitHub:
```r
source("https://raw.github.com/kevinsblake/PathotypeR/main/pathotype.R")
```

Alternatively, can download and then install using the filepath:
```r
source("dir1/dir2/pathotype.R")
```

## Pre-processing

_E. coli_ genomes of interest must first be run through [`AMRFinderPlus`](https://github.com/ncbi/amr). See their instructions for recommended usage.

The `AMRFinderPlus` output must be saved as a `.tsv` file. This can be done using the output flag: `-o ${outdir}/${sample}.tsv`. Copy all of these output files into one directory. The filepath of this directory will be the input for PathotypeR. 

## Functions

### `pathotypeR()`

#### Description
Function for assigning DEC pathotype to _E. coli_ genomes. First calls `amrfinder_process()`.

#### Usage

```r
library(dplyr)

pathotypeR(indir, output=c("patho_pred", "patho_prev", "vf_pres", "vf_prev"))
```

#### Arguments
`indir`		Filepath to directory containing `AMRFinderPlus` output files.

`output`	Specifies output. `patho_pred` = for each sample, outputs VF count and pathotype prediction; `patho_prev` = for each pathotype, outputs count (i.e. number of samples) and overall prevalence; `vf_pres` = for each sample, outputs VF presence/absence (1=present, 0=absent); `vf_prev` = for each VF, outputs count and overall prevalence. Default is `patho_pred`.


#### Examples

```r
# Outputs just sample names, VF count, and pathotype prediction
df <- pathotypeR("data/amrfinder")

# Outputs all VFs detected, count, and overall prevalence
df <- pathotypeR("data/amrfinder", output="vf_prev")
```

### `amrfinder_process()`

#### Description
Function for merging `AMRFinderPlus` output files into a single dataframe.

#### Usage

```r
amrfinder_process(indir)
```

#### Arguments

`indir`		Filepath to directory containing `AMRFinderPlus` output files.
`suffix`  Specifies the suffix added to the amrfinder output filename. The filename minus this suffix should be the same as the `Name` column in the amrfinder output. Default = ".tsv"

#### Examples

```r
df <- amrfinder_process("data/amrfinder")
```

## References

- Horesh _et al._ A comprehensive and high-quality collection of Escherichia coli genomes and their genes. _Microb Genom._ 2021 Feb;7(2):000499. doi: 10.1099/mgen.0.000499. PMID: 33417534.
- Jesser & Levy. Updates on defining and detecting diarrheagenic Escherichia coli pathotypes. _Curr Opin Infect Dis._ 2020 Oct; 33(5): 372–380. doi: 10.1097/QCO.0000000000000665. PMID: 32773499.
- Robins-Browne _et al._ Are _Escherichia coli_ pathotypes still relevant in the era of whole-genome sequencing? 2016 Nov 18;6:141. doi: 10.3389/fcimb.2016.00141. eCollection 2016.PMID: 27917373.
