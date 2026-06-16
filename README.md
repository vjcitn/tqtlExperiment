# tQTLExperiment is a package for xQTL data management

This package defines a data structure that
combines molecular phenotype data with genotype
data in PLINK bed format.

The functions in the package permit interoperation with
tensorqtl to collect statistics on association between
genotype and molecular phenotypes.

## Data structure review

```
> tqe
class: tQTLExperiment
features: 20  samples: 100 
assays( 1 ): pheno 
rowRanges: GRanges with 20 features
colData( 11 ) covariates: int, sex, age, expr_pc1 ... 
geno: 100 samples x 69638 variants [BEDMatrix - lazy]
plinkPrefix: /Library/Frameworks/R.framework/Versions/4.6/Resources/library/tQTLExperiment/extdata/chr22-n100 
use prepareTQTL() to write inputs and get CLI command
> assay(tqe[1:3,1:4])
                   691_692 693_694     688_689    683_684
ENSG00000100181 0.00000000       0 0.006378173 0.00000000
ENSG00000237438 0.00000000       0 0.000000000 0.00000000
ENSG00000177663 0.01798248       0 0.039776113 0.01103418
> tqtlGeno(tqe)[1:4,1:3]
          22:16849573A-G_G 22:16849971A-T_T 22:16850437G-A_A
0_691_692                1                1                1
0_693_694                1                1                1
0_688_689                1                1                1
0_683_684                1                1                1
> tqe = addGeneSymbols(tqe)
> rowRanges(tqe[1:3,])
GRanges object with 3 ranges and 2 metadata columns:
                  seqnames    ranges strand |    phenotype_id   gene_name
                     <Rle> <IRanges>  <Rle> |     <character> <character>
  ENSG00000100181       22  17082777      * | ENSG00000100181      TPTEP1
  ENSG00000237438       22  17517460      * | ENSG00000237438       CECR7
  ENSG00000177663       22  17565844      * | ENSG00000177663      IL17RA
  -------
  seqinfo: 1 sequence from hg38 genome; no seqlengths
```

## Running tensorqtl

### Preparing the data

At present, python/R interop is not well-established.  We export our data
to a specified, existing folder.  The `prepareTQTL` function does this
and provides a string that can be used to compute the association statistics.

```
> tt = prepareTQTL(tqe, outDir="/tmp/test1") 
> tt
[1] "python3 -m tensorqtl /Library/Frameworks/R.framework/Versions/4.6/Resources/library/tQTLExperiment/extdata/chr22-n100 /tmp/test1/pheno.bed /tmp/test1/tqtl_out --mode cis_nominal --maf_threshold 0.05 --window 1000000 -o /tmp/test1 --covariates /tmp/test1/covariates.tsv"
```

### Running tensorqtl

```
python3 -m tensorqtl /Library/Frameworks/R.framework/Versions/4.6/Resources/library/tQTLExperiment/extdata/chr22-n100 /tmp/test1/pheno.bed /tmp/test1/tqtl_out --mode cis_nominal --maf_threshold 0.05 --window 1000000 -o /tmp/test1 --covariates /tmp/test1/covariates.tsv
Warning: 'rfunc' cannot be imported. R with the 'qvalue' library and the 'rpy2' Python package are needed to compute q-values.
[Jun 16 07:16:37] Running TensorQTL v1.0.10: cis-QTL mapping
  * WARNING: using CPU!
  * reading phenotypes (/tmp/test1/pheno.bed)
  * cis-window detected as position ± 1,000,000
  * reading covariates (/tmp/test1/covariates.tsv)
  * loading genotypes
Mapping files: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████| 3/3 [00:00<00:00, 63.13it/s]
/Users/vincentcarey/Library/Python/3.12/lib/python/site-packages/pandera/_pandas_deprecated.py:149: FutureWarning: Importing pandas-specific classes and functions from the
top-level pandera module will be **removed in a future version of pandera**.
If you're using pandera to validate pandas objects, we highly recommend updating
your import:

```
# old import
import pandera as pa

# new import
import pandera.pandas as pa
```

If you're using pandera to validate objects from other compatible libraries
like pyspark or polars, see the supported libraries section of the documentation
for more information on how to import pandera:

https://pandera.readthedocs.io/en/stable/supported_libraries.html

To disable this warning, set the environment variable:

```
export DISABLE_PANDERA_IMPORT_WARNING=True
```

  warnings.warn(_future_warning, FutureWarning)
cis-QTL mapping: nominal associations for all variant-phenotype pairs
  * 100 samples
  * 20 phenotypes
  * 11 covariates
  * 69638 variants
  * applying in-sample 0.05 MAF filter
  * cis-window: ±1,000,000
  * checking phenotypes: 20/20
  * Computing associations
    Mapping chromosome 22
    processing phenotype 20/20
    time elapsed: 0.00 min
    * writing output
done.
[Jun 16 07:16:39] Finished mapping
```

### Retrieving results to R

Some improvements are needed ...

```
> tst = readTQTL("/tmp/test1")
> tst
$pairs
GRanges object with 71522 ranges and 9 metadata columns:
          seqnames    ranges strand |    phenotype_id     variant_id
             <Rle> <IRanges>  <Rle> |     <character>    <character>
      [1]       22  16849573      * | ENSG00000100181 22:16849573A-G
      [2]       22  16849971      * | ENSG00000100181 22:16849971A-T
      [3]       22  16850437      * | ENSG00000100181 22:16850437G-A
      [4]       22  16851225      * | ENSG00000100181 22:16851225C-T
      [5]       22  16851356      * | ENSG00000100181 22:16851356C-T
      ...      ...       ...    ... .             ...            ...
  [71518]       22  20309073      * | ENSG00000185065 22:20309073G-A
  [71519]       22  20309393      * | ENSG00000185065 22:20309393T-G
  [71520]       22  20311136      * | ENSG00000185065 22:20311136C-T
  [71521]       22  20311425      * | ENSG00000185065 22:20311425G-A
  [71522]       22  20318162      * | ENSG00000185065 22:20318162G-A
          start_distance        af ma_samples  ma_count pval_nominal
               <integer> <numeric>  <integer> <integer>    <numeric>
      [1]        -233204       0.4         64        80     0.295731
      [2]        -232806       0.4         64        80     0.295731
      [3]        -232340       0.4         64        80     0.295731
      [4]        -231552       0.4         64        80     0.295731
      [5]        -231421       0.4         64        80     0.295731
      ...            ...       ...        ...       ...          ...
  [71518]         873657     0.120         23        24     0.357219
  [71519]         873977     0.510         76        98     0.801957
  [71520]         875720     0.275         49        55     0.371372
  [71521]         876009     0.320         57        64     0.316733
  [71522]         882746     0.110         21        22     0.443948
                 slope    slope_se
             <numeric>   <numeric>
      [1]  -0.00204086  0.00194006
      [2]  -0.00204086  0.00194006
      [3]  -0.00204086  0.00194006
      [4]  -0.00204086  0.00194006
      [5]  -0.00204086  0.00194006
      ...          ...         ...
  [71518]  0.001083448 0.001170547
  [71519]  0.000194450 0.000772904
  [71520] -0.000770468 0.000857456
  [71521] -0.000879878 0.000873772
  [71522]  0.000927042 0.001205440
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths
```
