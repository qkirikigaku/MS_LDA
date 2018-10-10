# MS_LDA
MS_LDA extracts mutation signatures contained in samples from mutation catalogs of cancer genome with Bayesian model called latent Dirichlet allocation (LDA) and variational Bayes method (VB).

Implementation is done in Python 3.5 and C++.

## Requirements
* boost 1.58.0

## Preprocessing
To create a dataset from a mutation catalog of the COSMIC database (https://cancer.sanger.ac.uk/cosmic), execute the following command.

In advance, download mutation catalog (`CosmicMutantExport.tsv.gz`) from COSMIC official site and put `raw_data/CosmicMutatntExport.tsv`.

Here is an explanation when M1 and breast are adopted as a mutation dictionary and cancer type respectively.
For the details of the mutation dictionary, see the original paper.
```
sh Preprocessing/make_d1.sh 10
python Preprocessing/get_M1.py 400 breast
```

Then, an input file (`data1_o400_breast.txt`) is created under the `data/` directory.

The first row of the input file shows the number of samples and the number of mutational types, and the breakdown of each sample's mutation is shown after the second row.

## To compile
Execute `make compile`.

## To extract mutation signatures
Change the input destination and the output destination of source code as necessary.
Execute the following command.
```
scripts/MS_real.sh ${data_type} 400 ${cancer_type}
```
Here, data_type and cancer_type show mutation dictionary (1,2,3 or 4), and primary lesion respectively.

Example:
```
scripts/MS_real.sh 1 400 breast
```
Then, a result directory corresponding to the condition is created under the result directory (for example `result/data1_o400_breast_*`).

## Model selection and visualization of results

In order to select plausible number of signatures and visualization, execute the following command.
Since the scp command is included in `scripts/make_figure.sh` currently so as to refer to the result directory existing in the remote terminal, please comment the row containing `scp` out in `scripts/make_figure.sh`.

```
sh scripts/make_figure.sh ${data_type} 400 ${cancer_type}
```

For example, this will delete `result/data1_o400_breast_*` and `result/data1_o400_breast` will be created instead.
This result directory contains the distribution of mutations of each signature and graphs showing the transition of the value of the Variational Lower Bound for each signature number.
