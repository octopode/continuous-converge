# continuous-converge
scripts for detecting convergent adaptation of continuous traits

**"These are a work in progress, no guarantees!"**

As of 1/13/19, the scripts here implement two distinct methods: **PCOC** and **CPGLS**.

## PCOC
Original method: [Rey et al. 2018](http://dx.doi.org/10.1101/247296)

`pcoc_cont_scenarios.py` is a wrapper for the `pcoc_det.py` script published with 
the above paper. It attempts to identify amino acid sites that are convergently
adapted in association with a phenotypic trait. My script takes a continuous
trait; `pcoc_det` requires a binary one.

### Outline of Analysis
`pcoc_cont_scenarios` executes as follows:

1. Perform an ancestral trait reconstruction using the supplied species tree
and continuous trait table. At time of writing, this is a simple Brownian motion
(BM) reconstruction.

2. Generate a series of "convergent scenarios" for `pcoc_det` by proposing
"cutoffs" for discretizing the continuous trait into a binary one.
Cutoffs are established thus:

	1. List __all__ node trait values (terminal and internal) in order.
	
	2. Average trait values within a certain range of each other (termed trait 
	_precision_, set with `-p`).
	
	3. Place cutoffs at the midpoints between each pair of consecutive
	trait values.

3. Run `pcoc_det` for each unique scenario. (It is possible that two or more consecutive
cutoffs result in the same binary trait mapping.)

4. Collate all `pcoc_det` posterior probability (PP) values into a master dataframe
for downstream analysis.

5. Visualize these PP data in a heatmap or a Manhattan plot.

### Notes

* A series of command-line arguments (see `pcoc_cont_scenarios.py -h`) allow one
to specify which steps they want to run. 
e.g. Steps (4) and (5) above can be run on pre-generated PCOC results, but
trait values and precision must be specified the same as when the analysis was
first run, to ensure the same trait cutoff values.

* At time of writing, `pcoc_cont_scenarios` runs `pcoc_det` through a system call.
This is because `pcoc_det.py` as provided has no `main()`, and it means the user
should check the `subprocess.call()` call in `pcoc_cont_scenarios.py` to make
sure it's compatible with their shell environment. `pcoc_det` also requires
module `Bpp` (Bio++), which apparently changes a lot.

A straightforward solution to these issues is to use the Docker container 
distributed with PCOC:
```docker pull carinerey/pcoc```
and then install additional dependencies in it:
```
docker run -it carinerey/pcoc /bin/bash
pip install pandas biopython matplotlib
```
and then outside the container:
```docker commit <container ID> pcoc-matplot```

> The shell issue will go away when I refactor `pcoc_det.py` to use
a `main()` and implement a more efficient system for storing PCOC results indexed 
by _scenario_ rather than by _cutoff_.

### Example

To run `pcoc_cont_scenarios` in the Docker container, do something like this:
```
cd continuous-converge/
# or other dir containing both the repo and your data
# launch Docker container with $pwd mounted as /data/
docker run --mount type=bind,source="$(pwd)",target=/data/ -it pcoc-matplot /bin/bash
# run pcoc_cont_scenarios.py on your data. For the example dataset inside the repo:
xvfb-run python2 /data/pcoc_cont_scenarios.py -t data/example/FPs_data/tree/RAxML_bipartitions.FPs_62genes_MAFFT-gblocks-raxML-ultram.tree -o data/example/results -c example/FPs_data/trait/FPemWL_62genes_noref.tab -aa data/example/FPs_data/seq/FPs_62genes_MAFFT.fasta -p 2 -d -k _Avictoria -m master-table_2nm.tab -hm heatmap_2nm.pdf -mp manhattan_2nm.pdf
# ^xvfb-run is needed for graphics functionality in Docker
```

## CPGLS: Columnwise Phylogenetic Least-Squares
This is an original method based on phylogenetic regression _sensu_ Grafen 1989.

A separate least-squares regression is performed for each column of an amino acid
alignment. The predictor is AA state (a 20-D binary variable) and response is a
1-D continuous trait. 

At time of writing, the phylogenetic variance-covariance 
matrix is constructed under a simple BM model transformed with Pagel's lambda.
It is then adjusted to reflect the non-uniform substitution rates of amino acids.
Read on for how this is done.

### Outline of Analysis

1. Construct phylogenetic variance-covariance "**Q**" matrix for the entire dataset.

2. Load an empirical or database-derived transition rate matrix and normalize it
so that probabilities in each column sum to 1, rendering it asymmetrical.

Right now only BLOSUM62 is implemented.

3. Construct a biochemical variance-covariance "**B**" matrix for each column
in the alignment:

	1. Let **B** be a matrix with the same dimensions as **Q**, i.e. a row and 
	a column for each taxon.
	
	2. Let taxon _i_ encode amino acid _A_ at the column in question, 
	while taxon _j_ encodes amino acid _B_.
	
	3. Let **B**<sub>_i, j_</sub> = P(_A_->_B_) + P(_B_->_A_),
	i.e. the probability of transition from residue _A_ to residue _B_
	**or vice versa** per the normalized transition rate matrix in (2) above.
	
4. Calculate the element-wise non-exclusive sum of **Q** and **B**.
So for corresponding elements _p_ and _q_: _p_ + _q_ - _pq_.

5. Run the regression followed by a variety of hypothesis tests and return the
p-value for each, using any applicable multi-test correction.

### Notes

* CPGLS performs terribly in any simulation where either the "ancestral" or
"convergent" taxa contain several amino acid states. 2 possible solutions to this:

	1. _ad hoc_ "pooling" of AA states that correspond to sufficiently similar
	trait values. This imposes no explicit hypothesis about which residues group
	together functionally.
	
	2. Use a CAT-like profile mixture model 
	[Si Quang et al. 2008](http://dx.doi.org/10.1093/bioinformatics/btn445)
	to correlate continuous trait change with profile change along branches.
	
* CPGLS is roughly 10x faster than PCOC as presently implemented.

### Example

(no need to run from a Docker container or anything)

```
python2 cpgls.py -r example/FPs_data/trait/FPemWL_62genes_noref.tab -p example/FPs_data/seq/FPs_62genes_MAFFT.fasta -k _Avictoria -t example/FPs_data/tree/RAxML_bipartitions.FPs_62genes_MAFFT-gblocks-raxML-ultram.tree -m > FPs-pgls.tab
```