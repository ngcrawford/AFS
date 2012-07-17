
Miss estimating AFS has consequences

	- coverage, number of individuals (error and direction = per genotype quality threshold)

1.) Get simulations going
	
	- Simulate lineages of chromosomes, generate 'realistic' SF use MS to generate chromosome(s) and simulate reads from chromosomes
	- Or, draw SF, and impose them on sequences

2.) GATK frequency spectrum (AF in the INFO section)

3.) Samtools AF (map reads by indiv)

4.) Observed AF from different genotype quality thresholds

5.) Directional error (from true) for most categories

6.) 


Start from MS:

	- ten megabase genomes, theta of 1%, coverage as Gaussian with a coefficient of variation at 0.5

	- errors not associated with reduced quality scores
	- 10 and 100 samples
	- low coverage data
	- mean read counts 1x to 100x (intervals of 5)

To show it's needed:

	RMSD (root mean squared deviation) from the true spectrum (by category)
	- heat maps: rares, common, all. RMSD/Quatity threshold and read count
	- same for spectra derived from "our little trick" (string pulling)

