#########################################################
## Preprocessing of methylation and accessibility data ##
#########################################################

Steps:
(0) (If required) Remove sites that were not covered by any read: sometimes Bismark spits out a file were the non-observed sites have been removed, but sometimes it does not.
	raw_filter.sh
	
(1) (If required) Merge R1 and R2 files
	merge.R

(1) (Optional) Filtering non-CG methylation and non-GC accessibility: this is not required because this step is done by Bismark. If we find missmatches they are more likely to be disagreements with respect to the reference genome but true CG or GC dinucleotides.
	raw_filter.R

(2) (Optional) Remove strnad-specific information: convert all CpG/GpC sites to the positive strand
	assign_strand.R

(3) (Optional but recommended) Binarise CpG sites and calculate rates
	binarise.R

(4) Collapse single CpG sites into features using known annotations:
	annotate.R

(5) (Optional)  Merge loci with the same ensembl ID (for instance exons or introns in a gene)
	merge_sites.R

(6) (Optional) Second raw of filtering: by coverage and variance
	parsed_filter.R

Other:
(-) bulk.R: pseudobulk the single-cell data
(-) 