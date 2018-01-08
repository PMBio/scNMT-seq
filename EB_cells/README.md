scNMT-seq
=========

Source code of the manuscript ***Joint profiling of chromatin accessibility, DNA methylation and transcription in single cells*** ([bioRxiv](biorxiv.org/content/early/2017/05/17/138685)).

Abstract
--------
Parallel single-cell sequencing protocols represent powerful methods for investigating regulatory relationships, including epigenome-transcriptome interactions. Here, we report the first single-cell method for parallel chromatin accessibility, DNA methylation and transcriptome profiling. scNMT-seq (single-cell nucleosome, methylation and transcription sequencing) uses a GpC methyltransferase to label open chromatin followed by bisulfite and RNA sequencing. We validate scNMT-seq by applying it to mouse embryonic stem cells and embryoid body cells, finding links between all three molecular layers and revealing dynamic coupling between epigenomic layers during differentiation.

<p align="center"> 
<img src="https://github.com/rargelaguet/scNMT/raw/master/protocol.png" style="width: 50%; height: 50%"/>​
</p>

For more details you can read our preprint: https://www.biorxiv.org/content/early/2017/11/10/217554

Content
-------
* `/preprocessing/`: preprocessing scripts
* `/acc_conservation/`: accessibility conservation using clustering of single-cell accessibility profiles
* `/correlations/`: pairwise correlations between different omic layers
* `/differential/`: differential analysis between subpopulations
* `/dimensionality_reduction/`: dimensionality reduction (PCA, t-SNE, ZIFA, etc.)
* `/heatmap/`: heatmap plots of the molecular layers
* `/pseudobulk_profiles/`: methylation and accessibility profiles using a running window with pseudobulked data
* `/pseudotime/`: coupling of molecular layers along a pseudotime trajectory
* `/stats/`: calculation of methylation and accessibility statistics (mean, coverage, etc.)
* `/zoom/`: plot of zoomed methylation and accessibility values in single genes


Data
-------
The raw data is accessible at GEO (to-do: add link...).
The parsed data is provided upon request.

Contact
-------
* Computational analysis: Ricard Argelaguet (ricard@ebi.ac.uk) or Andreas Kapourani (c.a.kapourani@ed.ac.uk)
* Experimental protocol: Stephen Clark (stephen.clark@babraham.ac.uk)
