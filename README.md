# Migration-physiology

Data and code to go with the manuscript "Active migration is associated with specific and consistent changes to gut microbiota in shorebirds"
by Alice Risely, David Waite, Bea Ujvari, Bethany Hoye & Marcel Klaassen J. ANIMAL ECOLOGY 2017

Data analysed with Phyloseq and DESeq2 packages. Code to download this package from Bioconductor within this code.

File description

DATA FILES

shorebird.metadata.csv = metadata

shorebird.shared = OTU table

shorebird.taxonomy = Taxonomy table

shorebird.tree = Phylogenetic tree of OTUs

R CODE

migrant_resident_analysis.R = R code for full analysis with the above files with Phyloseq.

miseq.R = Source code with some helpful functions. Should only need if you want to scale reads (as sensitivity analysis), otherwise not needed to run published analysis.

R MARKDOWN FILE

R_markdown.pdf = R markdown file of analysis

BIOINFORMATICS - MOTHUR CODE

mothur_code.logfile = This is a copy of the mothur logfile outlining how we processed sequences



