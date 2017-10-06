# Migration-physiology

Data and code to go with the manuscript "Active migration is associated with specific and consistent changes to gut microbiota in shorebirds"
by Alice Risely, David Waite, Bea Ujvari, Bethany Hoye & Marcel Klaassen

Data analysed with Phyloseq and DESeq2 packages. Code to download this package from Bioconductor within this code.

File description

shorebird.metadata.csv = metadata

shorebird.shared = OTU table

shorebird.taxonomy = Taxonomy table

shorebird.tree = Phylogenetic tree of OTUs

migrant_resident_analysis.R = R code for full analysis with the above files with Phyloseq.

miseq.R = Source code with some helpful functions. Should only need if you want to scale reads (as sensitivity analysis), otherwise not needed to run published analysis.



