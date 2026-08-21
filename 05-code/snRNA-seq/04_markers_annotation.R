# ---- conserved markers per cluster across conditions ----
# Run FindConservedMarkers on clusters to help identify
markers_cluster0 <- FindConservedMarkers(seu_rpca, ident.1 = 0, grouping.var = "condition")
markers_cluster1 <- FindConservedMarkers(seu_rpca, ident.1 = 1, grouping.var = "condition") 
markers_cluster2 <- FindConservedMarkers(seu_rpca, ident.1 = 2, grouping.var = "condition")
markers_cluster3 <- FindConservedMarkers(seu_rpca, ident.1 = 3, grouping.var = "condition")
markers_cluster4 <- FindConservedMarkers(seu_rpca, ident.1 = 4, grouping.var = "condition")
markers_cluster5 <- FindConservedMarkers(seu_rpca, ident.1 = 5, grouping.var = "condition")
markers_cluster6 <- FindConservedMarkers(seu_rpca, ident.1 = 6, grouping.var = "condition")
markers_cluster7 <- FindConservedMarkers(seu_rpca, ident.1 = 7, grouping.var = "condition")
markers_cluster8 <- FindConservedMarkers(seu_rpca, ident.1 = 8, grouping.var = "condition")
markers_cluster9 <- FindConservedMarkers(seu_rpca, ident.1 = 9, grouping.var = "condition")
markers_cluster10 <- FindConservedMarkers(seu_rpca, ident.1 = 10, grouping.var = "condition")
markers_cluster11 <- FindConservedMarkers(seu_rpca, ident.1 = 11, grouping.var = "condition")
markers_cluster12 <- FindConservedMarkers(seu_rpca, ident.1 = 12, grouping.var = "condition")
markers_cluster13 <- FindConservedMarkers(seu_rpca, ident.1 = 13, grouping.var = "condition")
markers_cluster14 <- FindConservedMarkers(seu_rpca, ident.1 = 14, grouping.var = "condition")
markers_cluster15 <- FindConservedMarkers(seu_rpca, ident.1 = 15, grouping.var = "condition")
markers_cluster16 <- FindConservedMarkers(seu_rpca, ident.1 = 16, grouping.var = "condition")
markers_cluster17 <- FindConservedMarkers(seu_rpca, ident.1 = 17, grouping.var = "condition")
markers_cluster18 <- FindConservedMarkers(seu_rpca, ident.1 = 18, grouping.var = "condition")
markers_cluster19 <- FindConservedMarkers(seu_rpca, ident.1 = 19, grouping.var = "condition")

#
# Extract top 20 gene names from your dataframe
top20_rat <- head(rownames(markers_cluster0), 20)

top20_rat <- toupper(top20_rat)

mapped_orthologs <- orthologs(
  genes = top20_rat,
  species = "rat",
  human = TRUE
)
top20_human <- mapped_orthologs$human_symbol



# Export markers to tsv files
write.table(
  markers_cluster0,
  file = "03-analysis_scratch/markers_cluster0.tsv",
  sep = "\t",
  quote = FALSE,
  col.names = NA
)




# Renaming clusters 
print(Idents(seu_rpca))
seu_rpca <- RenameIdents(seu_rpca, '0' = 'Cardiomyocytes')