# Fix aliases in gene list -----

gene.list <- data.table::fread("MR_gene_list_final.txt", data.table = FALSE, header = FALSE)

colnames(gene.list) <- c("Symbol", "Table", "Pathway")



id.map <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv") %>%
  
  dplyr::select(Ensembl = ensembl_gene_id, Symbol = external_gene_name)



## Convert aliases to current symbols

aliases <- gene.list$Symbol[!gene.list$Symbol %in% id.map$Symbol]

alias.symbol <- tibble(
  
  Alias = aliases,
  
  Symbol = limma::alias2SymbolTable(aliases, species = "Hs")
  
) %>%
  
  filter(!is.na(Symbol)) %>%
  
  left_join(gene.list, by = c("Alias" = "Symbol"))



gene.list <- bind_rows(
  
  gene.list[!gene.list$Symbol %in% alias.symbol$Alias, ],
  
  dplyr::select(alias.symbol, -Alias)
  
) %>%
  
  left_join(id.map, by = "Symbol") %>%
  
  filter(!is.na(Ensembl)) %>%
  
  distinct(.keep_all = T)



rm(aliases, alias.symbol, id.map)
