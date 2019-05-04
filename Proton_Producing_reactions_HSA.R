# get the proton producing genes for each pathway

setwd("~/storage/Metabolic_reprogramming/Proton_producing_reactions")

# Get all reactions associated with Protons (C00080) -----------------------------

Proton_info <- keggGet("cpd:C00080")
rxns <- Proton_info[[1]]$REACTION

# Put reactions into a list of characters
rxns <- unlist(sapply(rxns, function(x) strsplit(x, " ")), use.names = F)
rxns <- paste0("rn:", rxns)

res <- matrix(nrow = length(rxns), ncol = 3)
colnames(res) <- c("reaction", "equation", "enzyme")
# Request 10 reactions at a time due to API limitations
for (i in seq(from=1, to=length(rxns), by=10)) {
  message(stringr::str_glue("Getting reactions {i} to {min(i + 9, length(rxns))}."))
  rxn <- keggGet(rxns[i:min(i + 9, length(rxns))])
  for (j in seq_along(rxn)) {
    res[i+ j - 1, ] <- c(rxn[[j]]$ENTRY,
                         rxn[[j]]$EQUATION,
                         paste(rxn[[j]]$ENZYME, collapse = ","))
  }
}
res <- as.data.table(res)

# 3 columns: reaction, equation, enzyme --------------------------------------

# reactions without EC number
rxn_wo_ec <- res[enzyme == ""]

# reactions with EC numbers like 2.7.2.-
rxn_w_ec_class <- res[grep("-", enzyme)]

# reactions with multiple EC numbers; ones with EC classes are excluded
rxn_w_multiple_ec <- res[grep(",", enzyme)][!reaction %in% rxn_w_ec_class$reaction]

# reactions with only one EC number
rxn_one_ec <- res[grep("^(\\d+\\.){3}\\d+$", enzyme)]

# Reactions with 1 EC number: get matching gene names in human ---------------
# The for loop is really slow!

# Function for getting enzyme - gene information
annotate.ec <- function(df) {
  ECres <- matrix(nrow = nrow(df), ncol = 3)
  colnames(ECres) <- c("EC", "name", "genes")
  df$enzyme <- paste0("ec:", df$enzyme)
  for (i in seq(from=1, to=nrow(df), by=10)) {
    message(stringr::str_glue("Getting enzymes {i} to {min(i + 9, nrow(df))}."))
    ec <- keggGet(df$enzyme[i:min(i + 9, nrow(df))])
    for (j in seq_along(ec)) {
      # Skip entries without matching gene names in human
      ecGene <- ec[[j]]$GENES[grep("^HSA: ", ec[[j]]$GENES)]
      if (length(ecGene) != 0) {
        ECres[i + j - 1, ] <- c(
          ec[[j]]$ENTRY,
          paste(ec[[j]]$NAME, collapse = " "),
          ec[[j]]$GENES[grep("^HSA: ", ec[[j]]$GENES)]
        )
      }
    }
  }
  df <- as.data.table(cbind(df, ECres))
  df <- df[!is.na(genes),.(reaction, equation, EC, name, genes)]
  return(df)
}

hsa_rxn_one_ec <- annotate.ec(rxn_one_ec)

# Reactions with multiple EC numbers: split them and get annotations ---------
#
# for each entry
rxn_w_multiple_ec <- rxn_w_multiple_ec %>%
  mutate(enzyme = strsplit(as.character(enzyme), ",")) %>%
  unnest() %>%
  filter(enzyme != "")


hsa_rxn_w_multiple_ec <- annotate.ec(rxn_w_multiple_ec)
# Save temp results as they take a while to run
save(res, hsa_rxn_one_ec, hsa_rxn_w_multiple_ec,
     file="./KEGG_Proton_reactions.RData")

# Convert entrez IDs to Ensembl gene IDs -------------------------------------


load("KEGG_Proton_reactions.RData")
annot <- data.table::fread("~/CSBL_shared/ID_mapping/Ensembl_symbol_entrez.csv") 

df <- rbind.data.frame(hsa_rxn_one_ec, hsa_rxn_w_multiple_ec) %>%
  mutate(genes = gsub("^HSA: ", "", genes)) %>%
  mutate(genes = strsplit(as.character(genes), " ")) %>%
  unnest() %>%
  filter(genes != "") %>%
  mutate(entrez = gsub("(^\\d+).*$", "\\1", genes),
         genes = gsub("^.*\\(([^(]+)\\)", "\\1", genes))

annot <- annot %>%
  dplyr::select(entrez = entrezgene, ensembl = ensembl_gene_id) %>%
  filter(!is.na(entrez)) %>%
  mutate_if(is.numeric, as.character) %>%
  filter(entrez %in% df$entrez)

# Prepare df equations to get number of Protons(s) consumed / produced -----------
parse_equation <- function(df) {
  eq_left <- sapply(df$equation, function(x) strsplit(x, " <=> ")[[1]][1])
  # Only keep the substrate / product C00002 (ATP)
  eq_left <- sapply(strsplit(eq_left, " \\+ "), function(x) ifelse(
    length(x[grep("C00080", x)]) == 0, 0, x[grep("C00080", x)]))
  eq_left[eq_left == "C00080"] <- 1
  eq_left[grep("[^0-9]", eq_left)] <- sapply(eq_left[grep("[^0-9]", eq_left)],
                                             function(x) strsplit(x, " ")[[1]][1])
  
  eq_right <- sapply(df$equation, function(x) strsplit(x, " <=> ")[[1]][2])
  eq_right <- sapply(strsplit(eq_right, " \\+ "), function(x) ifelse(
    length(x[grep("C00080", x)]) == 0, 0, x[grep("C00080", x)]))
  eq_right[eq_right == "C00080"] <- 1
  
  res <- data.frame(lhs = eq_left, rhs = eq_right)
  res <- res %>%
    mutate_if(is.factor, as.character) %>%
    mutate(net_Proton = as.numeric(rhs) - as.numeric(lhs))  # NA introduced by n ATP and n-1 ATP
  return(res$net_Proton)
}

df <- df %>%
  inner_join(annot, by = "entrez") %>%
  dplyr::select(-c(genes, entrez))
df$net_Proton <- parse_equation(df)

hsa_rxn_w_nATP <- df[is.na(df$net_Proton), ]
rxn_wo_ec <- res[enzyme == ""]
rxn_w_ec_class <- res[grep("-", enzyme)]
save(rxn_w_ec_class, rxn_wo_ec, hsa_rxn_w_nATP, file="error_rxns.RData")
df <- df[!is.na(df$net_Proton), ]
data.table::fwrite(df, file="./hsa_Proton_rxns.csv", row.names = F, col.names = T)

# proton producing genes. 
#dat <- data.table::fread("hsa_Proton_rxns.csv", data.table = FALSE)
