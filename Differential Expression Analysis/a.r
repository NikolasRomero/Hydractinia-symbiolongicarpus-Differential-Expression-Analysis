no_enriched_genes <- anti_join(annotated_degs, id2prot, by = "GeneID")
no_enriched_ids <- degs[degs$GeneID %in% no_enriched_genes$GeneID, ]
no_enriched_domains <- domains[domains$Protein_id %in% no_enriched_ids$protein_id, ]
results_no_enriched <- process_terms(terms_names, terms, no_enriched_domains, columns)

matching_indices_no_enriched <- grep(common_part, names(results_no_enriched))
matching_objects_no_enriched <- results_no_enriched[matching_indices_no_enriched]
all_terms_no_enriched <- unlist(matching_objects_no_enriched)
um_unique_terms_no_enriched <- length(unique(all_terms_no_enriched))

result_df_no_enriched <- data.frame()

for (i in seq_along(matching_objects_no_enriched)) {
  name <- names(matching_objects_no_enriched)[i]
  num_terms <- length(matching_objects_no_enriched[[i]])
  row <- data.frame(name = name, num_terms = num_terms)
  result_df_no_enriched <- rbind(result_df_no_enriched, row)
}

result_df_no_enriched$name <- gsub("unique_", "", result_df_no_enriched$name)
result_df_no_enriched$name <- gsub("_description", "", result_df_no_enriched$name)
result_df_no_enriched$name <- gsub("_domain_name", "", result_df_no_enriched$name)
sum_df_annot_no_enriched <- aggregate(num_terms ~ name, data = result_df_no_enriched, sum)

plot_domains_no_enriched <- result_df_no_enriched %>%
  filter(`num_terms` != 0)

r <- ggplot(plot_domains_no_enriched, aes(x = name, y = `num_terms`)) +
  geom_bar(stat = "identity", position = "dodge", fill = "#ee8037") +
  labs(title = "Not Enriched Genes IDCP Analysis", x = "Domain", y = "Number of IDCP") +
    theme_minimal() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
r

no_enriched_description_i <- grep(text_for_list1, names(matching_objects_no_enriched))
no_enriched_domain_name_i <- grep(text_for_list2, names(matching_objects_no_enriched))
no_enriched_description <- matching_objects_no_enriched[no_enriched_description_i]
no_enriched_domain_name <- matching_objects_no_enriched[no_enriched_domain_name_i]
names(no_enriched_description) <- gsub(text_for_list1, "", names(no_enriched_description))
names(no_enriched_domain_name) <- gsub(text_for_list2, "", names(no_enriched_domain_name))

no_enriched_list <- list()

for (name in names(no_enriched_description)) {
  no_enriched_list[[name]] <- c(no_enriched_description[[name]], no_enriched_domain_name[[name]])
}

names(no_enriched_list) <- gsub(text_for_list3, "", names(no_enriched_list))

original_list <- no_enriched_list

no_enriched_list <- Filter(function(x) any(x != 0), no_enriched_list)

extract_no_enriched_rows <- function(no_enriched_list, domains) {
  
  extracted_no_enriched_list <- list()
  
  for (name in names(no_enriched_list)) {
    extracted_no_enriched_rows <- domains[domains$Protein_id %in% no_enriched_list[[name]], ]
    extracted_no_enriched_list[[name]] <- extracted_no_enriched_rows
  }
  
  return(extracted_no_enriched_list)
}


extracted_no_enriched_list <- extract_no_enriched_rows(no_enriched_list, domains)
