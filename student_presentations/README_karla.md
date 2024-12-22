# ReDrugGraph
Term project for Molecular Informatics (CSCI 5580) 
## **Execute Scripts in this order**


## **1. get_disgenes.Rmd**
**Input**: Disease gene association data fromm DisGeNET and MONDO obo file.
**Output**: Propagated genes associations for all collected DisGeNET data.

**Files needed**: 
1. mondo_2024_03.obo
2. all_gene_disease_associations.tsv

**Functions**
1. map_mondo_to_ontology (DisGeNET df, obo_file): Maps MONDO disease identifiers to the disease identifiers contained in DisGeNET database(UMLS).
2. get_ances (list of disease identifiers, obo file): gets the ancestor term of all the terms found in the DisGeNET data set.
3. propagate_genes (list of disease identifiers from DisGeNET df, ancester df, DisGeNET data with mapped MONDO identifiers): function that propagates genes from children terms up to parent terms.

**Output**: 
1. not_propagated_disgenet_genes.tsv
2. tb_propagated_disgenet_genes.tsv

## **2. get_drug_data.RMD **
**Input**: Raw data downloaded from DGIdb website.
**Output**: Processed drug data containing CHEMBL drug identifiers.

**Files needed**: 
1. interactions.tsv

**Output**: 
1. dgidb_final_drugdata.tsv

## **3. karla_vela_gae_code.ipynb **
**Input**: BioGRID data, propagated genes obtained from DisGeNET, and processed DGIdb drug data
**Output**: performance metric plots for each GAE and checkpoint data files to provide for class presentation.

**Files needed**: 
1. propagated_disgenet_genes.tsv
2. biogrid_network.txt
3. dgidb_final_drugdata.tsv

**Functions**

1. def plot_roc_curve(title, model, data): plots ROC curve for GAE and VGAE models. 

2. plot_training_stats(title, losses, test_auc, test_ap, train_auc, train_ap): plots the AP and AUC after each epoch.

3. get_edge_dot_products( data, model, num_dz_nodes): computes the dot products between the encoded disease and gene nodes for the VGAE and GAE models.

4. get_ranked_edges(data_object, model, num_dz_ndoes): ranks the edges based on the dot products of the encoded nodes and returns the ranked edge list and dot products.

5. Class GCNEncoder(torch.nn.Module): defines and initializes the GAE model.

6. gae_train(train_data, gae_model, optimizer): trains the GAE model on the provided training data.

7. gae_test(test_data, gae_model): evaluates the GAE model on the provided test data.
 
 ## **3. karla_vela_class_version.ipynb**

 - This is a shortened version with added background and images for the class presentation.

**Input**: BioGRID data, propagated genes obtained from DisGeNET, and processed DGIdb drug data
**Output**: performance metric plots for each GAE and checkpoint data files to provide for class presentation.

**Files needed**: 
1. disgene_final_dataset.tsv
2. new_filtered_disgenet.tsv
3. filtered_biogrid_df.tsv

