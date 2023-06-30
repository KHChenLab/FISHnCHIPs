# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 14:07:10 2023

@author: chength
"""

import re

import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram

import igraph
from igraph import Graph
import leidenalg




def evaluate_gene_centric_panel(cell_bit_matrix_log1p, gcm_cellScaled, sig_cut, panel_info, nprobe = "nprobe"): 
  
  #Initialize dataframes
  signal_summary = pd.DataFrame()
  cell_counts = pd.DataFrame(np.zeros((0,3)))
  cell_counts.columns = ["program", "signal_cells", "signal_cutoff"]
  
  #Each unique value in the "Annotation" column is a program
  for program in np.unique(panel_info["Annotation"]):
    print(program)
    panel_subset = panel_info[panel_info["Annotation"] == program] #Extract data that for this iteration of program
    genes = panel_subset["gene"]
    print(genes)
    signal_cells = cell_bit_matrix_log1p[cell_bit_matrix_log1p[program] >= sig_cut] #filter cells with signal higher than the cutoff from the program(bit) in this iteration
    
    #label cells in cell_bit_matrix as signal cell or background noise
    tmp=[]
    for i in cell_bit_matrix_log1p[program]:
        if i >=sig_cut:
            tmp.append("signal")
        else:
            tmp.append("background")
    cell_bit_matrix_log1p[("labs_" + str(program))] = tmp
    
    cn = pd.DataFrame({
        "program": [program],
        "signal_cells": [signal_cells.shape[0]],
        "signal_cutoff": [sig_cut]
        })
    cell_counts=pd.concat((cell_counts,cn), axis=0)
    
    
    program_cell_counts = gcm_cellScaled.filter(list(genes), axis='index').T #extract rows that have the same name as genes from gcm_csllScaled
    program_cell_counts["cluster"] = cell_bit_matrix_log1p[("labs_" + str(program))] #create cluster column that indicate whether it is signal cell or background noise
    signal_program_genes = program_cell_counts.groupby("cluster").mean() #calculate mean of each column for signal cell and background cell separately
    signal_program_genes = signal_program_genes.T.reset_index(names = 'gene').drop(['background'], axis = 1)#converts column names into a column (gene) 
    
    if nprobe not in panel_info.columns:
      panel_info[nprobe] = np.ones(panel_info.shape[0])

    probe_info = panel_info[panel_info['gene'].isin(genes)].filter(items=['gene','nprobe']) #for gene columns that are in genes, extract gene and nprobe data into probe info column
    signal_program_genes = pd.merge(signal_program_genes, probe_info, on = "gene")
    #include columns with signal*nprobe values and cumulative sums
    signal_program_genes["signalxnprobe"] = signal_program_genes['signal'] * signal_program_genes['nprobe']
    signal_program_genes = signal_program_genes.sort_values(by = 'signalxnprobe', ascending=False)
    signal_program_genes.reset_index(drop = True, inplace = True)
    signal_program_genes["cum_signal"] = signal_program_genes['signal'].cumsum()
    signal_program_genes["cum_signalxnprobe"] = signal_program_genes['signalxnprobe'].cumsum()

    #program_report contains additional columns containing calculated values
    program_report = signal_program_genes.copy()
    program_report["program"] = np.repeat(program, program_report.shape[0])
    program_report["cum_signal_norm"] = program_report["cum_signal"] / float(program_report.loc[0, "signal"]) #normalize value by first value in the signal column
    program_report["cum_signal_nprobe_norm"] = program_report["cum_signalxnprobe"]/float(program_report.loc[0, "signalxnprobe"]) #normalize value by first value in the signalxnprobe column
    program_report = program_report.reset_index().rename(columns={"index": "gene_seq"}) #rename index column as gene_seq
    signal_summary = pd.concat((signal_summary, program_report), axis=0)

  return(signal_summary)





def get_cell_bit_matrix (gcm_cellScaled, panel_info, bitcol = "gene_group", nprobes = "nprobes"):
  subset_counts = gcm_cellScaled.filter(list(panel_info["gene"]), axis='index') #filter rows where the index is equivalent to values in the gene column of panel_info
  if nprobes in panel_info.columns:
    print("cell counts * nprobes ...")
    subset_counts = subset_counts * panel_info[nprobes]

  else:
    print("Number of probes not provided. Default: 1")
  
  tmp_panel = panel_info.copy()
  tmp_panel.set_index("gene", inplace = True)  
  subset_counts=pd.concat((subset_counts, tmp_panel[bitcol]), axis=1)
  subset_counts.rename(columns={bitcol : "bit_info"}, inplace=True)
  subset_counts = subset_counts.groupby("bit_info").sum() #sum values in each column for each gene group
  
  #change gene group to index column, drop it, then transpose the matrix
  cell_bit_matrix = subset_counts.T
  return(cell_bit_matrix)






def get_cnmf_panel (top_factors, correlation_matrix, min_corr = 0.02, max_corr = 0.3, max_n_neg = 1, max_n_high = 1, removedup = True, remove_smallgroup = True): 
  topn_genes = top_factors.drop(top_factors.filter(regex='_val').columns, axis=1) #drop columns with _val in its name
  panel = pd.DataFrame(np.zeros((0,2)))
  panel.columns = ["gene", "Annotation"]
  for f in topn_genes.columns:
    print(f)
    top50genes = topn_genes[f]
    genes = top50genes[top50genes.isin(correlation_matrix.columns)] #filter top50genes by values that are also present in correlation matrix rownames
    correlation_matrix.index = correlation_matrix.columns
    corrMat = correlation_matrix.loc[genes, genes]
    tmp = corrMat
    lowq = []
    while tmp.min().min() <= min_corr:
      subjects=[]
      for col in tmp.columns:
          subjects += tmp.index[tmp[col] == tmp.min().min()].tolist() #returns rows where tmp is minimum
      for g in subjects:
        n_neg = (corrMat.loc[g, :] < min_corr).sum()
        if (n_neg > max_n_neg):
          lowq.append(g)
          lowq = list(set(lowq))
          highq = list(set(genes) - set(lowq)) #returns a numpy array containing the elements in genes that are not in lowq.
          tmp = tmp.loc[highq, highq]

    #while the maximum value in tmp that is less than 0.999999999 is more than the max_cor, the loop will keep running
    while (tmp[tmp < 0.999999999].max(numeric_only=True).max(numeric_only=True)) > max_corr: 
      subjects=[]
      for col in tmp.columns:
          subjects += tmp.index[tmp[tmp < 0.999999999][col] == tmp[tmp < 0.999999999].max(numeric_only=True).max(numeric_only=True)].tolist() #returns rows where tmp is maximum and lower than 1
      n_highq = tmp.shape[0]
      for g in subjects:
        n_high = (corrMat.loc[g, :] > max_corr).sum()
        if (n_high > max_n_high):
          lowq.append(g)
          lowq = list(set(lowq))
          highq = list(set(genes) - set(lowq)) #returns a numpy array containing the elements in genes that are not in lowq.
          tmp = tmp.loc[highq, highq]

      n_highq_new = tmp.shape[0]
      if (n_highq == n_highq_new):
        break

    print(str(f) + "- new min_corr: " + str(tmp.min().min()))
    print(str(f) + "- new max_corr: " + str(tmp[tmp < 0.999999999].max(numeric_only=True).max(numeric_only=True)))
    if ((len(highq) > 1) and (len(lowq) > 1)):
      highqMat = corrMat.loc[highq, highq]
      lowqMat = correlation_matrix.loc[lowq, lowq]
      mats = list((highqMat, lowqMat)) #Putting both dataframes into a list 
      f_df = pd.DataFrame(np.zeros((0,2)))
      f_df.columns = ["gene", "Annotation"]
      for m in range(len(mats)):
        mat = mats[m]
        if m == 0:
            name = 'highq'
        else:
            name = 'lowq'
       
        #pdist computes the pairwise Euclidean distances between the rows of the matrix
        #linkage performs hierarchical clustering on the pairwise distances using the complete-linkage method
        dend = linkage(pdist(mat, metric='euclidean'), method='complete') 
        k = 1
        dend_cluster = fcluster(dend, t=k, criterion='maxclust') #assigns each observation to a cluster based on the dendrogram, t is the number of cluster
       
        # Color branches of the dendrogram
        plot=dendrogram(dend, labels = mat.columns, color_threshold=k, above_threshold_color='k')
        
        # Get labels of the dendrogram leaves and sort the cluster accordingly
        dend_labels = np.array(plot['ivl'])
        dend_cluster = dend_cluster[np.argsort(np.argsort(dend_labels))]

        m_df = pd.DataFrame({'gene': dend_labels, 
                             'Annotation': [str(f) + '_' + name + '_' + str(c) for c in dend_cluster]
                             })
        
        f_df = pd.concat((f_df,m_df), axis=0)

    panel = pd.concat((panel, f_df), axis=0)
    #Filter out rows where the Annotation column do not contain '_highq'
    panel = panel[panel['Annotation'].apply(lambda x: bool(re.search('_highq', x)))] 

  if removedup: 
    # Get the duplicate genes and filter them from the panel DataFrame
    dup_genes = panel.groupby('gene').size().reset_index(name='count').query('count >= 2')['gene']
    panel_clean = panel[~panel['gene'].isin(dup_genes)] #~ is used to invert the operation, so that only genes not found in dup_genes is kept


  if remove_smallgroup:
    #Get the Annotations that appears less than 3 times and filter them from the panel_clean DataFrame
    small_program = panel.groupby('Annotation').size().reset_index(name='count').query('count <= 2')['Annotation']
    panel_clean = panel_clean[~panel_clean['Annotation'].isin(small_program)]
  
  panel_clean.reset_index(drop = True, inplace = True)
  return(panel_clean)





#Calculate correlation matrix from data and save file
def get_correlation_matrix (data, out_rds_filename = "./correlation_matrix.csv"): 
  print("Computing GGC")
  t_data = data.T #assuming data is a dataframe
  #t_data <- Matrix::t(data) is this transpose supposed to be here?
  correlation_matrix = np.corrcoef(t_data, rowvar=0) #calculate correlation coefficient matrix from data
  correlation_matrix = pd.DataFrame(correlation_matrix, columns=data.index, index=data.index)

  print("Done.")
  if out_rds_filename != None:
    print("Saving the correlation matrix at: " + out_rds_filename)
    correlation_matrix.to_csv(out_rds_filename)
  
  return(correlation_matrix)







def get_cumulative_signals(gcm_cellScaled, panel_info, scRNA_celltype, multiply_nprobe = True):
    subset_counts = gcm_cellScaled.loc[np.unique(panel_info['gene']), :].T #extract rows that are indexed as genes found in panel info
    #create cluster column in subset_count with the same values as celltype in scRNA_celltype
    subset_counts = subset_counts.reset_index(drop=True)
    subset_counts = pd.concat((subset_counts, scRNA_celltype['Celltype']), axis=1) 
    subset_counts.rename(columns={'Celltype' : "cluster"}, inplace=True)
    
    if (multiply_nprobe):
      subset_counts = subset_counts * panel_info['nprobes']
    
    signal_in_clusters = subset_counts.groupby('cluster').mean() #calculate mean of each column for each cluster
    #signal_in_clusters.set_index("cluster", inplace=True) #python groupby function automatically sets the target column as index
    gene_sigs = signal_in_clusters.T
    gene_sigs = gene_sigs.reset_index(names="gene") #genes were the column name. After transposing, it is now the index names. Make it into a column
    
    cluster_names = np.unique(panel_info["cluster"])
    panel_sigs = pd.DataFrame()
    for c in cluster_names:
      print(c)
      panel_c = panel_info.loc[panel_info['cluster'] == c, :] #filter rows where cluster is the same as the current cluster in iteration
      panel_c_sigs = gene_sigs[gene_sigs['gene'].isin(panel_c['gene'])] #extract genes from gene_sigs that are present in panel_c
      panel_c_sigs.set_index('gene', inplace = True)
      panel_c_cumsigs = panel_c_sigs.cumsum() #iterates over rows and finds the cumulative sum in each column 
      panel_c_cumsigs.index = panel_c_sigs.index
      panel_c_cumsigs['cum_signal'] = panel_c_cumsigs.loc[:, c]
      panel_c_cumsigs['cum_signal_norm'] = panel_c_cumsigs['cum_signal']/float(panel_c_cumsigs["cum_signal"][0])
      tmp = panel_c_cumsigs.drop(['cum_signal','cum_signal_norm'], axis=1)
      panel_c_cumsigs['noise_2max'] = tmp.apply(lambda x: sorted(x, reverse=True)[1], axis=1)
      panel_c_cumsigs['noise_2max_cum'] = panel_c_cumsigs['noise_2max'].cumsum()
      panel_c_cumsigs['noise_sum'] = tmp.apply(lambda x: sum(x) - sorted(x, reverse=True)[1], axis=1)
      panel_c_cumsigs['ssr_2max'] = panel_c_cumsigs['cum_signal']/panel_c_cumsigs['noise_2max']
      panel_c_cumsigs['ssr_2max_cum'] = panel_c_cumsigs['cum_signal']/panel_c_cumsigs['noise_2max_cum']
      panel_c_cumsigs['ssr_sum'] = panel_c_cumsigs['cum_signal']/panel_c_cumsigs['noise_sum']
      #panel_c_cumsigs['gene'] = panel_c_cumsigs.index
      panel_c_cumsigs['gene_seq'] = np.arange(1, len(panel_c_cumsigs) + 1) #create column to sequence the genes in order
      panel_c_cumsigs = pd.merge(panel_c_cumsigs, panel_c, on = "gene") 
      panel_c_cumsigs = panel_c_cumsigs.sort_values(by='gene_seq')
      panel_sigs = pd.concat((panel_sigs, panel_c_cumsigs), axis=0)
    
    panel_sigs.reset_index(drop = True, inplace = True)
    return(panel_sigs)




def get_panel(correlation_matrix, ref_genes = None, min_corr = 0.7, min_ngenes = 3): 
  topCorr_mask = correlation_matrix
  topCorr_mask.index = topCorr_mask.columns
  topCorr_mask = (correlation_matrix >= min_corr) #create a boolean mask of correlations that are higher than min_corr
  n_high_corr_genes = topCorr_mask.sum(axis=1, skipna=True) #sum the number of genes with high correlation for each column

  n_high_corr_genes.name = "ngenes"
  if ref_genes is None:
    filt_corr_genes = n_high_corr_genes.loc[lambda x : x >= min_ngenes].index
    #cat(paste0(length(filt.corr.genes), " genes found correlated to at least ", min.ngenes, " genes with correlation > ", min.corr))
    #cat("\n")
    panel_genes = filt_corr_genes
    panel_genes.name = "gene"
  
  else:
    #initialize panel genes dataframe if ref_genes is present
    panel_genes = pd.DataFrame(np.zeros((0,4)))
    panel_genes.columns = ["gene", "corr", "ref", "cluster"]

    for c in np.unique(ref_genes["cluster"]):
      #extract genes that belong to the current iteration of cluster
      target = ref_genes.loc[ref_genes['cluster'] == c, :] 
      g = target['gene']
      #cat("\n")
      print(str(c), ": ", str(len(g)), " markers provided")
      missing_genes = g[~(g.isin(correlation_matrix.columns))] #identify genes that are not present in correlation matrix
      if len(missing_genes) > 0:
        print(missing_genes + " not found")
      
      corr_g = correlation_matrix.loc[:, g] #extract correlation matrix for selected genes
      g = list(set(g).difference(set(missing_genes))) #Remove missing_genes genes from g if present
      #select only rows that have at least one value in each column specified by g that is greater than min_corr
      topCorr_mat = corr_g.loc[corr_g[g].gt(min_corr).any(axis=1), g]
      print(str(len(topCorr_mat)) + " genes found (>" +str(min_corr) + ")")
      
      if min_ngenes != None:
          
        if (len(topCorr_mat) < (min_ngenes + len(g))):
          print("Not enough genes found highly correlated.")
          print("Getting top " + str(min_ngenes) + " genes highly correlated to the markers... ")
          print("The final panel may include super low correlated genes.")
          #Initialize dataframe
          topCorr_mat = pd.DataFrame(np.zeros((0, len(g) + 1)))
          topCorr_mat.columns = g + ["gene"]
          
          for i in range(len(g)):
            corr_g.rename(columns={corr_g.columns[i] : "corr"}, inplace=True) #rename each column to corr
            g_topN = corr_g.sort_values('corr', ascending=False).iloc[:min_ngenes + 1] #sort values in g_topN by column corr in descneding order and extract top n rows
            g_topN.rename(columns={g_topN.columns[i] : g[i]}, inplace=True) #rename the renamed column back to gene name
            g_topN = g_topN.reset_index().rename(columns={'index': 'gene'}) #make gene column into index
            topCorr_mat = pd.concat((topCorr_mat, g_topN), axis=0)
            topCorr_mat = topCorr_mat.drop_duplicates() #remove any duplicates from concatenated dataframe
            corr_g.rename(columns={corr_g.columns[i] : g[i]}, inplace=True)

        else:
          topCorr_mat.rename(columns={topCorr_mat.columns[0] : "corr"}, inplace=True) #rename first column of dataframe as corr
          topCorr_mat = topCorr_mat.sort_values('corr', ascending=False) 
          topCorr_mat.columns = g
          topCorr_mat = topCorr_mat.reset_index().rename(columns={'index': 'gene'})
        
        if (len(g) > 1): 
            #obtain maximum correlation from each column in topCorr_mat that has the same colname as g and assign it to a new column "corr"
            topCorr_mat["corr"] = topCorr_mat[g].apply(max, axis=1) 
            #obtain the index of the maximum correlation from each column in topCorr_mat that has the same colname as g and assign it to a new column "ref"
            topCorr_mat["ref"] = topCorr_mat[g].apply(lambda x: g[x.argmax()], axis=1)
        
        else:
          topCorr_mat["corr"] = topCorr_mat[g]
          topCorr_mat["ref"] = np.repeat(g, topCorr_mat.shape[0])
        
      
      else:
        topCorr_mat = topCorr_mat.reset_index().rename(columns={'index': 'gene'})
        topCorr_mat["ref"] = np.repeat(g, topCorr_mat.shape[0])
        topCorr_mat["corr"] = topCorr_mat[g]

      topCorr_mat['cluster'] = np.repeat(c, topCorr_mat.shape[0])
      topCorr_mat = topCorr_mat.sort_values('corr', ascending=False)[['gene', 'corr', 'ref', 'cluster']]
      panel_genes = pd.concat((panel_genes, topCorr_mat), axis=0)
    
    print(str(panel_genes.shape[0]) + " FISHnCHIPs panel genes found for all clusters")
    panel_genes.reset_index(drop = True, inplace = True)
  
  return(panel_genes)








def get_filtered_genes (data, min_cells = 0.001, max_cells = 0.5, min_expression = 0, filt_name_pattern = None, ignore_case = True): 
  #min_cells = 0.001 * data.shape[1]   
  #max_cells = 0.5 * data.shape[1]
  gene_in_ncells = (data > min_expression).sum(axis=1) #sum the number of cells where each gene's expression value is above min_expression
  #colnames(gene.in.ncells) <- c("ncells")
  print("Binarize on gene expression level - Done")
  #extract rownames where values in the 'ncells' column are larger than min_cells and smaller than max_cells
  filt_genes = gene_in_ncells[((gene_in_ncells>= min_cells) & (gene_in_ncells <= max_cells))].index.tolist()
  print("Filtering on gene distribution in cells - Done ")
  print(str(len(filt_genes)) + " genes qualified ...")
  if filt_name_pattern != None:
    pattern = re.compile(filt_name_pattern, flags=re.IGNORECASE) #establish pattern as filt_name_pattern using re.compile
    pattern_genes = [gene for gene in data.index if pattern.search(gene)] #search for pattern in data.index. Return True if present and False if absent
    print("Filtering genes by naming pattern...")
    print(str(len(pattern_genes)) + " genes in the dataset with the pattern.")
    filt_genes = set(filt_genes).difference(set(pattern_genes)) #identify genes that are present in filt_genes but absent in pattern_genes
    print(str(len(filt_genes)) + " genes passing the filters.")

  else:
    pattern_genes = None

  return(list(filt_genes), gene_in_ncells, pattern_genes)







def leiden_corr_matrix (correlation_matrix, edge_corr_threshold = 0.5,
                        partition_type = (leidenalg.RBConfigurationVertexPartition, leidenalg.ModularityVertexPartition, leidenalg.RBERVertexPartition, leidenalg.CPMVertexPartition,  
                                           leidenalg.SignificanceVertexPartition, leidenalg.SurpriseVertexPartition, leidenalg.CPMVertexPartition.Bipartite), 
                        seed = None, n_iterations = -1, max_comm_size = 0, degree_as_node_size = False):

    #Generates a graph from its weighted adjacency matrix  
    cor_g = Graph.Weighted_Adjacency(correlation_matrix, mode = "undirected", attr = "correlation",  loops = False)
    
    
    #Export edges from generated graph with attributes to DataFrame
    cor_edge_list = Graph.get_edge_dataframe(cor_g)
    cor_vert = Graph.get_vertex_dataframe(cor_g) #get vertex names
    #cor_edge_list['source'].replace(cor_vert['name'], inplace=True)
    #cor_edge_list['target'].replace(cor_vert['name'], inplace=True)
    #igraph.plot(cor_g, vertex_label = cor_vert.squeeze())
    print(cor_edge_list.shape)
    
    #Extract rows where correlation is larger than the threshold. Use a copy of cor_edge_list to prevent double indexing warning
    only_sig = cor_edge_list.loc[cor_edge_list['correlation'] > edge_corr_threshold, :].copy()
    print(only_sig.shape)
    
    #create new vertex dataframe from filtered dataset. Use a copy of cor_vert to prevent double indexing warning
    new_vert = cor_vert[cor_vert.index.isin(only_sig['target']) | cor_vert.index.isin(only_sig['source'])].copy()
    new_vert['enumerate'] = np.arange(0, len(new_vert))
    only_sig['source'].replace(new_vert['enumerate'], inplace=True)
    only_sig['target'].replace(new_vert['enumerate'], inplace=True)
    new_vert.set_index('enumerate', inplace = True)
    
    #Generates graph from dataframe containing edges above the threshold
    new_g = Graph.DataFrame(only_sig, directed = False, vertices = new_vert, use_vids = True)
    #igraph.plot(new_g, vertex_label = new_vert.squeeze())
    
    #Constructs a graph based on an adjacency matrix from the given file
    adjacency_matrix = Graph.get_adjacency(new_g)
    print(adjacency_matrix.shape)
    
    
    #Clustering via leiden algorithm using adjacency graph
    partition = leidenalg.find_partition(new_g, seed = seed, partition_type = partition_type, 
                                         n_iterations = n_iterations, max_comm_size = max_comm_size)
    
    # get the membership list
    membership = partition.membership
    leiden_partition = pd.concat((new_vert, pd.DataFrame(membership)), axis = 1)
    leiden_partition.columns = ["node", "leiden_cluster"]  
    leiden_partition = leiden_partition.sort_values(by = 'leiden_cluster')
    print(partition.summary())

    return(leiden_partition, new_g)



