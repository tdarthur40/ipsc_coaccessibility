#!/path/to/R-4.0.1

source("functions.R")
source("packages.R")
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(igraph))

args = commandArgs(trailingOnly=TRUE)
task = as.numeric(args[1])

edges_df="path_to_edges"
vertices_df="path_to_vertices"

### Two column dataframe of co-accessible ATAC-seq peaks
rnm_edges    = fread(edges_df, sep="\t",data.table=F)

### Dataframe with vertices in first column
rnm_vertices = fread(vertices_df,sep="\t",data.table=F)

### Make network
rnm_net = graph_from_data_frame(rnm_edges, directed = FALSE, vertices = rnm_vertices)

### To optimize the networks, I used a dataframe with permutations of parameter values (resolution, beta, n_iterations)
test_parameters = fread("analysis/networks/expression/optimize/scripts/test_parameters.txt",sep="\t",data.table=F)
tps = test_parameters[task,]

res = tps$resolution
beta = tps$beta
n_iterations = tps$n_iterations

### Perform Leiden Clustering on the co-accessibility network to find modules
leid_rnm = cluster_leiden(rnm_net, objective_function =  "modularity", 
           resolution_parameter = res, beta = beta, n_iterations = n_iterations)

## Calculate Modularity 
Q = modularity(rnm_net, membership(leid_rnm))

### Assign Vertices to Leiden modules (RNMs) 
membs = data.frame(peak_id = leid_rnm$names, original_module = leid_rnm$membership)
### modules with 500 or more elements
membs500 = membs[ membs$original_module %in% names(table(membs$original_module)) [ table(membs$original_module) > 500],]

### Create dataframe with Leiden performance for each parameter permutation to assess module quality
opt_df = data.frame(resolution = res , beta = beta, n_iterations =n_iterations, modularity = Q,
                    n_mods = length(unique(membs500$original_module)),
                   frac_ModPeaks = nrow(membs500)/nrow(membs))
opt_df$task_index = task

### Number RNMs by size
tabled = membs %>% group_by(original_module) %>% summarise(n_elements = length(original_module))
sorted = tabled[ order(-tabled$n_elements),]
sorted$new_module = seq(1,length(unique(sorted$original_module)),1)
sorted$task_index = task

verts1 = merge(membs,sorted)
mods2test = unique(names(table(verts1$new_module))[ table(verts1$new_module) > 500])

### Calculate intramodular degree
intramodular_degree = as.data.frame(rbindlist(lapply(mods2test, function(x){
    verts2test = verts1[ verts1$new_module == x,]
    edges2test = rnm_edges [ rnm_edges$peak1 %in% verts2test$peak_id & rnm_edges$peak2 %in% verts2test$peak_id ,]
    modtable = as.data.frame(table(c(edges2test$peak1, edges2test$peak2)))
    colnames(modtable) = c("peak_id","intramodular_degree")
    modtable$peak_id = as.character(modtable$peak_id)
    return(modtable)
})))

genome_degree = as.data.frame(table(c(rnm_edges$peak1, rnm_edges$peak2)))
colnames(genome_degree) = c("peak_id","genomewide_degree")
degree = merge(genome_degree,intramodular_degree, all=TRUE)

verts2 = merge(verts1, degree,all =TRUE)
verts2$original_module = NULL
verts2$n_elements = NULL
verts3 = verts2[ order(verts2$new_module, -verts2$intramodular_degree),]
colnames(verts3) = c("peak_id","module","task_index","genomewide_degree","intramodular_degree")

### Output the RNM assignments
fwrite(verts3,paste0("memberships/",task,".tsv"), sep="\t",row.names=F,quote=F)

### Output the model performance/optimization
fwrite(opt_df,paste0("parameter_optimization/",task,".tsv"), sep="\t",row.names=F,quote=F)
