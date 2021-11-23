##################################################################
#
#   R script for the data analysis as presented in
#   Medina Paz et al. (2021) in International Journal of Molecular Sciences
# 
###
#   for more information:
#   Francisco Medina Paz, francisco.medina@cinvestav.mx
#
#   11-13-2021
#
###################################################################

#################### Prepare R Environment ####################
rm(list=ls())

setwd("/Users/fmedi/Documents/CINVESTAV/Doctorado/articulo_maestria/Figuras_articulo_BMC/figuras_all_OTUs/Correcciones_13Nov_MPMI/")
options(max.print=999999, scipen=3, digits=5)

######### load R packages
library (phyloseq) 
library (ggplot2)
library (data.table)
library (plyr)
library (scales)
library (grid)
library (gplots)
library (VennDiagram)
library (biomformat)
library (vegan)
library (ggsignif)
library (edgeR)
if (Sys.getenv("JAVA_HOME")!="")
  Sys.setenv(JAVA_HOME="")
library (xlsx)

#################### Import Data ####################
mapFile <- "map.txt"
map <- import_qiime_sample_data(mapfilename=mapFile)
mapTab <- read.table(mapFile, header = TRUE, as.is=TRUE, sep="\t", comment.char="", row.names=1)
otus<-import_biom("new_otu_table_mod_final.biom", parseFunction= parse_taxonomy_greengenes)
otus
real_data<-merge_phyloseq(map,otus)
real_data<-filter_taxa(real_data, function(x) (sum(x) > 0), TRUE)
real_data

################### Defining parameters ############
topN=10
distance <- "bray"   ### It could be replace by unifrac, wunifrac dcpoa and jsd ##
Count_threshold<-10
tax_level <- c(2:6) ### 1 == "Kingdom", 2 == "Phylum", 3 == "Class", 4 == "Order", 5 == "Family", 6 == "Genus" 
FDR_thres <- 0.05
logFC_thres <- 2

################### Emule colors for ggplots2 ######
gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

################### Calculating Diversity Indices ###################
index_list<- c("Shannon", "Simpson", "Evenness", "Observed", "Chao1")
Diversity_Indices<- estimate_richness(real_data,split= TRUE, measures= index_list)
evenness<-c()
for (i in (1:nsamples(real_data))){
  evenness[i]<-c((Diversity_Indices[i,"Shannon"])/(log(Diversity_Indices[i,"Observed"],base = exp(1))))
}
Diversity_Indices[,"Evenness"]<-evenness
Diversity_Indices[,"SampleNames"]<-sample_data(real_data)[,4]
Diversity_Indices[,"Compartment"]<-sample_data(real_data)[,5]
Diversity_Indices[,"DevStage"]<-sample_data(real_data)[,6]
Diversity_Indices

################### Pairwise t-student test of Diversity indices ############# 
pairwise.t.test(Diversity_Indices$Shannon, Diversity_Indices$Compartment, p.adjust="bonferroni") #http://www.stat.columbia.edu/~martin/W2024/R3.pdf

################### Boxplot of Shannon Index per sample Compartment #################
ggplot(Diversity_Indices, aes(x=Compartment, y=Shannon, fill= Compartment)) + geom_boxplot() + theme_bw() + labs(title="Shannon Diversity Index", x="", y = "Shannon Index") + theme(legend.text = element_text (size=18), plot.title=element_text (size=22, hjust =0.5), axis.text = element_text(size=18, colour = "black"), axis.title=element_text(size=18), legend.title =element_text(size=18), axis.title.y = element_text(vjust = 2.5))+
  geom_signif(textsize = 8, comparisons = list(c("Endosphere", "Rhizosphere")), annotations = c ("***"), inherit.aes = T, map_signif_level=TRUE, y_position = c(7.3), test= t.test) + ylim(c(3,7.8)) + geom_jitter(alpha= 0.3, width = 0.1, size= 3) + scale_fill_manual (values=c("#ED8141", "#69b3a2")) +   scale_x_discrete(limits = rev(levels(factor(sample_data(real_data)$Compartment)))) + guides(fill = guide_legend(reverse=TRUE))

################### Principal Coordinate Analysis using Bray-Curtis distance #####
title= paste0 ("PCoA All Samples")
real_data.ord<-phyloseq::ordinate(real_data, method= "PCoA", distance= distance)
plot_ordination (real_data, real_data.ord,type = "samples", color = "Compartment", shape="DevStage", title = title, axes=c(1,2)) + theme_bw() +
  geom_point(size= 7) + theme(plot.title = element_text(size =18, hjust = 0.5),axis.title=element_text(size=18, vjust=-0.5), axis.text = element_text (size = 18, colour="black"), legend.title = element_text(size = 18), legend.text=element_text(size=18))+
  scale_color_manual(values = c( "#ED8141", "#69b3a2")) + guides(color = guide_legend(order = 1, reverse = TRUE), shape = guide_legend(order = 2))

############################ OTU Relative Abundance stacked barplot ###################
###################### Grouping samples ############
grouped_data<-merge_samples(real_data, group="Description")
mapFile <- "map_descriptions.txt"
map <- import_qiime_sample_data(mapfilename=mapFile)
rownames(map)<-rownames(sample_data(grouped_data))
sample_data(grouped_data)<-map
sample_data(grouped_data)

################### Calculating Relative Abundance ##################
grouped_data_rel<-transform_sample_counts(grouped_data, function(OTU) OTU / sum(OTU))
grouped_data_rel

################### Calculating Top 10 Most Abundant Taxa per Compartment ##################
Tax_level<- rank_names(grouped_data)[c(tax_level)]
samp_rel <- list()
i<-1
for (n in levels(factor(sample_data(grouped_data)$Compartment))){
  samp_rel[[n]] <- assign(n,subset_samples(grouped_data_rel,Compartment==(levels(factor(sample_data(grouped_data)$Compartment)))[i]))
  i<-(i+1)
}
top_samp <- paste0("top_",levels(factor(sample_data(grouped_data)$Compartment)))
top_samp_topN <- paste0(top_samp,"_",topN)
TopOTUs_Compartments <- paste0 ("TopOTUs_", levels(factor(sample_data(grouped_data)$Compartment)))
top_samp_rel <- list()
topN_samp <- list()
Compartments_TopOTUs <- list ()
i=1 ;j=1 ;k=1; b=1; c=1; h=1; d=1; q=1;
for (n in Tax_level){
  for (m in top_samp){
    top_samp_rel[[m]] <- assign(m,tax_glom(samp_rel[[i]],n, NArm=FALSE))
    i=i+1
  }
  for (c in top_samp_rel){
    for (d in (1:(length(tax_table(c)[,n])))){
      if (is.na (tax_table(c)[d,n]) == TRUE){
        for (j in c(1:6)){
          tax_table(c)[d,c(j:(grep(n,rank_names(c))))] <- ("Unclassified")
        }
      }
    }
    top_samp_rel[[top_samp[h]]] <- assign(top_samp[h],tax_glom(c, n, NArm=FALSE))
    h=h+1
  }
  print (top_samp_rel)
  for (l in (1:length(top_samp_rel))){
    if (any(names((sort(tapply(taxa_sums(top_samp_rel[[l]]), tax_table(top_samp_rel[[l]])[, n], sum),TRUE)[seq_len(topN)])/nsamples(grouped_data)) == "Unclassified") == TRUE)
    {
      Compartments_TopOTUs[[l]] <- (sort(taxa_sums(top_samp_rel[[l]]),TRUE)[seq_len(topN+1)])/(nsamples(samp_rel[[l]]))
      topN_samp[[l]] <- assign(top_samp_topN[l],prune_taxa(names(Compartments_TopOTUs[[l]]),top_samp_rel[[l]]))
      color_var <- sort(names ((sort(tapply(taxa_sums(top_samp_rel[[l]]), tax_table(top_samp_rel[[l]])[, n], sum),TRUE)[seq_len(topN+1)])/nsamples(grouped_data)))
    } else if (any(names((sort(tapply(taxa_sums(top_samp_rel[[l]]), tax_table(top_samp_rel[[l]])[, n], sum),TRUE)[seq_len(topN)])/nsamples(grouped_data)) == "Unclassified") == FALSE)
    {
      Compartments_TopOTUs[[l]] <- (sort(taxa_sums(top_samp_rel[[l]]),TRUE)[seq_len(topN)])/(nsamples(samp_rel[[l]]))
      topN_samp[[l]] <- assign(top_samp_topN[l],prune_taxa(names(Compartments_TopOTUs[[l]]),top_samp_rel[[l]]))
      color_var <- sort(names ((sort(tapply(taxa_sums(top_samp_rel[[l]]), tax_table(top_samp_rel[[l]])[, n], sum),TRUE)[seq_len(topN)])/nsamples(grouped_data)))
    }
  }
  for (a in (1:length(top_samp_rel))){
    for (k in (1:length(tax_table(top_samp_rel[[a]])[,n]))){
      if (any((as.vector(tax_table(top_samp_rel[[a]])[k,n]) == as.vector(tax_table(topN_samp[[a]])[,n]))) == F){
        if (as.vector(tax_table(top_samp_rel[[a]])[k,n]) == "Unclassified"){} else{
          for (o in (1:length(tax_table(topN_samp[[a]])[,n]))){
            if (as.vector(tax_table(top_samp_rel[[a]])[k,n]) != as.vector(tax_table(topN_samp[[a]])[o,n])){
              for (p in c(1:6)){
                tax_table(top_samp_rel[[a]])[k,c(p:grep(n,rank_names(grouped_data)))] <- ("Other")
              }
              break }
          }
        }
      }
    }
  }
  #  }
  for (m in top_samp){
    top_samp_rel[[m]] <- assign(m,tax_glom(top_samp_rel[[q]],n, NArm=FALSE))
    q=q+1
  }
  for (l in (1:length(top_samp_rel))){
    if (any(names((sort(tapply(taxa_sums(top_samp_rel[[l]]), tax_table(top_samp_rel[[l]])[, n], sum),TRUE)[seq_len(topN)])/nsamples(grouped_data)) == "Other", na.rm=T) == TRUE)
    {
      Compartments_TopOTUs[[l]] <- (sort(taxa_sums(top_samp_rel[[l]]),TRUE)[seq_len(topN+2)])/(nsamples(samp_rel[[l]]))
      print ((sort(tapply(taxa_sums(top_samp_rel[[l]]), tax_table(top_samp_rel[[l]])[, n], sum),TRUE)[seq_len(topN+1)])/(nsamples(samp_rel[[l]])))
      topN_samp[[l]] <- assign(top_samp_topN[l],prune_taxa(names(Compartments_TopOTUs[[l]]),top_samp_rel[[l]]))
      color_var <- sort(names ((sort(tapply(taxa_sums(top_samp_rel[[l]]), tax_table(top_samp_rel[[l]])[, n], sum),TRUE)[seq_len(topN+2)])/nsamples(grouped_data)))
    } else if (any(names((sort(tapply(taxa_sums(top_samp_rel[[l]]), tax_table(top_samp_rel[[l]])[, n], sum),TRUE)[seq_len(topN+1)])/nsamples(grouped_data)) == "Other", na.rm=T) == FALSE)
    {
      Compartments_TopOTUs[[l]] <- (sort(taxa_sums(top_samp_rel[[l]]),TRUE)[seq_len(topN)])/(nsamples(samp_rel[[l]]))
      print ((sort(tapply(taxa_sums(top_samp_rel[[l]]), tax_table(top_samp_rel[[l]])[, n], sum),TRUE)[seq_len(topN)])/(nsamples(samp_rel[[l]])))
      topN_samp[[l]] <- assign(top_samp_topN[l],prune_taxa(names(Compartments_TopOTUs[[l]]),top_samp_rel[[l]]))
      color_var <- sort(names ((sort(tapply(taxa_sums(top_samp_rel[[l]]), tax_table(top_samp_rel[[l]])[, n], sum),TRUE)[seq_len(topN)])/nsamples(grouped_data)))
      print (top_samp_rel[[l]])
    }
  }
  for (o in (1:length(topN_samp))){
    if (n == "Kingdom"){title <- paste0("Top ", topN," ","Kingdoms ",levels(factor(sample_data(grouped_data)$Compartment))[o], " Samples")} else if (n == "Phylum"){title <- paste0("Top ", topN," ","Phyla ",levels(factor(sample_data(grouped_data)$Compartment))[o], " Samples")}
    else if (n == "Class"){title <- paste0("Top ", topN," ","Classes ",levels(factor(sample_data(grouped_data)$Compartment))[o], " Samples")} else if (n == "Order"){title <- paste0("Top ", topN," ","Orders ",levels(factor(sample_data(grouped_data)$Compartment))[o], " Samples")}
    else if (n == "Family"){title <- paste0("Top ", topN," ","Families ",levels(factor(sample_data(grouped_data)$Compartment))[o], " Samples")} else if (n == "Genus"){title <- paste0("Top ", topN," ","Genera ",levels(factor(sample_data(grouped_data)$Compartment))[o], " Samples")} 
    else {title <- paste0("Top ", topN," ","Species ",levels(factor(sample_data(grouped_data)$Compartment))[o], " Samples")}
    color_pal<- gg_color_hue(length(color_var))
    names(color_pal)<-c(color_var)
    color_pal["Unclassified"] <- "#999999"; color_pal["RF3"]<- "olivedrab2"; color_pal["Thaumarchaeota"] <- "#B385FF";color_pal["TM7"] <- "gold"; color_pal["Planctomycetes"] <- "darkslateblue"
    color_pal["Aurantimonadaceae"] <- "#F8766D";color_pal["Chitinophagaceae"] <- "#DB8E00";color_pal["Cytophagales"] <- "#AEA200"; color_pal["Burkholderiales"]<- "#64B200";
    color_pal["Cytophagaceae"] <- "#AEA200"; color_pal["Nitrososphaeraceae"] <- "coral" ;color_pal["Paenibacillaceae"] <- "#00BD5C";color_pal["Myxococcales"]<- "gold";
    color_pal["Rickettsiales"]<- "goldenrod4";color_pal["Bacillales"]<- "darkkhaki";color_pal["Enterobacteriales"]<- "darkolivegreen3";color_pal["Pseudomonadales"]<- "cornflowerblue";
    color_pal["Rhizobiaceae"] <- "#00C1A7"; color_pal["Sinobacteraceae"] <- "#00A6FF";color_pal["Sphingomonadaceae"] <- "#B385FF";color_pal["Comamonadaceae"] <- "tomato2";
    color_pal["Bacillaceae"] <- "skyblue1";color_pal["Enterobacteriaceae"] <- "moccasin";color_pal["Planococcaceae"] <- "mistyrose3";color_pal["Pseudomonadaceae"] <- "orange4";
    color_pal["Staphylococcaceae"] <- "red4";color_pal["mitochondria"] <- "seagreen4";  color_pal["Streptomycetaceae"] <- "sienna1";
    color_pal["Caulobacteraceae"] <- "#00BADE";color_pal["Hyphomonadaceae"] <- "#EF67EB";color_pal["Syntrophobacteraceae"] <- "#F8766D";color_pal ["Pseudomonas"]<- "darkkhaki";
    color_pal["Rhodospirillaceae"] <- "gold2"; color_pal["ML615J-28"]<- "darkkhaki"; color_pal["Nitrososphaerales"]<- "coral"; color_pal["[Saprospirales]"]<-"skyblue" ;
    color_pal["Betaproteobacteria"] <- "darkkhaki"; color_pal["Deltaproteobacteria"]<- "coral"; color_pal ["[Saprospirae]"]<-"skyblue";color_pal["Sphingomonas"]<-"firebrick1"
    color_pal["Proteobacteria"] <- "#B385FF"; color_pal["Nitrospirae"] <- "#FF6A98"; color_pal["Firmicutes"] <- "deepskyblue"; color_pal["Chloroflexi"] <- "darkkhaki"; color_pal["Verrucomicrobia"]<- "skyblue";
    color_pal["Other"] <- "lemonchiffon2";
    color_pal["Paenibacillus"] <- "#00BFC4";
    color_pal["Asticcacaulis"] <- "gold2" ;  color_pal["Gemmatimonadetes"] <- "#00A6FF";color_pal["Gammaproteobacteria"]<-"cadetblue"
    color_pal["Bradyrhizobiaceae"]<- "#00B0F6"; color_pal["Syntrophobacteraceae"]<- "gold"; color_pal["Skermanella"]<- "darkslateblue";
    color_pal["Bacilli"]<- "#AEA200"; color_pal["S085"]<- "firebrick3"; color_pal["Flavobacteriia"]<- "darkslateblue"; color_pal["TM7-3"]<- "gold"; color_pal["Tenericutes"]<- "darkslateblue";
    color_pal["Saccharibacillus"]<- "tan3";color_pal["Sporosarcina"]<- "gold";color_pal["Staphylococcus"]<- "turquoise4"; color_pal["Acidobacteria"]<- "gold";
    color_pal["Brucellaceae"]<-"forestgreen";color_pal["Lupinus"]<- "darkolivegreen1";color_pal["Azospirillum"]<- "darkseagreen1";color_pal["Alcaligenaceae"]<- "dodgerblue";color_pal["Burkholderiaceae"]<- "darkorchid1";
    color_pal["Methylobacteriaceae"]<- "darkslateblue";color_pal["Methylobacterium"]<- "goldenrod2";color_pal["Burkholderia"]<- "darkslateblue";color_pal["Moraxellaceae"]<- "lightsalmon";color_pal["Dickeya"]<- "mediumpurple4";
    color_pal["Acinetobacter"]<- "dodgerblue3";color_pal["Streptococcus"]<- "goldenrod2";
    color_pal["Rhizobium"] <- "#9590FF";color_pal["Agrobacterium"] <- "lightslateblue";color_pal["Ochrobactrum"] <- "blueviolet";color_pal["Kaistobacter"] <- "mediumpurple4";color_pal["Steroidobacter"] <- "purple";
    color_pal["Cellvibrio"]<- "darkviolet";color_pal["Sinorhizobium"] <- "mediumorchid2";color_pal["Didymeles"]<- "violet";color_pal ["Pseudomonas"]<- "darkslateblue";color_pal["Balneimonas"] <- "purple4";
    color_pal["Bacillus"] <- "deepskyblue";
    color_pal["Nitrospira"]<- "pink";
    color_pal["Chryseobacterium"]<- "darkgoldenrod";color_pal["Flavisolibacter"] <- "darkgoldenrod3";
    color_pal["Candidatus Nitrososphaera"] <- "mediumseagreen";
    color_pal["Streptomyces"]<- "yellow4";color_pal["Rubrobacter"] <- "gold";
    mdf <- psmelt(topN_samp[[o]])
    order <- as.character(unique(mdf[,n]))
    if (any(order == "Unclassified") == TRUE)
    {
      order<- order[-grep("Unclassified", order)]
    }
    if (any(order == "Other") == TRUE)
    {
      order <- order [-grep ("Other", order)]
      order <- sort(order,F)
      order <- append(order, "Other")
      order <- append(order, "Unclassified")
    } else if (any (order == "Other") == FALSE)
    {
      order <- sort(order,F)
      order <- append(order, "Unclassified")
    }
    mdf$Names <- factor(mdf[,n], levels = c(order))
    addline_format <- function(x,...){
      gsub('\\s','\n',x)
    }
    labels=addline_format(mdf$SampleName)
    color_pal<-color_pal[order]
    p<-ggplot(mdf, aes_string(x = 'DevStage', y = "Abundance")) + theme_bw()+
      theme(axis.text.x = element_text(size = 18,angle = 0, hjust = 0.5, colour = "black"), axis.title.x = element_text(size = 18,hjust=0.5,vjust=-0.5),axis.title.y = element_text(size = 18,vjust=2, hjust=0.5)) + ggtitle(title) +
      geom_bar (aes(fill = factor(mdf[,n],levels=order)), stat = 'identity', position = 'stack') +
      ylab("Relative Abundance") + theme(plot.title=element_text(size=18, hjust=0.5), title=element_text(size=18, hjust= 0.38),axis.text.y=element_text(size=18, colour = "black"), legend.text=element_text(size=18), legend.title=element_text(size=18) + scale_y_continuous(limit = c(0:1)))+ 
      scale_colour_manual(values = color_pal) + scale_fill_manual(values = color_pal, name = n)
    print (p)
  }
  i=1;j=1 ;k=1; b=1; c=1; h=1; d=1; q=1;
}

######################## Venn Diagrams ##################
samp <- list()
for (n in unique(sample_data(real_data)$Compartment)){
  samp[[n]] <- assign(n,subset_samples(real_data,Compartment==n))
}
Tax_level<- rank_names(real_data)[c(tax_level)]
taxa_Compartments <- paste0(unique(sample_data(real_data)$Compartment),"_taxa")
taxa_Compartments_pruned <- paste0(taxa_Compartments,"_pruned")
taxa_Compartments_names <- paste0(taxa_Compartments_pruned,"_names")
Compartments_venn<-paste0(unique(sample_data(real_data)$Compartment),"_venn")
Compartments_taxa<- list()
pruned_Compartments_taxa<- list()
names_Compartments_taxa<- list()
venn_list <- list()
venn_Compartments<-list()
universe <- list()
i=1; j=1; k=1; a=1; b=1; c=1
for (n in Tax_level){
  if (n == "Kingdom"){t <- 1} else if (n == "Phylum"){t <- 2} else if (n == "Class"){t <- 3} else if (n == "Order"){t <- 4}
  else if (n == "Family"){t <- 5} else if (n == "Genus"){t <- 6} else {t <- 7}
  for (m in taxa_Compartments){
    Compartments_taxa[[m]] <- assign(m,tax_glom(samp[[i]],rank_names(samp[[i]])[t]))
    i=i+1
  }
  for (l in taxa_Compartments_pruned){
    pruned_Compartments_taxa[[l]] <- assign(l,prune_taxa(taxa_sums(Compartments_taxa[[j]]) > Count_threshold, Compartments_taxa[[j]]))
    j=j+1
  }
  for (o in taxa_Compartments_names){
    names_Compartments_taxa[[o]] <- assign(o,names(sort(tapply(taxa_sums(pruned_Compartments_taxa[[k]]), tax_table(pruned_Compartments_taxa[[k]])[, n], sum),TRUE)))
    k=k+1
  }
  for (p in names_Compartments_taxa){
    venn_list[[unique(sample_data(real_data)$Compartment)[a]]]<-(p)
    if (a == 1){
      universe <- names_Compartments_taxa[[a]]
    }
    else{
      universe <- append(universe, names_Compartments_taxa[[a]])
    }
    a=a+1
  }
  universe <- unique(universe)
  for (r in Compartments_venn){
    venn_Compartments[[r]]<-assign(r, universe %in% names_Compartments_taxa[[b]])
    b=b+1
  }
  if (n == "Kingdom"){title<- "Kingdoms All Samples"} else if (n == "Phylum"){title<- "Phyla All Samples"}
  else if (n == "Class"){title<- "Classses All Samples"} else if (n == "Order"){title<- "Orders All Samples"}
  else if (n == "Family"){title<- "Families All Samples"} else if (n == "Genus"){title<- "Genera All Samples"} 
  else{title<- "Species All Samples"}
  p <- venn.diagram (venn_list, height = 100, width = 100, filename=NULL,cex = 1.5 ,cat.cex = 1.3, fill=rainbow(2),main.cex = 1.5,main = title, category.names = c("Endosphere", "Rhizosphere"))  
  grid.newpage()
  grid.draw(p)
  Endosphere_and_Rhizosphere_shared <- universe [Rhizosphere_venn & Endosphere_venn]
  Only_Rhizosphere <- universe  [Rhizosphere_venn & !Endosphere_venn]
  Only_Endosphere <- universe  [!Rhizosphere_venn & Endosphere_venn]
  i=1;j=1;k=1;a=1;b=1;c=1
  max.len= max(length( Endosphere_and_Rhizosphere_shared ), length(Only_Rhizosphere), length(Only_Endosphere))
  Endosphere_and_Rhizosphere_shared = c(Endosphere_and_Rhizosphere_shared, rep(NA, max.len - length(Endosphere_and_Rhizosphere_shared)))
  Only_Rhizosphere = c(Only_Rhizosphere, rep(NA, max.len - length(Only_Rhizosphere)))
  Only_Endosphere = c(Only_Endosphere, rep(NA, max.len - length(Only_Endosphere)))
#  matrix <- list( All_shared=  All_shared, Endosphere_and_soil_shared= Endosphere_and_Rhizosphere_shared , Only_Rhizosphere= Only_Rhizosphere, Only_Endosphere=Only_Endosphere)
#  matrix_filename <- paste0("matrix_",n,".xlsx")
#  write.xlsx(matrix, matrix_filename)
}

####################### PERMANOVA  using adonis and  ANOSIM (vegan package) ##################
set.seed (42)
BC <- phyloseq::distance(real_data, distance)
map <- data.frame(sample_data(real_data))
adonis(formula=BC~Compartment,data=map, permutations = 999, method = distance)
BC.ano<-anosim(BC, data.frame(map)$Compartment, permutations = 999, distance , strata = NULL, parallel = getOption("mc.cores"))
BC.ano

BC <- phyloseq::distance (subset_samples(real_data,Compartment != c("Endosphere")),distance)
map <- data.frame(sample_data(subset_samples(real_data,Compartment != c("Endosphere"))))
adonis(formula=BC~DevStage,data=map,permutations = 999, method = distance)
BC.ano<-anosim(BC, data.frame(map)$DevStage, permutations = 999, distance , strata = NULL, parallel = getOption("mc.cores"))
BC.ano

BC <- phyloseq::distance (subset_samples(real_data,Compartment != c("Rhizosphere")),distance)
map <- data.frame(sample_data(subset_samples(real_data,Compartment != c("Rhizosphere"))))
adonis(formula=BC~DevStage,data=map,permutations = 999, method = distance)
BC.ano<-anosim(BC, data.frame(map)$DevStage, permutations = 999, distance , strata = NULL, parallel = getOption("mc.cores"))
BC.ano

############## Differential Abundance Analysis Developmental Stage #####################
real_data<-merge_phyloseq(map,otus)
otutable <- data.frame(otu_table(subset_samples(real_data,Compartment == ("Rhizosphere"))))
head(otutable)
otutable <- otutable[rowSums(cpm(otutable) > 0) >= 2,]
dim(otutable)
head(otutable)
grp<- as.factor(paste(sample_data(subset_samples(real_data,Compartment == ("Rhizosphere")))$Compartment,sample_data(subset_samples(real_data,Compartment == ("Rhizosphere")))$DevStage,sep="."))
grp
dge <- DGEList(counts = otutable, group = grp)
dge <- calcNormFactors(dge)
dge <- estimateCommonDisp(dge)
dge
de12 <- exactTest(dge, pair=c(1,2) , dispersion = dge$common.dispersion)
de13 <- exactTest(dge, pair=c(1,3) , dispersion = dge$common.dispersion)
de23 <- exactTest(dge, pair=c(2,3) , dispersion = dge$common.dispersion)
de <-list(de12, de13, de23)
for (a in de) {
  b <- a$comparison
  de_tab <- topTags(a, n=Inf)$table
  de_tab <- signif(de_tab, digits = 3)
  de_tab$FDR <- p.adjust(de_tab$PValue, method="BH") 
  head(de_tab)
  deOTUs <- rownames(de_tab)[de_tab$FDR < FDR_thres & abs(de_tab$logFC) > logFC_thres]
  pos_de_tab <- de_tab[(de_tab$FDR < FDR_thres & abs(de_tab$logFC) >= logFC_thres & de_tab$logFC > 0),]
  if (nrow(pos_de_tab) != 0){
    v<-rownames(tax_table(prune_taxa(rownames(pos_de_tab), real_data)))
    pos_de_tab<-as.data.frame(pos_de_tab)[v,]
    pos_de_tab <- cbind (pos_de_tab,tax_table(prune_taxa(rownames(pos_de_tab), real_data)), stringsAsFactors = FALSE)
    pos_de_tab <- cbind (pos_de_tab,rowSums(otu_table(prune_taxa(rownames(pos_de_tab), subset_samples(subset_samples(real_data,Compartment == c("Rhizosphere")),DevStage == c(substr(b[1],start=13,stop=nchar(b[1]))))))), stringsAsFactors = FALSE)
    pos_de_tab <- cbind (pos_de_tab,rowSums(otu_table(prune_taxa(rownames(pos_de_tab), subset_samples(subset_samples(real_data,Compartment == c("Rhizosphere")),DevStage == c(substr(b[2],start=13,stop=nchar(b[2]))))))), stringsAsFactors = FALSE)
    pos_de_tab <- cbind (pos_de_tab,rowSums(otu_table(prune_taxa(rownames(pos_de_tab), subset_samples(subset_samples(real_data,Compartment == c("Rhizosphere")),DevStage == c(substr(b[1],start=13,stop=nchar(b[1])))))))/sum(otu_table((subset_samples(subset_samples(real_data,Compartment == c("Rhizosphere")),DevStage == c(substr(b[1],start=13,stop=nchar(b[1]))))))), stringsAsFactors = FALSE)
    pos_de_tab <- cbind (pos_de_tab,rowSums(otu_table(prune_taxa(rownames(pos_de_tab), subset_samples(subset_samples(real_data,Compartment == c("Rhizosphere")),DevStage == c(substr(b[2],start=13,stop=nchar(b[2])))))))/sum(otu_table((subset_samples(subset_samples(real_data,Compartment == c("Rhizosphere")),DevStage == c(substr(b[2],start=13,stop=nchar(b[2]))))))), stringsAsFactors = FALSE)
    colnames(pos_de_tab)[12]<-paste0("Abundance_Rhiz_",c(substr(b[1],start=13,stop=nchar(b[1]))))
    colnames(pos_de_tab)[13]<-paste0("Abundance_Rhiz_",c(substr(b[2],start=13,stop=nchar(b[2]))))
    colnames(pos_de_tab)[14]<-paste0("Abundance_Rhiz_","Pod-filling stage")
    colnames(pos_de_tab)[15]<-paste0("Abundance_Rhiz_rel_",c(substr(b[1],start=13,stop=nchar(b[1]))))
    pos_de_tab<- pos_de_tab[with(pos_de_tab, order(-pos_de_tab$logFC)), ]
    print (pos_de_tab)} else {pos_de_tab<- "None"
    print(pos_de_tab)}
  neg_de_tab <- de_tab[(de_tab$FDR < FDR_thres & abs(de_tab$logFC) >= logFC_thres & de_tab$logFC < 0),]
  if (nrow(neg_de_tab) != 0){
    v<-rownames(tax_table(prune_taxa(rownames(neg_de_tab), real_data)))
    neg_de_tab<-as.data.frame(neg_de_tab)[v,]
    neg_de_tab <- cbind (neg_de_tab,tax_table(prune_taxa(rownames(neg_de_tab), real_data)), stringsAsFactors = FALSE)
    neg_de_tab <- cbind (neg_de_tab,rowSums(otu_table(prune_taxa(rownames(neg_de_tab), subset_samples(subset_samples(real_data,Compartment == c("Rhizosphere")),DevStage == c(substr(b[1],start=13,stop=nchar(b[1]))))))), stringsAsFactors = FALSE)
    neg_de_tab <- cbind (neg_de_tab,rowSums(otu_table(prune_taxa(rownames(neg_de_tab), subset_samples(subset_samples(real_data,Compartment == c("Rhizosphere")),DevStage == c(substr(b[2],start=13,stop=nchar(b[2]))))))), stringsAsFactors = FALSE)
    neg_de_tab <- cbind (neg_de_tab,rowSums(otu_table(prune_taxa(rownames(neg_de_tab), subset_samples(subset_samples(real_data,Compartment == c("Rhizosphere")),DevStage == c(substr(b[1],start=13,stop=nchar(b[1])))))))/sum(otu_table((subset_samples(subset_samples(real_data,Compartment == c("Rhizosphere")),DevStage == c(substr(b[1],start=13,stop=nchar(b[1]))))))), stringsAsFactors = FALSE)
    neg_de_tab <- cbind (neg_de_tab,rowSums(otu_table(prune_taxa(rownames(neg_de_tab), subset_samples(subset_samples(real_data,Compartment == c("Rhizosphere")),DevStage == c(substr(b[2],start=13,stop=nchar(b[2])))))))/sum(otu_table((subset_samples(subset_samples(real_data,Compartment == c("Rhizosphere")),DevStage == c(substr(b[2],start=13,stop=nchar(b[2]))))))), stringsAsFactors = FALSE)
    colnames(neg_de_tab)[12]<-paste0("Abundance_Rhiz_",c(substr(b[1],start=13,stop=nchar(b[1]))))
    colnames(neg_de_tab)[13]<-paste0("Abundance_Rhiz_",c(substr(b[2],start=13,stop=nchar(b[2]))))
    colnames(neg_de_tab)[14]<-paste0("Abundance_Rhiz_","Pod-filling stage")
    colnames(neg_de_tab)[15]<-paste0("Abundance_Rhiz_rel_",c(substr(b[1],start=13,stop=nchar(b[1]))))
    neg_de_tab<- neg_de_tab[with(neg_de_tab, order(neg_de_tab$logFC)), ]
    print (neg_de_tab)} else {neg_de_tab<- "None"
    print(neg_de_tab)}
  #  matrix <- list(pos_de_tab)
  #  matrix_filename <- paste0("Rhiz_",c(substr(b[1],start=13,stop=nchar(b[1]))),"_vs_","Rhiz_",c(substr(b[2],start=13,stop=nchar(b[2]))),".xlsx")
  #  write.xlsx(matrix, matrix_filename, sheetName = paste0("Rhiz_",c(substr(b[1],start=13,stop=nchar(b[1]))),"_enriched"))
  #  matrix<- list(neg_de_tab)
  #  write.xlsx (matrix, matrix_filename, sheetName = paste0("Rhiz_",c(substr(b[2],start=13,stop=nchar(b[2]))),"_enriched"), append = T)
  plotSmear(dge, pair=c(b[[1]],b[[2]]), de.tags=deOTUs, cex = 0.65, ylab= "logFC")
  title <- paste0("LogFC:Rhiz", c(substr(b[1],start=13,stop=nchar(b[1]))), "-Rhiz.", c(substr(b[2],start=13,stop=nchar(b[2]))))
  plotSmear(a, de.tags=deOTUs, ylab = title, cex = 0.65)
}

otutable <- data.frame(otu_table(subset_samples(real_data,Compartment == ("Endosphere"))))
head(otutable)
otutable <- otutable[rowSums(cpm(otutable) > 0) >= 2,]
dim(otutable)
head(otutable)
grp<- as.factor(paste(sample_data(subset_samples(real_data,Compartment == ("Endosphere")))$Compartment,sample_data(subset_samples(real_data,Compartment == ("Endosphere")))$DevStage,sep="."))
grp
dge <- DGEList(counts = otutable, group = grp)
dge <- calcNormFactors(dge)
dge <- estimateCommonDisp(dge)
dge
de12 <- exactTest(dge, pair=c(1,2) , dispersion = dge$common.dispersion)
de13 <- exactTest(dge, pair=c(1,3) , dispersion = dge$common.dispersion)
de23 <- exactTest(dge, pair=c(2,3) , dispersion = dge$common.dispersion)
de <-list(de12, de13, de23)
for (a in de) {
  b <- a$comparison
  de_tab <- topTags(a, n=Inf)$table
  de_tab <- signif(de_tab, digits = 3)
  de_tab$FDR <- p.adjust(de_tab$PValue, method="BH") 
  head(de_tab)
  deOTUs <- rownames(de_tab)[de_tab$FDR < FDR_thres & abs(de_tab$logFC) > logFC_thres]
  pos_de_tab <- de_tab[(de_tab$FDR < FDR_thres & abs(de_tab$logFC) >= logFC_thres & de_tab$logFC > 0),]
  if (nrow(pos_de_tab) != 0){
    v<-rownames(tax_table(prune_taxa(rownames(pos_de_tab), real_data)))
    pos_de_tab<-as.data.frame(pos_de_tab)[v,]
    pos_de_tab <- cbind (pos_de_tab,tax_table(prune_taxa(rownames(pos_de_tab), real_data)), stringsAsFactors = FALSE)
    pos_de_tab <- cbind (pos_de_tab,rowSums(otu_table(prune_taxa(rownames(pos_de_tab), subset_samples(subset_samples(real_data,Compartment == c("Endosphere")),DevStage == c(substr(b[1],start=12,stop=nchar(b[1]))))))), stringsAsFactors = FALSE)
    pos_de_tab <- cbind (pos_de_tab,rowSums(otu_table(prune_taxa(rownames(pos_de_tab), subset_samples(subset_samples(real_data,Compartment == c("Endosphere")),DevStage == c(substr(b[2],start=12,stop=nchar(b[2]))))))), stringsAsFactors = FALSE)
    pos_de_tab <- cbind (pos_de_tab,rowSums(otu_table(prune_taxa(rownames(pos_de_tab), subset_samples(subset_samples(real_data,Compartment == c("Endosphere")),DevStage == c(substr(b[1],start=12,stop=nchar(b[1])))))))/sum(otu_table((subset_samples(subset_samples(real_data,Compartment == c("Endosphere")),DevStage == c(substr(b[1],start=12,stop=nchar(b[1]))))))), stringsAsFactors = FALSE)
    pos_de_tab <- cbind (pos_de_tab,rowSums(otu_table(prune_taxa(rownames(pos_de_tab), subset_samples(subset_samples(real_data,Compartment == c("Endosphere")),DevStage == c(substr(b[2],start=12,stop=nchar(b[2])))))))/sum(otu_table((subset_samples(subset_samples(real_data,Compartment == c("Endosphere")),DevStage == c(substr(b[2],start=12,stop=nchar(b[2]))))))), stringsAsFactors = FALSE)
    colnames(pos_de_tab)[12]<-paste0("Abundance_Endo_",c(substr(b[1],start=12,stop=nchar(b[1]))))
    colnames(pos_de_tab)[13]<-paste0("Abundance_Endo_",c(substr(b[2],start=12,stop=nchar(b[2]))))
    colnames(pos_de_tab)[14]<-paste0("Abundance_Endo_rel_",c(substr(b[1],start=12,stop=nchar(b[1]))))
    colnames(pos_de_tab)[15]<-paste0("Abundance_Endo_rel_",c(substr(b[2],start=12,stop=nchar(b[2]))))
    pos_de_tab<- pos_de_tab[with(pos_de_tab, order(-pos_de_tab$logFC)), ]
    print (pos_de_tab)} else {pos_de_tab<- "None"
    print(pos_de_tab)}
  neg_de_tab <- de_tab[(de_tab$FDR < FDR_thres & abs(de_tab$logFC) >= logFC_thres & de_tab$logFC < 0),]
  if (nrow(neg_de_tab) != 0){
    v<-rownames(tax_table(prune_taxa(rownames(neg_de_tab), real_data)))
    neg_de_tab<-as.data.frame(neg_de_tab)[v,]
    neg_de_tab <- cbind (neg_de_tab,tax_table(prune_taxa(rownames(neg_de_tab), real_data)), stringsAsFactors = FALSE)
    neg_de_tab <- cbind (neg_de_tab,rowSums(otu_table(prune_taxa(rownames(neg_de_tab), subset_samples(subset_samples(real_data,Compartment == c("Endosphere")),DevStage == c(substr(b[1],start=12,stop=nchar(b[1]))))))), stringsAsFactors = FALSE)
    neg_de_tab <- cbind (neg_de_tab,rowSums(otu_table(prune_taxa(rownames(neg_de_tab), subset_samples(subset_samples(real_data,Compartment == c("Endosphere")),DevStage == c(substr(b[2],start=12,stop=nchar(b[2]))))))), stringsAsFactors = FALSE)
    neg_de_tab <- cbind (neg_de_tab,rowSums(otu_table(prune_taxa(rownames(neg_de_tab), subset_samples(subset_samples(real_data,Compartment == c("Endosphere")),DevStage == c(substr(b[1],start=12,stop=nchar(b[1])))))))/sum(otu_table((subset_samples(subset_samples(real_data,Compartment == c("Endosphere")),DevStage == c(substr(b[1],start=12,stop=nchar(b[1]))))))), stringsAsFactors = FALSE)
    neg_de_tab <- cbind (neg_de_tab,rowSums(otu_table(prune_taxa(rownames(neg_de_tab), subset_samples(subset_samples(real_data,Compartment == c("Endosphere")),DevStage == c(substr(b[2],start=12,stop=nchar(b[2])))))))/sum(otu_table((subset_samples(subset_samples(real_data,Compartment == c("Endosphere")),DevStage == c(substr(b[2],start=12,stop=nchar(b[2]))))))), stringsAsFactors = FALSE)
    colnames(neg_de_tab)[12]<-paste0("Abundance_Endo_",c(substr(b[1],start=12,stop=nchar(b[1]))))
    colnames(neg_de_tab)[13]<-paste0("Abundance_Endo_",c(substr(b[2],start=12,stop=nchar(b[2]))))
    colnames(neg_de_tab)[14]<-paste0("Abundance_Endo_rel_",c(substr(b[1],start=12,stop=nchar(b[1]))))
    colnames(neg_de_tab)[15]<-paste0("Abundance_Endo_rel_",c(substr(b[2],start=12,stop=nchar(b[2]))))
    neg_de_tab<- neg_de_tab[with(neg_de_tab, order(neg_de_tab$logFC)), ]
    print (neg_de_tab)} else {neg_de_tab<- "None"
    print(neg_de_tab)}
  #  matrix <- list(pos_de_tab)
  #  matrix_filename <- paste0("Endo_",c(substr(b[1],start=12,stop=nchar(b[1]))),"_vs_","Endo_",c(substr(b[2],start=12,stop=nchar(b[2]))),".xlsx")
  #  write.xlsx (matrix, matrix_filename, sheetName = paste0("Endo_",c(substr(b[1],start=12,stop=nchar(b[1]))),"_enriched"))
  #  matrix<- list(neg_de_tab)
  #  write.xlsx (matrix, matrix_filename, sheetName = paste0("Endo_",c(substr(b[2],start=12,stop=nchar(b[2]))),"_enriched"), append = T)
  plotSmear(dge, pair=c(b[[1]],b[[2]]), de.tags=deOTUs, cex = 0.65, ylab= "logFC")
  title <- paste0("LogFC:Endo.", c(substr(b[1],start=12,stop=nchar(b[1]))), "-Endo.", c(substr(b[2],start=12,stop=nchar(b[2]))))
  plotSmear(a, de.tags=deOTUs, ylab = title, cex = 0.65)
}

######################### Rarefaction curves ######################
sample_sums(real_data)
set.seed(42)

calculate_rarefaction_curves <- function(real_data, measures, depths) {
  require('plyr') # ldply
  require('reshape2') # melt
  
  estimate_rarified_richness <- function(real_data, measures, depth) {
    if(max(sample_sums(real_data)) < depth) return()
    real_data <- prune_samples(sample_sums(real_data) >= depth, real_data)
    rarified_real_data <- rarefy_even_depth(real_data, depth, verbose = FALSE)
    alpha_diversity <- estimate_richness(rarified_real_data, measures = measures)
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, real_data = real_data, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  #original  rarefaction_curve_data$Depth <- as.numeric(unique(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  rarefaction_curve_data
}

rarefaction_curve_data <- calculate_rarefaction_curves(real_data,c('Observed'),rep(c(1000,1:110*10000), each=3))
summary(rarefaction_curve_data)
rarefaction_curve_data_summary <- ddply(rarefaction_curve_data, c('Depth', 'Sample', 'Measure'), summarise, Alpha_diversity_mean = mean(Alpha_diversity), Alpha_diversity_sd = sd(Alpha_diversity))
head(rarefaction_curve_data_summary$Sample)
head(rarefaction_curve_data_summary)
sample_data(real_data)$X.SampleID
rarefaction_curve_data_summary$Sample<-sub("\\.","-",data.frame(rarefaction_curve_data_summary)$Sample)
rarefaction_curve_data_summary_verbose <- merge(rarefaction_curve_data_summary, data.frame(sample_data(real_data)), by.x= "Sample", by.y = "X.SampleID",all=T)
colnames(rarefaction_curve_data_summary_verbose)[4]<-"Richness"

head(rarefaction_curve_data_summary_verbose) 

ggplot(data = rarefaction_curve_data_summary_verbose,mapping = aes(x = Depth,y = Richness, ymin = Richness - Alpha_diversity_sd, ymax = Richness + Alpha_diversity_sd, colour = Compartment, group = Sample)) + theme_bw() +
  geom_line() + theme(legend.text = element_text(size=18),axis.title.x=element_text(size=18),axis.text.x=element_text(size=18, colour = "black"),axis.title.y=element_text(size=18),axis.text.y=element_text(size=18, colour = "black"),plot.title = element_text(hjust=0.475), title = element_text(size=18), legend.title = element_text(size=18)) +
  ggtitle ("Rarefaction Curve") + geom_pointrange() + facet_wrap(facets = ~ Measure, scales = 'free_y')+ scale_color_manual(values=c("#ED8141","#69b3a2")) + theme(strip.background = element_blank(),strip.text= element_blank())+ guides(colour = guide_legend(reverse=TRUE))
