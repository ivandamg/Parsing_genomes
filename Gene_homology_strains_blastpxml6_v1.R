########################################################################################
###############
############### GENE homology between strains using blast .xml (output6)
########## 12 oct 2017

########################################################################################

# Libraries
#install.packages('genoPlotR')
library('genoPlotR')
library(RColorBrewer)
library(pheatmap)
#######################################
# Set working directory

setwd('~/Documents/Noemie/WGC_N16961_A1552_Sa5Y_v2/Comp_A1552_N16961_Sa5Y/NonSharedGenes/db_prot_genomes/')

###############################
# Import xml files
# xml best blast hit per gene.

filesToProcess <- dir(pattern = "*\\.xml$")  #files to pro# if event 3 merged
listOfFiles <- lapply(filesToProcess, function(x) tryCatch(read.table(x, header = F),
                                                           error= function (e) cbind.data.frame(V1="NA",V2="NA",V3="NA",
                                                                                                V4="NA",V5="NA",V6="NA",
                                                                                                V7=0,V8=0,V9="NA",
                                                                                                V10="NA",V11="NA",V12="NA")))

names(listOfFiles)<-gsub(".xml","",gsub("blast","",filesToProcess))

listOfFiles<-listOfFiles[grep("Rec",names(listOfFiles),invert=T)]

# Create column header

colnam<-strsplit(names(listOfFiles), "_" )
colnam<-lapply(colnam, function (x) x[c(2,4)])
colnam<-lapply(colnam, function (x) gsub(".xml","",x))
   
for (i in 1:length(filesToProcess)) {
  NAME<-c(unlist(colnam[[i]]),"Identity","Length","Mnn","Humm","Start_strain1","End_strain1","Start_strain2","End_strain2","Evalue","Bitscore")
  
  colnames(listOfFiles[[i]]) <- NAME }

listOfFiles[[1]]
names(listOfFiles)

###############################################################################################3
########### FILTERING THE DATA
# filter 100% identity
listOfFiles<-lapply(listOfFiles, function(x) x[x$Identity>79.9,] )

# length no higher than 20% 

listOfFiles<-lapply(listOfFiles, function(x) x[(x$End_strain1- x$Start_strain1) / (x$End_strain2- x$Start_strain2) >0.8 ,] )

listOfFiles<-lapply(listOfFiles, function(x) x[(x$End_strain1- x$Start_strain1) / (x$End_strain2- x$Start_strain2) <1.2 ,] )

# Reciprocal blast hits only



###############################################################################################

# not reciprocal hits
Comp1<-listOfFiles




merged<-merge(Comp1[[1]][,1:2],Comp1[[2]][,1:2],by= "N16961",all = T)
#merged<-merge(merged,Comp1[[3]][,1:2],by= "A1552",all = T)

merged[merged$N16961=="VC0515",]



###
# Merge the homologs of N16961 (A1552 and Sa5y ) To all the genes in N16961.

# import protein names of N16961

labN16961<-read.table("../VC-N16961_ASM674v1_protein.labels",h=F)

labN16961$V1<-gsub(">","",labN16961[,1])
 colnames(labN16961)<-"N16961"
head(labN16961)

labN16961[labN16961=="VC0512"]

merged<-merge(merged,labN16961,by="N16961",all=T)


merged$N16961<-as.character(merged$N16961)
merged$A1552<-as.character(merged$A1552)
merged$Sa5Y<-as.character(merged$Sa5Y)


merged<-merged[order(merged$N16961),]


Colorcitos<-c(rep("gray",dim(merged[grep("VCA",merged$N16961),])[1]),
rep("black" ,length(merged$N16961) - dim(merged[grep("VCA",merged$N16961),])[1])     )



# info 
length(unique(merged$N16961))
sum(is.na(merged$A1552))
sum(is.na(merged$Sa5Y))


write.table(merged,"~/Documents/Noemie/WGC_N16961_A1552_Sa5Y_v2/Comp_A1552_N16961_Sa5Y/NonSharedGenes/Gene_orthology_N16961-A1552-Sa5Y.txt",sep = "\t",quote = F,row.names = F)

# info about merged data
tail(merged,100)
merged[300:600,]
merged[merged$N16961=="VC0515",]
table(merged$N16961)[table(merged$N16961)==2]


# plot 
merged[!is.na(merged)] <- 1
merged[is.na(merged)] <- 0


legenda<-rownames(merged)

legenda[merged$A1552!=0 | merged$Sa5Y!=0]<-""

merged$N16961<-as.numeric(merged$N16961)
merged$A1552<-as.numeric(merged$A1552)
merged$Sa5Y<-as.numeric(merged$Sa5Y)
dim(merged)
length(legenda)
head(merged)

# genes not present in other strains.
rownames(merged[merged$Sa5Y==0,])
rownames(merged[merged$A1552==0,])

rownames(merged[merged$N16961==0,])


pdf("~/Documents/Noemie/WGC_N16961_A1552_Sa5Y_v2/Comp_A1552_N16961_Sa5Y/NonSharedGenes/Gene_orthology_N16961-A1552-Sa5Y.pdf", height=2, width = 8)

pheatmap(t(merged),cluster_rows = F,cluster_cols = F,color = c("white", "black"),cellwidth =0.1 ,cellheight = 30,fontsize_row = 10,legend = F, show_colnames = F)

dev.off()


####################################################
####################################################
####################################################
####################################################














