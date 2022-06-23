library(data.table)
library(ggplot2)
library(ggpubr)
library(stringr)

BathMet = fread("/Users/smedina/Desktop/Data/MethylationData/BathMet2Merged.txt")
Bismark = fread("/Users/smedina/Desktop/Data/MethylationData/BismarkMerged.txt")
DS = fread("/Users/smedina/Desktop/Data/MethylationData/DeepSignalMerged.txt")
liste=list(Bismark,DS,BathMet)
names(liste)=c( "Bismark", "DeepSignalPlant", "BathMet2")
tab = rbindlist(listeT, idcol = 'Type', fill=TRUE, use.names=TRUE)

DSTranspos = fread(DSMerged)
BismarkTranspos = fread(BmMerged)


plotCorreDSBmparCont <- function(filelist) {
  tabT <- rbindlist(filelist, idcol = 'Type', fill=TRUE, use.names=TRUE)
  missing_contexte = nrow(tabT[Contexte == ""])
  warning("Missing contexte in ", missing_contexte)
  tabT[Contexte == "", Contexte := 'CHH']
  #print (tabT[, .N, by = Contexte])
  print(tabT[Contexte == "CHG" & is.na(Pourcentage_de_methylation), .N])
  aT = tabT[, list("%m"=mean(Pourcentage_de_methylation, na.rm  = T)), by = c('Type', 'ID', 'Chr', "Contexte")]
  bT = dcast(aT, Chr+ID+Contexte ~ Type)
  bT[, Annotation := "TE"]
  bT[str_detect(ID, "G"), Annotation := "Gene"]
  table(bT$Contexte)
  table(bT$Annotation)
  
  plotT <- ggplot(bT, aes(x = Bismark, y = DeepSignalPlant )) + geom_point(alpha = 0.2) + 
    geom_smooth(col ='deeppink')+ 
    stat_cor(p.accuracy = 0.001, r.accuracy = 0.01) + facet_grid(Contexte~Annotation)+
    labs(title="Diagramme de corrélation", 
         subtitle="DeepSignalPlant / Bismark par élément et contexte de cytosine")
  return (plotT)
  }

plotCorreDSBmparCont(listeT)
plotCorreDSBmparCont(mylistbd)

DSTranspos = fread("/Users/smedina/Desktop/Stage/DSMergedTranspos.txt")
BismarkTranspos = fread("/Users/smedina/Desktop/Stage/BismarkMergedTranspos.txt")
listeT=list(DSTranspos,BismarkTranspos)
names(listeT)=c("DeepSignalPlant","Bismark")

mylist <- list("DSPlant" = DSTranspos, "Bismark" = BismarkTranspos)

#names(listeT)=c("DeepSignalPlant","Bismark")
#listeT=list(DSMerged[ID %in% geneid],BmMerged[ID %in% geneid])

barPlotparGene <- function(filelist, geneid, debug = F){
  tab <- rbindlist(filelist, idcol = 'Type', fill=TRUE, use.names=TRUE)
  print(dim(tab))
  missing_contexte = nrow(tab[Contexte == ""])
  warning("Missing contexte in ", missing_contexte)
  tab[Contexte == "", Contexte := 'CHH']
  sub <- tab[ID %in% geneid, c("ID","Pourcentage_de_methylation","Contexte","Pos_Locale","Type", "Couverture")]
  subNA <- sub[is.na(Pourcentage_de_methylation) ]
  if (debug){
    print(tab)
    print(sub)
    print(subNA)
    }
  A <- ggplot(sub, aes(x=Pos_Locale, y=Pourcentage_de_methylation)) +
    geom_bar(aes(fill = Contexte), stat="identity",position=position_dodge()) + theme_bw() +
    labs(x = "Position locale", y = "% methylation") + geom_line(aes(y = Couverture), size = 0.1) + 
    geom_point(aes(x = Pos_Locale, y = -10,color = "Cytosine"), size = 0.2) + coord_cartesian(ylim = c(-11, 100))
  # if NA in table
  if (nrow(subNA) > 0) {
    A <- A + geom_col(data = subNA, aes(Pos_Locale, y = 100), width = 0.3, fill = "grey") +
      geom_point(data = subNA, aes(Pos_Locale, y = -5, color = "Pas de couverture"), width = 3,  size = 0.2, alpha = 0.5) +
      scale_color_manual(values = c("Pas de couverture" = "black", "Cytosine"="orange"))
    }
  # if multiple genes
  if(length(geneid) > 1) {
    A <- A + facet_wrap(ID ~ Type, ncol = length(filelist), scales = "free")
  } else {
    A <- A + facet_wrap(~Type, ncol = 1, scales = "free_x") +scale_color_manual(values = c("Cytosine"="orange"))
  } 
  A <- A + labs(color = "Position")
  A
}

barPlotparGene(filelist = mylist, c("AT1G01010", "AT1G38065", "AT5TE75660"))

getProfile <- function(filelist, geneid){
  tab <- rbindlist(filelist, idcol = 'Type', fill=TRUE, use.names=TRUE)
  missing_contexte = nrow(tab[Contexte == ""])
  warning("Missing contexte : ", missing_contexte)
  tab[Contexte == "", Contexte := 'CHH']
  sub <- tab[ID %in% geneid, c("ID","Pourcentage_de_methylation","Contexte","Pos_Locale","Type", "Couverture")]
  p <- ggplot(sub, aes(x=Contexte, y=Pourcentage_de_methylation, color=Contexte)) + theme_bw() +
    labs(title="BoxPlot de la méthylation par contexte", subtitle = unique(sub$ID), y = "% methylation") +
    geom_boxplot(outlier.colour="grey", outlier.shape=8, outlier.size=2) +
    stat_summary(fun=mean, geom="point", shape=13, size=4) + ylim(c(0,100))
  return (p)
}
 
getProfile(mylist,"AT5G64550")

barPlotparGene(mylistbd, "BdiBd21-3.1G0000200")

mylistbd <- list("Bismark" = BdBismarkMerged, "DeepSignalPlant" =bddpmerged)

nrow(subNA) #donne  le nombre de lignes d'une dt
length(subNA) #donne le nombre de colonnes d'une dt


ggplot(tab[ID == "AT1G38065"], aes(Pos_Locale, Pourcentage_de_methylation, fill = Contexte)) + 
  geom_col(position = "dodge", width = 2) + theme_bw() + facet_wrap(~Type, ncol = 1) 

ggplot(tab[ID == "AT1G38065"], aes(Pos_Locale, Pourcentage_de_methylation, color = Contexte)) + 
  geom_point() + theme_bw() + facet_wrap(~Type, ncol = 1) 





Test= fread("/Users/smedina/Desktop/Stage/BismarkMTTest.txt")
DSTranspos[DSTranspos$Contexte==""]
BismarkTranspos[!BismarkTranspos$Contexte_x == Test$Contexte_y]


sub = tab[, list("Total_cytosine"=.N), by = c('Type', 'ID', 'Chr')] #list permet de nommer la colonne
missing = tab[is.na(Couverture) | Couverture == 0, list("Missing_cytosine" = .N), by = c('Type', 'ID')]
missingloc = tab[is.na(Couverture) | Couverture == 0, list("Missing_cytosine" = .N), by = c('Type', 'ID','Chr','Pourcentage_de_methylation')]


c = dcast(missing,ID~Type)
fwrite(c, file="/Users/smedina/Desktop/Data/MissingparType.txt", sep='\t', na =NA)


sub2 = tab[!is.na(Couverture) & Couverture != 0 , list("Total_C_outil"=.N), by = c('Type', 'ID', 'Chr')]
sub2[ID=='AT1G38065']
missing[ID=='AT1G38065']

BdBismarkMerged <- fread("/Volumes/DATA/Bd05_1_val_1_bismark_hisat2_pe.deduplicated.CX_reportMerged.txt")
bddpmerged <- fread("/Volumes/DATA/GitHub/DPtoBmMerged.txt")

bddpmerged[Contexte == "CHG" & !is.na(Pourcentage_de_methylation)]
#a = tab[, mean(Pourcentage_de_methylation, na.rm = T), by = c("Type", "ID")]#[is.na(V1)]

joined_dt <- merge(sub, missing, by = c("Type", 'ID'),all.x = TRUE)
joined_dt[is.na(Missing_cytosine), Missing_cytosine := 0]
joined_dt[, perc := Missing_cytosine / Total_cytosine *100]
#joined_dt[, perc := round(perc, digits = 2)]

tricyto <- setorder(joined_dt, Total_cytosine)
fwrite(tricyto, file="TableTriC.csv")
tripercmissing <- setorder(joined_dt, -perc)
fwrite(tripercmissing, file="TableTriMissing.csv")


# Filter


dcast.data.table(joined_dt, ID+Chr~Type, value.var = "perc")  #change les données en wide
dcast.data.table(joined_dt, ID+Chr~Type, value.var = "Missing_cytosine")  #change les données en wide
dcast.data.table(joined_dt, ID+Chr~Type, value.var = "Total_cytosine")








