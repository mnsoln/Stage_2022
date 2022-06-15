library(data.table)
library(ggplot2)
library(ggpubr)
library(stringr)

BathMet = fread("/Users/smedina/Desktop/Data/MethylationData/BathMet2Merged.txt")
Bismark = fread("/Users/smedina/Desktop/Data/MethylationData/BismarkMerged.txt")
DS = fread("/Users/smedina/Desktop/Data/MethylationData/DeepSignalMerged.txt")
#Bed = fread("/Users/smedina/Desktop/Data/MethylationData/BedMerged.txt")
liste=list(Bismark,DS,BathMet)
names(liste)=c( "Bismark", "DeepSignalPlant", "BathMet2")
tab = rbindlist(liste, idcol = 'Type', fill=TRUE, use.names=TRUE)

DSTranspos = fread(DSMerged)
BismarkTranspos = fread(BmMerged)
listeT=list(DSTranspos,BismarkTranspos)
names(listeT)=c("DeepSignalPlant","Bismark")

plotCorreDSBmparCont <- function(filelist) {
  tabT <- rbindlist(filelist, idcol = 'Type', fill=TRUE, use.names=TRUE)
  missing_contexte = nrow(tabT[Contexte == ""])
  warning("Missing contexte in ", missing_contexte)
  tabT[Contexte == "", Contexte := 'CHH']
  aT = tabT[, list("%m"=mean(Pourcentage_de_methylation, na.rm  = T)), by = c('Type', 'ID', 'Chr', "Contexte")]
  bT = dcast(aT, Chr+ID+Contexte ~ Type)
  bT[, Annotation := "TE"]
  bT[str_detect(ID, "G"), Annotation := "Gene"]
  table(bT$Contexte)
  DSTranspos[Contexte == ""]
  table(bT$Annotation)
  
  plotT <- ggplot(bT, aes(x = Bismark, y = DeepSignalPlant )) + geom_point(alpha = 0.2) + 
    geom_smooth(col ='deeppink')+ 
    stat_cor(p.accuracy = 0.001, r.accuracy = 0.01) + facet_grid(Contexte~Annotation)+
    labs(title="Diagramme de corrélation", 
         subtitle="DeepSignalPlant / Bismark par élément et contexte de cytosine")
  return (plotT)
  }

plotCorreDSBmparCont(listeT)

DSTranspos = fread("/Users/smedina/Desktop/Stage/DSMergedTranspos.txt")
BismarkTranspos = fread("/Users/smedina/Desktop/Stage/BismarkMergedTranspos.txt")

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

barPlotparGene(filelist = mylist, c("AT1G01010", "AT1G38065"))#, "AT5TE75660"))

getProfile <- function(filelist, geneid){
  tab <- rbindlist(filelist, idcol = 'Type', fill=TRUE, use.names=TRUE)
  missing_contexte = nrow(tab[Contexte == ""])
  warning("Missing contexte : ", missing_contexte)
  tab[Contexte == "", Contexte := 'CHH']
  sub <- tab[ID %in% geneid, c("ID","Pourcentage_de_methylation","Contexte","Pos_Locale","Type", "Couverture")]
  p <- ggplot(sub, aes(x=Contexte, y=Pourcentage_de_methylation, color=Contexte)) + theme_bw() +
    labs(title="BoxPlot de la méthylation par contexte", subtitle = unique(sub$ID), y = "% methylation") +
    geom_boxplot(outlier.colour="grey", outlier.shape=8, outlier.size=2) +
    stat_summary(fun=mean, geom="point", shape=13, size=4)
  return (p)
}

getProfile(mylist,"AT1G40090")


nrow(subNA) #donne  le nombre de lignes d'une dt
length(subNA) #donne le nombre de colonnes d'une dt


ggplot(tab[ID == "AT1G38065"], aes(Pos_Locale, Pourcentage_de_methylation, fill = Contexte)) + 
  geom_col(position = "dodge", width = 2) + theme_bw() + facet_wrap(~Type, ncol = 1) 

ggplot(tab[ID == "AT1G38065"], aes(Pos_Locale, Pourcentage_de_methylation, color = Contexte)) + 
  geom_point() + theme_bw() + facet_wrap(~Type, ncol = 1) 






    
p <- ggplot(df3, aes(x=dose, y=len, fill=supp)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=len-sd, ymax=len+sd), width=.2,
                position=position_dodge(.9))



Test= fread("/Users/smedina/Desktop/Stage/BismarkMTTest.txt")
DSTranspos[DSTranspos$Contexte==""]
BismarkTranspos[!BismarkTranspos$Contexte_x == Test$Contexte_y]


Bismarktraite[Chr == "Chr5" & Pos_Globale < 15]
Test[Chr == "Chr5" & Pos_Globale < 15]




b[Contexte == "CG" & DeepSignalPlant >= 50 & Bismark <= 30]
b[ID == 'AT1G38065']



#elements avec l + et le - de c
data <- tab[with(tab,order(-)),]

data <- data[1:10,]








sub = tab[, list("Total_cytosine"=.N), by = c('Type', 'ID', 'Chr','Pourcentage_de_methylation')] #list permet de nommer la colonne
missing = tab[is.na(Couverture) | Couverture == 0, list("Missing_cytosine" = .N), by = c('Type', 'ID')]
missingloc = tab[is.na(Couverture) | Couverture == 0, list("Missing_cytosine" = .N), by = c('Type', 'ID','Chr','Pourcentage_de_methylation')]


ggplot(b, aes(x = Bismark, y = DeepSignalPlant)) + geom_point(alpha = 0.2) + 
  #geom_smooth(method = lm, se = FALSE)+ 
  geom_smooth()+ 
  stat_cor(p.accuracy = 0.001, r.accuracy = 0.01) 


c = dcast(missing,ID~Type)
fwrite(c, file="/Users/smedina/Desktop/Data/MissingparType.txt", sep='\t', na =NA)


sub2 = tab[!is.na(Couverture) & Couverture != 0 , list("Total_C_outil"=.N), by = c('Type', 'ID', 'Chr')]
sub2[ID=='AT1G38065']
missing[ID=='AT1G38065']



#a = tab[, mean(Pourcentage_de_methylation, na.rm = T), by = c("Type", "ID")]#[is.na(V1)]

joined_dt <- merge(sub, missing, by = c("Type", 'ID'),all.x = TRUE)
joined_dt[is.na(Missing_cytosine), Missing_cytosine := 0]
joined_dt[, perc := Missing_cytosine / Total_cytosine *100]
joined_dt[, perc := round(perc, digits = 2)]
joined_dt
setorder(joined_dt, -perc)

# Filter


mat = dcast(joined_dt, ID ~ Type, value.var = "perc", fun.aggregate = mean)



dcast.data.table(joined_dt, ID+Chr~Type, value.var = "perc")  #change les données en wide
dcast.data.table(joined_dt, ID+Chr~Type, value.var = "Missing_cytosine")  #change les données en wide
dcast.data.table(joined_dt, ID+Chr~Type, value.var = "Total_cytosine")

ggplot(bT, aes(Chr, perc, fill = Contexte)) + geom_boxplot() + facet_wrap(~Type)


joined_pos <- merge(sub, missingloc, by = c("Type", 'ID','Pos_Locale','Chr','Pourcentage_de_methylation'),all.x = TRUE)
joined_pos[is.na(Missing_cytosine), Missing_cytosine := 0]
ggplot(joined_pos, aes(Pos_Locale,Missing_cytosine, fill = Type)) + geom_count() + facet_wrap(~Chr)

ggplot(joined_pos, aes(Pos_Locale,as.integer(Pourcentage_de_methylation), fill = Type)) + geom_count() + facet_wrap(~Chr)



joined_dt[ID == "AT1G38065"]
tab[ID == "AT1G38065" & Type == "Bismark"]







