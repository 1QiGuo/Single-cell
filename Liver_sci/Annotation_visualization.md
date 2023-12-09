# Heatmap

```{r}
library("readxl")
library(reshape2) 
library("pheatmap")
library("ggsci")
library("RColorBrewer")
library(qs)
library("Seurat")
setwd("/fs/ess/PAS1475/guoqi/Grant/NIDDK/")
marker <- read_excel("./markers.xlsx", sheet = 1)
draw_colnames_45 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 0.75, hjust = 1, rot = 42, gp = grid::gpar(...)
  )
  return(res)
}
assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_45",
  ns = asNamespace("pheatmap")
)
#two genes are not exist - Serpina3k, F2

annotation_heatmap<-function(marker_f,object,clusterorder){
  library("Polychrome")
  marker_f<-melt(marker_f,measure.vars = colnames(marker_f),variable_name = "Celltype",value.name="marker")
  colnames(marker_f)<-c("Celltype","value")
  marker_f<-marker_f[!is.na(marker_f$value),]
  #compute average expression 
  avg_data<-data.frame(rep(0,nrow(marker_f)))
  setdiff(marker_f$value,rownames(liver_integrate))
  DefaultAssay(object)<-"RNA"
  for(i in sort(unique(Idents(object)))){
    object_subset<-subset(object,idents=i)
    df<-AverageExpression(object_subset, assays = "RNA",features = marker_f$value, slot = "data")
    df<-as.data.frame(df$RNA)
    avg_data<-cbind(avg_data,df)
  }
  avg_data<-avg_data[,-1]
  colnames(avg_data)<-sort(unique(Idents(object)))
  #create pheatmap data
  marker_f$Celltype<-factor(marker_f$Celltype,levels = unique(marker_f$Celltype))
  marker = data.frame(marker = marker_f$Celltype)
  color = marker
  levels(color) <- Polychrome::dark.colors(15)
  color <- list(marker = levels(color))
  names(color$marker)<- levels(marker$marker)
  separation_sequence <- cumsum(table(marker_f$Celltype))
  gaps_row = separation_sequence
  p <- pheatmap(avg_data[,clusterorder],
                color = colorRampPalette(c("blue","white","red"))(100),
                cluster_rows = F,
                annotation_row = marker,
                annotation_colors = color,
                cluster_cols = F,
                scale = "row",border_color = "NA",
                gaps_row = separation_sequence,fontsize = 15,
                annotation_names_row = F
  )
  
  return(p)
}

FeaturePlot(liver_integrate,features=c("Alb","Apoc2"))
FeaturePlot(liver_integrate,features=c("Clec4f","Apoc2"))
Idents(liver_integrate)<-liver_integrate$seurat_clusters
order<-c("0","1","5","6","12","13","14","2","16","17","15","9","10","3","18","11","4","7","8")
annotation_heatmap(marker,liver_integrate,clusterorder = order)

```
# Dot plot
```{r}
DotPlot(liver_integrate,features=marker_f$value)+
  theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1),axis.text = element_text(size=17))+
  scale_colour_gradient(low = "#00008b", high = "#ffa500")+
  geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5) 
annotation_heatmap(marker,liver_integrate,clusterorder = order)
```
