## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=F-----------------------------------------------------------
library(devtools)
library(lionessR)
library(igraph)
library(reshape2)
library(limma)

## ------------------------------------------------------------------------
data(OSdata)

## ------------------------------------------------------------------------
nsel=500
cvar <- apply(as.array(as.matrix(exp)), 1, sd)
dat <- cbind(cvar, exp)
dat <- dat[order(dat[,1], decreasing=T),]
dat <- dat[1:nsel, -1]
dat <- as.matrix(dat)

## ------------------------------------------------------------------------
groupyes <- which(targets$mets=="yes")
groupno <- which(targets$mets=="no")
netyes <- cor(t(dat[,groupyes]))
netno <- cor(t(dat[,groupno]))
netdiff <- netyes-netno

## ------------------------------------------------------------------------
cormat2 <- rep(1:nsel, each=nsel)
cormat1 <- rep(1:nsel,nsel)
el <- cbind(cormat1, cormat2, c(netdiff))
melted <- melt(upper.tri(netdiff))
melted <- melted[which(melted$value),]
values <- netdiff[which(upper.tri(netdiff))]
melted <- cbind(melted[,1:2], values)
genes <- row.names(netdiff)
melted[,1] <- genes[melted[,1]]
melted[,2] <- genes[melted[,2]]
row.names(melted) <- paste(melted[,1], melted[,2], sep="_")
tosub <- melted
tosel <- row.names(tosub[which(abs(tosub[,3])>0.5),])

## ------------------------------------------------------------------------
cormat <- lioness(dat, netFun)
row.names(cormat) <- paste(cormat[,1], cormat[,2], sep="_")
corsub <- cormat[which(row.names(cormat) %in% tosel),3:ncol(cormat)]
corsub <- as.matrix(corsub)

## ------------------------------------------------------------------------
group <- factor(targets[,2])
design <- model.matrix(~0+group)
cont.matrix <- makeContrasts(yesvsno = (groupyes - groupno), levels = design)  
fit <- lmFit(corsub, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2e <- eBayes(fit2)
toptable <- topTable(fit2e, number=nrow(corsub), adjust="fdr")

## ------------------------------------------------------------------------
toptable_edges <- t(matrix(unlist(c(strsplit(row.names(toptable), "_"))),2))
z <- cbind(toptable_edges[1:50,], toptable$logFC[1:50])
g <- graph.data.frame(z, directed=F)
E(g)$weight <- as.numeric(z[,3])
E(g)$color[E(g)$weight<0] <- "blue"
E(g)$color[E(g)$weight>0] <- "red"
E(g)$weight <- 1

## ------------------------------------------------------------------------
topgeneslist <- unique(c(toptable_edges[1:50,]))
fit <- lmFit(exp, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2e <- eBayes(fit2)
topDE <- topTable(fit2e, number=nrow(exp), adjust="fdr")
topDE <- topDE[which(row.names(topDE) %in% topgeneslist),]
topgenesDE <- cbind(row.names(topDE), topDE$t)

## ------------------------------------------------------------------------
# add t-statistic to network nodes
nodeorder <- cbind(V(g)$name, 1:length(V(g)))
nodes <- merge(nodeorder, topgenesDE, by.x=1, by.y=1)
nodes <- nodes[order(as.numeric(as.character(nodes[,2]))),]
nodes[,3] <- as.numeric(as.character(nodes[,3]))
nodes <- nodes[,-2]
V(g)$weight <- nodes[,2]

# make a color palette
mypalette4 <- colorRampPalette(c("blue","white","white","red"), space="Lab")(256) 
breaks2a <- seq(min(V(g)$weight), 0, length.out=128)
breaks2b <- seq(0.00001, max(V(g)$weight)+0.1,length.out=128)
breaks4 <- c(breaks2a,breaks2b)

# select bins for colors
bincol <- rep(NA, length(V(g)))
for(i in 1:length(V(g))){
	bincol[i] <- min(which(breaks4>V(g)$weight[i]))
}
bincol <- mypalette4[bincol]
	
# add colors to nodes
V(g)$color <- bincol

## ---- dpi=70, fig.width=10, fig.height=10--------------------------------
par(mar=c(0,0,0,0))
plot(g, vertex.label.cex=0.7, vertex.size=10, vertex.label.color = "black", vertex.label.font=3, edge.width=10*(abs(as.numeric(z[,3]))-0.7), vertex.color=V(g)$color)

