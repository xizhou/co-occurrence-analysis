

MeSH term co-occurrence analysis
==========
It is a probabilistic model based on the Hypergeometric distribution  for detecting statistically significant co-occurred MeSH items.

## Authors

[周晓北] (Zhou Xiaobei)  
[周淼]  (Zhou Miao)  
[黄德生] (Huang Desheng)    
[崔雷] (Cui Lei)  

## Download data and code
```
 $ git clone https://github.com/xizhou/co-occurrence-analysis

## Install pubMR package
```r
install.packages("devtools")
devtools::install_github("xizhou/pubMR")
```

The detailed methods are introduced in our paper.

## Parametric simulation 
We set up a simulation to reflect the reality of MeSH term co-occurrence data and we simulate different variables. Detailed information are elaborated in our article.
The code is following to realize this process:

```r
dir <- "~/co-occurrence-analysis"
setwd(dir)
source(".code/code.R")
set.seed(100)
nt <- 100
nf <- 200
pf <- runif(nf,min=0.05,max=0.2)
pt1 <- runif(nt,min=0.05,max=0.2)
pt2 <- runif(nt,min=0.05,max=0.2)

rho01 <- rep(0.1,nt)
r200_01 <- simu(rho01,pf,pt1,pt2,nt,nf,200)
r500_01 <- simu(rho01,pf,pt1,pt2,nt,nf,500)
r1000_01 <- simu(rho01,pf,pt1,pt2,nt,nf,1000)
r2000_01 <- simu(rho01,pf,pt1,pt2,nt,nf,2000)

rho02 <- rep(0.2,nt)
r200_02 <- simu(rho02,pf,pt1,pt2,nt,nf,200)
r500_02 <- simu(rho02,pf,pt1,pt2,nt,nf,500)
r1000_02 <- simu(rho02,pf,pt1,pt2,nt,nf,1000)
r2000_02 <- simu(rho02,pf,pt1,pt2,nt,nf,2000)
```
The function "simu()" is in "code.R".

- We calculated type I and II errors to evaluate the performance of our model:

```r
res1 <- rbind(r200_01,r500_01,r1000_01,r2000_01)
res1 <- cbind(c(200,500,1000,2000),res1)
res2 <- rbind(r200_02,r500_02,r1000_02,r2000_02)
res <- cbind(res1,res2)
```
The result is shown in Table 1 in the manuscript.

![Image text](https://raw.githubusercontent.com/xizhou/co-occurrence-analysis/main/table1.png)

## Nonparametric simulation 

Firstly you need download genome.xml file using pubMR:
```r
#dir <- "~/co-occurrence-analysis"
#setwd(dir)
#m <- '"Genome Biol"[Journal]("2001/01/01"[PDAT]:"2020/01/01"[PDAT])'
#library(pubMR)
#library(XML)
#library(data.table)
#obj1 <- txtList(input=m,outputType='xml')
#saveXML(obj1,file=".data/genome.xml")
```

or else you can unzip genome.xml file in data directory, then run the following Rcode: 

```r
dir <- "~/co-occurrence-analysis"
setwd(dir)
source("./code/code.R")
library(pubMR)
library(data.table)
library(tidyr)
obj <- txtList(input="./data/genome.xml",inputType="xml")
obj1=data.table(PMID=obj@PMID,MS=obj@MH)
MS <- obj1[,MS]
idx <- sapply(MS,is.null)
obj1 <- obj1[!idx,]
obj1 = obj1 %>% unnest(MS) %>%as.data.table
obj1[,N:=.N,by=MS]
obj1 <- obj1[N>50,]
v0 <- v <- table(obj1[,c("MS","PMID")])
nc <- ncol(v)
v <- crossprod(t(v))
p=hyp(v,length(obj@PMID))
diag(p) <- 1
p[is.na(p)] <- 0
s <- 1-p

cosine <- function(x)
{
   x <- as.matrix(x)
   sim <- x / sqrt(rowSums(x*x))
   sim <- sim %*% t(sim)
   sim
}

n <- obj1[,length(unique(PMID))]
tt <- (1+sqrt(1+8*n))/2
id <- diag(v)<tt
s1 <- v
s1 <- cosine(s1)
diag(s1) <- 0

n <- 100
v1 <- matrix(0,nrow=n,nc=nc)
for( i in seq(nrow(v1)))
{
   id <-sample(nc,500)
   v1[i,id] <- 1
}
rownames(v1) <- paste0("word_",seq(n))
v1 <- rbind(v0,v1)


v1 <- crossprod(t(v1))
p1=hyp(v1,length(obj@PMID))
diag(p1) <- 1
s2 <- 1-p1

n <- obj1[,length(unique(PMID))]
tt <- (1+sqrt(1+8*n))/2
id <- diag(v1)<tt
s3 <- v1
s3 <- cosine(s3)
diag(s3) <- 0

library(igraph)
png("fig5.png",w=4000,h=4200)
par(mfrow=c(2,2), mar=c(3,5,5,4),oma=c(3,3,3,3))
set.seed(100)
g <- graph.adjacency(s, mode = "undirected", weighted =T, diag = F)
g <- as_data_frame(g)
g <- g[order(g[,3],decreasing=T),]
g <- g[1:100,]
g <- graph_from_data_frame(g)
E(g)$color <- "blue"
V(g)$color <- "black"
V(g)$color[grep("word",V(g)$name)] <- "red"
plot(g,vertex.size=2.5,vertex.label.cex=3)
title("h-test (without noise)",cex.main=7)
mtext("(a)", side=3, padj=0,adj=0,cex=8)
set.seed(100)
g <- graph.adjacency(s1, mode = "undirected", weighted =T, diag = F)
g <- as_data_frame(g)
g <- g[order(g[,3],decreasing=T),]
g <- g[1:100,]
g <- graph_from_data_frame(g)
E(g)$color <- "blue"
V(g)$color <- "black"
V(g)$color[grep("word",V(g)$name)] <- "red"
plot(g,vertex.size=2.5,vertex.label.cex=3,main="d-method")
title("d-method (without noise)", cex.main=7)
mtext("(b)", side=3, padj=0,adj=0,cex=8)
set.seed(100)
g <- graph.adjacency(s2, mode = "undirected", weighted =T, diag = F)
g <- as_data_frame(g)
g <- g[order(g[,3],decreasing=T),]
g <- g[1:100,]
g <- graph_from_data_frame(g)
E(g)$color <- "blue"
V(g)$color <- "black"
V(g)$color[grep("word",V(g)$name)] <- "red"
plot(g,vertex.size=2.5,vertex.label.cex=3,main="d-method")
title("h-test (with noise)", cex.main=7)
mtext("(c)", side=3, padj=0,adj=0,cex=8)
set.seed(100)
g <- graph.adjacency(s3, mode = "undirected", weighted =T, diag = F)
g <- as_data_frame(g)
g <- g[order(g[,3],decreasing=T),]
g <- g[1:100,]
g <- graph_from_data_frame(g)
E(g)$color <- "blue"
V(g)$color <- "black"
V(g)$color[grep("word",V(g)$name)] <- "red"
plot(g,vertex.size=2.5,vertex.label.cex=3,main="d-method")
title("d-method (with noise)", cex.main=7)
mtext("(d)", side=3, padj=0,adj=0,cex=8)
dev.off()
```
The result is shown in Figure 5 in the manuscript.
![Image text](https://raw.githubusercontent.com/xizhou/co-occurrence-analysis/main/fig5.png)

## Applications
We chose a bibliometric study about pelvic organ prolapse (POP) as our real data application to illustrate how the probabilistic model can be used for MeSH word co-occurrence analysis and compared our new result with  the original one.

### Visualizing the result 

**Visualizing p-value matrix**

Firstly, you can unzip zuo.xml file in data directory.The code is in the following:

```r
library(pubMR)
dir <- "~/co-occurrence-analysis"
setwd(dir)
source("./code/code.R")
obj <- txtList(input="./data/zuo.xml",inputType="xml")
#obj <- AB(input="zuo.xml")
obj1=data.table(PMID=obj@PMID,MS=obj@MAJR)
MS <- obj1[,MS]
idx <- sapply(MS,is.null)
obj1 <- obj1[!idx,]
obj1 = obj1 %>% unnest(MS) 
v <- table(obj1[,c("MS","PMID")])
v <- crossprod(t(v))
v1 <- v
diag(v1) <- NA
p=hyp(v,length(obj@PMID))
diag(p) <- 1
#write.csv(p,file="Supplementary_Table_1.csv")
s <- 1-p
g <- graph.adjacency(s, mode = "undirected", weighted =T, diag = F)
g <- as_data_frame(g)
#g <- g[order(g[,3],decreasing=T),]
#g <- g[1:100,]
g <- g[g[,3]>0.95,]
co <- reshape2:::dcast(data=g,to~from)
rownames(co) <- co[,1]
co <- co[,-1]
co[is.na(co)] <- 0
co <- as.matrix(co)

png("fig1.png",w=2000,h=2000)
Corrplot(co,p.mat=1-co,type="full",method="circle",tl.srt=45,tl.col=1,tl.cex=1.5,cl.length=5,is.corr=FALSE,cl.cex=2,pch.cex=4,pch.col="green",insig="label_sig",sig.level=1e-14)
dev.off()

```
ALL function are defined in the file "code.R".  The result in form of a high dimensional p-value matrix (3192×3192) produced by the model is shown in Figure 1 in the manuscript.

![Image text](https://raw.githubusercontent.com/xizhou/co-occurrence-analysis/main/fig1.png)


**Network of MeSH terms**
- Then, we built a network structure of MeSH terms (0 < pval < 0.01) of POP dataset  to explore the relationship among these MeSH terms.

The code are following:

```r
png("fig2.png",w=2000,h=2000)
set.seed(100)
g <- graph.adjacency(s, mode = "undirected", weighted =T, diag = F)
g <- as_data_frame(g)
g <- g[g[,3]>0.99,]
f <- unique(c(g[,1],g[,2]))
write.csv(f,file="g.csv")
g <- graph_from_data_frame(g,directed=F)
E(g)$width <- 1
E(g)$width[E(g)$weight==1] <- 2.5
E(g)$color <- "blue"
V(g)$color <- "black"
plot(g,vertex.size=2.5,vertex.label.cex=2)
dev.off()
```

The result is presented in Figure 2 in the manucript.

![Image text](https://raw.githubusercontent.com/xizhou/co-occurrence-analysis/main/fig2.png)


**Comparison of results**

We compare our results of our method with the results of original method. Then, we visualize the results of comparison.

```r
library(data.table)
f <- fread("./data/venn.csv",header=T)
library(VennDiagram)
library(scales)
library(VennDiagram)
venn.plot <- draw.pairwise.venn(
area1 = length(f[!htest=="",htest]), 
area2 = length(f[!zuo=="",zuo]),
cross.area = sum(f[!htest=="",htest]%in%f[!zuo=="",zuo]),
category = c("h-test", "d-method"),
fill = c("blue", "yellow"), 
lty = "blank",
cex = 2, 
cat.cex = 3,
cat.pos = c(230, 10), 
cat.dist = c(0.07,0.03), 
cat.just = list(c(0, 0), c(0, 0)),
ext.pos = 0, 
ext.dist = -0.13,
ext.length =0, 
ext.line.lwd = 2,
ext.line.lty = "dashed",
alpha=0.7,
euler.d=T)
png(filename="venn.png",width=550,height=300)
grid.draw(venn.plot);
dev.off()
```

We can get the result (Figure 3):


"![Image text](https://raw.githubusercontent.com/xizhou/co-occurrence-analysis/main/venn.png)"



