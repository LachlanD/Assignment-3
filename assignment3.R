library(edgeR)
library(locfit)
library(gridExtra)
library(grid)

counts = read.delim("assignment3data.tsv")

d=DGEList(counts=counts[,1:30], group=factor(c(rep("F", 15),rep("M",15))))

d=calcNormFactors(d)

d=estimateDisp(d)

et=exactTest(d)

tt = topTags(et, n=10, adjust.method="BH", sort.by="P")
tt$table=round(tt$table,4)
tt$table["Gene"]=counts[rownames(tt),][31]
tt$table["Chromosome"]=counts[rownames(tt),][32]
#tt$table[1:3]=round(tt$table[1:3],4)
grid.table(tt$table)



de = topTags(et, n=15, adjust.method="none", sort.by="P")

heatmap(d$counts[rownames(d) %in% rownames(de),])

ua =decideTestsDGE(et, adjust.method = "none",p.value = 0.05)
bh = decideTestsDGE(et, adjust.method = "BH", p.value = 0.05)
summary(ua)
sum(ua!=0)
sum(ua!=0)/length(ua)
summary(bh)
sum(bh!=0)
sum(bh!=0)/length(bh)

sum(bh[counts[32]=='Y'])
sum(counts[32]=='Y')


et$table[counts[31]=="XIST",]

set.seed(1)
random.group = sample(c(rep("F",15),rep("M",15)), size = 30, replace = F)
random.group

random.dge=DGEList(counts=counts[,1:30], group=factor(random.group))
random.dge=calcNormFactors(random.dge)
random.dge=estimateDisp(random.dge)
r.et=exactTest(random.dge)

r.tt = topTags(r.et, n=10, adjust.method="none", sort.by="P")
r.tt$table["Gene"]=counts[rownames(r.tt),][31]
r.tt$table["Chromosome"]=counts[rownames(r.tt),][32]
grid.table(r.tt$table)


r.de = topTags(r.et, n=15, adjust.method="BH", sort.by="P")

heatmap(random.dge$counts[rownames(random.dge) %in% rownames(r.de),])

r.ua =decideTestsDGE(r.et, adjust.method = "none",p.value = 0.05)
summary(r.ua)
sum(r.ua!=0)
sum(r.ua!=0)/length(r.ua)

r.bh = decideTestsDGE(r.et, adjust.method = "BH", p.value = 0.05)
summary(r.bh)
sum(r.bh!=0)
sum(r.bh!=0)/length(r.bh)

sum(r.bh[counts[32]=='Y'])
sum(counts[32]=='Y')

r.et$table[counts[31]=="XIST",]
