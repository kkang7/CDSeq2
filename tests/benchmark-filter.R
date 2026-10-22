library(CDSeq2)
# A shared DEG in A/B must not become a C-specific positive label.
a <- data.frame(gene = c('shared','null'), logFC = c(1,0), adj.P.Val = c(.001,.9))
b <- a
c <- a
c$adj.P.Val <- .9
set.seed(11)
before <- .Random.seed
x <- filter_genes_benchmark(list(A=a,B=b,C=c),c('A','B','C'),
    n_deg_per_ct=0,n_shared_deg=1,n_non_deg=1,seed=2)
stopifnot(identical(before,.Random.seed))
z <- x$labels[x$labels$gene=='shared',]
stopifnot(z$DE_A,z$DE_B,!z$DE_C)
legacy <- lapply(list(A=a,B=b,C=c),function(d) {
    names(d)[2:3] <- c('log2FoldChange','padj');d
})
y <- filter_genes_benchmark(legacy,c('A','B','C'),n_deg_per_ct=0,
    n_shared_deg=1,n_non_deg=1,seed=2)
stopifnot(identical(x,y))
one <- filter_genes_benchmark(list(A=a), 'A', n_deg_per_ct=1,
    n_shared_deg=0, n_non_deg=1, seed=2)
stopifnot(nrow(one$labels)==2, sum(one$labels$DE_A)==1)
