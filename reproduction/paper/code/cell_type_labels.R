# Display labels only; archived data identifiers and numerical inputs are unchanged.
ct_name<-function(x,wrapped=FALSE){
 d<-c(Astro='Astrocytes',Excit='Excitatory neurons',Inhib='Inhibitory neurons',Oligo='Oligodendrocytes',OPC='Oligodendrocyte precursor cells',Immune='Immune cells',Vascular='Vascular cells')
 if(wrapped)d[c('Excit','Inhib','OPC')]<-c('Excitatory\nneurons','Inhibitory\nneurons','Oligodendrocyte\nprecursor cells')
 y<-unname(d[as.character(x)]);y[is.na(y)]<-as.character(x)[is.na(y)];y
}
