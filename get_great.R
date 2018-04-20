library(rGREAT)

bed_list = as.matrix(read.table('bed_list.txt', header = F))
print(bed_list)

go_mf = c()
go_bp = c()
go_cc = c()

for (i in seq(1, dim(bed_list)[1])){
	bed_file = bed_list[i]
	print(bed_file)
	### read bed file
	bed = read.table(paste('index_set_bed/' ,bed_file, sep=''), header = F)[,c(1,2,3)]

	### submit to GREAT
	job = submitGreatJob(bed, species = 'mm10', request_interval = 300)

	pdf(paste('index_set_go_dist/', toString(bed_file), '.dist.pdf', sep=''), width=15, height=5)
	par(mfrow = c(1, 3))
	res = plotRegionGeneAssociationGraphs(job)
	dev.off()

	### get GO analysis result
	tb = getEnrichmentTables(job, category = c("GO"))

	### GO Molecular Function
	tb_x = tb[['GO Molecular Function']]
	tb_x_padj = p.adjust(as.matrix(tb_x[order(tb_x[,1]),][8]), method='fdr')
	### cbind 
	go_mf = cbind(go_mf, tb_x_padj)

	### GO Biological Process
	tb_x = tb[['GO Biological Process']]
	tb_x_padj = p.adjust(as.matrix(tb_x[order(tb_x[,1]),][8]), method='fdr')
	### cbind 
	go_bp = cbind(go_bp, tb_x_padj)

	### GO Cellular Component
	tb_x = tb[['GO Cellular Component']]
	tb_x_padj = p.adjust(as.matrix(tb_x[order(tb_x[,1]),][8]), method='fdr')
	### cbind 
	go_cc = cbind(go_cc, tb_x_padj)
}

tb_x = tb[['GO Molecular Function']]
id = tb_x[order(tb_x[,1]),1]
names = tb_x[order(tb_x[,1]),2]
total_num = tb_x[order(tb_x[,1]),9]
go_mf = cbind(id, names, total_num, go_mf)

tb_x = tb[['GO Biological Process']]
id = tb_x[order(tb_x[,1]),1]
names = tb_x[order(tb_x[,1]),2]
total_num = tb_x[order(tb_x[,1]),9]
go_bp = cbind(id, names, go_bp)

tb_x = tb[['GO Cellular Component']]
id = tb_x[order(tb_x[,1]),1]
names = tb_x[order(tb_x[,1]),2]
total_num = tb_x[order(tb_x[,1]),9]
go_cc = cbind(id, names, go_cc)


write.table(go_mf, file = 'index_set_go_padj/go_mf.txt', sep='\t', quote=F, col.names=F, row.names=F)
write.table(go_bp, file = 'index_set_go_padj/go_bp.txt', sep='\t', quote=F, col.names=F, row.names=F)
write.table(go_cc, file = 'index_set_go_padj/go_cc.txt', sep='\t', quote=F, col.names=F, row.names=F)

write.table(total_num, file = 'index_set_go_padj/go_bp_totalnum.txt', sep='\t', quote=F, col.names=F, row.names=F)



data=read.table('index_set_go_padj/go_bp.txt', header = F, sep='\t', quote='')
bp_p = -log10(as.matrix(data[, c(-1,-2,-3)]))
go_total_num = data[,3]
rownames(bp_p) = data[,2]
colnames(bp_p) = seq(0,dim(bp_p)[2]-1)
bp_p = bp_p[,-dim(bp_p)[2]]
rowmax = apply(bp_p,1,max)
bp_p_sig = bp_p[as.logical((rowmax >= 5) * (go_total_num<=200)),]

bp_p_sig_lim = bp_p_sig
bp_p_sig_lim[bp_p_sig_lim>16] =16

library(pheatmap)
my_colorbar = colorRampPalette(c('white', 'red'))(n = 128)
png('bp_p.png', width=9000,height=9000)
pheatmap(bp_p_sig_lim, color=my_colorbar, cluster_cols = F,cluster_rows=T,annotation_names_row=FALSE,annotation_names_col=F,show_rownames=T,show_colnames=T, clustering_method='ward.D2')
dev.off()



