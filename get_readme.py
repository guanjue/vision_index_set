data0 = open('index_set_list.txt', 'r')
data1 = []

for records in data0:
	data1.append(records.split())
data0.close()

data2 = open('README.md', 'w')

data2.write('# Snapshot\n')
data2.write('\n')
data2.write('### Hematopoietic cell differentiation in VISION (ValIdated Systematic IntegratiON of hematopoietic epigenomes) project\n')
data2.write('\n')
data2.write('### Index-set information:\n')
data2.write('\n')
data2.write('\n')

data2.write('#### all index-set info: ')
data2.write('(https://github.com/guanjue/vision_index_set/blob/master/index_set_bed/atac_20cell.fun.txt)' + '\n')
data2.write('\n')
data2.write('\n')

for info in data1:
	### write index-set header
	data2.write('###### The cell differentiation tree for index-set: ' + info[0] + '\n')
	index_set_id = info[0].split('.')[0]
	index_set_pattern = info[0].split('.')[1]
	data2.write('<img src="https://github.com/guanjue/vision_index_set/blob/master/signal_tree/'+str(index_set_id)+'.signal_list.txt'+str(index_set_pattern)+'.tree.png" width="300"/> <img src="https://github.com/guanjue/vision_index_set/blob/master/fun_tree/'+str(index_set_id)+'.ideas_list.txt'+str(index_set_pattern)+'.tree.png" width="300"/> ' + '\n')
	data2.write('<img src="https://github.com/guanjue/vision_index_set/blob/master/signal_violin/'+info[0]+'.violin.png" height="100" width="300"/> <img src="https://github.com/guanjue/vision_index_set/blob/master/fun_bar/'+info[0]+'.bar.png" height="100" width="300"/> ' + '\n')
	data2.write('\n')
	data2.write('bed file for index-set-' + info[0] + ': ')
	data2.write('(https://github.com/guanjue/vision_index_set/blob/master/index_set_bed/'+info[0]+'.index_sed.bed)' + '\n')
	data2.write('\n')

data2.close()