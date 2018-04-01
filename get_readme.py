data0 = open('index_set_list.txt', 'r')
data1 = []

for records in data0:
	data1.append(records.split())
data0.close()

data2 = open('README.md', 'w')

data2.write('# Snapshot\n')
data2.write('\n')
data2.write('### Epigenetic modification of chromatin plays a pivotal role in regulating gene expression during cell differentiation. The scale and complexity of epigenetic data pose a significant challenge for biologists to identify the regulatory events controlling each stage of cell differentiation. Here, we present a model-free method, called Snapshot, that uses epigenetic data to generate a hierarchical visualization for the DNA regions segregating with respect to chromatin state along any given cell differentiation hierarchy of interest. Different cell type hierarchies may be used to highlight the epigenetic history specific to particular lineages of cell differentiation. We demonstrate the utility of Snapshot using data from the VISION project, an international project for ValIdated Systematic IntegratiON of epigenomic data in mouse and human hematopoiesis.\n')
data2.write('\n')
data2.write('### Hematopoietic cell differentiation in VISION (ValIdated Systematic IntegratiON of hematopoietic epigenomes) project\n')
data2.write('\n')
data2.write('\n')
for info in data1:
	### write index-set header
	data2.write('###### The cell differentiation tree for index-set: ' + info[0] + '\n')
	index_set_id = info[0].split('.')[0]
	index_set_pattern = info[0].split('.')[1]
	data2.write('<img src="https://github.com/guanjue/vision_index_set/blob/master/signal_tree/'+str(index_set_id)+'.signal_list.txt'+str(index_set_pattern)+'.tree.png" width="300"/> <img src="https://github.com/guanjue/vision_index_set/blob/master/fun_tree/'+str(index_set_id)+'.function_list.txt'+str(index_set_pattern)+'.tree.png" width="300"/> ' + '\n')
	data2.write('<img src="https://github.com/guanjue/vision_index_set/blob/master/signal_violin/'+info[0]+'.violin.png" width="300"/> <img src="https://github.com/guanjue/vision_index_set/blob/master/fun_bar/'+info[0]+'.bar.png" width="300"/> ' + '\n')
	data2.write('\n')

data2.close()