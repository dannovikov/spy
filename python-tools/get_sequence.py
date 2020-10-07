import sys

try:
	seq_ids = names = [sys.argv[1]]
	fasta_path = sys.argv[2]
	output_path = sys.argv[3]
except:
	fasta_path = 'J:\\Workspace\\coast-to-coast\\168_nextstrain_masked_n.fasta'
	output_path = 'J:\\Workspace\\coast-to-coast\\viz_seqs\\prim\\input.fasta'

	# find each seq_id in fasta path and
	# write names and seqs to output_path

	names = ['0','1','2','4','5','11','12','13','19','22','23','24','29','39']
	# dijkstra seq_ids
	# seq_ids = ['Wuhan-Hu-1/2019','USA/CA8/2020','USA/MA1/2020',
	# 			'USA/WA-S4/2020','USA/IL2/2020','USA/WA-UW25/2020',
	# 			'USA/WA-UW75/2020','Switzerland/GE5373/2020',
	# 			'USA/CA1/2020','Wuhan/WIV02/2019']

	#prim seq_ids
	seq_ids = ['Wuhan-Hu-1/2019','USA/CA8/2020','USA/MA1/2020','USA/WA-S4/2020',
				'USA/IL2/2020','USA/WA-UW25/2020','USA/WA-UW75/2020','USA/WA-UW69/2020',
				'USA/WA-UW93/2020','USA/MN3-MDH3/2020','USA/CA1/2020',
				'Wuhan/IPBCAMS-WH-02/2019','Wuhan/WIV02/2019','USA/WA-UW30/2020']

assert len(names) == len(seq_ids)
name_index = 0
for seq_id in seq_ids:
	print('seq_id:', seq_id)
	with open(fasta_path, 'r') as f:
		with open(output_path, 'a') as out:
			found = False
			for i, val in enumerate(f):
				if i % 2 == 0: #is label
					if val[1:].strip() == seq_id:
						print(val[1:].strip(), seq_id,  val[1:].strip() == seq_id)
						found = True
						out.write(val)
						# out.write('>')
						# out.write(names[name_index])
						# out.write('\n')
						out.write(next(f))
						name_index += 1
						break
			if not found:
				print(seq_id, 'Sequence not found in fasta.')

