import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

labels = []
final_seqs = []

print("Reading input file.")
with open(input_file, "r") as f:
	print("Cleaning input...")
	new_seq = ''
	for index, seq in enumerate(f.readlines()):
		if seq == '\n':
			print("-------")
			continue
		if seq[0] != '>':   # is seq
			for i in seq:
				if i != 'A' and i!='G' and i!='C' and i!='G':
					new_seq += 'N'
				else:
					new_seq += i
		else:   # is label
			if new_seq != '':
				print("adding seq:", new_seq[0:5], "len: ", len(new_seq))
				final_seqs.append(new_seq)
				new_seq = ''
			else:
				print("seq is empty")
			labels.append(seq)
	final_seqs.append(new_seq)
			
			
	print("Done.")
	# for i,x in enumerate(labels):
	# 	print(x)
	# 	print(final_seqs[i][:5])
	print(len(labels), len(final_seqs))
	print(labels[0], final_seqs[0][:5])

print("Writing output file.")
with open(output_file, "w") as o:
	assert len(labels) == len(final_seqs)
	for i in range(len(labels)):
		o.write(labels[i])
		o.write(final_seqs[i]+'\n')
		
print("Successfully cleaned input.")


