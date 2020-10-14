import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

labels = []
final_seqs = []

print("Reading input file.")
with open(input_file, "r") as f:
	print("Cleaning input...")
	for index, seq in enumerate(f.readlines()):
		print(index)
		if index % 2 != 0:   # Odd lines are seqs
			new_seq = ''
			for i in seq:
				if i != 'A' and i!='T' and i!='C' and i!='G' and i!='\n':
					new_seq += 'N'
				else:
					new_seq += i
			final_seqs.append(new_seq)
		else:   # even lines are labels
			labels.append(seq)
	print("Done.")

print("Writing output file.")
with open(output_file, "w") as o:
	assert len(labels) == len(final_seqs)
	for i in range(len(labels)):
		o.write(labels[i])
		o.write(final_seqs[i]+'\n')
		
print("Successfully cleaned input.")


