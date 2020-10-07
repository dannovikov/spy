import sys

node_input_file = sys.argv[1]
node_output_file = sys.argv[2]
edge_input_file = sys.argv[3]
edge_output_file = sys.argv[4]

def transformNodes(node_input_file, node_output_file):
	with open(node_input_file, 'r') as in_file:
		with open(node_output_file, 'w') as out_file:

			node_dict = {}
			changed_nodes= {}
			header = in_file.readline()
			lines = in_file.readlines()

			#build nodecount dict
			for i, line in enumerate(lines):
				lines[i] = line.strip().split(',')
				lines[i][1] = int(lines[i][1])
				if lines[i][1] in node_dict:
					node_dict[lines[i][1]].append(lines[i][0])
				else:
					node_dict[lines[i][1]] = [lines[i][0]]

			#change node numbering
			node = 0
			final_dict = {}
			while bool(node_dict):
				min_key = min(node_dict)
				final_dict[node] = node_dict.pop(min_key)
				changed_nodes[min_key]=node
				node+=1

			#write outfile
			out_file.write(header)
			for key in final_dict:
				for strain in final_dict[key]:
					out_file.write(f"{strain},{key}\n")
	return changed_nodes


def transformEdges(edge_input_file, edge_output_file, changed_nodes):
	with open(edge_input_file, 'r') as in_file:
		with open(edge_output_file, 'w') as out_file:
			header = in_file.readline()
			lines = in_file.readlines()
			for i, line in enumerate(lines):
				lines[i] = [int(j) for j in line.strip().split(',')]
				lines[i][0], lines[i][1] = changed_nodes[lines[i][0]], changed_nodes[lines[i][1]]
			out_file.write(header)
			for line in lines:
				out_file.write("{},{},{}\n".format(line[0], line[1], line[2]))


changed_nodes = transformNodes(node_input_file, node_output_file)
transformEdges(edge_input_file, edge_output_file, changed_nodes)

