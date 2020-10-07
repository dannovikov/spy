import subprocess
import sys
import os
import visualization_adapter as adapter

print(sys.argv)
if "dijkstra" in sys.argv:
	programs = ['dijkstra']
elif "prim1" in sys.argv:
	programs = ['prim1']
elif "prim2" in sys.argv:
	programs = ['prim2']
else:
	programs = ['dijkstra', 'prim1', 'prim2']

input_source = ''

def main():
	global programs
	global input_source

	if'compile' in sys.argv and 'only' in sys.argv:
		for program in programs:
			compile(program)
			exit()

	if 'compile' in sys.argv:
		for program in programs:
			compile(program)


	if "showimages" in sys.argv:
		show_images = True
	else:
		show_images = False

	if "collapseSiblings" in sys.argv:
		collapse_siblings = True
	else:
		collapse_siblings = False

	print('Running programs:', programs)
	if "tests" in sys.argv:
		tests = [i+1 for i in range(12)]
	else:
		tests = []

	if "tests" in sys.argv:
		for test in tests:
			for program in programs:
				input_source = f'C:\\Users\\Dan\\Desktop\\projects\\spy\\tests\\Oct1\\test{test}\\{program}'
				runJar(program)
				runVisualizationAdapter(program)
				drawTree(program, True, show_images, input_source, collapse_siblings)
				drawTree(program, False, show_images, input_source, collapse_siblings)
	else:
		for program in programs:
			input_source = f'C:\\Users\\Dan\\Desktop\\projects\\spy\\run\\full_seqs\\{program}'

			runJar(program)
			runVisualizationAdapter(program)
			drawTree(program, True, show_images,input_source, collapse_siblings)
			drawTree(program, False, show_images,input_source, collapse_siblings)

def compile(program):
	previous_directory = os.getcwd()
	os.chdir("C:\\Users\\Dan\\Desktop\\projects\\spy")
	subprocess.call(['mvn', 'clean','install'], shell=True)
	os.replace('.\\spy\\spy.jar', f'C:\\Users\\Dan\\Desktop\\spy_jars\\{program}.jar')
	os.chdir(previous_directory)

def runJar(program):
	global input_source		
	jar = f'C:\\Users\\Dan\\Desktop\\spy_jars\\{program}.jar'

	command = f"java -jar {jar}\
    -i {input_source}/input.fasta\
    -r {input_source}/ref.fas\
    -e {input_source}/edges.txt\
    -v {input_source}/nodes.txt\
    -s {input_source}/final_seqs.txt\
    -mp"
	print([i.strip() for i in command.split(' ') if i != ''])
	subprocess.call([i.strip() for i in command.split(' ') if i != ''])
	print('Jar executed successfully')


def runVisualizationAdapter(program):
	global input_source

	adapter.adapt(f'{input_source}\\nodes.txt',
				f'{input_source}\\viz_nodes.txt',
				f'{input_source}\\edges.txt',
				f'{input_source}\\viz_edges.txt',
				f'{input_source}\\final_seqs.txt',
				f'{input_source}\\viz_final_seqs.txt')
	print('Jar output adapted successfully')


def drawTree(program, showLabels, show_images, input_source, collapse_siblings):


	if show_images:
		s = "showimages"
	else:
		s=''

	if collapse_siblings:
		if showLabels:
			command = f"python drawTree.py {program} showlabels {s} collapseSiblings -source {input_source}"
		else:
			command = f'python drawTree.py {program} {s} collapseSiblings -source {input_source}'
	else:
		if showLabels:
			command = f"python drawTree.py {program} showlabels {s} -source {input_source}"
		else:
			command = f'python drawTree.py {program} {s} -source {input_source}'
	print([i.strip() for i in command.split(' ') if i != ''])
	subprocess.call([i.strip() for i in command.split(' ') if i != ''])


if __name__ == '__main__':
	main()