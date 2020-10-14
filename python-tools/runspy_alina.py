import subprocess
import sys
import os

input_source = "/alina-data1/dan/spy/andrew_data/"
programs = ['dijkstra', 'prim1', 'prim2']

def runJar(program):
	global input_source		
	jar = f'{input_source}/jars/{program}.jar'

	command = f"java -jar {jar}\
    -i {input_source}/input_30k.fasta\
    -r {input_source}/ref.fas\
    -e {input_source}/{program}_edges.txt\
    -v {input_source}/{program}_nodes.txt\
    -s {input_source}/{program}_final_seqs.txt\
    -mp"
	print([i.strip() for i in command.split(' ') if i != ''])
	subprocess.call([i.strip() for i in command.split(' ') if i != ''])
	print('Jar executed successfully')

for program in programs:
	runJar(program)