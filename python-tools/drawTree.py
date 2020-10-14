import networkx as nx
import pandas as pd
import pylab as plb
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pygraphviz
import sys
from PIL import Image

import os
os.getcwd()

if "dijkstra" in sys.argv:
    program = 'dijkstra'
elif "prim1" in sys.argv:
    program = 'prim1'
elif "prim2" in sys.argv:
    program = 'prim2'
else:
    program = 'dijkstra'

if "showlabels" in sys.argv:
    show_labels = True
else:
    show_labels = False

if "collapseSiblings" in sys.argv:
    user_collapse = True
else:
    user_collapse = False


#user_collapse = True
if show_labels:
    viz_out_file = f'PhyloTree_visualization_{program}_labeled.png'
else:
    viz_out_file = f'PhyloTree_visualization_{program}.png'

if '-source' in sys.argv:
    source_dir = sys.argv[-1]
else:
    source_dir = f'C:\\Users\\Dan\\Desktop\\projects\\spy\\run\\paint_tests\\{program}'
#source_dir = r'../spy/tests/cprime/peek/mp/subset/500/final/mp'
#read nodes, edges, metadata
nodes = nodes_file = pd.read_csv(fr'{source_dir}/viz_nodes.txt', header=0, delimiter=',')
edges = edges_file = pd.read_csv(fr'{source_dir}/viz_edges.txt', header=0, delimiter=',')
metadata = metadata_file = pd.read_csv(fr'{source_dir}/metadata.csv', delimiter=',')

final_edges = edges

#Count node frequency
nodeCount = {}
for index, row in nodes.iterrows():
    if row['Vertex'] not in nodeCount:
        nodeCount.update({row['Vertex']: 1})
    else:
        nodeCount[row['Vertex']] += 1



nodeCount = pd.DataFrame.from_dict(nodeCount, orient='index')



#Create Graph
G = nx.DiGraph()

#Add nodes to G 
G.add_nodes_from(i for i in range(len(nodeCount)))
 


G.nodes()



#Add edges to G
for index, row in edges.iterrows():
    s, t = row['Source'], row['Target']
    G.add_edge(int(s), int(t))



#Node size should be proportional to number of sequences associated 
size_multiplier = 1
node_size_list = nodeCount.values * size_multiplier


from networkx.drawing.nx_agraph import graphviz_layout, pygraphviz_layout



#Collapse Leaf nodes
#Count how many leaves are removed per node



#Collapsing Leaf Nodes

#Select first node in group of leaf nodes and collapse all others
leafNodes = []
leafNodes_kept = []
parents = []
for node in G.nodes:
    if G.out_degree[node] == 0:
        #print('node:', node, 'parent:', list(G.predecessors(node))[0])
        try:
            parent = list(G.predecessors(node))[0] #there is only one predecessor
            if parent not in parents:
                parents.append(parent)
                leafNodes_kept.append(node)
                for child in list(G.successors(parent)):
                    if G.out_degree[child] == 0:
                        if child != node:
                            leafNodes.append(child)
        except:
            print(node)
    

#Create list of nodes without the leaves
G_noleaves = G.copy()


#Create list of edges that were removed during leaf-collapse
if user_collapse:
    G_noleaves.remove_nodes_from(leafNodes)
    removed_edges = [edge for edge in G.edges if edge not in list(G_noleaves.edges)]


    #Create dictionary of collapsed_edges per node
    collapsed_edges = {}
    for node in G.nodes():
        e = G.edges(node)
        for edge in e:
            if edge in removed_edges:
                if node in collapsed_edges:
                    collapsed_edges[node].append(edge)
                else:
                    collapsed_edges.update({node:[edge]})


    #Create dict of counts of collapsed leaves per parent
    collapsed_count = {}
    for node, edges in collapsed_edges.items():
        for edge in edges:
            if node in collapsed_count:
                collapsed_count[node] += 1
            else:
                collapsed_count.update({node : 1})
    collapsed_count = pd.DataFrame.from_dict(collapsed_count, orient='index')
    collapsed_count



    #Create new node_size_list for G_noleaves  
    #removing leafnodes that are no longer in G_noleaves

    noleaves_nodeCount = nodeCount.copy() #contains sizes for nodes not in G
    for node in leafNodes:
        try:
            noleaves_nodeCount = noleaves_nodeCount.drop(node)   
        except:
            print(node)


    collapsed_node_size = noleaves_nodeCount.copy()

    #add to collapsed_node_size[node] the sizes of its siblings
    parents=[]
    for leaf in leafNodes:
        parent = list(G.predecessors(leaf))[0]
        if parent not in parents:
            parents.append(parent)
            for i in list(G_noleaves.successors(parent)):
                if G_noleaves.out_degree(i) == 0:   
                    sibling = i
                    break      #choose one sibling
            for child in list(G.successors(parent)):
                if G.out_degree(child) == 0:
                    if child != sibling: #for all other siblings add size to sibling
                        #perhaps here we can also add the locations? or perhaps strains?
                        collapsed_node_size.loc[sibling]+=nodeCount.loc[child] 



#join on Strain
metadata = metadata.rename(columns={"strain":"Strain"})
try:
    nodes_meta = pd.merge(nodes, metadata, on='Strain', how='left')
except:
    nodes_meta = pd.concat(nodes, metadata, on='Strain', how='left')
nodes_meta[0:1]



#replace the location columns with one unified location column:
#for coast to coast we need a to check if 
def buildLocationColumn(df):
    regions = ['Asia', 'Oceania', 'Europe','South America']
    states = ['Minnesota','Illinois', 'Wisconsin',
             'Washington', 'California', 'Texas',
              'Arizona', 'Massachusetts','Connecticut']
    df['location'] = ''
    for i, row in df.iterrows():
        if df.at[i,'division'] in states:
            df.at[i,'location'] = df.at[i,'division']
        elif df.at[i,'region'] in regions:
            df.at[i,'location'] = df.at[i,'region']
        else:
            print(i,'location not listed')
            df.at[i, 'color'] = df.at[i,'region']
    return df[['Strain', 'Vertex', 'location']]



node_locations = buildLocationColumn(nodes_meta)

color_dict = {
    'Asia': '#965A27',
    'Europe': '#C77B23',
    'Oceania': '#E4B646',
    'South America': '#eee066',
    'Minnesota': '#420db7',
    'Illinois': '#4337d7',
    'Wisconsin': '#5a68de',
    'Washington': '#5800b0',
    'California': '#782fc6',
    'Texas': '#985cd9',
    'Arizona': '#ba8aee',
    'Massachusetts': '#dfbcfe',
    'Connecticut': '#b4f8e2',
    'North America': '#000000',
    'default': '#000000',
}

color_dict_rgb = {
    'Asia': [i/255 for i in [150, 90, 39]],
    'Europe': [i/255 for i in [199, 123, 35]],
    'Oceania': [i/255 for i in [228, 182, 70]],
    'South America': [i/255 for i in [238, 224, 102]],
    'Minnesota': [i/255 for i in [66, 13, 183]],
    'Illinois': [i/255 for i in [67, 55, 215]],
    'Wisconsin': [i/255 for i in [90, 104, 222]],
    'Washington': [i/255 for i in [88, 0, 176]],
    'California': [i/255 for i in [120, 47, 198]],
    'Texas': [i/255 for i in [152, 92, 217]],
    'Arizona': [i/255 for i in [186, 138, 238]],
    'Massachusetts': [i/255 for i in [223, 188, 254]],
    'Connecticut': [i/255 for i in [180, 248, 226]],
    'North America': [i/255 for i in [120, 47, 198]],
    'default': [0,0,0]
}



def AppendColors(df, color_dict):
    df['color'] = ''
    for i, row in df.iterrows():
            df.at[i,'color'] = color_dict[row.location]
    return df[['Strain', 'Vertex', 'location','color']]



colored_locations = AppendColors(node_locations, color_dict)



def groupVerticesAndCountLocations(df, loc):
    #Return df counting the number of times each location occurs for each node.
    data = df.loc[:, ['Vertex', loc]]
    data['Count'] = 0
    result = data.groupby(['Vertex', loc]).count()
    return(result)



locations = groupVerticesAndCountLocations(colored_locations, 'location')
ld = locations.to_dict(orient='index') #as in location dictionary



regions = [
'Asia',
'Europe',
'Oceania',
'South',
'Minnesota',
'Illinois',
'Wisconsin',
'Washington',
'California',
'Texas',
'Arizona',
'Massachusetts',
'Connecticut',
]
if user_collapse:
    #Add sibling's regions to collapsed leaf node
    #key of ld is a touple (node, region)
    for leaf in leafNodes:
        for region in regions:
            try: 
                region_count = ld[(leaf, region)]
                parent = list(G.predecessors(leaf))[0] #parent of the leafnode
                for i in list(G_noleaves.successors(parent)): #for each of leaf's parent's children
                    if G_noleaves.out_degree(i) == 0: #find kept leaf node sibling
                        sibling = i
                        print('sibling:', i)
                        break
                try:
                    print('before:', region, ld[(sibling, region)]['Count'])
                    ld[(sibling, region)]['Count'] += region_count['Count'] #add region count
                    print('after:', region, ld[(sibling, region)]['Count'])
                except:
                    ld[(sibling, region)] = region_count
            except:
                continue

    #labeling edges with mutations count
    collapsed_mutations = {}
    for kept_edge in G_noleaves.edges:
        collapsed_mutations[kept_edge] = [final_edges[final_edges.Source == kept_edge[0]][final_edges.Target == kept_edge[1]].Number_mutations.values[0]]
        if kept_edge[1] in leafNodes_kept:
            print(kept_edge, collapsed_mutations[kept_edge])
            for source in collapsed_edges:
                if source == kept_edge[0]:
                    for collapsed_edge in collapsed_edges[source]:
                        collapsed_mutations[kept_edge].append(final_edges[final_edges.Source == collapsed_edge[0]][final_edges.Target == collapsed_edge[1]].Number_mutations.values[0])
else:
    collapsed_node_size = node_size_list

#piechart_locations contains all collapsed node's data
piechart_locations = pd.DataFrame.from_dict(ld, orient='index')
piechart_locations.reset_index(inplace=True)
piechart_locations.columns = ['Vertex', 'location', 'Count']


#define color map
from matplotlib import cm
from matplotlib.colors import ListedColormap

def buildColorList(locations):
    n = len(locations.location.unique())
    print('n=', n)
    cmap = cm.get_cmap('terrain',n+2)
    color_map = {}
    for i,x in enumerate(locations.location.unique()):
        color_map[x] = list(cmap(i))
        
    return color_map


def createColorString(n, locations, color_dict):
        rows = locations[locations['Vertex'] == n] #get rows where node_id is target node n
        #print(n, rows)
        total_regions = sum(rows.Count) #
        s = ''
        for i,c in enumerate(list(rows.location)):
            #get country's color from colordict
            rgb = color_dict[c]
            color = '{}'.format(matplotlib.colors.rgb2hex(rgb))
            #Get country frequency out of total represented regions in node n
            frequency = rows[['location', 'Count']][rows['location'] == c].Count.values[0] / total_regions
            if len(list(rows.location)) == 1:
                #s+='{};{}:{};{}:'.format(color, 1, color, 0)
                s+='{};{}:'.format(color, 1)
            else:
                s+='{};{}:'.format(color, frequency)
        return s

#show_labels = False

gvg = pygraphviz.agraph.AGraph(directed=True) #graphviz graph
collapsed_sibling_color = '#0094ff'
#add nodes
for node in G_noleaves:
    if user_collapse:
        node_label = str(collapsed_node_size.loc[node][0])
        size = np.log(collapsed_node_size.loc[node]) * 2
    else:
        node_label = str(node_size_list[node][0])
        size = np.log(node_size_list[node]) * 2
    
    if size[0] < 3.5:
        size[0] = 3.5
    fill = createColorString(node, piechart_locations, color_dict)
    if node in leafNodes_kept:
        gvg.add_node(node, width=size[0], height=size[0], label='',
                     style='wedged', fillcolor=fill, penwidth=15, color='black')
    else:
        gvg.add_node(node, width=size[0], height=size[0], label='',
                     style='wedged', fillcolor=fill)
    if show_labels:
        gvg.get_node(node).attr['label']= str(node)+'|' + node_label
        
        
#add edges
for edge in list(G_noleaves.edges):
    if user_collapse:
        if edge[1] in leafNodes_kept and len(collapsed_mutations[edge]) > 1:
            edge_label = f'{max(collapsed_mutations[edge])}*'
        else:
            assert len(collapsed_mutations[edge]) == 1
            edge_label = f'{collapsed_mutations[edge][0]}'
    else:
        edge_label = str(final_edges[final_edges.Source == edge[0]][final_edges.Target == edge[1]].Number_mutations.values[0])
    if edge[1] in leafNodes_kept:
        gvg.add_edge(edge, color = 'black', style = 'dotted')
    else:
        gvg.add_edge(edge)
    if show_labels:
        gvg.get_edge(edge[0], edge[1]).attr['label']=edge_label

#draw graph
prog_layout = 'dot'
gvg.layout(prog=prog_layout, args=' -Gratio=0.5 -Gsize=20,8 -Gnodesep=2')
gvg.draw(f'{source_dir}/{viz_out_file}', prog=prog_layout, args=f'-Epenwidth=10                                         -Nfixedsize=True                                         -Gpad=10                                         -Efontsize=150                                         -Nfontsize=125                                        -Nfontcolor=white                                         ')
#                                         -Glabel={program} \
#                                         -Gfontsize=200 \
#                                         -Glabelloc=t \

if "showimages" in sys.argv:
    Image.open((f'{source_dir}/{viz_out_file}')).show()

legend_graph = pygraphviz.agraph.AGraph(directed=True)
legendcounter = 1
for c in color_dict:
    if c == 'Africa':
        fillcolor='black'
    legend_graph.add_node(legendcounter, style='wedged', shape='circle', width = 1.5, height=1.5,
                          fillcolor=f'{matplotlib.colors.rgb2hex(color_dict[c])}:', label=c,
                          fontsize=15,fontcolor='white')
    legendcounter+=1
legend_graph.add_node(legendcounter+1, style='wedged', shape='circle', width=1.5, height=1.5, 
                      fillcolor='black:', label = 'Missing', fontcolor='white', fontsize=15)

legend_graph.add_node(legendcounter+2, style='wedged', shape='circle', width=1.5, height=1.5, 
                      fillcolor='white:', label = 'Collapsed', fontcolor='black', fontsize=15, penwidth=5)
legend_graph.draw('PhyloTree_legend.png', prog='dot')

#Image(filename='PhyloTree_legend.png')

if user_collapse:
    print(program)
    print('edge','\t', 'mutations\t', 'count')
    for edge in G_noleaves.edges:
        print(edge)
        mutations_set = set(collapsed_mutations[edge])
        for i in mutations_set:
            print('\t',i,'\t\t', collapsed_mutations[edge].count(i))

    for i in list(G_noleaves.nodes):
        print(f'({i})',nodes_meta[nodes_meta.Vertex==i][['Strain']].iloc[0])


