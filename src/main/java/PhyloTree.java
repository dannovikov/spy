import javafx.util.Pair;
import org.biojava.nbio.core.sequence.DNASequence;
import org.jgrapht.Graph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.traverse.BreadthFirstIterator;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;

import java.io.File;
import java.io.FileWriter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

public class PhyloTree {
    Graph<PhyloNode, PhyloEdge> g;
    public static LinkedList<Integer> var_pos;

    LinkedList<PhyloNode> nodes = new LinkedList<PhyloNode>();
    public static PriorityQueue<PhyloNode> pq = new PriorityQueue<PhyloNode>();
    static LinkedList<PhyloNode> collapsed_nodes = new LinkedList<PhyloNode>();
    static PhyloNode last;

    public PhyloTree(LinkedList<Integer> var_pos) {
        g = new DefaultDirectedGraph<>(PhyloEdge.class);
        this.var_pos = var_pos;
    }

    public void buildPhylo(Map<Integer, Set<String>> h_dist_map,
                           LinkedHashMap<String, DNASequence> seqs,
                           String ref_id, DNASequence ref_seq, boolean multiple_parents) {

        //Create root node from reference sequence
        PhyloNode root = new PhyloNode(ref_id, ref_seq, 0);
        nodes.add(root);

        //Create nodes for all sequences in order of distance to root
        Object[] distance_keys = h_dist_map.keySet().toArray();
        Arrays.sort(distance_keys);

        int i = 1;
        for (Object distance_key : distance_keys) 
        {
            System.out.println(String.format("Creating nodes with dist=%d to root. There are %d such sequences.", distance_key, h_dist_map.get(distance_key).size()));
            Set<String> seq_ids = h_dist_map.get(distance_key);
            for (String id : seq_ids) {
                nodes.add(makeNode(id, seqs.get(id), Integer.parseInt(distance_key.toString())));
            }
        }

        last = root;
        //Add all nodes to priority queue. 
        for(PhyloNode n : nodes)
        {
            n.parents.add(root);
            pq.add(n);
        }
        
        //Build Tree
        PhyloNode v;
        while (!pq.isEmpty()) {

            System.out.println(String.format("\nInserting sequence #%d/%d into graph with %d nodes", i, seqs.size() + 1, g.vertexSet().size())); ++i;
            
            v = pq.poll(); 
            System.out.println("Node: " + v);

            if(v.asInt() == root.asInt()) {
                g.addVertex(v);
                continue;
            } 

            // if chooseParents finds a parent with hamming distance 0, then it will collapse the node to that parent
            // otherwise it will find the best parent, set the parent of v, and return false.
            boolean node_is_collapsed = chooseParents(v, g, multiple_parents);

            //if v is not collapsed, add v to tree, fill gaps from parent, update parents, add edges
            if(!node_is_collapsed) 
            { 
                g.addVertex(v);
                last = v;

                System.out.println("Adding " + v.parents.size() + "edges.");

                //add edges from all parents to v
                boolean fill = true;  //fill gaps only from first parent
                for (PhyloNode parent: v.parents) 
                {
                    addEdgeAndFillGaps(parent, v, fill);
                    fill = false;
                }

            }
        }
    }



    private PhyloNode makeNode(String seq_name, DNASequence seq, int h_dist_ref) 
    {
        PhyloNode n = new PhyloNode(seq_name, seq, h_dist_ref);
        return n;
    }


    private static boolean chooseParents(PhyloNode x,
                                         Graph<PhyloNode, PhyloEdge> g,
                                         boolean multiple_parents )
    {
        System.out.println("Finding parents for node: " + x);
        int max_dist_ref = -99;
        for(int i = g.vertexSet().size()-1; i >= 0; i--) {

            PhyloNode candidate = (PhyloNode) g.vertexSet().toArray()[i];

            if (candidate.asInt() == x.asInt()) continue;

            int distance_to_candidate = SeqAlgs.hamDist(var_pos, x.seq, candidate.seq);
            int distance_to_current_parent = SeqAlgs.hamDist(var_pos, x.parents.getFirst().seq, x.seq);

            if(distance_to_candidate == 0) {
                collapseNodeToParent(x, candidate);
                return true;
            }

            if (max_dist_ref != -99 && candidate.dist_ref != max_dist_ref) {
                return false;
            }
            

            if (candidate.dist_ref + distance_to_candidate == x.dist_ref && candidate.asInt() != 0) {
                if (distance_to_candidate < distance_to_current_parent)
                    x.parents.clear();
                x.parents.add(candidate); 
                
                if (!multiple_parents) 
                    return false;
                else {
                    max_dist_ref = candidate.dist_ref;
                }
            }
        }
        System.out.println("No parent candidate found");
        return false;
    }


    private void addEdgeAndFillGaps(PhyloNode source, PhyloNode target, boolean fill) 
    {
        if(fill) {
            System.out.println("Filling gaps in node: " + target + " from " + source +  "...");
            LinkedList<Integer> changedPositions = fillGapsFromParent(target, source);
        }
        Set<Pair<Integer, Character>> mutations = SeqAlgs.findMutations(var_pos, source.seq, target.seq);
        System.out.println("Adding edge from " + source + " to " + target);
        System.out.println(mutations.size() + " mutations.");
        g.addEdge(source, target, new PhyloEdge(mutations));
    }


    private static void collapseNodeToParent(PhyloNode v, PhyloNode p) {
        // Instead of adding v to tree, we add its sequence ID to its parent p
        System.out.println("Collapsing node: " + v + " to " + p);
        String seq_id = v.seq_ids.stream().findFirst().get();
        DNASequence seq = v.getSeq();
        p.updateSeq(seq_id, seq);
        v.parents.clear();
        v.parents.add(p); 
        collapsed_nodes.add(v);
    }


    private LinkedList<Integer> fillGapsFromParent(PhyloNode node, PhyloNode parent)
    {
        char[] node_seq = node.seq.getSequenceAsString().toUpperCase().toCharArray();
        LinkedList<Integer> changedPositions = new LinkedList<Integer>();
        for(int i = 0; i < node.seq.getLength(); i++)
        {
            if (node.seq.getCompoundAt(i+1).getShortName().equals("N"))
            {
                if (!parent.seq.getCompoundAt(i+1).getShortName().equals("N"))
                {
                    // Fill position i from parent
                    
                    node_seq[i] = parent.seq.getCompoundAt(i+1).getShortName().charAt(0);
                    changedPositions.add(i);
                }
            }
        }
        try {
            node.seq = new DNASequence(new String(node_seq));
        } catch (CompoundNotFoundException e) {e.printStackTrace();}
        
        System.out.println(changedPositions.size() + " changed positions.");
        return changedPositions;
    }


    public void exportEdgesCsv(File out_file) throws FileNotFoundException {
        PrintWriter f = new PrintWriter(out_file);
        f.println("Source,Target,Number_mutations");
        for (PhyloEdge e: g.edgeSet()) {
            PhyloNode source = g.getEdgeSource(e);
            PhyloNode target = g.getEdgeTarget(e);
            f.println(String.format("%s,%s,%d", source, target, e.mutations.size()));
        }
        f.close();
    }
    public void exportNodesCsv(File out_file) throws FileNotFoundException {
        PrintWriter f = new PrintWriter(out_file);
        f.println("Strain,Vertex");
        for (PhyloNode v: g.vertexSet())
            for(String s: v.seq_ids) f.println(String.format("%s,%s", s, v));
        f.close();
    }
    public void exportSeqsCsv(File out_file) throws FileNotFoundException {
        //Save node_seqs
        PrintWriter f = new PrintWriter(out_file);
        f.println("Node,Seq");
        for(PhyloNode node : g.vertexSet())
            f.println( node + "," + node.seq);
        f.close();
    }  
}

