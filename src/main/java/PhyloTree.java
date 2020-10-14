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
        last = root;

        //Create nodes from input fasta sequences in order of distance to root
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

        // Set initial last inserted node as root
        // this is important for intial node ordering in pq
        last = root;

        //Add all nodes to priority queue. 
        for(PhyloNode n : nodes)
        {
            n.parents.add(root);
            pq.add(n);
        }
        
        //Build Tree
        PhyloNode v;
        while(!pq.isEmpty()) {
            //System.out.println(Arrays.toString(pq.toArray()));
            System.out.println(String.format("\nInserting sequence #%d/%d into graph with %d nodes", i, seqs.size() + 1, g.vertexSet().size())); ++i;
            
            v = pq.poll(); 
            
            System.out.println("Node: " + v);

            boolean node_is_collapsed = false;
            if(v.asInt() == root.asInt())
            {
                g.addVertex(v);
            } 
            else {
                //check if any of v.parents are hamming distance 0. If they are, collapse v.
                for (PhyloNode p: v.parents)
                {
                    if(SeqAlgs.hamDist(var_pos, v.seq, p.seq) == 0)
                    {
                        //Collapse p if it has hamming distance of 0 to its parent
                        collapseNodeToParent(v, p);
                        node_is_collapsed = true;
                        break;
                    }
                }
                //if not collapsed, add v to tree, fill gaps from parent, update parents, add edges
                if(!node_is_collapsed) 
                { 
                    g.addVertex(v);
                    last = v;

                    System.out.println("Adding edges...");

                    boolean parentIsCollapsed = false;

                    //Identify parents of v with minimal hamming distance to root.
                    for (PhyloNode parent : v.parents)
                    { 
                        if(collapsed_nodes.contains(parent))
                        {   
                            //if any parent is collapsed, add edge from that parent's parent, break
                            addEdgeAndFillGaps(parent.parents.getFirst(), v, true);
                            parentIsCollapsed = true;
                            break;
                        }
                    }
                    if (!parentIsCollapsed)
                    {
                        //if there are no collapsed parents, then add edges from all parents with minimial hamming distance to root
                        boolean fill = true;
                        for (PhyloNode parent: v.parents) 
                        {
                            addEdgeAndFillGaps(parent, v, fill);
                            fill = false;
                        }
                    }
                }
            }
        }
    }



    private void addEdgeAndFillGaps(PhyloNode source, PhyloNode target, boolean fill) 
    {
        if(fill) {
            System.out.println("Filling gaps in node: " + target + "...");
            LinkedList<Integer> changedPositions = fillGapsFromParent(target, source);
            System.out.println("Updating parents...");
            for (Object x: pq.toArray())
                updateParents((PhyloNode) x);
        }
        Set<Pair<Integer, Character>> mutations = SeqAlgs.findMutations(var_pos, source.seq, target.seq);
        g.addEdge(source, target, new PhyloEdge(mutations));
    }

    private PhyloNode makeNode(String seq_name, DNASequence seq, int h_dist_ref) 
    {
        PhyloNode n = new PhyloNode(seq_name, seq, h_dist_ref);
        return n;
    }

    private void collapseNodeToParent(PhyloNode v, PhyloNode p) {
        //Instead of adding v to tree, we add its sequence ID to its parent
        System.out.println("Collapsing node: " + v);

        String seq_id = v.seq_ids.stream().findFirst().get();
        DNASequence seq = v.getSeq();

        p.updateSeq(seq_id, seq);

        v.parents.clear();
        v.parents.add(p); 

        collapsed_nodes.add(v);
    }

    // private void updateParents(PhyloNode x) 
    // { 
    //     //dijkstra
    //     int distance_last_inserted = SeqAlgs.hamDist(var_pos, x.seq, last.seq);
    //     int distance_first_parent = SeqAlgs.hamDist(var_pos, x.parents.getFirst().seq, x.seq);

    //     if(last.dist_ref + distance_last_inserted == x.dist_ref)
    //     {
    //         if (distance_last_inserted < distance_first_parent)
    //         {
    //             x.parents.clear();
    //         }
    //         x.parents.add(last); 
    //     }
    //     PhyloNode closest_to_ref = x.parents.getFirst();
    //     LinkedList<PhyloNode> parent_candidates = new LinkedList<PhyloNode>();
    //     parent_candidates.add(closest_to_ref);

    //     //Keep only the parents that are closest to root
    //     for (PhyloNode parent : x.parents) {
    //         if(parent.dist_ref < closest_to_ref.dist_ref) 
    //         {   //Otherwise, find closest parent to root
    //             closest_to_ref = parent;
    //             parent_candidates.clear();
    //             parent_candidates.add(closest_to_ref);
    //         }
    //         else if (parent.dist_ref == closest_to_ref.dist_ref)
    //         {   //if there are multiple closest nodes to root, append them all to parent candidates
    //             parent_candidates.add(parent);
    //         }
    //     }
    //     x.parents = (LinkedList<PhyloNode>) parent_candidates.clone();
    // }

    private void updateParents(PhyloNode x) 
    { 
        //prim
        int distance_last_inserted = SeqAlgs.hamDist(var_pos, x.seq, last.seq);
        int distance_first_parent = SeqAlgs.hamDist(var_pos, x.parents.getFirst().seq, x.seq);

        if(distance_last_inserted < x.dist_closest)
        {
            decreaseKey(x, distance_last_inserted);
            x.parents.clear();
            x.parents.add(last); 
        }
        else if(distance_last_inserted == x.dist_closest)
        {
            x.parents.add(last);
        }

        PhyloNode closest_to_ref = x.parents.getFirst();
        LinkedList<PhyloNode> parent_candidates = new LinkedList<PhyloNode>();
        parent_candidates.add(closest_to_ref);

        //Keep only the parents that are closest to root
        for (PhyloNode parent : x.parents) {
            if(parent.dist_ref < closest_to_ref.dist_ref) 
            {   //Otherwise, find closest parent to root
                closest_to_ref = parent;
                parent_candidates.clear();
                parent_candidates.add(closest_to_ref);
            }
            else if (parent.dist_ref == closest_to_ref.dist_ref)
            {   //if there are multiple closest nodes to root, append them all to parent candidates
                parent_candidates.add(parent);
            }
        }
        x.parents = (LinkedList<PhyloNode>) parent_candidates.clone();
    }

    private void decreaseKey(PhyloNode x, int distance_closest)
    {
        pq.remove(x);
        System.out.println("Updating key of " + x + " from " + x.dist_closest + " to " + distance_closest);
        x.dist_closest = distance_closest;
        pq.add(x);
    }


    private LinkedList<Integer> fillGapsFromParent(PhyloNode node, PhyloNode parent)
    {
        //System.out.println("Filling gaps for node " + node + " from " + parent);
        LinkedList<Integer> changedPositions = new LinkedList<Integer>();
        for(int i : var_pos)
        {
            if (node.seq.getCompoundAt(i+1).getShortName().equals("N"))
            {
                if (!parent.seq.getCompoundAt(i+1).getShortName().equals("N"))
                {
                    //System.out.println("Filling gaps for node " + node + " from " + parent);
                    char[] node_seq = node.seq.getSequenceAsString().toUpperCase().toCharArray();
                    node_seq[i] = parent.seq.getCompoundAt(i+1).getShortName().charAt(0);
                    try {
                        node.seq = new DNASequence(new String(node_seq));
                    } catch (CompoundNotFoundException e) {
                        e.printStackTrace();
                    }
                    changedPositions.add(i);
                }
            }
        }
        System.out.println( changedPositions.size() + " changed positions.");
        return changedPositions;
    }


    private List<PhyloNode> findParentCandidates(DNASequence seq) {
        int min_h_dist = Integer.MAX_VALUE;
        List<PhyloNode> parent_candidates = new LinkedList<>();
        for (PhyloNode v:  g.vertexSet()) {
            int h_dist = SeqAlgs.hamDist(var_pos, v.seq, seq);
            if (h_dist == min_h_dist) parent_candidates.add(v);
            if (h_dist < min_h_dist) {
                min_h_dist = h_dist;
                parent_candidates = new LinkedList<>();
                parent_candidates.add(v);
            }
        }
        return parent_candidates;
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

