import javafx.util.Pair;
import org.biojava.nbio.core.sequence.DNASequence;
import org.jgrapht.Graph;
import org.jgrapht.graph.DefaultDirectedGraph;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.*;

public class PhyloTree {
    Graph<PhyloNode, PhyloEdge> g;
    LinkedList<Integer> var_pos;
    public PhyloTree(LinkedList<Integer> var_pos) {
        g = new DefaultDirectedGraph<>(PhyloEdge.class);
        this.var_pos = var_pos;
    }
    public void buildPhylo(Map<Integer, Set<String>> h_dist_map,
                           LinkedHashMap<String, DNASequence> seqs,
                           String ref_id, DNASequence ref_seq) {
        PhyloNode root = new PhyloNode(ref_id, ref_seq);
        g.addVertex(root);
        Object[] keys = h_dist_map.keySet().toArray();
        Arrays.sort(keys);
        int i = 1;
        for (Object k : keys) {
            System.out.println(String.format("Inserting sequences with dist=%d to root. There are %d such sequences.", k, h_dist_map.get(k).size()));
            for (String id : h_dist_map.get(k)) {
                System.out.print(String.format("Inserting sequence #%d/%d into graph with %d nodes\r", i, seqs.size(), g.vertexSet().size()));
                ++i;
                insertNode(id, seqs.get(id));
            }
        }
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
    private void insertNode(String seq_name, DNASequence seq) {
        List<PhyloNode> parent_candidates = findParentCandidates(seq);
        if (parent_candidates.size() != 1)
            System.out.println(String.format("Node %s has %d parent candidater. Nodes in graph: %d",
                    seq_name, parent_candidates.size(), g.vertexSet().size()));

        int h_dist = SeqAlgs.hamDist(this.var_pos, parent_candidates.get(0).seq, seq);
        if (h_dist == 0) {
            parent_candidates.get(0).updateSeq(seq_name, seq);
        }
        else{
            PhyloNode n = new PhyloNode(seq_name, seq);
            g.addVertex(n);
            for (PhyloNode parent: parent_candidates) {
                Set<Pair<Integer, Character>> mutations = SeqAlgs.findMutations(this.var_pos, parent.seq, seq);
                System.out.print(mutations);
                g.addEdge(parent, n, new PhyloEdge(mutations));
            }
        }
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
}
