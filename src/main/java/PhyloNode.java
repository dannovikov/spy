import javafx.util.Pair;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;

import java.util.HashSet;
import java.util.Set;
import java.util.LinkedList;
import java.util.PriorityQueue;

public class PhyloNode implements Comparable<PhyloNode>{
    static int max_id = -1;
    public int id; 
    Set<String> seq_ids;
    DNASequence seq;
    public LinkedList<PhyloNode> parents;
    public int dist_ref;
    public int dist_closest;

    public PhyloNode(String id, DNASequence seq, int dist_ref) {
        seq_ids = new HashSet<>();
        seq_ids.add(id);
        this.id = ++max_id;
        this.seq = seq;
        this.dist_ref = dist_ref;
        this.dist_closest = dist_ref;
        this.parents = new LinkedList<PhyloNode>();
    }

    public PhyloNode(PhyloNode n) {
        this.seq_ids = n.seq_ids;
        this.id = n.id;
        this.seq = n.seq;
        this.dist_ref = n.dist_ref;
        this.dist_closest = n.dist_closest;
        this.parents=n.parents;
    }

    @Override
    public int compareTo(PhyloNode o)
    {
        //Dijkstra
        //Compare distance to ref only
        if (this.dist_ref > o.dist_ref) return 1;
        if (this.dist_ref < o.dist_ref) return -1;
        return 0;
    }


    public void updateSeq(String id, DNASequence seq)
    {
        /* Handles hamming distance == 0 case. Add sequence to parent. */
        seq_ids.add(id);
        char[] seq_arr = this.seq.getSequenceAsString().toCharArray();

        // Not sure why this is necessary. It seems like it won't change any position in seq_arr.
        for (int i=0; i<seq_arr.length; ++i)
            if (seq_arr[i]!='A' && seq_arr[i]!='C' && seq_arr[i]!= 'T' && seq_arr[i]!='G')
                seq_arr[i]=seq.getCompoundAt(i+1).getShortName().charAt(0);
        try {
            this.seq = new DNASequence(new String(seq_arr));
        } catch (CompoundNotFoundException e) {
            e.printStackTrace();
        }
    }


    public DNASequence getSeq() {
        return this.seq;
    }
    public int getHDC() {
        return this.dist_closest;
    }
    public int getHDR() {
        return this.dist_ref;
    }

    public String toString() {
        return Integer.toString(id);
    }

    public int asInt() {
        return id;
    }

    public int hashCode() {
        return toString().hashCode();
    }

    public boolean equals(Object o) {
        return (o instanceof PhyloNode) && (toString().equals(o.toString()));
    }

}

