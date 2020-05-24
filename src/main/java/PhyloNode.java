import javafx.util.Pair;
import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;

import java.util.HashSet;
import java.util.Set;

public class PhyloNode {
    static int max_id = -1;
    int id;
    Set<String> seq_ids;
    DNASequence seq;
    public PhyloNode(String id, DNASequence seq) {
        seq_ids = new HashSet<>();
        seq_ids.add(id);
        this.id = ++max_id;
        this.seq = seq;
    }
    public String toString() {
        return Integer.toString(id);
    }

    public int hashCode() {
        return toString().hashCode();
    }

    public boolean equals(Object o) {
        return (o instanceof PhyloNode) && (toString().equals(o.toString()));
    }
    public void updateSeq(String id, DNASequence seq) {
        seq_ids.add(id);
        char[] seq_ar = this.seq.getSequenceAsString().toCharArray();
        for (int i=0; i<seq_ar.length; ++i)
            if (seq_ar[i]!='A' && seq_ar[i]!='C' && seq_ar[i]!= 'G' && seq_ar[i]!='T')
                seq_ar[i]=seq.getCompoundAt(i+1).getShortName().charAt(0);
        try {
            this.seq = new DNASequence(new String(seq_ar));
        } catch (CompoundNotFoundException e) {
            e.printStackTrace();
        }
    }
}
