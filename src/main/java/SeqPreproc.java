import org.biojava.nbio.core.exceptions.CompoundNotFoundException;
import org.biojava.nbio.core.sequence.DNASequence;

import java.util.*;

class ReferenceWrongAlphabet extends Exception {
    public String toString() {
        return "No gaps or ambiguities are allowed in the reference!";
    }
}

public class SeqPreproc {
    public static DNASequence preprocessReference(DNASequence ref) {
        DNASequence new_ref = null;
        try {
            char[] seq = ref.getSequenceAsString().toUpperCase().toCharArray();
            for (char c : seq)
                if (c != 'A' && c != 'C' && c != 'G' && c != 'T')
                    throw new ReferenceWrongAlphabet();
            new_ref = new DNASequence(new String(seq));
        } catch (Exception e) {
            e.printStackTrace();
        }
        return new_ref;
    }
    public static void preprocessSeqs(LinkedHashMap<String, DNASequence> seqs) {
        for (Map.Entry<String, DNASequence> entry : seqs.entrySet()) {
            char[] seq = entry.getValue().getSequenceAsString().toUpperCase().toCharArray();
            for (int i=0; i<seq.length; ++i)
                if (seq[i]!='A' && seq[i]!='C' && seq[i]!= 'G' && seq[i]!='T') seq[i]='-';
            try {
                seqs.put(entry.getKey(), new DNASequence(new String(seq)));
            } catch (CompoundNotFoundException e) {
                e.printStackTrace();
            }
        }
    }
    public static void fillGapsRef(LinkedHashMap<String, DNASequence> seqs,
                                                                 DNASequence ref) {
        String ref_str = ref.getSequenceAsString();
        for (Map.Entry<String, DNASequence> entry : seqs.entrySet()) {
            char[] seq = entry.getValue().getSequenceAsString().toCharArray();
            for (int i=0; i<ref_str.length(); ++i) {
                if (seq[i]=='-') seq[i]=ref_str.charAt(i);
            }
            try {
                seqs.put(entry.getKey(), new DNASequence(new String(seq)));
            } catch (CompoundNotFoundException e) {
                e.printStackTrace();
            }
        }
    }
    public static LinkedList<Integer> getVariablePositions(LinkedHashMap<String, DNASequence> seqs) {
        LinkedList<Integer> variablePositions = new LinkedList<>();
        int seq_length = seqs.entrySet().iterator().next().getValue().getLength();
        for(int i=0; i<seq_length; ++i) {
            Set<String> nucls = new HashSet<>();
            for (Map.Entry<String, DNASequence> entry : seqs.entrySet()) {
                String nucl = entry.getValue().getCompoundAt(i+1).getShortName();
                if (nucl.equals("-")) {
                    continue;
                }
                nucls.add(nucl);
                if (nucls.size() >= 2) {
                    variablePositions.add(i);
                    break;
                }
            }
        }
        return variablePositions;
    }
    public static Map<Integer, Set<String>> getHamDistToRef(DNASequence ref, LinkedHashMap<String,DNASequence> seqs,
                                                            LinkedList<Integer> var_pos) {
        Map<Integer, Set<String>> h_dist_map = new HashMap<>();
        for (Map.Entry<String, DNASequence> entry : seqs.entrySet()) {
            int h_dist = SeqAlgs.hamDist(var_pos, ref, entry.getValue());
            if (!h_dist_map.containsKey(h_dist)) {
                h_dist_map.put(h_dist, new HashSet<>());
            }
            h_dist_map.get(h_dist).add(entry.getKey());
        }
        return h_dist_map;
    }
}
