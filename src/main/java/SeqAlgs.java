import javafx.util.Pair;
import org.biojava.nbio.core.sequence.DNASequence;

import java.util.HashSet;
import java.util.LinkedList;
import java.util.Set;

public class SeqAlgs {
    public static int hamDist(LinkedList<Integer> var_pos, DNASequence s1, DNASequence s2) {
        int h_dist = 0;
        for (Integer i : var_pos) {
            if (s1.getCompoundAt(i).getShortName().toString().equals("-") || s2.getCompoundAt(i).getShortName().toString().equals("-")) continue;
            if (!s1.getCompoundAt(i + 1).equals(s2.getCompoundAt(i + 1))) ++h_dist;
        }
        return h_dist;
    }

    public static Set<Pair<Integer, Character>> findMutations(LinkedList<Integer> var_pos,
                                                              DNASequence parent, DNASequence child) {
        Set<Pair<Integer, Character>> mutations = new HashSet<>();
        for (Integer i : var_pos) {
            if (parent.getCompoundAt(i).getShortName().toString().equals("-") || child.getCompoundAt(i).getShortName().toString().equals("-")) continue;
            if (!parent.getCompoundAt(i+1).equals(child.getCompoundAt(i+1))) {
                mutations.add(new Pair<>(i, child.getCompoundAt(i+1).getShortName().charAt(0)));
            }
        }
        return mutations;
    }
}
