import javafx.util.Pair;
import org.jgrapht.graph.DefaultEdge;

import java.util.Set;

public class PhyloEdge extends DefaultEdge {
    Set<Pair<Integer,Character>> mutations;
    public PhyloEdge(Set<Pair<Integer,Character>> mutations) {
        this.mutations = mutations;
    }

    @Override
    public String toString() {
        return "(" + getSource() + " : " + getTarget() + ")";
    }
}
