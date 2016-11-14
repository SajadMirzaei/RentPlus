package object;

import java.util.List;

public class DistanceObject {
	private Node[] nodes;
	private Float distance;
	private Short compatibleRange;
	public Node[] getNodes() {
		return nodes;
	}
	public void setNodes(Node[] nodes) {
		this.nodes = nodes;
	}
	public Float getDistance() {
		return distance;
	}
	public void setDistance(Float distance) {
		this.distance = distance;
	}
	public Short getCompatibleRange() {
		return compatibleRange;
	}
	public void setCompatibleRange(Short compatibleRange) {
		this.compatibleRange = compatibleRange;
	}
	
	@Override
	public int hashCode() {
		return nodes.hashCode() ^ distance.hashCode();
	}
	
	@Override
	public String toString() {
		return "<[" + nodes[0] + "," + nodes[1] + "], " + distance + ", " + compatibleRange + ">";
	}
}
