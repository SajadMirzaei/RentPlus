package object;

import java.util.List;

public class DistanceObject {
	private List<Node> nodes;
	private Float distance;
	private Integer compatibleRange;
	public List<Node> getNodes() {
		return nodes;
	}
	public void setNodes(List<Node> nodes) {
		this.nodes = nodes;
	}
	public Float getDistance() {
		return distance;
	}
	public void setDistance(Float distance) {
		this.distance = distance;
	}
	public Integer getCompatibleRange() {
		return compatibleRange;
	}
	public void setCompatibleRange(Integer compatibleRange) {
		this.compatibleRange = compatibleRange;
	}
	
	@Override
	public int hashCode() {
		return nodes.hashCode() ^ distance.hashCode();
	}
	
	@Override
	public String toString() {
		return "<[" + nodes.get(0) + "," + nodes.get(1) + "], " + distance + ", " + compatibleRange + ">";
	}
}
