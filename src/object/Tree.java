package object;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import tools.Util;

public class Tree {
	private Node root;
	private List<Node> nodes = new ArrayList<Node>();
	private Set<Set<Short>> splits = new HashSet<Set<Short>>();
	private String string = "";
	private boolean updateFlag = true;
	private HashMap<Node, Double> nodeLengthMap;

	public Node getRoot() {
		return root;
	}

	public void setRoot(Node root) {
		this.root = root;
	}

	public List<Node> getNodes() {
		return nodes;
	}

	public void setNodes(List<Node> nodes) {
		this.nodes = nodes;
	}

	public Set<Set<Short>> getSplits() {
		return splits;
	}

	public void setSplits(Set<Set<Short>> splits) {
		this.splits = splits;
	}
	
	public void updateSplits(){
		updateSplit(root);
		splits.clear();
		for (Node node : nodes) {
			splits.add(node.getSplit());
		}
	}

	private void updateSplit(Node node) {
		if (node.getChildren().isEmpty()) {
			node.getSplit().add(Short.valueOf(node.getId()));
			return;
		}
		node.getSplit().clear();
		for (Node child : node.getChildren()) {
			updateSplit(child);
			node.getSplit().addAll(child.getSplit());
		}
	}
	
	public List<Node> getTaxa(){
		List<Node> taxa = new ArrayList<Node>();
		for (Node node : nodes) {
			if (node.isLeaf()) {
				taxa.add(node);
			}
		}
		return taxa;
	}

	public Node getExactNodeForSplit(Set<Short> bestSplit) {
		for (Node node : nodes) {
			if (node.getSplit().equals(bestSplit)) {
				if (node.getChildren().size() != 1) {
					return node;
				}
			}
		}
		return null;
	}
	@Override
	public String toString() {
		if (Util.BRANCH_LENGTH) {
			return toStringWithBranchLength();
		}
//		if (updateFlag) {
			string = Util.getString(root);
			updateFlag = false;
//		}
		return string;
	}
	
	public String toString2() {
		return Util.getString(root, new HashMap<Node, String>());
	}
	
	public String toStringWithBranchLength() {
		string = Util.getStringWithBrachLengths(root);
		return string;
	}

	public void update() {
		updateNodes();
		updateSplits();
	}

	public void updateNodes() {
		nodes.clear();
		updateNode(root);
	}

	private void updateNode(Node node) {
		if (!nodes.contains(node)) {
			nodes.add(node);
		}else{
			return;
		}
		if (node.isLeaf()) {
			return;
		}
		for (Node child : node.getChildren()) {
			updateNode(child);
		}		
	}
	
	public Node getBestSplitNode(Set<Short> split) {
		Node bestNode = root;
		for (Node node : nodes) {
			Set<Short> nodeSplit = node.getSplit();
			if (node.getSplit().equals(split))
			{
				return node;
			}
			if (nodeSplit.containsAll(split)) {
				if (nodeSplit.size() < bestNode.getSplit().size() && node.getChildren().size() > 1) {
					bestNode = node;
				}
			}
		}
		return bestNode;
	}
	
	public void cleanUp(){
		boolean done = false;
		while(!done){
			List<Node> nodesToRemove = new ArrayList<Node>();
			for (Node node : nodes) {
				if (node.getId().equals("5318")) {
					System.err.println();
				}
				if (node.isLeaf() && !node.getId().contains("n")) {
					nodesToRemove.add(node);
					node.getParent().getChildren().remove(node);
					if (node.getParent2()!= null) {
						node.getParent2().getChildren().remove(node);
					}
					break;
				}
				if (node.getChildren().size() == 1){
					if (node.getParent() == null) {
						nodesToRemove.add(node);
						root = node.getChildren().get(0);
						break;
					}else{
						Node child = node.getChildren().get(0);
						child.setParent(node.getParent());
						node.getParent().getChildren().remove(node);
						node.getParent().getChildren().add(child);
						nodesToRemove.add(node);
						break;
					}
				}
			}
			if (nodesToRemove.size() == 0) {
				done = true;
			}else{
				nodes.removeAll(nodesToRemove);
			}
		}
	}
	
	public double getSplitLength(Set<Integer> split){
		if (nodeLengthMap == null) {
			nodeLengthMap = new HashMap<Node, Double>();
		}
		Node bestNode = root;
		for (Node node : nodes) {
			if (node.getSplit().containsAll(split)) {
				if (node.getSplit().size() < bestNode.getSplit().size()) {
					bestNode = node;
				}
			}
		}
		double length = getLengthOfNode(bestNode);
		return length;
	}

	private double getLengthOfNode(Node node) {
		if (node.isLeaf()) {
			return 0;
		}
		if (nodeLengthMap.get(node) != null) {
			return nodeLengthMap.get(node);
		}else{
			double longest = 0;
			for (Node child : node.getChildren()) {
				double childLength = Double.valueOf(child.getBranchLength()) + getLengthOfNode(child);
				if (childLength > longest) {
					longest = childLength;
				}
			}
			nodeLengthMap.put(node, longest);
			return longest;
		}
	}
	
	public void setUpdateFlag(boolean updateFlag) {
		this.updateFlag = updateFlag;
	}
	
	public List<Set<Short>> getClades(){
		List<Set<Short>> list = new ArrayList<Set<Short>>();
		getClade(root, list);
		return list;
	}
	
	public Set<Short> getClade(Node node, List<Set<Short>> list){
		if (node.getChildren().size() == 0) {
			list.add(node.getSplit());
			return node.getSplit();
		}
		Set<Short> clade = new HashSet<Short>();
		for (Node child : node.getChildren()) {
			clade.addAll(getClade(child, list));
		}
		list.add(clade);
		return clade;
	}
}
