package object;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import tools.Util;

public class Node {
	private String id;
	private String branchLength;
	private String info;
	private Node parent;
	private Node parent2;
	private List<Node> children;
	private Set<Short> split;
	
	public Node() {
		children = new ArrayList<Node>();
		split = new HashSet<Short>();
	}
	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id;
	}
	public Node getParent() {
		return parent;
	}
	public void setParent(Node parent) {
		this.parent = parent;
	}
	public void setParent2(Node parent2) {
		this.parent2 = parent2;
	}
	public Node getParent2() {
		return parent2;
	}
	public List<Node> getChildren() {
		return children;
	}
	public void setChildren(List<Node> children) {
		this.children = children;
	}
	public void setInfo(String info) {
		this.info = info;
	}
	public String getInfo() {
		return info;
	}
	public String getBranchLength() {
		if (branchLength == null || branchLength.isEmpty()) {
			return "0";
		}
		return branchLength;
	}
	
	public Node copy(Node parent, Node parent2, HashMap<String, Node> map){
		Node node = map.get(id);
		if (node != null) {
			if (parent == null && parent2 != null) {
				node.setParent2(parent2);
			}else if (parent != null && parent2 == null){
				node.setParent(parent);
			}
			return node;
		}
		node = new Node();
		node.setId(id);
		node.setBranchLength(branchLength);
		node.setInfo(info);
		map.put(id, node);
		for (Node child : children) {
			Node newChild = new Node();
			if (child != null) {
				if (child.getParent().equals(this)) {
					newChild = child.copy(node, null, map);
				}else{
					newChild = child.copy(null, node, map);
				}
			}else{
				newChild = null;
			}
			node.getChildren().add(newChild);
		}
		if (parent == null && parent2 != null) {
			node.setParent2(parent2);
		}else if (parent != null && parent2 == null){
			node.setParent(parent);
		}
		return node;
	}
	
	public void setBranchLength(String time) {
		this.branchLength = time;
	}
	public Set<Short> getSplit() {
		return split;
	}
	public void setSplit(Set<Short> split) {
		this.split = split;
	}
	@Override
	public String toString() {
		return id;
//		return Util.getString(this);
	}
	public boolean isLeaf() {
		if (children.size() == 0) {
			return true;
		}
		return false;
	}
}
