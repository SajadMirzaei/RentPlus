package object;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import main.Main;
import tools.Util;

public class Node {
	private String id;
	private String branchLength;
	private String info;
	private Node parent;
	private List<Node> children;
	private Set<Short> split;
	
	public static Map<Set<Short>, Set<Short>> generalSplitMap = new HashMap<Set<Short>, Set<Short>>();
	
	public Node() {
		children = new ArrayList<Node>(0);
		split = new HashSet<Short>();
	}
	public Node(int splitSize) {
		children = new ArrayList<Node>(0);
		split = new HashSet<Short>(splitSize);
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
	
	public void setBranchLength(String time) {
		this.branchLength = time;
	}
	public Set<Short> getSplit() {
		return split;
	}
//	public void addToSplit(Set<Short> split){
//		this.split.addAll(split);
//	}
//	public void addToSplit(Short taxa){
//		split.add(taxa);
//	}
	
	public void addToSplit(Set<Short> splitToAdd){
		Set<Short> newSplit = new HashSet<Short>();
		if (generalSplitMap.get(split) != null) {
			newSplit.addAll(generalSplitMap.get(split));
		}
		newSplit.addAll(splitToAdd);
		Set<Short> existingSplit = generalSplitMap.get(newSplit);
		if (existingSplit == null) {
			generalSplitMap.put(newSplit, newSplit);
			this.split = newSplit;
		}else{
			this.split = existingSplit;
		}
	}
	public void addToSplit(Short taxa){
		Set<Short> newSplit = new HashSet<Short>();
		newSplit.add(taxa);
		addToSplit(newSplit);
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
	public void setSplit(Set<Short> split) {
		Set<Short> existingSplit = generalSplitMap.get(split);
		if (existingSplit == null) {
			existingSplit = new HashSet<Short>(split);
			generalSplitMap.put(existingSplit, existingSplit);
		}
		this.split = existingSplit;
	}
}
