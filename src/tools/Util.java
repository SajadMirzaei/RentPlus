package tools;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import object.Node;
import object.Tree;

public class Util {
	public static Tree getPerfectPhylogeny(int[][] matrix, int[] root){
		Tree result = new Tree();
		// Getting related matrix
		int[][] sortedMatrix = sortMatrix(matrix, root);
		int[][] leavesEdges = new int[matrix.length][matrix[0].length];
		for (int i = 0; i < leavesEdges.length; i++) {
			int counter = 0;
			for (int j = 0; j < leavesEdges[0].length; j++) {
				if (sortedMatrix[i][j] != root[j]) {
					leavesEdges[i][counter] = j+1;
					counter++;
				}
			}
		}
		// Building the tree
		Node rootNode = new Node();
		result.setRoot(rootNode);
		for (int i = 0; i < leavesEdges.length; i++) {
			Node temp = rootNode;
			for (int j = 0; j < leavesEdges[i].length; j++) {
				boolean reach = false;
				if (leavesEdges[i][j] != 0) {
					while(!reach){
						for (Node node : temp.getChildren()) {
							if (node.getId().equals(String.valueOf(leavesEdges[i][j]))) {
								temp = node;
								reach = true;
							}
						}
						if (!reach) {
							Node node = new Node();
							node.setId(String.valueOf(leavesEdges[i][j]));
							temp.getChildren().add(node);
							node.setParent(temp);
							result.getNodes().add(node);
							temp = node;
							reach = true;
						}
					}
				}else {
					Node leaf = new Node();
					leaf.setId("(" + (i+1) + ")");
					temp.getChildren().add(leaf);
					leaf.setParent(temp);
					result.getNodes().add(leaf);
					break;
				}
				if (j == leavesEdges[i].length - 1) {
					Node leaf = new Node();
					leaf.setId("(" + (i+1) + ")");
					temp.getChildren().add(leaf);
					leaf.setParent(temp);
					result.getNodes().add(leaf);
				}
			}
		}
		// Cleaning up the tree
		Set<Node> nodesToRemove = new HashSet<Node>();
		for (Node node : result.getNodes()) {
			if (node.getChildren().size() == 1) {
				Node child = node.getChildren().get(0);
				Node parent = node.getParent();
				parent.getChildren().add(child);
				parent.getChildren().remove(node);
				child.setParent(parent);
				nodesToRemove.add(node);
			}
		}
		result.getNodes().removeAll(nodesToRemove);
		// fixing leaf names
		for (Node node : result.getNodes()) {
			if (node.isLeaf()) {
				String id = node.getId().replaceAll("\\(", "").replaceAll("\\)", "");
				node.setId(id);
			}else{
				node.setId(null);
			}
		}
		return result;
	}
	
	public static String getString(Node node){
		String result = "";
		if (node.isLeaf()) {
			result = node.getId();
		}else if (node.getChildren().size() > 1){
			result = "(";
			List<String> childStrings = new ArrayList<String>();
			for (Node child : node.getChildren()) {
				String childString = getString(child);
				childStrings.add(childString);
			}
			Collections.sort(childStrings);
			for (String childString : childStrings) {
				result = result + childString + ",";
			}
			result = result.substring(0,result.length()-1);
			result += ")";
		}else{
			result = getString(node.getChildren().get(0));
		}
		return result;
	}
	
	public static String getString(Node node, HashMap<Node, String> map){
//		if (node.getId().equals("5397")) {
//			System.err.println();
//		}
		if (map.get(node) != null) {
			return map.get(node);
		}
		String result = "";
		if (node.isLeaf()) {
			result = node.getId();
		}else if (node.getChildren().size() > 1){
			result = "(";
			List<String> childStrings = new ArrayList<String>();
			for (Node child : node.getChildren()) {
				String childString = getString(child, map);
				childStrings.add(childString);
			}
			Collections.sort(childStrings);
			for (String childString : childStrings) {
				result = result + childString + ",";
			}
			result = result.substring(0,result.length()-1);
			result += ")";
		}else{
			result = getString(node.getChildren().get(0), map);
		}
		map.put(node, result);
		return result;
	}
	
	public static String getStringWithBrachLengths(Node node){
		String result = "";
		if (node.isLeaf()) {
			if (node.getBranchLength() != null && !node.getBranchLength().isEmpty()) {
				result = node.getId() + ":" + node.getBranchLength();
			}else{
				result = node.getId();
			}
		}else{
			result = "(";
			for (Node child : node.getChildren()) {
				result = result + getStringWithBrachLengths(child) + ",";
			}
			result = result.substring(0,result.length()-1);
			result += ")";
			if (node.getBranchLength() != null && !node.getBranchLength().isEmpty()) {
				result += ":" + node.getBranchLength();
			}
		}
		return result;
	}

	private static int[][] sortMatrix(int[][] matrix, int[] root) {
		int[] countingArray = new int[matrix[0].length];
		for (int i = 0; i < countingArray.length; i++) {
			for (int j = 0; j < matrix.length; j++) {
				if (matrix[j][i] != root[i]) {
					countingArray[i]++;
				}
			}
		}
		int[] sortedArray = new int[countingArray.length];
		for (int i = 0; i < countingArray.length; i++) {
			int max = -1;
			int maxCol = 0;
			for (int j = 0; j < countingArray.length; j++) {
				if (countingArray[j] > max) {
					max = countingArray[j];
					maxCol = j;
				}
			}
			countingArray[maxCol] = -1;
			sortedArray[i] = maxCol;
		}
		int[][] result = new int[matrix.length][matrix[0].length];
		for (int i = 0; i < sortedArray.length; i++) {
			for (int j = 0; j < result.length; j++) {
				result[j][i] = matrix[j][sortedArray[i]];
			}
		}
		return result;
	}
	
	public static Tree getTreeFromNewickWithTime(String s){
		Tree tree = new Tree();
		List<Node> stack = new ArrayList<Node>(); 
		int i = 0;
		String[] array = s.split("\\)");
		s = s.replace(array[array.length-1],"");
		s = s.replaceAll("\\;", "");
		char previousChar = ',';
		while(i<s.length()){
			int j=i;
			while(',' != s.charAt(j) && ')' != s.charAt(j)){
				if ('(' == s.charAt(j)) {
					stack.add(null);
				}
				j++;
			}
			String temp = s.substring(i,j).replaceAll("\\(","").replaceAll("\\)", "").replaceAll(",", "");
//			if (temp.length() > 0 && !"".equals(temp.split(":")[0])){
			if (temp.length() > 0 && ',' == previousChar){
				Node node = new Node();
				node.setId(temp.split(":")[0]);
				node.setBranchLength(temp.split(":")[1]);
				stack.add(node);
				tree.getNodes().add(node);
			}
//			if ("".equals(temp.split(":")[0])) {
			if (')' == previousChar) {
				stack.get(stack.size()-1).setBranchLength(temp.split(":")[1]);
			}
			if (')' == s.charAt(j)) {
				List<Node> children = new ArrayList<Node>();
				while(stack.get(stack.size()-1) != null){
					Node child = stack.remove(stack.size()-1);
					children.add(child);
				}
				Node parent = new Node();
				parent.getChildren().addAll(children);
				for (Node node : children) {
					node.setParent(parent);
				}
				stack.remove(stack.size()-1);
				stack.add(parent);
				tree.getNodes().add(parent);
			}
			previousChar = s.charAt(j);
			i=j+1;
		}
		Node parent = stack.remove(stack.size()-1);
		tree.getNodes().add(parent);
		tree.setRoot(parent);
		return tree;
	}
	
	public static Tree getTreeFromNewick(String s){
		Tree tree = new Tree();
		ArrayList<Node> stack = new ArrayList<Node>(); 
		int i = 0;
		s = s.replaceAll("\\;", "");
		while(i<s.length()){
			int j=i;
			while(',' != s.charAt(j) && ')' != s.charAt(j)){
				if ('(' == s.charAt(j)) {
					stack.add(null);
				}
				j++;
			}
			String temp = s.substring(i,j).replaceAll("\\(","").replaceAll("\\)", "").replaceAll(",", "");
			if (temp.length() > 0){
				Node node = new Node();
				node.setId(temp);
				stack.add(node);
				tree.getNodes().add(node);
			}
			if (')' == s.charAt(j)) {
				List<Node> children = new ArrayList<Node>();
				while(stack.get(stack.size()-1) != null){
					Node child = stack.remove(stack.size()-1);
					children.add(child);
				}
				Node parent = new Node();
				parent.getChildren().addAll(children);
				for (Node node : children) {
					node.setParent(parent);
				}
				stack.remove(stack.size()-1);
				stack.add(parent);
				tree.getNodes().add(parent);
			}
			i=j+1;
		}
		Node parent = stack.remove(stack.size()-1);
		tree.getNodes().add(parent);
		tree.setRoot(parent);
		return tree;
	}
	
}
