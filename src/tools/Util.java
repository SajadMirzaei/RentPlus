package tools;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import object.Node;
import object.Tree;

public class Util {
	public static boolean BRANCH_LENGTH = false;
	
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
				Node node = new Node(1);
				node.setId(temp.split(":")[0]);
				node.setBranchLength(temp.split(":")[1]);
				stack.add(node);
				tree.addNode(node);
			}
//			if ("".equals(temp.split(":")[0])) {
			if (')' == previousChar) {
				stack.get(stack.size()-1).setBranchLength(temp.split(":")[1]);
			}
			if (')' == s.charAt(j)) {
				List<Node> children = new ArrayList<Node>();
				int capacity = 0;
				while(stack.get(stack.size()-1) != null){
					Node child = stack.remove(stack.size()-1);
					children.add(child);
					capacity += child.getSplit().size();
				}
				Node parent = new Node(capacity);
				parent.getChildren().addAll(children);
				for (Node node : children) {
					node.setParent(parent);
				}
				stack.remove(stack.size()-1);
				stack.add(parent);
				tree.addNode(parent);
			}
			previousChar = s.charAt(j);
			i=j+1;
		}
		Node parent = stack.remove(stack.size()-1);
		tree.addNode(parent);
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
				Node node = new Node(1);
				node.setId(temp);
				stack.add(node);
				tree.addNode(node);
			}
			if (')' == s.charAt(j)) {
				List<Node> children = new ArrayList<Node>();
				int capacity = 0;
				while(stack.get(stack.size()-1) != null){
					Node child = stack.remove(stack.size()-1);
					children.add(child);
					capacity += child.getSplit().size();
				}
				Node parent = new Node(capacity);
				parent.getChildren().addAll(children);
				for (Node node : children) {
					node.setParent(parent);
				}
				stack.remove(stack.size()-1);
				stack.add(parent);
				tree.addNode(parent);
			}
			i=j+1;
		}
		Node parent = stack.remove(stack.size()-1);
		tree.addNode(parent);
		tree.setRoot(parent);
		return tree;
	}
	
	public static InputStream load(String path){
		InputStream input = Util.class.getResourceAsStream(path);
		if (input == null){
			input = Util.class.getResourceAsStream("/"+path);
		}
		return input;
	}
}
