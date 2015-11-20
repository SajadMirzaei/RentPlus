package main;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import object.Node;
import object.Tree;
import tools.Util;

public class Main {
	private int[][] matrix;
	private int[] root;
	private boolean[][] compatibleMatrix;
	private Tree[] localTrees;
	private Tree[] trueTrees;
	private Tree[] timeTrees;

	private List<double[][]> distanceMatrices = new ArrayList<double[][]>();
	private double[] timeTreeLengths;
	private double[] trueTreeLengths;
	private double[] localTreeLengths;
	private List<Integer> positions;
	private Set<Integer> indicesToRemove;
	private Set<Integer> overalSplit;
	private double TOTAL_SEQ_LENGTH = 0;
	private double UPGMA_THRESHOD = 0.2;
	
	private String outputName;
	
	private int[][] bestRegions;
	private List<Set<Integer>> originalSplits;
	
	public Main(String args[]) {
		outputName = args[0];
		Date startTime = new Date();
		System.out.println("Start Time: " + startTime);
		overalSplit = new HashSet<Integer>();
		originalSplits = new ArrayList<Set<Integer>>();
		buildDataFromHapMatrix(args[0]);
		for (int i = 0; i < matrix.length; i++) {
			overalSplit.add(i+1);
		}
		timeTreeLengths = new double[matrix[0].length];
		trueTreeLengths = new double[matrix[0].length];
		localTreeLengths = new double[matrix[0].length];
		if (args.length > 1) {
			getTrueTreesFromMsFile(args[1]);
			getDistanceMatricesForTrueTrees();
		}
		positions = getPositionArrayFromInputFile(args[0]);
		makeDistanceMatricesFast2();
		makeUPGMATimeTrees2();
		buildNewRoot();
		localTrees = initializeTrees();
		double[] sortedTreeLengths = new double[localTrees.length];
		for (int i = 0; i < sortedTreeLengths.length; i++) {
			sortedTreeLengths[i] = timeTreeLengths[i];
		}
		bestRegions = new int[localTrees.length][2];
		Arrays.sort(sortedTreeLengths);
		outputTreesSeperate(timeTrees, outputName, "Upgma");
		buildLocalTrees();
	    Date endTime = new Date();
	    scaleTimes(timeTreeLengths);
	    if (args.length > 1) {
	    	compareAllTreesRooted(positions);
	    	compareAllTrees(positions);
	    	System.out.println("-----------------------");
		}else{
			outputInferredTrees(positions);
		}
	    outputTreesSeperate(localTrees,args[0], "");
		System.out.println("Running Time: "
				+ (double) (endTime.getTime() - startTime.getTime()) / 1000
				+ " Seconds");
	}



	
	private void makeUPGMATimeTrees2() {
		System.out.print("Making UPGMA Trees");
		Date startDate = new Date();
		timeTrees = new Tree[positions.size()];
		for (int i = 0; i < timeTrees.length; i++) {
//			System.out.println("Index " + i);
			Tree upgmaTree = new Tree();
			double[][] distanceMatrix = distanceMatrices.get(i);
			HashMap<Node, List<Integer>> clusters = new LinkedHashMap<Node, List<Integer>>();
			//initiating clusters
			for (int j = 0; j < matrix.length; j++) {
				Node node = new Node();
				node.setId(String.valueOf(j+1));
				node.setInfo("0.0");
				List<Integer> nodesCovering = new ArrayList<Integer>();
				nodesCovering.add(j+1);
				clusters.put(node, nodesCovering);
			}
			//
			while(clusters.keySet().size() > 1){
				List<Node> clusterNodes = new ArrayList<Node>();
				clusterNodes.addAll(clusters.keySet());
				double minDist = Double.MAX_VALUE;
				Set<Node> nodesToConcat = new HashSet<Node>();
				Map<Set<Node>, Double> validNodesToContactMap = new HashMap<Set<Node>, Double>();
				if (clusterNodes.size() == 2) {
					nodesToConcat.addAll(clusterNodes);
					minDist = getDistanceOfClusters(clusters.get(clusterNodes.get(0)),clusters.get(clusterNodes.get(1)), distanceMatrix);
				}else{
					for (int j = 0; j < clusterNodes.size()-1; j++) {
						for (int j2 = j+1; j2 < clusterNodes.size(); j2++) {
							Node node1 = clusterNodes.get(j);
							Node node2 = clusterNodes.get(j2);
							Set<Integer> combinedSplit = new HashSet<Integer>();
							combinedSplit.addAll(clusters.get(node1));
							combinedSplit.addAll(clusters.get(node2));
							double distance = getDistanceOfClusters(clusters.get(node1),clusters.get(node2), distanceMatrix);
							if (distance <= minDist*(1+UPGMA_THRESHOD) && isCompatible(originalSplits.get(i), combinedSplit)) {
								if (distance < minDist) {
									minDist = distance;
									List<Set<Node>> iterationList = new ArrayList<Set<Node>>();
									iterationList.addAll(validNodesToContactMap.keySet());
									for (Set<Node> set : iterationList) {
										if (validNodesToContactMap.get(set) > minDist*(1+UPGMA_THRESHOD)) {
											validNodesToContactMap.remove(set);
										}
									}
									nodesToConcat.clear();
									nodesToConcat.add(node1);
									nodesToConcat.add(node2);
								}
								Set<Node> validContact = new HashSet<Node>();
								validContact.add(node1);
								validContact.add(node2);
								validNodesToContactMap.put(validContact, distance);
							}
						}
					}
				}
				findBestPairs(validNodesToContactMap, clusters, i);
//				int index = 1;
//				while(validNodesToContactMap.size() > 1 && (i+index < positions.size() || i-index >= 0)){
//					List<Set<Node>> iterationList = new ArrayList<Set<Node>>();
//					iterationList.addAll(validNodesToContactMap.keySet());
//					for (Set<Node> set : iterationList) {
//						Set<Integer> combinedSplit = new HashSet<Integer>();
//						for (Node node : set) {
//							combinedSplit.addAll(clusters.get(node));
//						}
//						List<Node> nodes = new ArrayList<Node>();
//						nodes.addAll(set);
//						if (i-index >= 0 && !isCompatible(originalSplits.get(i-index), combinedSplit)) {
//							validNodesToContactMap.remove(set);
//						}else if (i+index <= positions.size()-1 && !isCompatible(originalSplits.get(i+index), combinedSplit)) {
//							validNodesToContactMap.remove(set);
//						}
//					}
//					index++;
//				}
				if (validNodesToContactMap.size() == 1) {
					nodesToConcat = validNodesToContactMap.keySet().iterator().next();
				}
				Node parent = new Node();
				parent.getChildren().addAll(nodesToConcat);
				List<Integer> coveringTaxa = new ArrayList<Integer>();
				parent.setInfo(String.valueOf(minDist/2));
				for (Node node : nodesToConcat) {
					node.setParent(parent);
					double branchLength = (minDist/2)-Double.valueOf(node.getInfo());
					node.setBranchLength(String.valueOf(branchLength));
					coveringTaxa.addAll(clusters.get(node));
					clusters.remove(node);
				}
				clusters.put(parent, coveringTaxa);
			}
			upgmaTree.setRoot(clusters.keySet().iterator().next());
			upgmaTree.updateNodes();
			upgmaTree.updateSplits();
			timeTrees[i] = upgmaTree;
			Map<Node, Double> distanceToRootMap = buildDistanceToRootMap(upgmaTree);
			timeTreeLengths[i] = distanceToRootMap.get(upgmaTree.getTaxa().get(0));
		}
		System.out.println(" [" + (double) (new Date().getTime() - startDate.getTime()) / 1000
				+ " Seconds] ");
	}

	private void findBestPairs(Map<Set<Node>, Double> validNodesToContactMap, HashMap<Node, List<Integer>> clusters, int i){
		int index = 1;
		while(validNodesToContactMap.size() > 1 && (i+index < positions.size() || i-index >= 0)){
			List<Set<Node>> iterationList = new ArrayList<Set<Node>>();
			iterationList.addAll(validNodesToContactMap.keySet());
			for (Set<Node> set : iterationList) {
				Set<Integer> combinedSplit = new HashSet<Integer>();
				for (Node node : set) {
					combinedSplit.addAll(clusters.get(node));
				}
				List<Node> nodes = new ArrayList<Node>();
				nodes.addAll(set);
				if (i-index >= 0 && !isCompatible(originalSplits.get(i-index), combinedSplit)) {
					validNodesToContactMap.remove(set);
				}else if (i+index <= positions.size()-1 && !isCompatible(originalSplits.get(i+index), combinedSplit)) {
					validNodesToContactMap.remove(set);
				}
			}
			index++;
		}
	}

	private double getDistanceOfClusters(List<Integer> list, List<Integer> list2, double[][] distanceMatrix) {
		double distance = 0.0;
		for (Integer i : list) {
			for (Integer j : list2) {
				distance += distanceMatrix[i-1][j-1];
			}
		}
		return distance/(list.size()*list2.size());
	}

	private void makeDistanceMatricesFast2() {
		System.out.print("Making Distance Matrices Fast 2");
		Date startDate = new Date();
		for (int i = 0; i < positions.size(); i++) {
			double[][] distanceMatrix = new double[matrix.length][matrix.length];
			distanceMatrices.add(distanceMatrix);
		}
		for (int j = 0; j < matrix.length-1; j++) {
			for (int j2 = j+1; j2 < matrix.length; j2++) {
				Set<Integer> split = new HashSet<Integer>();
				split.add(j+1);
				split.add(j2+1);
				List<int[]> normalRegions = new ArrayList<int[]>();
				Map<Integer, int[]> specialRegions = new HashMap<Integer, int[]>();
				int startingIndex = 0;
				int previousIndex = -1;
				for (int k = 0; k < positions.size(); k++) {
					Set<Integer> originalSplit = originalSplits.get(k);
					if (originalSplit.contains(j+1) ^ originalSplit.contains(j2+1)) {
						if (startingIndex != 0 || k!= 0) {
							if (k > startingIndex+1) {
								int[] region = new int[]{startingIndex,k};
								normalRegions.add(region);
							}
							if (previousIndex != -1) {
								int[] specialRegion = new int[]{previousIndex,k};
								specialRegions.put(startingIndex, specialRegion);
							}
							if (k == positions.size()-1) {
								int[] specialRegion = new int[]{startingIndex,k+1};
								specialRegions.put(k, specialRegion);
							}
							previousIndex = startingIndex;
							startingIndex = k;
						}else{
							previousIndex = 0;
						}
					}
				}
				if (startingIndex < positions.size()-1) {
					int[] region = new int[]{startingIndex,positions.size()-1};
					normalRegions.add(region);
				}
				for (int[] region1 : normalRegions) {
					double score = getSplitScoreForRegion2(split, region1);
					for (int i = region1[0]; i <= region1[1]; i++) {
						distanceMatrices.get(i)[j][j2] = score;
						distanceMatrices.get(i)[j2][j] = score;
					}
				}
				for (int i : specialRegions.keySet()) {
					double score = getSplitScoreForRegion2(split, specialRegions.get(i));
					distanceMatrices.get(i)[j][j2] = score;
					distanceMatrices.get(i)[j2][j] = score;
				}
			}
		}
		System.out.println(" [" + (double) (new Date().getTime() - startDate.getTime()) / 1000
				+ " Seconds] ");
	}
	
	private List<Integer> getPositionArrayFromMSFile(String url) {
		List<Integer> positions = new ArrayList<Integer>();
		try {
			FileInputStream fis = new FileInputStream(url);
			InputStreamReader reader = new InputStreamReader(fis);
			BufferedReader br = new BufferedReader(reader);
			String line = br.readLine();
			String recombLengthString = "1";
			if (line.contains("-r")) {
				recombLengthString = line.split("-r")[1].split(" ")[2];
			}
//			else{
//				recombLengthString = line;
//			}
			double recombLength = Double.valueOf(recombLengthString);
			TOTAL_SEQ_LENGTH = recombLength;
//			line = br.readLine();
			while (line != null) {
				if (line.contains("positions")) {
					String[] positionStrings = line.split(":")[1].substring(1)
							.split(" ");
					for (int i = 0; i < positionStrings.length; i++) {
						double position = Double.valueOf(positionStrings[i])*recombLength;
						positions.add((int) position);
					}
				}
				line = br.readLine();
			}
			if (positions.isEmpty()) {
				int divider = (int)recombLength/(localTrees.length+1);
				for (int i = 0; i < localTrees.length; i++) {
					positions.add((i+1)*divider);
				}
			}else{
				List<Integer> temPositions = new ArrayList<Integer>();
				for (int i = 0; i < positions.size(); i++) {
					if (!indicesToRemove.contains(i)) {
						temPositions.add(positions.get(i));
					}
				}
				positions.clear();
				positions.addAll(temPositions);
			}
			br.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return positions;
	}
	
	private List<Integer> getPositionArrayFromInputFile(String url) {
		List<Integer> positions = new ArrayList<Integer>();
		try {
			FileInputStream fis = new FileInputStream(url);
			InputStreamReader reader = new InputStreamReader(fis);
			BufferedReader br = new BufferedReader(reader);
			String line = br.readLine();
			if (line.matches(".*[2-9].*") || line.contains(".")) {
				String[] positionStrings = line.split(" ");
				int maxLength = Integer.MIN_VALUE;
				if (positionStrings[0].contains(".")) {
					for (int i = 0; i < positionStrings.length; i++) {
						if (positionStrings[i].length() > maxLength) {
							maxLength = positionStrings[i].length();
						}
					}
					for (int i = 0; i < positionStrings.length; i++) {
						double pos = Double.valueOf(positionStrings[i])*Math.pow(10,(maxLength-2));
						positions.add((int) pos);
					}
				}else{
					for (int i = 0; i < positionStrings.length; i++) {
						positions.add(Integer.valueOf(positionStrings[i]));
					}
				}
				TOTAL_SEQ_LENGTH = positions.get(positions.size()-1)-positions.get(0)+1;
				List<Integer> temPositions = new ArrayList<Integer>();
				for (int i = 0; i < positions.size(); i++) {
					if (!indicesToRemove.contains(i)) {
						temPositions.add(positions.get(i));
					}
				}
				positions.clear();
				positions.addAll(temPositions);
			}else {
				for (int i = 0; i < localTrees.length; i++) {
					positions.add((i+1)*10);
				}
				TOTAL_SEQ_LENGTH = positions.get(localTrees.length*10);
			}
			br.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return positions;
	}

	private void scaleTimes(double[] times) {
		double avg = 0;
		for (double d : times) {
			avg += d;
		}
		avg = avg/times.length;
		double multiplier = 2/avg;
		for (int i = 0; i < times.length; i++) {
			times[i] = times[i]*multiplier;
		}
	}

	private void compareAllTrees(List<Integer> positions) {
		int numTotalSplits = 0;
		int numTotalInfered = 0;
		int numInferedSplits = 0;
		int numTimeSplits = 0;
		int numArgSplits = 0;
		int numMargSplits = 0;
		int numSharedSplits = 0;
		int counter = 0;
		int notShareCounter = 0;
		for (int i = 0; i < localTrees.length; i++) {
			Tree localTree = localTrees[i];
			Tree timeTree = timeTrees[i];
			Tree trueTree = trueTrees[i];
//			Tree margTree = margTrees[i];
//			Tree argTree = argTrees[i];
			Set<Set<Integer>> realTrueSplits = new HashSet<Set<Integer>>();
			Set<Set<Integer>> realTimeSplits = new HashSet<Set<Integer>>();
			Set<Set<Integer>> realLocalSplits = new HashSet<Set<Integer>>();
			Set<Set<Integer>> realArgSplits = new HashSet<Set<Integer>>();
			Set<Set<Integer>> realMargSplits = new HashSet<Set<Integer>>();
			for (Set<Integer> split : trueTree.getSplits()) {
				if (1 < split.size()
						&& split.size() < trueTree.getRoot().getSplit().size() - 1
						&& !realTrueSplits.contains(reverseSplit(split))) {
					realTrueSplits.add(split);
				}
			}
			numTotalSplits += realTrueSplits.size();
			for (Set<Integer> split : localTree.getSplits()) {
				if (1 < split.size()
						&& split.size() < overalSplit.size() - 1
						&& !realLocalSplits.contains(reverseSplit(split))) {
					realLocalSplits.add(split);
				}
			}
			numTotalInfered += realLocalSplits.size();
			
			for (Set<Integer> split : timeTree.getSplits()) {
				if (1 < split.size()
						&& split.size() < timeTree.getRoot().getSplit().size() - 1
						&& !realTimeSplits.contains(reverseSplit(split))) {
					realTimeSplits.add(split);
				}
			}
			
//			for (Set<Integer> split : argTree.getSplits()) {
//				if (1 < split.size()
//						&& split.size() < argTree.getRoot().getSplit().size() - 1
//						&& !realArgSplits.contains(reverseSplit(split))) {
//					realArgSplits.add(split);
//				}
//			}
			
//			for (Set<Integer> split : margTree.getSplits()) {
//				if (1 < split.size()
//						&& split.size() < margTree.getRoot().getSplit().size() - 1
//						&& !realMargSplits.contains(reverseSplit(split))) {
//					realMargSplits.add(split);
//				}
//			}
			
			Set<Set<Integer>> tmpLocalSplit = new HashSet<Set<Integer>>();
			Set<Set<Integer>> tmpTimeSplit = new HashSet<Set<Integer>>();
			Set<Set<Integer>> tmpArgSplit = new HashSet<Set<Integer>>();
			Set<Set<Integer>> tmpMargSplit = new HashSet<Set<Integer>>();
			for (Set<Integer> set : realLocalSplits) {
				if (realTrueSplits.contains(set)
						|| realTrueSplits.contains(reverseSplit(set))) {
					tmpLocalSplit.add(set);
				}
			}
			for (Set<Integer> set : realTimeSplits) {
				if (realTrueSplits.contains(set)
						|| realTrueSplits.contains(reverseSplit(set))) {
					tmpTimeSplit.add(set);
				}
			}
			for (Set<Integer> set : realArgSplits) {
				if (realTrueSplits.contains(set)
						|| realTrueSplits.contains(reverseSplit(set))) {
					tmpArgSplit.add(set);
				}
			}
			for (Set<Integer> set : realMargSplits) {
				if (realTrueSplits.contains(set)
						|| realTrueSplits.contains(reverseSplit(set))) {
					tmpMargSplit.add(set);
				}
			}
			Set<Set<Integer>> sharedSplits = new HashSet<Set<Integer>>();
			for (Set<Integer> set : tmpTimeSplit) {
				if (tmpLocalSplit.contains(set)
						|| tmpLocalSplit.contains(reverseSplit(set))) {
					sharedSplits.add(set);
				}
			}
			numInferedSplits += tmpLocalSplit.size();
			numTimeSplits += tmpTimeSplit.size();
			numArgSplits += tmpArgSplit.size();
			numMargSplits += tmpMargSplit.size();
			numSharedSplits += sharedSplits.size();
//			String split = originalSplits.get(i).isEmpty()?getSingleMutation(i):originalSplits.get(i).toString();
//			System.out.println("At index " + (i) + ", " + tmpLocalSplit.size()
//					+ " (Inf), and " + tmpTimeSplit.size() + " (time)" + " of "
//					+ realTrueSplits.size()
//					+ " total splits are inferred.  Shared: "+ sharedSplits.size() 
//					+ " ---  Inferred Tree: " + localTree.toString()+ " [" + localTreeLengths[i] + "]"
//					+ " UPGMA Tree: " + timeTree.toString() + " [" + timeTreeLengths[i]+ "]"
////					+ " ArgWeaver Tree: " + argTree.toString() + " [" + argTreeLengths[i]+ "]"
//					+ "  True Tree: " + trueTree + " ["	+ trueTreeLengths[i] + "]"
//					 + " Best Compatible Region: [" + bestRegions[i][0] + ","
//					 + bestRegions[i][1] + "]"
//					+ " {" + positions.get(i) + "}"
//					+ " Original Split: " + split
//					);
			if (tmpLocalSplit.size() < tmpTimeSplit.size()) {
				counter += tmpTimeSplit.size() - tmpLocalSplit.size();
			}
			if (sharedSplits.size() < tmpTimeSplit.size()
					&& sharedSplits.size() < tmpLocalSplit.size()) {
				notShareCounter += tmpTimeSplit.size() <= sharedSplits.size() ? tmpTimeSplit
						.size() - sharedSplits.size()
						: tmpLocalSplit.size() - sharedSplits.size();
			}
		}
		System.out.println("Inferred: " + numInferedSplits + " out of "
				+ numTotalSplits + " splits = %" + (double) numInferedSplits
				* 100.0 / (double) numTotalSplits);
		System.out.println("Time: " + numTimeSplits + " out of "
				+ numTotalSplits + " splits = %" + (double) numTimeSplits
				* 100.0 / (double) numTotalSplits);
//		System.out.println("ArgWeaver: " + numArgSplits + " out of "
//				+ numTotalSplits + " splits = %" + (double) numArgSplits
//				* 100.0 / (double) numTotalSplits);
//		System.out.println("Margarita: " + numMargSplits + " out of "
//				+ numTotalSplits + " splits = %" + (double) numMargSplits
//				* 100.0 / (double) numTotalSplits);
//		System.out.println("Shared: " + numSharedSplits + " out of "
//				+ numTotalInfered + " splits = %" + (double) numSharedSplits
//				* 100.0 / (double) numTotalSplits);
//		System.out.println("Total Potential Fixable Indices: " + counter
//				+ " = %" + (double) counter * 100.0 / (double) numTotalSplits);
//		System.out.println("Total Potential Extra Fixable Indices: "
//				+ notShareCounter + " = %" + (double) notShareCounter * 100.0
//				/ (double) numTotalSplits);
	}
	
	private void compareAllTreesRooted(List<Integer> positions) {
		int numTotalClades = 0;
		int numTotalInfered = 0;
		int numInferedSplits = 0;
		int numTimeSplits = 0;
		int numSharedSplits = 0;
		int counter = 0;
		int notShareCounter = 0;
		for (int i = 0; i < localTrees.length; i++) {
			Tree localTree = localTrees[i];
			Tree timeTree = timeTrees[i];
			Tree trueTree = trueTrees[i];
			Set<Set<Integer>> realTrueClades = new HashSet<Set<Integer>>();
			Set<Set<Integer>> realTimeClades = new HashSet<Set<Integer>>();
			Set<Set<Integer>> realLocalClades = new HashSet<Set<Integer>>();
			for (Set<Integer> clade : trueTree.getClades()) {
				if (1 < clade.size()) {
					realTrueClades.add(clade);
				}
			}
			numTotalClades += realTrueClades.size();
			for (Set<Integer> clade : localTree.getClades()) {
				if (1 < clade.size()) {
					realLocalClades.add(clade);
				}
			}
			numTotalInfered += realLocalClades.size();
			
			for (Set<Integer> clade : timeTree.getClades()) {
				if (1 < clade.size()) {
					realTimeClades.add(clade);
				}
			}
			
			Set<Set<Integer>> tmpLocalClades = new HashSet<Set<Integer>>();
			Set<Set<Integer>> tmpTimeClades = new HashSet<Set<Integer>>();
			
			for (Set<Integer> set : realLocalClades) {
				if (realTrueClades.contains(set)) {
					tmpLocalClades.add(set);
				}
			}
			for (Set<Integer> set : realTimeClades) {
				if (realTrueClades.contains(set)) {
					tmpTimeClades.add(set);
				}
			}
			Set<Set<Integer>> sharedClades = new HashSet<Set<Integer>>();
			for (Set<Integer> set : tmpTimeClades) {
				if (tmpLocalClades.contains(set)) {
					sharedClades.add(set);
				}
			}
			numInferedSplits += tmpLocalClades.size();
			numTimeSplits += tmpTimeClades.size();
			numSharedSplits += sharedClades.size();
			String split = originalSplits.get(i).isEmpty()?getSingleMutation(i):originalSplits.get(i).toString();
			System.out.println("At index " + (i) + ", " + tmpLocalClades.size()
					+ " (Inf), and " + tmpTimeClades.size() + " (time)" + " of "
					+ realTrueClades.size()
					+ " total clades are inferred.  Shared: "+ sharedClades.size() 
					+ " ---  Inferred Tree: " + localTree.toString()+ " [" + localTreeLengths[i] + "]"
					+ " UPGMA Tree: " + timeTree.toString() + " [" + timeTreeLengths[i]+ "]"
					+ "  True Tree: " + trueTree + " ["	+ trueTreeLengths[i] + "]"
					+ " {" + positions.get(i) + "}"
					+ " Original Split: " + split
					);
			if (tmpLocalClades.size() < tmpTimeClades.size()) {
				counter += tmpTimeClades.size() - tmpLocalClades.size();
			}
			if (sharedClades.size() < tmpTimeClades.size()
					&& sharedClades.size() < tmpLocalClades.size()) {
				notShareCounter += tmpTimeClades.size() <= sharedClades.size() ? tmpTimeClades
						.size() - sharedClades.size()
						: tmpLocalClades.size() - sharedClades.size();
			}
		}
		System.out.println("Inf: " + numInferedSplits + " out of "
				+ numTotalClades + " clades = %" + (double) numInferedSplits
				* 100.0 / (double) numTotalClades);
		System.out.println("Time: " + numTimeSplits + " out of "
				+ numTotalClades + " clades = %" + (double) numTimeSplits
				* 100.0 / (double) numTotalClades);
//		System.out.println("Shared: " + numSharedSplits + " out of "
//				+ numTotalInfered + " clades = %" + (double) numSharedSplits
//				* 100.0 / (double) numTotalClades);
//		System.out.println("Total Potential Fixable clades: " + counter
//				+ " = %" + (double) counter * 100.0 / (double) numTotalClades);
//		System.out.println("Total Potential Extra Fixable clades: "
//				+ notShareCounter + " = %" + (double) notShareCounter * 100.0
//				/ (double) numTotalClades);
	}
	
	private void outputInferredTrees(List<Integer> positions) {
		for (int i = 0; i < localTrees.length; i++) {
			Tree localTree = localTrees[i];
			Tree timeTree = timeTrees[i];
			String split = originalSplits.get(i).isEmpty()?getSingleMutation(i):originalSplits.get(i).toString();
			System.out.println("At index " + (i) + ", " + " Inferred Tree: " + localTree.toString()+ " [" + localTreeLengths[i] + "]"
					+ " Time Tree: " + timeTree.toString() + " [" + timeTreeLengths[i]+ "]"
					 + " Best Compatible Region: [" + bestRegions[i][0] + ","
					 + bestRegions[i][1] + "]"
					+ " {" + positions.get(i) + "}"
					+ " Original Split: " + split
					);
		}
	}
	
	private void outputTreesSeperate(Tree[] trees, String inputUrl, String level) {
		try {
			BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(
			          new FileOutputStream(inputUrl + ".trees" + level)));
			for (int i = 0; i < trees.length; i++) {
				Tree localTree = trees[i];
				String position = String.valueOf(positions.get(i));
				writer.write(position + "\t" + localTree.toString());
				writer.newLine();
				localTree.setUpdateFlag(true);
			}
			writer.flush();
			writer.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (UnsupportedEncodingException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private String getSingleMutation(int i) {
		for (int j = 0; j < matrix.length; j++) {
			if (matrix[j][i] != root[i]) {
				return String.valueOf(j+1);
			}
		}
		return null;
	}

	private Set<Integer> reverseSplit(Set<Integer> split) {
		Set<Integer> reverseSplit = new HashSet<Integer>();
		for (Integer integer : overalSplit) {
			if (!split.contains(integer)) {
				reverseSplit.add(integer);
			}
		}
		return reverseSplit;
	}

	private void buildDataFromHapMatrix(String url) {
		try {
			// reading input file saving as samples
			FileInputStream fis = new FileInputStream(url);
			InputStreamReader reader = new InputStreamReader(fis);
			BufferedReader br = new BufferedReader(reader);
			String line = br.readLine();
			List<String> rows = new ArrayList<String>();
			line = br.readLine();
			while (line != null && !line.isEmpty()) {
				rows.add(line);
				line = br.readLine();
			}
			br.close();
			Set<Integer> indicesToRemove = new HashSet<Integer>();
			for (int i = 0; i < rows.get(0).length(); i++) {
				char first = rows.get(0).charAt(i);
				boolean hasInfo = false;
				for (int j = 1; j < rows.size(); j++) {
					if (first != rows.get(j).charAt(i)) {
						hasInfo = true;
					}
				}
				if (!hasInfo) {
					indicesToRemove.add(i);
				}
			}
			this.indicesToRemove = indicesToRemove;
			List<StringBuilder> builders = new ArrayList<StringBuilder>();
			for (int i = 0; i < rows.size(); i++) {
				builders.add(new StringBuilder());
			}
			for (int i = 0; i < rows.get(0).length(); i++) {
				if (!indicesToRemove.contains(i)) {
					for (int j = 0; j < builders.size(); j++) {
						builders.get(j).append(rows.get(j).charAt(i));
					}
				}
			}
			rows.clear();
			for (int i = 0; i < builders.size(); i++) {
				rows.add(builders.get(i).toString());
			}
			matrix = new int[rows.size()][rows.get(0).length()];
			for (int i1 = 0; i1 < matrix.length; i1++) {
				for (int j1 = 0; j1 < matrix[0].length; j1++) {
					matrix[i1][j1] = Integer.valueOf(rows.get(i1).substring(j1,
							j1 + 1));
				}
			}
			root = buildRoot(matrix);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void getTrueTreesFromMsFile(String url) {
		List<Integer> positions = new ArrayList<Integer>();
		List<Tree> trees = new ArrayList<Tree>();
		try {
			FileInputStream fis = new FileInputStream(url);
			InputStreamReader reader = new InputStreamReader(fis);
			BufferedReader br = new BufferedReader(reader);
			String line = br.readLine();
			String recombLengthString = line.split("-r")[1].split(" ")[2];
			double recombLength = Double.valueOf(recombLengthString);
			line = br.readLine();
			while (line != null) {
				if (line.contains(";")) {
					String positionString = line.split("]")[0].substring(1);
					if (positions.isEmpty()) {
						positions.add(Integer.valueOf(positionString));
					} else {
						positions.add(positions.get(positions.size() - 1)
								+ Integer.valueOf(positionString));
					}
					String treeString = line.split("]")[1];
					Tree tree = Util.getTreeFromNewickWithTime(treeString);
					trees.add(tree);
				}
				if (line.contains("positions")) {
					String[] positionArray = line.split(":")[1].substring(1)
							.split(" ");
					trueTrees = new Tree[positionArray.length];
					for (int i = 0; i < positionArray.length; i++) {
						double truePosition = Double.valueOf(positionArray[i])
								* recombLength;
						int index = 0;
						for (int j = 0; j < positions.size(); j++) {
							if (truePosition < positions.get(j)) {
								index = j;
								break;
							}
						}
						trueTrees[i] = trees.get(index);
					}
				}
				line = br.readLine();
			}
			for (Tree tree : trueTrees) {
				tree.updateSplits();
			}
			br.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	private void timeSplitRule() {
		System.out.print("Time Tree Rule");
		Date startDate = new Date();
		int counter = 0;
		for (int i = 0; i < timeTrees.length; i++) {
			Tree tree = timeTrees[i];
			for (Set<Integer> split : tree.getSplits()) {
				counter += addSplitTooTree(i, localTrees[i], split, false, false, true);
			}
		}
		System.out.print(" [" + (double) (new Date().getTime() - startDate.getTime()) / 1000
				+ " Seconds] ");
		System.out.println("(" + counter + ")");
	}

	private void oneSPREventRule(boolean time) {
		System.out.print("One SPR Rule ");
		Date startDate = new Date();
		int counter = 0;
		for (int i = 0; i < localTrees.length - 1; i++) {
			Tree tree1 = localTrees[i];
			Tree tree2 = localTrees[i + 1];
			Set<Set<Integer>> oldSplits1 = new HashSet<Set<Integer>>();
			oldSplits1.addAll(tree1.getSplits());
			Set<Set<Integer>> oldSplits2 = new HashSet<Set<Integer>>();
			oldSplits2.addAll(tree2.getSplits());
			Set<Set<Integer>> toRemoveSplits = new HashSet<Set<Integer>>();

			Set<Set<Integer>> tree1Splits = new HashSet<Set<Integer>>();
			tree1Splits.addAll(tree1.getSplits());
			Set<Set<Integer>> tree2Splits = new HashSet<Set<Integer>>();
			tree2Splits.addAll(tree2.getSplits());
			tree1Splits.removeAll(tree2.getSplits());
			tree2Splits.removeAll(tree1.getSplits());
			for (Set<Integer> split1 : tree1Splits) {
				for (Set<Integer> split2 : tree2Splits) {
					if (!isCompatible(split1, split2)) {
						Set<Integer> tmpSplit = new HashSet<Integer>();
						tmpSplit.addAll(split1);
						tmpSplit.retainAll(split2);
						toRemoveSplits.add(tmpSplit);
					}
				}
			}
			for (Set<Integer> split : toRemoveSplits) {
				Set<Integer> reverseSplit = reverseSplit(split);
				Node nodeInTree1 = null;
				Node nodeInTree2 = null;
				if (tree1.getSplits().contains(split)
						&& tree2.getSplits().contains(split)) {
					nodeInTree1 = tree1.getExactNodeForSplit(split);
					nodeInTree2 = tree2.getExactNodeForSplit(split);
				} else if (tree1.getSplits().contains(reverseSplit)
						&& tree2.getSplits().contains(reverseSplit)) {
					nodeInTree1 = tree1.getExactNodeForSplit(reverseSplit);
					nodeInTree2 = tree2.getExactNodeForSplit(reverseSplit);
				}
				if (nodeInTree1 != null && nodeInTree2 != null) {
					// separating the subtree
					Tree subTree1 = new Tree();
					subTree1.setRoot(nodeInTree1);
					subTree1.update();
					Tree subTree2 = new Tree();
					subTree2.setRoot(nodeInTree2);
					subTree2.update();

					Node parent1 = nodeInTree1.getParent();
					Node parent2 = nodeInTree2.getParent();
					parent1.getChildren().remove(nodeInTree1);
					parent2.getChildren().remove(nodeInTree2);

					nodeInTree1.setParent(null);
					nodeInTree2.setParent(null);

					tree1.update();
					tree2.update();
					subTree1.update();
					subTree2.update();
					boolean breakBool = false;
					// combining splits
					if (isCompatible(tree1, tree2)
							&& isCompatible(subTree1, subTree2)) {
						combineSplits(i, tree1, tree2, time);
						combineSplits(i, subTree1, subTree2, time);
						breakBool = true;
					}

					// building the original tree
					nodeInTree1.setParent(parent1);
					nodeInTree2.setParent(parent2);

					parent1.getChildren().add(nodeInTree1);
					parent2.getChildren().add(nodeInTree2);

					tree1.update();
					tree2.update();
					if (breakBool) {
						Set<Set<Integer>> newSplits1 = new HashSet<Set<Integer>>();
						newSplits1.addAll(tree1.getSplits());
						Set<Set<Integer>> newSplits2 = new HashSet<Set<Integer>>();
						newSplits2.addAll(tree2.getSplits());
						newSplits1.removeAll(oldSplits1);
						newSplits2.removeAll(oldSplits2);
						for (Set<Integer> set : newSplits1) {
							counter++;
//							counter += propagateSplit(i, set);
						}
						for (Set<Integer> set : newSplits2) {
							counter++;
//							counter += propagateSplit(i + 1, set);
						}
						break;
					}
				}
			}
		}
		System.out.print(" [" + (double) (new Date().getTime() - startDate.getTime()) / 1000
				+ " Seconds] ");
		System.out.println("(" + counter + ")");
	}

	private void twoSPREventRule(boolean time) {
		System.out.print("Two SPR Rule ");
		Date startDate = new Date();
		int counter = 0;
		for (int i = 0; i < localTrees.length - 1; i++) {
			Tree tree1 = localTrees[i];
			Tree tree2 = localTrees[i + 1];
			Set<Set<Integer>> oldSplits1 = new HashSet<Set<Integer>>();
			oldSplits1.addAll(tree1.getSplits());
			Set<Set<Integer>> oldSplits2 = new HashSet<Set<Integer>>();
			oldSplits2.addAll(tree2.getSplits());
			Set<Set<Integer>> toRemoveSplits = new HashSet<Set<Integer>>();

			Set<Set<Integer>> tree1Splits = new HashSet<Set<Integer>>();
			tree1Splits.addAll(tree1.getSplits());
			Set<Set<Integer>> tree2Splits = new HashSet<Set<Integer>>();
			tree2Splits.addAll(tree2.getSplits());
			tree1Splits.removeAll(tree2.getSplits());
			tree2Splits.removeAll(tree1.getSplits());
			for (Set<Integer> split1 : tree1Splits) {
				for (Set<Integer> split2 : tree2Splits) {
					if (!isCompatible(split1, split2)) {
						Set<Integer> tmpSplit = new HashSet<Integer>();
						tmpSplit.addAll(split1);
						tmpSplit.retainAll(split2);
						toRemoveSplits.add(tmpSplit);
					}
				}
			}
			List<Set<Integer>> toRemoveSplitList = new ArrayList<Set<Integer>>(
					toRemoveSplits);
			for (int j = 0; j < toRemoveSplitList.size() - 1; j++) {
				for (int j2 = j + 1; j2 < toRemoveSplitList.size(); j2++) {
					Set<Integer> split1 = toRemoveSplitList.get(j);
					Set<Integer> split2 = toRemoveSplitList.get(j2);

					Set<Integer> reverseSplit1 = reverseSplit(split1);
					Set<Integer> reverseSplit2 = reverseSplit(split2);

					Node node1InTree1 = null;
					Node node1InTree2 = null;
					Node node2InTree1 = null;
					Node node2InTree2 = null;

					if (tree1.getSplits().contains(split1)
							&& tree2.getSplits().contains(split1)
							&& tree1.getSplits().contains(split2)
							&& tree2.getSplits().contains(split2)) {
						node1InTree1 = tree1.getExactNodeForSplit(split1);
						node1InTree2 = tree2.getExactNodeForSplit(split1);
						node2InTree1 = tree1.getExactNodeForSplit(split2);
						node2InTree2 = tree2.getExactNodeForSplit(split2);
					} else if (tree1.getSplits().contains(reverseSplit1)
							&& tree2.getSplits().contains(reverseSplit1)
							&& tree1.getSplits().contains(reverseSplit2)
							&& tree2.getSplits().contains(reverseSplit2)) {
						node1InTree1 = tree1.getExactNodeForSplit(reverseSplit1);
						node1InTree2 = tree2.getExactNodeForSplit(reverseSplit1);
						node2InTree1 = tree1.getExactNodeForSplit(reverseSplit2);
						node2InTree2 = tree2.getExactNodeForSplit(reverseSplit2);
					}
					if (node1InTree1 != null && node1InTree2 != null
							&& node2InTree1 != null && node2InTree2 != null) {
						// separating the subtree
						Tree subTree11 = new Tree();
						subTree11.setRoot(node1InTree1);
						subTree11.update();
						Tree subTree12 = new Tree();
						subTree12.setRoot(node1InTree2);
						subTree12.update();
						Tree subTree21 = new Tree();
						subTree21.setRoot(node2InTree1);
						subTree21.update();
						Tree subTree22 = new Tree();
						subTree22.setRoot(node2InTree2);
						subTree22.update();

						Node parent11 = node1InTree1.getParent();
						Node parent12 = node1InTree2.getParent();
						parent11.getChildren().remove(node1InTree1);
						parent12.getChildren().remove(node1InTree2);
						Node parent21 = node2InTree1.getParent();
						Node parent22 = node2InTree2.getParent();
						parent21.getChildren().remove(node2InTree1);
						parent22.getChildren().remove(node2InTree2);

						Node grandParent1 = null;
						boolean parentRemoval1 = false;
						if (parent11.equals(parent21)) {
							grandParent1 = parent11.getParent();
							if (parent11.getChildren().isEmpty()) {
								grandParent1.getChildren().remove(parent11);
								parentRemoval1 = true;
							}
						}
						Node grandParent2 = null;
						boolean parentRemoval2 = false;
						if (parent12.equals(parent22)) {
							grandParent2 = parent12.getParent();
							if (parent12.getChildren().isEmpty()) {
								grandParent2.getChildren().remove(parent12);
								parentRemoval2 = true;
							}
						}

						node1InTree1.setParent(null);
						node1InTree2.setParent(null);
						node2InTree1.setParent(null);
						node2InTree2.setParent(null);

						tree1.update();
						tree2.update();

						subTree11.update();
						subTree12.update();
						subTree21.update();
						subTree22.update();

						boolean breakBool = false;
						// combining splits
						if (isCompatible(tree1, tree2)
								&& isCompatible(subTree11, subTree12)
								&& isCompatible(subTree21, subTree22)) {
							combineSplits(i, tree1, tree2, time);
							combineSplits(i, subTree11, subTree12, time);
							combineSplits(i, subTree21, subTree22, time);
							breakBool = true;
						}

						if (parentRemoval1) {
							grandParent1.getChildren().add(parent11);
						}
						if (parentRemoval2) {
							grandParent2.getChildren().add(parent12);
						}
						// building the original tree
						node1InTree1.setParent(parent11);
						node1InTree2.setParent(parent12);
						node2InTree1.setParent(parent21);
						node2InTree2.setParent(parent22);

						parent11.getChildren().add(node1InTree1);
						parent12.getChildren().add(node1InTree2);
						parent21.getChildren().add(node2InTree1);
						parent22.getChildren().add(node2InTree2);

						tree1.update();
						tree2.update();
						if (breakBool) {
							Set<Set<Integer>> newSplits1 = new HashSet<Set<Integer>>();
							newSplits1.addAll(tree1.getSplits());
							Set<Set<Integer>> newSplits2 = new HashSet<Set<Integer>>();
							newSplits2.addAll(tree2.getSplits());
							newSplits1.removeAll(oldSplits1);
							newSplits2.removeAll(oldSplits2);
							for (Set<Integer> set : newSplits1) {
								counter++;
//								counter += propagateSplit(i, set);
							}
							for (Set<Integer> set : newSplits2) {
								counter++;
//								counter += propagateSplit(i + 1, set);
							}
							break;
						}
					}
				}
			}
		}
		System.out.print(" [" + (double) (new Date().getTime() - startDate.getTime()) / 1000
				+ " Seconds] ");
		System.out.println("(" + counter + ")");
	}

	private void combineSplits(int index, Tree tree1, Tree tree2, boolean time) {
		for (Set<Integer> split : tree1.getSplits()) {
			if (!time || isCompatible(split, timeTrees[index+1])) {
				addSplitTooTree(0, tree2, split, false, true, false);
			}
		}
		for (Set<Integer> split : tree2.getSplits()) {
			if (!time || isCompatible(split, timeTrees[index+1])) {
				addSplitTooTree(0, tree1, split, false, true, false);
			}
		}
	}

	private boolean isCompatible(Set<Integer> split1, Set<Integer> split2) {
		if (split1.size() <= 1 || split2.size() <= 1
				|| split1.size() >= overalSplit.size() - 1
				|| split2.size() >= overalSplit.size() - 1
				|| split1.equals(split2)) {
			return true;
		}
		boolean[] split1Bool = new boolean[overalSplit.size()];
		boolean[] split2Bool = new boolean[overalSplit.size()];
		for (Integer i : split1) {
			split1Bool[i-1] = true;
		}
		for (Integer i : split2) {
			split2Bool[i-1] = true;
		}
		boolean gammit00 = false;
		boolean gammit01 = false;
		boolean gammit10 = false;
		boolean gammit11 = false;
		for (int i = 0; i < split1Bool.length; i++) {
			if (!split1Bool[i] && !split2Bool[i]) {
				gammit00 = true;
			}else if (!split1Bool[i] && split2Bool[i]) {
				gammit01 = true;
			}else if (split1Bool[i] && !split2Bool[i]) {
				gammit10 = true;
			}else if (split1Bool[i] && split2Bool[i]) {
				gammit11 = true;
			}
			if (gammit00 && gammit01 && gammit10 && gammit11) {
				return false;
			}
		}
		return true;
	}

	private boolean isCompatible(Tree tree1, Tree tree2) {
		Set<Set<Integer>> tree1Splits = new HashSet<Set<Integer>>();
		tree1Splits.addAll(tree1.getSplits());
		Set<Set<Integer>> tree2Splits = new HashSet<Set<Integer>>();
		tree2Splits.addAll(tree2.getSplits());
		tree1Splits.removeAll(tree2.getSplits());
		tree2Splits.removeAll(tree1.getSplits());
		for (Set<Integer> split1 : tree1Splits) {
			for (Set<Integer> split2 : tree2Splits) {
				if (!isCompatible(split1, split2)) {
					return false;
				}
			}
		}
		return true;
	}
	
	private void propagationRule() {
		System.out.print("Propagation Rule ");
		Date startDate = new Date();
		int counter = 0;
		for (int i = 0; i < localTrees.length; i++) {
			Tree localTree = localTrees[i];
			for (Set<Integer> split : localTree.getSplits()) {
				if (split.size() > 1
						&& split.size() < localTree.getRoot().getSplit().size() - 1) {
					counter += propagateSplit(i, split);
				}
			}
		}
		System.out.print(" [" + (double) (new Date().getTime() - startDate.getTime()) / 1000
				+ " Seconds] ");
		System.out.println("(" + counter + ")");
	}

	private void propagationRuleNew() {
		System.out.print("Propagation Rule Check");
		Date startDate = new Date();
		for (int i = 0; i < localTrees.length; i++) {
			Tree localTree = localTrees[i];
			for (Set<Integer> split : localTree.getSplits()) {
				if (split.size() > 1
						&& split.size() < localTree.getRoot().getSplit().size() - 1) {
					propagateSplit3(i, split);
				}
			}
		}
		System.out.println(" [" + (double) (new Date().getTime() - startDate.getTime()) / 1000
				+ " Seconds] ");
	}
	
	private void propagationHeightRule() {
		System.out.print("Propagation Rule Height");
		Date startDate = new Date();
		int counter = 0;
		for (int i = 0; i < localTrees.length; i++) {
			Tree localTree = localTrees[i];
			for (Set<Integer> split : localTree.getSplits()) {
				if (split.size() > 1
						&& split.size() < overalSplit.size() - 1) {
					counter += propagateSplitHeight(i, split);
				}
			}
		}
		System.out.print(" [" + (double) (new Date().getTime() - startDate.getTime()) / 1000
				+ " Seconds] ");
		System.out.println("(" + counter + ")");
	}
	
	private class SubsetIterator<E> {
	    private final List<E> set;
	    private final int max;
	    private int index;

	    public SubsetIterator(List<E> originalList) {
	        set = originalList;
	        max = (1 << set.size());
	        index = 0;
	    }

	    public boolean hasNext() {
	        return index < max;
	    }

	    public List<E> next() {
	        List<E> newSet = new ArrayList<E>();
	        int flag = 1;      
	        for (E element : set) {
	            if ((index & flag) != 0) {
	                newSet.add(element);
	            }
	            flag <<= 1;
	        }
	        ++index;
	        return newSet;
	    }
	}

	private void uniqueRefinementRule(boolean time) {
		System.out.print("Unique Refinement Rule");
		Date startDate = new Date();
		int counter = 0;
		for (int i = 0; i < localTrees.length; i++) {
			Tree tree = localTrees[i];
			for (Node node : tree.getNodes()) {
				if (node.getChildren().size() > 2) {
					Set<Set<Integer>> validSplits = new HashSet<Set<Integer>>();
					List<Set<Integer>> splits = new ArrayList<Set<Integer>>();
					for (Node child : node.getChildren()) {
						splits.add(child.getSplit());
					}
					SubsetIterator<Set<Integer>> it = new SubsetIterator<Set<Integer>>(splits);
					while (it.hasNext()) {
						List<Set<Integer>> list = it.next();
						if (1 < list.size() && list.size() < splits.size()) {
							Set<Integer> concatenatedSplit = new HashSet<Integer>();
							for (Set<Integer> subset : list) {
								concatenatedSplit.addAll(subset);
							}
							if (!tree.getSplits().contains(
									reverseSplit(concatenatedSplit))) {
								if (0 < i) {
									if (i < root.length - 1) {
										if (isCompatible(concatenatedSplit,
												localTrees[i - 1])
												&& isCompatible(
														concatenatedSplit,
														localTrees[i + 1])) {
											validSplits.add(concatenatedSplit);
										}
									} else {
										if (isCompatible(concatenatedSplit,
												localTrees[i - 1])) {
											validSplits.add(concatenatedSplit);
										}
									}
								} else {
									if (isCompatible(concatenatedSplit,
											localTrees[i + 1])) {
										validSplits.add(concatenatedSplit);
									}
								}
							}
						}
						if (validSplits.size() > 1) {
							break;
						}
			        }
//					for (Set<Set<Integer>> set : powerSet) {
//						if (1 < set.size() && set.size() < splits.size()) {
//							Set<Integer> concatenatedSplit = new HashSet<Integer>();
//							for (Set<Integer> subset : set) {
//								concatenatedSplit.addAll(subset);
//							}
//							if (!tree.getSplits().contains(
//									reverseSplit(concatenatedSplit))) {
//								if (0 < i) {
//									if (i < root.length - 1) {
//										if (isCompatible(concatenatedSplit,
//												localTrees[i - 1])
//												&& isCompatible(
//														concatenatedSplit,
//														localTrees[i + 1])) {
//											validSplits.add(concatenatedSplit);
//										}
//									} else {
//										if (isCompatible(concatenatedSplit,
//												localTrees[i - 1])) {
//											validSplits.add(concatenatedSplit);
//										}
//									}
//								} else {
//									if (isCompatible(concatenatedSplit,
//											localTrees[i + 1])) {
//										validSplits.add(concatenatedSplit);
//									}
//								}
//							}
//						}
//						if (validSplits.size() > 1) {
//							break;
//						}
//					}
					if (validSplits.size() == 1) 
					{
						if (!time || isCompatible(validSplits.iterator().next(), timeTrees[i])) {
							Set<Integer> validSplit = validSplits.iterator().next();
							List<Node> nodesToJoin = new ArrayList<Node>();
							for (Node child : node.getChildren()) {
								if (validSplit.containsAll(child.getSplit())) {
									nodesToJoin.add(child);
								}
							}
							Node newNode = new Node();
							newNode.setParent(node);
							newNode.setSplit(validSplit);
							newNode.getChildren().addAll(nodesToJoin);
							
							node.getChildren().removeAll(nodesToJoin);
							node.getChildren().add(newNode);
							for (Node child : nodesToJoin) {
								child.setParent(newNode);
							}
							tree.getSplits().add(validSplit);
							counter ++;
							// using propagation
//							counter += propagateSplit(i, validSplit);
						}
					}
				}
			}
			tree.update();
		}
		System.out.print(" [" + (double) (new Date().getTime() - startDate.getTime()) / 1000
				+ " Seconds] ");
		System.out.println("(" + counter + ")");
	}

	private boolean isCompatible(Set<Integer> concatenatedSplit, Tree tree) {
		boolean compatible = true;
		// checking compatibility of the split
		for (Set<Integer> split : tree.getSplits()) {
			if (!isCompatible(concatenatedSplit, split)) {
				compatible = false;
				break;
			}
		}
		return compatible;
	}

	public Set<Set<Set<Integer>>> powerSet(Set<Set<Integer>> originalSet) {
		Set<Set<Set<Integer>>> sets = new HashSet<Set<Set<Integer>>>();
		if (originalSet.isEmpty()) {
			sets.add(new HashSet<Set<Integer>>());
			return sets;
		}
		List<Set<Integer>> list = new ArrayList<Set<Integer>>(originalSet);
		Set<Integer> head = list.get(0);
		Set<Set<Integer>> rest = new HashSet<Set<Integer>>(list.subList(1,
				list.size()));
		for (Set<Set<Integer>> set : powerSet(rest)) {
			Set<Set<Integer>> newSet = new HashSet<Set<Integer>>();
			newSet.add(head);
			newSet.addAll(set);
			sets.add(newSet);
			sets.add(set);
		}
		return sets;
	}

	private void fullyCompatibleRegionRuleNew() {
		System.out.print("Fully Compatible 2");
		Date startDate = new Date();
//		initiateCompatibleMatrix();
//		List<int[]> compatibleRegions = getAllCompatibleRegions();
		List<Set<Set<Integer>>> splitsToAdd = new ArrayList<Set<Set<Integer>>>();
		for (int i = 0; i < localTrees.length; i++) {
			// if (i==186) {
			// System.err.println();
			// }
			// finding the best region that fits the index
			int[] bestRegion = new int[2];
			bestRegion[0] = i;
			bestRegion[1] = i;
			int r = 1;
			boolean compatible = true;
			while(compatible){
				if(i+r < positions.size()){
					for (int j = bestRegion[0]; j <= bestRegion[1]; j++) {
						if (!compatibleMatrix[j][i+r]) {
							compatible = false;
							break;
						}
					}
					if (compatible) {
						bestRegion[1]++;
					}
				}else{
					compatible = false;
				}
				if (i-r >= 0) {
					for (int j = bestRegion[0]; j <= bestRegion[1]; j++) {
						if (!compatibleMatrix[i-r][j]) {
							compatible = false;
							break;
						}
					}
					if (compatible) {
						bestRegion[0]--;
					}
				}else{
					compatible = false;
				}
				r++;
			}
			// finding splits to add
			Set<Set<Integer>> splitSet = new LinkedHashSet<Set<Integer>>();
			int index = 0;
			if (i - bestRegion[0] > bestRegion[1] - i) {
				index = i - bestRegion[0];
			} else {
				index = bestRegion[1] - i;
			}
			for (int j = 1; j <= index; j++) {
				Set<Set<Integer>> prevSplitSet = new LinkedHashSet<Set<Integer>>();
				if (bestRegion[0] <= i - j && i - j >= 0) {
					Set<Set<Integer>> splits = localTrees[i - j].getSplits();
					prevSplitSet.addAll(splits);
				}
				if (i + j <= bestRegion[1] && i - j < root.length) {
					Set<Set<Integer>> splits = localTrees[i + j].getSplits();
					Set<Set<Integer>> splitsToRemove = new LinkedHashSet<Set<Integer>>();
					Set<Set<Integer>> addedSplits = new LinkedHashSet<Set<Integer>>();
					for (Set<Integer> split : splits) {
						if (prevSplitSet.isEmpty()) {
							splitSet.addAll(splits);
							break;
						}
						for (Set<Integer> prevSplit : prevSplitSet) {
							if (!isCompatible(split, prevSplit)) {
								if (isCompatible(split, localTrees[i])
										&& isCompatible(split, timeTrees[i])) {
									addedSplits.add(split);
									splitsToRemove.add(prevSplit);
								}
							}
						}
					}
					prevSplitSet.removeAll(splitsToRemove);
					prevSplitSet.addAll(addedSplits);
					prevSplitSet.addAll(splits);
				}
				splitSet.addAll(prevSplitSet);
			}
			splitsToAdd.add(splitSet);
			bestRegions[i] = bestRegion;
		}
		int counter = 0;
		for (int i = 0; i < splitsToAdd.size(); i++) {
			Set<Set<Integer>> splitSet = splitsToAdd.get(i);
			for (Set<Integer> split : splitSet) {
				counter += addSplitTooTree(i, localTrees[i], split, false, true, false);
			}
		}
		System.out.print(" [" + (double) (new Date().getTime() - startDate.getTime()) / 1000
				+ " Seconds] ");
		System.out.println("(" + counter + ")");
	}
	
	private void fullyCompatibleWithCheckNew(){
		System.out.print("Fully Compatible 1");
		Date startDate = new Date();
		initiateCompatibleMatrix();
		List<Set<Set<Integer>>> splitsToAdd = new ArrayList<Set<Set<Integer>>>();
		for (int i = 0; i < localTrees.length; i++) {
			// finding the best region that fits the index
			int[] bestRegion = new int[2];
			bestRegion[0] = i;
			bestRegion[1] = i;
			int r = 1;
			boolean compatible = true;
			while(compatible){
				if(i+r < positions.size()){
					for (int j = bestRegion[0]; j <= bestRegion[1]; j++) {
						if (!compatibleMatrix[j][i+r]) {
							compatible = false;
							break;
						}
					}
					if (compatible) {
						bestRegion[1]++;
					}
				}else{
					compatible = false;
				}
				if (i-r >= 0) {
					for (int j = bestRegion[0]; j <= bestRegion[1]; j++) {
						if (!compatibleMatrix[i-r][j]) {
							compatible = false;
							break;
						}
					}
					if (compatible) {
						bestRegion[0]--;
					}
				}else{
					compatible = false;
				}
				r++;
			}
			// finding splits to add
			Set<Set<Integer>> splitSet = new LinkedHashSet<Set<Integer>>();
			int index = 0;
			if (i - bestRegion[0] > bestRegion[1] - i) {
				index = i - bestRegion[0];
			} else {
				index = bestRegion[1] - i;
			}
			for (int j = 1; j <= index; j++) {
				Set<Set<Integer>> prevSplitSet = new LinkedHashSet<Set<Integer>>();
				if (bestRegion[0] <= i - j && i - j >= 0) {
					Set<Set<Integer>> splits = localTrees[i - j].getSplits();
					prevSplitSet.addAll(splits);
				}
				if (i + j <= bestRegion[1] && i - j < root.length) {
					Set<Set<Integer>> splits = localTrees[i + j].getSplits();
					Set<Set<Integer>> splitsToRemove = new LinkedHashSet<Set<Integer>>();
					Set<Set<Integer>> addedSplits = new LinkedHashSet<Set<Integer>>();
					for (Set<Integer> split : splits) {
						if (prevSplitSet.isEmpty()) {
							splitSet.addAll(splits);
							break;
						}
						for (Set<Integer> prevSplit : prevSplitSet) {
							if (!isCompatible(split, prevSplit)) {
								if (isCompatible(split, localTrees[i])
										&& isCompatible(split, timeTrees[i])) {
									addedSplits.add(split);
									splitsToRemove.add(prevSplit);
								}
							}
						}
					}
					prevSplitSet.removeAll(splitsToRemove);
					prevSplitSet.addAll(addedSplits);
					prevSplitSet.addAll(splits);
				}
				splitSet.addAll(prevSplitSet);
			}
			splitsToAdd.add(splitSet);
			bestRegions[i] = bestRegion;
		}
		int counter = 0;
		for (int i = 0; i < splitsToAdd.size(); i++) {
			Set<Set<Integer>> splitSet = splitsToAdd.get(i);
			for (Set<Integer> split : splitSet) {
//				if (isCompatible(split, timeTrees[i])) {
				if (timeTrees[i].getSplits().contains(split)){
					counter += addSplitTooTree(i, localTrees[i], split, false, true, false);
				}
			}
		}
//		for (int i = 0; i < splitsToAdd.size(); i++) {
//			Set<Set<Integer>> splitSet = splitsToAdd.get(i);
//			for (Set<Integer> split : splitSet) {
//				counter += addSplitTooTree(i, localTrees[i], split, false, true, false);
//			}
//		}
		System.out.print(" [" + (double) (new Date().getTime() - startDate.getTime()) / 1000
				+ " Seconds] ");
		System.out.println("(" + counter + ")");
	}

	private int addSplitTooTree(int index, Tree tree, Set<Integer> addedSplit,
			boolean propagate, boolean eligibilityCheck, boolean normalPropagation) {
		int counter = 0;
		if (!tree.getSplits().contains(addedSplit)
				&& !tree.getSplits().contains(reverseSplit(addedSplit))) {
			Set<Integer> bestSplit = new HashSet<Integer>();
			bestSplit = tree.getRoot().getSplit();
			boolean compatible = true;
			// checking compatibility of the split
			for (Set<Integer> split : tree.getSplits()) {
				if (split.size() < addedSplit.size()) {
					List<Integer> copyList = new ArrayList<Integer>();
					copyList.addAll(split);
					copyList.retainAll(addedSplit);
					if (copyList.size() > 0 && copyList.size() < split.size()) {
						compatible = false;
						break;
					}
				} else {
					List<Integer> copyList = new ArrayList<Integer>();
					copyList.addAll(split);
					copyList.retainAll(addedSplit);
					if (0 < copyList.size()
							&& copyList.size() < addedSplit.size()) {
						compatible = false;
						break;
					} else if (copyList.size() == addedSplit.size()) {
						if (split.size() < bestSplit.size()) {
							bestSplit = split;
						}
					}
				}
			}
			if (compatible) {
				boolean valid = true;
				if (eligibilityCheck) {
					valid = isTimeEligible(index, tree, addedSplit);
				}
				if (valid) {
					Node splitNode = tree.getExactNodeForSplit(bestSplit);
					Node newNode = new Node();
					newNode.setParent(splitNode);
					List<Node> nodesToRemove = new ArrayList<Node>();
					for (Node child : splitNode.getChildren()) {
						if (addedSplit.containsAll(child.getSplit()))
							nodesToRemove.add(child);
					}
					splitNode.getChildren().removeAll(nodesToRemove);
					splitNode.getChildren().add(newNode);
					newNode.getChildren().addAll(nodesToRemove);
					for (Node node : nodesToRemove) {
						node.setParent(newNode);
						newNode.getSplit().addAll(node.getSplit());
					}
					tree.getNodes().add(newNode);
					tree.getSplits().add(newNode.getSplit());
					counter++;
					if (propagate) {
						if (normalPropagation) {
							counter += propagateSplit(index, addedSplit);
						}else{
							counter += propagateSplitHeight(index, addedSplit);
						}
					}
				}
			}
		}
		return counter;
	}
	
	private int addSplitTooTreeOnly(int index, Tree tree, Set<Integer> addedSplit) {
		int counter = 0;
		if (!tree.getSplits().contains(addedSplit)
				&& !tree.getSplits().contains(reverseSplit(addedSplit))) {
			Set<Integer> bestSplit = new HashSet<Integer>();
			bestSplit = tree.getRoot().getSplit();
			boolean compatible = true;
			// checking compatibility of the split
			for (Set<Integer> split : tree.getSplits()) {
				if (split.size() < addedSplit.size()) {
					List<Integer> copyList = new ArrayList<Integer>();
					copyList.addAll(split);
					copyList.retainAll(addedSplit);
					if (copyList.size() > 0 && copyList.size() < split.size()) {
						compatible = false;
						break;
					}
				} else {
					List<Integer> copyList = new ArrayList<Integer>();
					copyList.addAll(split);
					copyList.retainAll(addedSplit);
					if (0 < copyList.size()
							&& copyList.size() < addedSplit.size()) {
						compatible = false;
						break;
					} else if (copyList.size() == addedSplit.size()) {
						if (split.size() < bestSplit.size()) {
							bestSplit = split;
						}
					}
				}
			}
			if (compatible) {
				Node splitNode = tree.getExactNodeForSplit(bestSplit);
				Node newNode = new Node();
				newNode.setParent(splitNode);
				List<Node> nodesToRemove = new ArrayList<Node>();
				for (Node child : splitNode.getChildren()) {
					if (addedSplit.containsAll(child.getSplit()))
						nodesToRemove.add(child);
				}
				splitNode.getChildren().removeAll(nodesToRemove);
				splitNode.getChildren().add(newNode);
				newNode.getChildren().addAll(nodesToRemove);
				for (Node node : nodesToRemove) {
					node.setParent(newNode);
					newNode.getSplit().addAll(node.getSplit());
				}
				tree.getNodes().add(newNode);
				tree.getSplits().add(newNode.getSplit());
				counter++;
			}else{
				return -1;
			}
		}else{
			return -1;
		}
		return counter;
	}

	private boolean isTimeEligible(int index, Tree tree, Set<Integer> addedSplit) {
		// double mrcaTime = 0;
		// List<Integer> split = new ArrayList<>();
		// split.addAll(addedSplit);
		// for (int i = 0; i < split.size()-1; i++) {
		// for (int j = i+1; j < split.size(); j++) {
		// double time = distanceMatrices.get(index)[i][j];
		// if (time > mrcaTime) {
		// mrcaTime = time;
		// }
		// }
		// }
		// if (mrcaTime/2 > treeLengths[index]*0.8) {
		// return false;
		// }
		return true;
	}

	private void getDistanceMatricesForTrueTrees() {
		for (int i = 0; i < trueTrees.length; i++) {
			Node taxa = new Node();
			for (Node node : trueTrees[i].getNodes()) {
				if (node.isLeaf()) {
					taxa = node;
				}
			}
			Map<Node, Double> distanceToRootMap = new HashMap<Node, Double>();
			calculatedDistance(distanceToRootMap, taxa);
			trueTreeLengths[i] = distanceToRootMap.get(taxa);
		}
	}

	public Node LCA(Node root, Node a, Node b) {

		if (root == null) {
			return null;
		}

		// If the root is one of a or b, then it is the LCA
		if (root == a || root == b) {
			return root;
		}

		List<Node> childrenCheck = new ArrayList<Node>();
		for (Node child : root.getChildren()) {
			childrenCheck.add(LCA(child, a, b));
		}

		List<Node> validNodes = new ArrayList<Node>();
		for (Node node : childrenCheck) {
			if (node != null) {
				validNodes.add(node);
			}
		}
		if (validNodes.size() > 1) {
			return root;
		} else if (validNodes.size() == 0) {
			return null;
		} else {

			return validNodes.get(0);
		}
	}

	private Map<Node, Double> buildDistanceToRootMap(Tree tree) {
		Map<Node, Double> distanceToRootMap = new HashMap<Node, Double>();
		for (Node node : tree.getNodes()) {
			calculatedDistance(distanceToRootMap, node);
		}
		return distanceToRootMap;
	}

	private void calculatedDistance(Map<Node, Double> distanceToRootMap,
			Node node) {
		if (distanceToRootMap.get(node) != null)
			return;
		if (distanceToRootMap.get(node.getParent()) == null) {
			if (node.getParent() == null) {
				// node is root
				distanceToRootMap.put(node, 0.0);
				return;
			} else {
				calculatedDistance(distanceToRootMap, node.getParent());
			}
		}
		distanceToRootMap.put(node, Double.valueOf(node.getBranchLength())
				+ distanceToRootMap.get(node.getParent()));
	}

	private int propagateSplitHeight(int index, Set<Integer> addedSplit) {
		int i = index;
		double ratio = 0;
		if (index < root.length - 1) {
			ratio = Double.valueOf(timeTrees[index + 1].getBestSplitNode(addedSplit).getInfo())/Double.valueOf(timeTrees[index].getBestSplitNode(addedSplit).getInfo());
			index++;
			while(ratio < 1.5 && ratio > 0.75 && addSplitTooTreeOnly(index, localTrees[index], addedSplit) >= 0){
				if (index == root.length-1) {
					break;
				}
				ratio = Double.valueOf(timeTrees[index + 1].getBestSplitNode(addedSplit).getInfo())/Double.valueOf(timeTrees[index].getBestSplitNode(addedSplit).getInfo());
				index++;
			}
		}
		index = i;
		if (index > 0) {
			ratio = Double.valueOf(timeTrees[index - 1].getBestSplitNode(addedSplit).getInfo())/Double.valueOf(timeTrees[index].getBestSplitNode(addedSplit).getInfo());
			index--;
			while(ratio < 1.5 && ratio > 0.75 && addSplitTooTreeOnly(index, localTrees[index], addedSplit) >= 0){
				if (index == 0) {
					break;
				}
				ratio = Double.valueOf(timeTrees[index - 1].getBestSplitNode(addedSplit).getInfo())/Double.valueOf(timeTrees[index].getBestSplitNode(addedSplit).getInfo());
				index--;
			}
		}
		return 0;
	}
	
	private int propagateSplit(int index, Set<Integer> addedSplit) {
		int i = index;
		if (index < root.length - 1) {
			index++;
			while(addSplitTooTreeOnly(index, localTrees[index], addedSplit) >= 0){
				if (index == root.length-1) {
					break;
				}
				index++;
			}
		}
		index = i;
		if (index > 0) {
			index--;
			while(addSplitTooTreeOnly(index, localTrees[index], addedSplit) >= 0){
				if (index == 0) {
					break;
				}
				index--;
			}
		}
		return 0;
	}
	
	private int propagateSplit3(int index, Set<Integer> addedSplit) {
		int i = index;
		if (index < root.length - 1) {
			index++;
			while(timeTrees[index].getSplits().contains(addedSplit) && addSplitTooTreeOnly(index, localTrees[index], addedSplit) >= 0){
				if (index == root.length-1) {
					break;
				}
				index++;
			}
		}
		index = i;
		if (index > 0) {
			index--;
			while(timeTrees[index].getSplits().contains(addedSplit) && addSplitTooTreeOnly(index, localTrees[index], addedSplit) >= 0){
				if (index == 0) {
					break;
				}
				index--;
			}
		}
		return 0;
	}

	private void initiateCompatibleMatrix() {
		compatibleMatrix = new boolean[root.length][root.length];
		for (int i = 0; i < root.length - 1; i++) {
			for (int j = i + 1; j < root.length; j++) {
				boolean compatible = isCompatible(originalSplits.get(i), originalSplits.get(j));
				compatibleMatrix[i][j] = compatible;
			}
		}
	}
	
	private Tree[] initializeTrees() {
		System.out.print("Initializing Trees ");
		Tree[] localTrees = new Tree[matrix[0].length];
		int counter = 0;
		originalSplits = new ArrayList<Set<Integer>>();
		for (int i = 0; i < matrix[0].length; i++) {
 			Tree tree = new Tree();
			Node rootNode = new Node();
			rootNode.setId("root");
			tree.setRoot(rootNode);
			tree.getNodes().add(rootNode);
			
			if (isInformative(i)) {
				counter ++;
				Node inputSplitNode = new Node();
				inputSplitNode.setParent(rootNode);
				tree.getNodes().add(inputSplitNode);
				for (int j = 0; j < matrix.length; j++) {
					Node node = new Node();
					node.setId(String.valueOf(j + 1));
					node.getSplit().add(j + 1);
					rootNode.getSplit().add(j + 1);
					tree.getNodes().add(node);
					if (matrix[j][i] == root[i]) {
						rootNode.getChildren().add(node);
						node.setParent(rootNode);
					} else {
						inputSplitNode.getChildren().add(node);
						node.setParent(inputSplitNode);
						inputSplitNode.getSplit().add(j + 1);
					}
					tree.getSplits().add(node.getSplit());
				}
				rootNode.getChildren().add(inputSplitNode);
				tree.getSplits().add(inputSplitNode.getSplit());
				Set<Integer> originalSplit = new HashSet<Integer>();
				originalSplit.addAll(inputSplitNode.getSplit());
				originalSplits.add(originalSplit);
			} else {
				for (int j = 0; j < matrix.length; j++) {
					Node node = new Node();
					node.setId(String.valueOf(j + 1));
					node.getSplit().add(j + 1);
					rootNode.getSplit().add(j + 1);
					rootNode.getChildren().add(node);
					node.setParent(rootNode);
					tree.getNodes().add(node);
					tree.getSplits().add(node.getSplit());
				}
				originalSplits.add(new HashSet<Integer>());
			}
			localTrees[i] = tree;
			tree.getSplits().add(rootNode.getSplit());

		}
		System.out.println("(" + counter + ")");
		return localTrees;
	}
	
	private boolean isInformative(int i) {
		int mutationCount = 0;
		for (int j = 0; j < matrix.length; j++) {
			if (matrix[j][i] != root[i]) {
				mutationCount++;
			}
		}
		if (mutationCount == 1 || matrix.length - mutationCount == 1) {
			return false;
		}
		return true;
	}

	private int[] buildRoot(int[][] matrix) {
		int[] root = new int[matrix[0].length];
		for (int j = 0; j < matrix[0].length; j++) {
			int zeroCounter = 0;
			int oneCounter = 0;
			Set<Integer> originalSplit = new HashSet<Integer>();
			for (int i = 0; i < matrix.length; i++) {
				if (matrix[i][j] == 0) {
					zeroCounter++;
				} else {
					originalSplit.add(i+1);
					oneCounter++;
				}
			}
			if (zeroCounter == 1 || oneCounter == 1) {
				originalSplits.add(new HashSet<Integer>());
			}else{
				originalSplits.add(originalSplit);
			}
			if (oneCounter > zeroCounter) {
				root[j] = 1;
			}
		}
		return root;
	}
	
	private void buildNewRoot() {
		for (int j = 0; j < matrix[0].length; j++) {
			Tree timeTree = timeTrees[j];
			Set<Integer> zeroSplit = new HashSet<Integer>();
			Set<Integer> oneSplit = new HashSet<Integer>();
			zeroSplit.addAll(reverseSplit(originalSplits.get(j)));
			oneSplit.addAll(originalSplits.get(j));
			Node zeroNode = timeTree.getBestSplitNode(zeroSplit);
			Node oneNode = timeTree.getBestSplitNode(oneSplit);
			double zeroHeight = zeroNode.getChildren().isEmpty() ? Double.valueOf(zeroNode.getBranchLength()):Double.valueOf(zeroNode.getInfo());
			double oneHeight = oneNode.getChildren().isEmpty() ? Double.valueOf(oneNode.getBranchLength()):Double.valueOf(oneNode.getInfo());
			if (zeroHeight > oneHeight) {
				root[j] = 0;
			}else{
				root[j] = 1;
			}
		}
	}

	private void timeSplitOverrideRule() {
		System.out.print("Incompatible Region Rule");
		Date startDate = new Date();
		int counter = 0;
		for (int i = 0; i < timeTrees.length - 1; i++) {
			// if (i > 1000) {
			// System.out.println();
			// }
			Set<Integer> split1 = originalSplits.get(i);
			if (!split1.isEmpty()) {
				HashMap<Set<Integer>, Integer> passedSplits = new LinkedHashMap<Set<Integer>, Integer>();
				for (int j = i + 1; j < timeTrees.length; j++) {
					Set<Integer> split2 = originalSplits.get(j);
					if (!split2.isEmpty()) {
						if (split1.equals(split2)
								|| split1.equals(reverseSplit(split2))) {
							break;
						} else if (!isCompatible(originalSplits.get(i), originalSplits.get(j))) {
							Set<Set<Integer>> keySet = passedSplits.keySet();
							List<Set<Integer>> keyList = new ArrayList<Set<Integer>>();
							keyList.addAll(keySet);
							boolean takenCareOf = false;
							for (int k = keyList.size() - 1; k >= 0; k--) {
								Set<Integer> split = keyList.get(k);
								if (!isCompatible(split, split2)) {
									counter += resolveConflict(passedSplits.get(split), j);
									takenCareOf = true;
									break;
								}
							}
							if (!takenCareOf) {
								counter += resolveConflict(i, j);
							}
							break;
						} else {
							Set<Set<Integer>> keySet = passedSplits.keySet();
							List<Set<Integer>> keyList = new ArrayList<Set<Integer>>();
							keyList.addAll(keySet);
							for (int k = keyList.size() - 1; k >= 0; k--) {
								Set<Integer> split = keyList.get(k);
								if (!isCompatible(split, split2)) {
									counter += resolveConflict(passedSplits.get(split), j);
									passedSplits.remove(split);
								}
							}
							passedSplits.remove(reverseSplit(split2));
							passedSplits.put(split2, j);
						}
					}
				}
			}
		}
		System.out.print(" [" + (double) (new Date().getTime() - startDate.getTime()) / 1000
				+ " Seconds] ");
		System.out.println("(" + counter + ")");
	}

	private int resolveConflict(int i, int j) {
		Set<Integer> split1 = originalSplits.get(i);
		Set<Integer> split2 = originalSplits.get(j);
		int counter = 0;
		for (int k = i+1; k < j; k++) {
			if (timeTrees[k].getSplits().contains(split1)) {
				counter += addSplitTooTree(k, localTrees[k], split1, false, false, true);
			}else if (timeTrees[k].getSplits().contains(split2)){
				counter += addSplitTooTree(k, localTrees[k], split2, false, false, true);
			}
		}
		return counter;
	}

	private double getSplitScoreForRegion2(Set<Integer> split, int[] region) {
		double score = 0.0;
		int count = 0;
		for (int i = region[0]; i <= region[1]; i++) {
			if (i == positions.size()) {
				break;
			}
			count = 0;
			for (Integer integer : split) {
				count += matrix[integer-1][i];
			}
			if (count != split.size() && count != 0) {
				score+= 1.0;
			}
		}
		int lowestPosition = 0;
		int highestPosition = 0;
		lowestPosition = positions.get(region[0]);
		if (region[1] == positions.size()) {
			highestPosition = (int) TOTAL_SEQ_LENGTH;
		}else{
			highestPosition = positions.get(region[1]);
		}
		if (highestPosition == lowestPosition){
			return 1;
		}else{
			return score/(highestPosition - lowestPosition + 1);
		}
	}

	private void RefineNonBinaryNodes() {
		for (int i = 0; i < localTrees.length; i++) {
			Tree tree = localTrees[i];
			Tree upgmaTree = timeTrees[i];
			List<Node> nodesToBreak = new ArrayList<Node>();
			for (Node node : tree.getNodes()) {
				if (node.getChildren().size() > 2){
					nodesToBreak.add(node);
				}
			}
			while (!nodesToBreak.isEmpty()){
				Node node = nodesToBreak.remove(0);
				Node upgmaNode = upgmaTree.getBestSplitNode(node.getSplit());
				Set<Integer> split1 = upgmaNode.getChildren().get(0).getSplit();
				Set<Integer> split2 = upgmaNode.getChildren().get(1).getSplit();
				Set<Integer> tmp;
				// Assuming split1 is the bigger one
				if (split1.size() < split2.size()) {
					tmp = split2;
					split2 = split1;
					split1 = tmp;
				}
				Set<Node> set1 = new HashSet<Node>();
				Set<Node> set2 = new HashSet<Node>();
				// Check if splits fits better for split1 or split2
				Node worstNode = null;
				double worstrate = 0;
				for (Node child : node.getChildren()) {
					Set<Integer> childSplit = child.getSplit();
					int count = 0;
					for (Integer integer : childSplit) {
						if (split1.contains(integer)) {
							count++;
						}
					}
					if (count*2 >= childSplit.size()) {
						if (worstNode == null) {
							worstNode = child;
							worstrate = (double)count/(double)childSplit.size();
						}else{
							if ((double)count/(double)childSplit.size() < worstrate) {
								worstrate = (double)count/(double)childSplit.size();
								worstNode = child;
							}
						}
						set1.add(child);
					}else{
						if (worstNode == null) {
							worstNode = child;
							worstrate = (double)count/(double)childSplit.size();
						}else{
							if ((double)count/(double)childSplit.size() < worstrate) {
								worstrate = (double)count/(double)childSplit.size();
								worstNode = child;
							}
						}
						set2.add(child);
					}
				}
				if (set2.size() == 0) {
					set2.add(worstNode);
					set1.remove(worstNode);
				}
				if (set1.size() == 0) {
					set1.add(worstNode);
					set2.remove(worstNode);
				}
				//Dividing the node;
				Node node1 = new Node();
				Node node2 = new Node();
				node1.getChildren().addAll(set1);
				node2.getChildren().addAll(set2);
				for (Node n : set1) {
					n.setParent(node1);
					node1.getSplit().addAll(n.getSplit());
				}
				for (Node n : set2) {
					n.setParent(node2);
					node2.getSplit().addAll(n.getSplit());
				}
				node.getChildren().clear();
				node.getChildren().add(node1);
				node.getChildren().add(node2);
				tree.getNodes().add(node1);
				tree.getNodes().add(node2);
				tree.getSplits().add(node1.getSplit());
				tree.getSplits().add(node2.getSplit());
				if (node1.getChildren().size() > 2) {
					nodesToBreak.add(node1);
				}
				if (node2.getChildren().size() > 2) {
					nodesToBreak.add(node2);
				}
			}
		}
	}
	
	private Tree[] buildLocalTrees() {
		fullyCompatibleWithCheckNew();
		fullyCompatibleRegionRuleNew();
		
		timeSplitOverrideRule();
		
		propagationRuleNew();
		propagationHeightRule();
		propagationRule();
		
		timeSplitRule();
		
//		uniqueRefinementRule(false);
		
		oneSPREventRule(false);
		twoSPREventRule(false);
		propagationRule();
		
		RefineNonBinaryNodes();
		
		return localTrees;
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		new Main(args);
	}
}
