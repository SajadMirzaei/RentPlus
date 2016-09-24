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
import java.util.Comparator;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import java.util.TreeSet;

import object.DistanceObject;
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

	private List<float[][]> distanceMatrices = new ArrayList<float[][]>();
	private double[] timeTreeLengths;
	private double[] trueTreeLengths;
	private double[] localTreeLengths;
	private List<Integer> positions;
	private List<Integer> singletons;
	private Set<Integer> indicesToRemove;
	private Set<Short> overalSplit;
	private float TOTAL_SEQ_LENGTH = 0;
	private double UPGMA_THRESHOD = 0.2;
	private double UPGMA_RANGE_CHECK = 1;
	private int WINDOW_SIZE = 5;
	private float theta;
	public static int NUM_TAXA;
	
	private String outputName;
	
	private int[][] bestRegions;
	private List<Set<Short>> originalSplits;
	
	private void logMemory(){
		System.out.println("-------------------------");
		System.out.println("Used Memory: " + (Runtime.getRuntime().totalMemory()-Runtime.getRuntime().freeMemory())/ 1024 /1024 + "MB");
		System.out.println("Total Memory: " + Runtime.getRuntime().totalMemory()/ 1024 /1024 + "MB");
		System.out.println("-------------------------");
	}
	
	public Main(String args[]) {
		boolean needHelp = false;
		if (args.length == 0) {
			needHelp = true;
		}else{
			for (int i = 0; i < args.length; i++) {
				if (args[i]=="-h" || args[i]=="-help"){
					needHelp = true;
					break;
				}
			}
		}
		if (needHelp) {
			try {
				InputStreamReader in = new InputStreamReader(Util.load("/help.txt"));
				BufferedReader br = new BufferedReader(in);
				String line = br.readLine();
				while (line != null) {
					System.out.println(line);
					line = br.readLine();
				}
				br.close();
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			} catch (IOException e) {
				e.printStackTrace();
			}
			return;
		}
		if (args[0].equals("-t")) {
			Util.BRANCH_LENGTH = true;
			args[0] = args[1];
			if (args.length > 2) {
				args[1] = args[2];
			}else{
				args[1] = null;
			}
		}
		outputName = args[0];
		Date startTime = new Date();
		System.out.println("Start Time: " + startTime);
		overalSplit = new HashSet<Short>();
		originalSplits = new ArrayList<Set<Short>>();
		singletons = new ArrayList<Integer>();
		logMemory();
		buildDataFromHapMatrix(args[0]);
		logMemory();
		NUM_TAXA = matrix.length;
		System.out.println("Initiating TMRCAs");
		for (int i = 0; i < matrix.length; i++) {
			overalSplit.add((short) (i+1));
		}
		timeTreeLengths = new double[matrix[0].length];
		trueTreeLengths = new double[matrix[0].length];
		localTreeLengths = new double[matrix[0].length];
		logMemory();
		if (args.length > 1 && args[1] != null) {
			System.out.println("Reading MS files");
			getTrueTreesFromMsFile(args[1]);
			System.gc();
			getDistanceMatricesForTrueTrees();
			positions = getPositionArrayFromMSFile(args[1]);
		}else{
			System.out.println("Reading positions");
			positions = getPositionArrayFromInputFile(args[0]);
		}
		// Calculating theta
		float sum = 0;
		for (int i = 1; i <= matrix.length-1; i++) {
			sum+= 1.0/i;
		}
		theta = matrix[0].length/sum;
		
		logMemory();
		makeDistanceMatrices();
		System.gc();
		logMemory();
//		makeUPGMATimeTrees();
		makeGuideTrees();
		System.gc();
		logMemory();
		buildNewRoot();
		localTrees = initializeTrees();
		logMemory();
		double[] sortedTreeLengths = new double[localTrees.length];
		for (int i = 0; i < sortedTreeLengths.length; i++) {
			sortedTreeLengths[i] = timeTreeLengths[i];
		}
		bestRegions = new int[localTrees.length][2];
		Arrays.sort(sortedTreeLengths);
		buildLocalTrees();
		if (Util.BRANCH_LENGTH) {
			inferBranchLengths();
		}
		logMemory();
	    Date endTime = new Date();
//	    outputTreesSeperate(timeTrees, outputName, "Upgma");
//	    scaleTimes(timeTreeLengths);
	    if (args.length > 1 && args[1] != null) {
//	    	compareAllTreesRooted(positions);
	    	compareAllTrees(positions);
	    	System.out.println("-----------------------");
		}else{
			outputInferredTrees(positions);
//			outputTimes(args[0]);
		}
	    outputTreesSeperate(localTrees,outputName, "");
	    outputTMRCAs(outputName);
		System.out.println("Running Time: "
				+ (double) (endTime.getTime() - startTime.getTime()) / 1000
				+ " Seconds");
		logMemory();
		
	}


	private void makeDistanceMatrices(){
		// initiating all distance matrices for all sites
		System.out.print("Making Distance Matrices");
		Date startDate = new Date();
		for (int i = 1; i <= positions.size(); i++) {
			float[][] distanceMatrix = new float[matrix.length][matrix.length];
			distanceMatrices.add(distanceMatrix);
		}
		for (int j = 0; j < matrix.length-1; j++) {
			for (int j2 = j+1; j2 < matrix.length; j2++) {
				// indices for next or previous non informative site with different values for j and j2
				// number of mutations from i to indexP and indexM
				int indexP = 0;
				int indexM = 0;
				int countFP = 0;
				int countFM = 0;
				// number of informative mutations within the window size 
				int countWIP = 0;
				int countWIM = 0;
				// number of mutations within the window size including singletons
				int countWP = 0;
				int countWM = 0;
				
				boolean selfInformative = false;
				boolean selfNonInformative = false;
				
				// initializing indexP, countFP, countWIP, and count WP
				int i2 = 1;
				while (i2<=positions.size()-1){
					if (matrix[j][i2] + matrix[j2][i2] == 1) {
						if (isInformative(i2)) {
							indexP = i2;
							if (i2 <= WINDOW_SIZE) {
								countWIP++;
							}else{
								countFP++;
								break;
							}
						}
						if (i2 <= WINDOW_SIZE) {
							countWP++;
						}
						countFP++;
					}
					i2++;
				}
				
				if (matrix[j][0] + matrix[j2][0] == 1) {
					if (isInformative(0)){
						selfInformative = true;
					}else{
						selfNonInformative = true;
					}
				}
				//calculating distance for index 0
				float distance = 0;
				if (countWIP > 0 || selfInformative){
					if (selfInformative || selfNonInformative) {
						distance = (countWP+1)/(float)positions.get(WINDOW_SIZE);
					}else{
						distance = (countWP)/(float)positions.get(WINDOW_SIZE);
					}
				}else{
					if (selfNonInformative) {
						distance = (countFP+1)/(float)positions.get(indexP);
					}else{
						distance = countFP/(float)positions.get(indexP);
					}
					distance = countFP/(float)positions.get(indexP);
				}
				distanceMatrices.get(0)[j][j2] = distance;
				distanceMatrices.get(0)[j2][j] = distance;
				for (int i = 1; i < positions.size(); i++) {
					selfInformative = false;
					selfNonInformative = false;
					if (matrix[j][i] + matrix[j2][i] == 1) {
						if (isInformative(i)){
							selfInformative = true;
						}else{
							selfNonInformative = true;
						}
					}
					
					if (matrix[j][i-1] + matrix[j2][i-1] == 1) {
						if (isInformative(i-1)) {
							countWIM++;
						}
						countWM++;
						countFM++;
					}
					
					if (i-WINDOW_SIZE-1 >= 0) {
						if (matrix[j][i-WINDOW_SIZE-1] + matrix[j2][i-WINDOW_SIZE-1] == 1) {
							countWM--;
							if (isInformative(i-WINDOW_SIZE-1)) {
								countWIM--;
								indexM = i-WINDOW_SIZE-1;
								countFM = countWM+1;
							}
						}
					}
					if (matrix[j][i] + matrix[j2][i] == 1) {
						if (isInformative(i)) {
							countWIP--;
						}
						countWP--;
						countFP--;
					}
					if (i+WINDOW_SIZE <= positions.size()-1) {
						if (matrix[j][i+WINDOW_SIZE] + matrix[j2][i+WINDOW_SIZE] == 1) {
							if (isInformative(i+WINDOW_SIZE)) {
								countWIP++;
							}
							countWP++;
						}
					}
					if(indexP <= i+WINDOW_SIZE){
						i2 = 1;
						while (i+WINDOW_SIZE+i2 <= positions.size()-1){
							if (matrix[j][i+WINDOW_SIZE+i2] + matrix[j2][i+WINDOW_SIZE+i2] == 1) {
								countFP++;
								if (isInformative(i+WINDOW_SIZE+i2)) {
									indexP = i+WINDOW_SIZE+i2;
									break;
								}
							}
							i2++;
							if (i+WINDOW_SIZE+i2 == positions.size()-1){
								indexP = positions.size()-1;
							}
						}
					}
					
					
					// Calculating the distances
					int minTravel = 0;
					int totalCount = 0;
					if (countWIM >=1 || selfInformative) {
						minTravel = Math.max(0, i-WINDOW_SIZE);
						totalCount = countWM;
					}else{
						minTravel = indexM;
						totalCount = countFM;
					}
					
					int maxTravel = 0;
					if (countWIP >= 1 || selfInformative) {
						maxTravel = Math.min(i+WINDOW_SIZE, positions.size()-1);
						totalCount += countWP;
					}else{
						maxTravel = indexP;
						totalCount += countFP;
					}
					if (selfNonInformative || selfInformative) {
						totalCount++;
					}
					
					int lowestPosition = 0;
					int highestPosition = 0;
					float score;
					if (minTravel == 0) {
						lowestPosition = 0;
					}else{
						lowestPosition = (positions.get(minTravel)+positions.get(minTravel-1))/2;
					}
					if (maxTravel >= positions.size()-1) {
						highestPosition = (int) TOTAL_SEQ_LENGTH;
					}else{
						highestPosition = (positions.get(maxTravel) + positions.get(maxTravel+1))/2;
					}
					if (highestPosition == lowestPosition){
						score = 1;
					}else{
//						score = (double)totalCount/(double)(highestPosition - lowestPosition);
						score = (float)totalCount/(theta*(float)(highestPosition - lowestPosition)/TOTAL_SEQ_LENGTH);
					}
					distanceMatrices.get(i)[j][j2] = score;
					distanceMatrices.get(i)[j2][j] = score;
				}
			}
		}
		System.out.println(" [" + (double) (new Date().getTime() - startDate.getTime()) / 1000
				+ " Seconds] ");
	}
	
	private void makeUPGMATimeTrees() {
		System.out.print("Making Guide Trees");
		Date startDate = new Date();
		timeTrees = new Tree[positions.size()];
		for (int i = 0; i < timeTrees.length; i++) {
			Tree upgmaTree = new Tree();
			float[][] distanceMatrix = distanceMatrices.get(i);
			HashMap<Node, List<Short>> clusters = new LinkedHashMap<Node, List<Short>>();
			HashMap<Set<Node>, Float> distanceMap = new HashMap<Set<Node>, Float>();
			//initiating clusters
			for (int j = 0; j < matrix.length; j++) {
				Node node = new Node();
				node.setId(String.valueOf(j+1));
				node.setInfo("0.0");
				List<Short> taxaCovering = new ArrayList<Short>();
				taxaCovering.add((short) (j+1));
				clusters.put(node, taxaCovering);
			}
			//
			while(clusters.keySet().size() > 1){
				List<Node> clusterNodes = new ArrayList<Node>();
				clusterNodes.addAll(clusters.keySet());
				float minDist = Float.MAX_VALUE;
				Set<Node> nodesToConcat = new HashSet<Node>();
				Map<Set<Node>, Float> validNodesToContactMap = new HashMap<Set<Node>, Float>();
				if (clusterNodes.size() == 2) {
					nodesToConcat.addAll(clusterNodes);
					minDist = getDistanceOfClusters(clusters.get(clusterNodes.get(0)),clusters.get(clusterNodes.get(1)), distanceMatrix);
					validNodesToContactMap.put(nodesToConcat, minDist);
				}else{
					for (int j = 0; j < clusterNodes.size()-1; j++) {
						for (int j2 = j+1; j2 < clusterNodes.size(); j2++) {
							Node node1 = clusterNodes.get(j);
							Node node2 = clusterNodes.get(j2);
							Set<Short> combinedSplit = new HashSet<Short>();
							combinedSplit.addAll(clusters.get(node1));
							combinedSplit.addAll(clusters.get(node2));
							if (isCompatible(originalSplits.get(i), combinedSplit)) {
								Set<Node> pair = new HashSet<Node>();
								pair.add(node1);
								pair.add(node2);
								float distance;
//								if (combinedSplit.size() > 2 && distanceMap.containsKey(pair)) {
//									distance = distanceMap.get(pair);
//								} else {
									distance = getDistanceOfClusters(clusters.get(node1), clusters.get(node2),
											distanceMatrix);
//									if (combinedSplit.size() > 2) {
//										distanceMap.put(pair, distance);
//									}
//								}
								// System.out.println(distanceMap.size());
								if (distance <= minDist * (1 + UPGMA_THRESHOD)) {
									if (distance < minDist) {
										minDist = distance;
										List<Set<Node>> iterationList = new ArrayList<Set<Node>>();
										iterationList.addAll(validNodesToContactMap.keySet());
										for (Set<Node> set : iterationList) {
											if (validNodesToContactMap.get(set) > minDist * (1 + UPGMA_THRESHOD)) {
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
				}
				LinkedHashMap<Set<Node>, Float> sortedMap = findBestPairs(validNodesToContactMap, clusters, i);
				while(sortedMap.size() > 0){
					List<Set<Node>> reversedlist = new ArrayList<Set<Node>>(sortedMap.keySet());
					nodesToConcat = reversedlist.get(reversedlist.size()-1);
					if (clusters.keySet().containsAll(nodesToConcat)) {
						Node parent = new Node();
						parent.getChildren().addAll(nodesToConcat);
						List<Short> coveringTaxa = new ArrayList<Short>();
						parent.setInfo(String.valueOf(sortedMap.get(nodesToConcat)/2));
						for (Node node : nodesToConcat) {
							node.setParent(parent);
							double branchLength = (sortedMap.get(nodesToConcat)/2)-Double.valueOf(node.getInfo());
							node.setBranchLength(String.valueOf(branchLength));
							coveringTaxa.addAll(clusters.get(node));
							clusters.remove(node);
//							distanceMap.remove(nodesToConcat);
						}
						clusters.put(parent, coveringTaxa);
					}
					sortedMap.remove(nodesToConcat);
				}
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


	private void makeGuideTrees() {
		System.out.print("Making Guide Trees");
		Date startDate = new Date();
		timeTrees = new Tree[positions.size()];
		Map<Set<Short>, short[]> compatibleRegionMap = new HashMap<Set<Short>, short[]>();
		for (int i = 0; i < timeTrees.length; i++) {
			Tree upgmaTree = new Tree();
			float[][] distanceMatrix = distanceMatrices.get(i);
			HashMap<Node, List<Short>> clusters = new LinkedHashMap<Node, List<Short>>();
//			HashMap<Node, List<DistanceObject>> nodeDistanceObjectMap = new HashMap<Node, List<DistanceObject>>();
			
			TreeSet<DistanceObject> distanceSet = new TreeSet<DistanceObject>(new Comparator<DistanceObject>() {
				@Override
				public int compare(DistanceObject o1, DistanceObject o2) {
					return o1.getDistance() < o2.getDistance() ? -1 : 1;
				}
			});
			
			TreeSet<DistanceObject> compatibleDistanceSet = new TreeSet<DistanceObject>(new Comparator<DistanceObject>() {
				@Override
				public int compare(DistanceObject o1, DistanceObject o2) {
					return o1.getCompatibleRange() <= o2.getCompatibleRange() ? 1 : -1;
				}
			});
			
			Set<Node> removedNodes = new HashSet<Node>();
			// initiating clusters
			for (int j = 0; j < matrix.length; j++) {
				Node node = new Node(1);
				node.setId(String.valueOf(j + 1));
				node.setInfo("0.0");
				List<Short> taxaCovering = new ArrayList<Short>();
				taxaCovering.add((short) (j + 1));
				clusters.put(node, taxaCovering);
//				List<DistanceObject> distanceObjects = new ArrayList<DistanceObject>();
//				nodeDistanceObjectMap.put(node, distanceObjects);
			}
			//
			List<Node> clusterNodes = new ArrayList<Node>();
			clusterNodes.addAll(clusters.keySet());
			Node[] nodesToConcat = new Node[2];
			// Making the initial TreeMap with the distances sorted
			for (int j = 0; j < clusterNodes.size() - 1; j++) {
				for (int j2 = j + 1; j2 < clusterNodes.size(); j2++) {
					Node node1 = clusterNodes.get(j);
					Node node2 = clusterNodes.get(j2);
					Set<Short> combinedSplit = new HashSet<Short>();
					combinedSplit.addAll(clusters.get(node1));
					combinedSplit.addAll(clusters.get(node2));
					if (isCompatible(originalSplits.get(i), combinedSplit)) {
						Node[] pair = new Node[2];
						pair[0] = node1;
						pair[1] = node2;
						Float distance;
						distance = getDistanceOfClusters(clusters.get(node1), clusters.get(node2), distanceMatrix);
						DistanceObject distanceObject = new DistanceObject();
						distanceObject.setNodes(pair);
						distanceObject.setDistance(distance);
						distanceObject.setCompatibleRange((short)0);
						distanceSet.add(distanceObject);
//						nodeDistanceObjectMap.get(node1).add(distanceObject);
//						nodeDistanceObjectMap.get(node2).add(distanceObject);
					}
				}
			}
			while (clusters.keySet().size() > 1) {
				// Making the second Map sorted by the compatible distance
				DistanceObject firstElement = distanceSet.first();
				while (removedNodes.contains(firstElement.getNodes()[0])
						|| removedNodes.contains(firstElement.getNodes()[1])) {
					distanceSet.pollFirst();
					firstElement = distanceSet.first();
				}
				for (DistanceObject key : distanceSet) {
					if (!removedNodes.contains(key.getNodes()[0])
							&& !removedNodes.contains(key.getNodes()[1])) {
						if (key.getDistance() == 0.0) {
							key.setCompatibleRange(Short.MAX_VALUE);
							compatibleDistanceSet.add(key);
						}else if (key.getDistance() <= firstElement.getDistance() * (1 + UPGMA_THRESHOD)) {
							if (key.getCompatibleRange() == 0) {
								Set<Short> coveringSet = new HashSet<Short>(clusters.get(key.getNodes()[0]));
								coveringSet.addAll(clusters.get(key.getNodes()[1]));
								short[] range = compatibleRegionMap.get(coveringSet);
								short compatipleDistance = 0; 
								if (range != null && range[0] <= i && i <= range[1]) {
									if (range[0] == 0) {
										compatipleDistance = (short) (range[1] - i);
									}else if (range[1] == positions.size()-1){
										compatipleDistance = (short) (i - range[0]);
									}else{
										compatipleDistance = (short) Math.min(range[1] - i, i - range[0]);
									}
//									compatipleDistance = (short) (range[1] - range[0] + 1);
								} else {
									range = findCompatibleRegion(key.getNodes(), coveringSet, i);
									if (range[0] == 0) {
										compatipleDistance = (short) (range[1] - i);
									}else if (range[1] == positions.size()-1){
										compatipleDistance = (short) (i - range[0]);
									}else{
										compatipleDistance = (short) Math.min(range[1] - i, i - range[0]);
									}
//									compatipleDistance = (short) (range[1] - range[0] + 1);
									compatibleRegionMap.put(coveringSet, range);
								}
//								short compatipleDistance = findCompatibleDistance(key.getNodes(), clusters, i);
								key.setCompatibleRange(compatipleDistance);
								compatibleDistanceSet.add(key);
							}
						} else {
							break;
						}
					}
				}
				DistanceObject firstConcatObject = compatibleDistanceSet.first();
				while (removedNodes.contains(firstConcatObject.getNodes()[0])
						|| removedNodes.contains(firstConcatObject.getNodes()[1])) {
					compatibleDistanceSet.pollFirst();
					firstConcatObject = compatibleDistanceSet.first();
				}
				nodesToConcat = firstConcatObject.getNodes();
				int capacity = nodesToConcat[0].getSplit().size() + nodesToConcat[1].getSplit().size();
				Node newNode = new Node(capacity);
				newNode.getChildren().add(nodesToConcat[0]);
				newNode.getChildren().add(nodesToConcat[1]);
				
				List<Short> coveringTaxa = new ArrayList<Short>();
				newNode.setInfo(String.valueOf(firstConcatObject.getDistance() / 2));
				for (Node node : nodesToConcat) {
					node.setParent(newNode);
					double branchLength = (firstConcatObject.getDistance() / 2) - Double.valueOf(node.getInfo());
					node.setBranchLength(String.valueOf(branchLength));
					coveringTaxa.addAll(clusters.get(node));
					clusters.remove(node);
					removedNodes.add(node);
//					for (DistanceObject dObject : nodeDistanceObjectMap.get(node)) {
//						distanceSet.remove(dObject);
//						compatibleDistanceSet.remove(dObject);
//					}
				}
//				nodeDistanceObjectMap.put(newNode, new ArrayList<DistanceObject>());
				for (Node node : clusters.keySet()) {
//					Set<Short> newCovering = new HashSet<Short>();
//					if (isInformative(i)) {
//						newCovering.addAll(coveringTaxa);
//						newCovering.addAll(clusters.get(node));
//					}
//					if (isCompatible(newCovering, originalSplits.get(i))) {
						Float newDistance = getDistanceOfClusters(coveringTaxa, clusters.get(node), distanceMatrix);
						Node[] newPair = new Node[2];
						newPair[0] = node;
						newPair[1] = newNode;
						DistanceObject newObject = new DistanceObject();
						newObject.setNodes(newPair);
						newObject.setDistance(newDistance);
						newObject.setCompatibleRange((short)0);
						distanceSet.add(newObject);
//						nodeDistanceObjectMap.get(node).add(newObject);
//						nodeDistanceObjectMap.get(newNode).add(newObject);
//					}
				}
				clusters.put(newNode, coveringTaxa);
			}
			upgmaTree.setRoot(clusters.keySet().iterator().next());
			upgmaTree.updateNodes();
			upgmaTree.updateSplits();
			timeTrees[i] = upgmaTree;
			Map<Node, Double> distanceToRootMap = buildDistanceToRootMap(upgmaTree);
			timeTreeLengths[i] = distanceToRootMap.get(upgmaTree.getTaxa().get(0));
		}
		System.out.println(" [" + (double) (new Date().getTime() - startDate.getTime()) / 1000 + " Seconds] ");
	}

	private float getDistanceOfClusters(List<Short> list, List<Short> list2, float[][] distanceMatrix) {
		float distance = 0.0F;
		for (Short i : list) {
			for (Short j : list2) {
				distance += distanceMatrix[i-1][j-1];
			}
		}
		return distance/(list.size()*list2.size());
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
			float recombLength = Float.valueOf(recombLengthString);
			TOTAL_SEQ_LENGTH = recombLength;
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
	
	private  Short findCompatibleDistance(Node[] nodes, HashMap<Node, List<Short>> clusters, int i){
		int index = 1;
		short counter = 0;
		Set<Short> combinedSplit = new HashSet<Short>();
		for (Node node : nodes) {
			combinedSplit.addAll(clusters.get(node));
		}
		while((i+index < positions.size() || i-index >= 0) && index < UPGMA_RANGE_CHECK*positions.size()){
				if (i-index >= 0 && isCompatible(originalSplits.get(i-index), combinedSplit)) {
					counter++;
				}
//				else{
//					break;
//				}
				else if (i-index >= 0){
					break;
				}
				if (i+index <= positions.size()-1 && isCompatible(originalSplits.get(i+index), combinedSplit)) {
					counter++;
				}
//				else{
//					break;
//				}
				else if(i+index <= positions.size()-1){
					break;
				}
			index++;
		}
		return counter;
	}
	
	private  short[] findCompatibleRegion(Node[] nodes, Set<Short> combinedSplit, int i){
		short index = 1;
		short minIndex = (short) i;
		short maxIndex = (short) i;
		while(i-index >= 0 && index < UPGMA_RANGE_CHECK*positions.size()){
				if (i-index >= 0 && isCompatible(originalSplits.get(i-index), combinedSplit)) {
					minIndex = (short) (i-index);
				}
				else{
					break;
				}
			index++;
		}
		index=1;
		while(i+index < positions.size() && index < UPGMA_RANGE_CHECK*positions.size()){
			if (i+index <= positions.size()-1 && isCompatible(originalSplits.get(i+index), combinedSplit)) {
				maxIndex = (short) (i+index);
			}
			else{
				break;
			}
			index++;
		}
		return new short[]{minIndex,maxIndex};
	}
	
	private LinkedHashMap<Set<Node>, Float> findBestPairs(Map<Set<Node>, Float> validNodesToContactMap, HashMap<Node, List<Short>> clusters, int i){
		int index = 1;
		LinkedHashMap<Set<Node>, Float> sortedMap = new LinkedHashMap<Set<Node>, Float>();
		while(validNodesToContactMap.size() > 1 && (i+index < positions.size() || i-index >= 0) && index < UPGMA_RANGE_CHECK*positions.size()){
			List<Set<Node>> iterationList = new ArrayList<Set<Node>>();
			iterationList.addAll(validNodesToContactMap.keySet());
			for (Set<Node> set : iterationList) {
				Set<Short> combinedSplit = new HashSet<Short>();
				for (Node node : set) {
					combinedSplit.addAll(clusters.get(node));
				}
				List<Node> nodes = new ArrayList<Node>();
				nodes.addAll(set);
				if (i-index >= 0 && !isCompatible(originalSplits.get(i-index), combinedSplit)) {
					sortedMap.put(set, validNodesToContactMap.get(set));
					validNodesToContactMap.remove(set);
					if (validNodesToContactMap.size() == 1) {
						break;
					}
				}
				else if (i+index <= positions.size()-1 && !isCompatible(originalSplits.get(i+index), combinedSplit)) {
					sortedMap.put(set, validNodesToContactMap.get(set));
					validNodesToContactMap.remove(set);
					if (validNodesToContactMap.size() == 1) {
						break;
					}
				}
			}
			index++;
		}
		if (validNodesToContactMap.size() > 1) {
			Set<Entry<Set<Node>, Float>> entryList = new TreeSet<Entry<Set<Node>,Float>>(new Comparator<Entry<Set<Node>, Float>>() {
				@Override
				public int compare(Entry<Set<Node>, Float> o1, Entry<Set<Node>, Float> o2) {
					if (o1.getValue().equals(o2.getValue())) {
						return 1;
					}else{
						return o2.getValue().compareTo(o1.getValue());
					}
				}
			});
			entryList.addAll(validNodesToContactMap.entrySet());
			for (Entry<Set<Node>, Float> entry : entryList) {
				sortedMap.put(entry.getKey(), entry.getValue());
			}
		}else{
			sortedMap.putAll(validNodesToContactMap);
		}
		return sortedMap;
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

//	private void scaleTimes(double[] times) {
//		double avg = 0;
//		for (double d : times) {
//			avg += d;
//		}
//		avg = avg/times.length;
//		double multiplier = 2/avg;
//		for (int i = 0; i < times.length; i++) {
//			times[i] = times[i]*multiplier;
//		}
//	}
	
	private void outputTimes(String inputUrl) {
		try {
			BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(
			          new FileOutputStream(inputUrl + ".times")));
			for (double d : timeTreeLengths) {
				writer.write(String.valueOf(d));
				writer.newLine();
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
			Set<Set<Short>> realTrueSplits = new HashSet<Set<Short>>();
			Set<Set<Short>> realTimeSplits = new HashSet<Set<Short>>();
			Set<Set<Short>> realLocalSplits = new HashSet<Set<Short>>();
			Set<Set<Short>> realArgSplits = new HashSet<Set<Short>>();
			Set<Set<Short>> realMargSplits = new HashSet<Set<Short>>();
			for (Set<Short> split : trueTree.getSplits()) {
				if (1 < split.size()
						&& split.size() < trueTree.getRoot().getSplit().size() - 1
						&& !realTrueSplits.contains(reverseSplit(split))) {
					realTrueSplits.add(split);
				}
			}
			numTotalSplits += realTrueSplits.size();
			for (Set<Short> split : localTree.getSplits()) {
				if (1 < split.size()
						&& split.size() < overalSplit.size() - 1
						&& !realLocalSplits.contains(reverseSplit(split))) {
					realLocalSplits.add(split);
				}
			}
			numTotalInfered += realLocalSplits.size();
			
			for (Set<Short> split : timeTree.getSplits()) {
				if (1 < split.size()
						&& split.size() < timeTree.getRoot().getSplit().size() - 1
						&& !realTimeSplits.contains(reverseSplit(split))) {
					realTimeSplits.add(split);
				}
			}
			
			Set<Set<Short>> tmpLocalSplit = new HashSet<Set<Short>>();
			Set<Set<Short>> tmpTimeSplit = new HashSet<Set<Short>>();
			Set<Set<Short>> tmpArgSplit = new HashSet<Set<Short>>();
			Set<Set<Short>> tmpMargSplit = new HashSet<Set<Short>>();
			for (Set<Short> set : realLocalSplits) {
				if (realTrueSplits.contains(set)
						|| realTrueSplits.contains(reverseSplit(set))) {
					tmpLocalSplit.add(set);
				}
			}
			for (Set<Short> set : realTimeSplits) {
				if (realTrueSplits.contains(set)
						|| realTrueSplits.contains(reverseSplit(set))) {
					tmpTimeSplit.add(set);
				}
			}
			for (Set<Short> set : realArgSplits) {
				if (realTrueSplits.contains(set)
						|| realTrueSplits.contains(reverseSplit(set))) {
					tmpArgSplit.add(set);
				}
			}
			for (Set<Short> set : realMargSplits) {
				if (realTrueSplits.contains(set)
						|| realTrueSplits.contains(reverseSplit(set))) {
					tmpMargSplit.add(set);
				}
			}
			Set<Set<Short>> sharedSplits = new HashSet<Set<Short>>();
			for (Set<Short> set : tmpTimeSplit) {
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
			String split = originalSplits.get(i).isEmpty()?getSingleMutation(i):originalSplits.get(i).toString();
			System.out.println("At index " + (i) + ", " + tmpLocalSplit.size()
					+ " (Inf), and " + tmpTimeSplit.size() + " (time)" + " of "
					+ realTrueSplits.size()
					+ " total splits are inferred.  Shared: "+ sharedSplits.size() 
					+ " ---  Inferred Tree: " + localTree.toString()+ " [" + localTreeLengths[i] + "]"
//					+ " UPGMA Tree: " + timeTree.toString() + " [" + timeTreeLengths[i]+ "]"
					+ "  True Tree: " + trueTree + " ["	+ trueTreeLengths[i] + "]"
//					 + " Best Compatible Region: [" + bestRegions[i][0] + ","
//					 + bestRegions[i][1] + "]"
					+ " {" + positions.get(i) + "}"
					+ " Original Split: " + split
					);
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
//		System.out.println("Time: " + numTimeSplits + " out of "
//				+ numTotalSplits + " splits = %" + (double) numTimeSplits
//				* 100.0 / (double) numTotalSplits);
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
			Set<Set<Short>> realTrueClades = new HashSet<Set<Short>>();
			Set<Set<Short>> realTimeClades = new HashSet<Set<Short>>();
			Set<Set<Short>> realLocalClades = new HashSet<Set<Short>>();
			for (Set<Short> clade : trueTree.getClades()) {
				if (1 < clade.size()) {
					realTrueClades.add(clade);
				}
			}
			numTotalClades += realTrueClades.size();
			for (Set<Short> clade : localTree.getClades()) {
				if (1 < clade.size()) {
					realLocalClades.add(clade);
				}
			}
			numTotalInfered += realLocalClades.size();
			
			for (Set<Short> clade : timeTree.getClades()) {
				if (1 < clade.size()) {
					realTimeClades.add(clade);
				}
			}
			
			Set<Set<Short>> tmpLocalClades = new HashSet<Set<Short>>();
			Set<Set<Short>> tmpTimeClades = new HashSet<Set<Short>>();
			
			for (Set<Short> set : realLocalClades) {
				if (realTrueClades.contains(set)) {
					tmpLocalClades.add(set);
				}
			}
			for (Set<Short> set : realTimeClades) {
				if (realTrueClades.contains(set)) {
					tmpTimeClades.add(set);
				}
			}
			Set<Set<Short>> sharedClades = new HashSet<Set<Short>>();
			for (Set<Short> set : tmpTimeClades) {
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
//					+ " UPGMA Tree: " + timeTree.toString() + " [" + timeTreeLengths[i]+ "]"
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
//		System.out.println("Time: " + numTimeSplits + " out of "
//				+ numTotalClades + " clades = %" + (double) numTimeSplits
//				* 100.0 / (double) numTotalClades);
	}
	
	private void outputInferredTrees(List<Integer> positions) {
		for (int i = 0; i < localTrees.length; i++) {
			Tree localTree = localTrees[i];
//			Tree timeTree = timeTrees[i];
			String split = originalSplits.get(i).isEmpty()?getSingleMutation(i):originalSplits.get(i).toString();
			System.out.println("At index " + (i) + ", " + " Inferred Tree: " + localTree.toString()+ " [" + localTreeLengths[i] + "]"
//					+ " UPGMA Tree: " + timeTree.toString() + " [" + timeTreeLengths[i]+ "]"
//					 + " Best Compatible Region: [" + bestRegions[i][0] + ","
//					 + bestRegions[i][1] + "]"
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
		int zero = -1;
		int one = -1;
		for (int j = 0; j < matrix.length; j++) {
			if (matrix[j][i] == 0) {
				if (zero == -1){
					zero = j;
				}else{
					if (one != -1) {
						return String.valueOf(one+1);
					}
				}
			}else{
				if (one == -1){
					one = j;
				}else{
					if (zero != -1) {
						return String.valueOf(zero+1);
					}
				}
			}
		}
		return String.valueOf(matrix.length);
	}

	private Set<Short> reverseSplit(Set<Short> split) {
		Set<Short> reverseSplit = new HashSet<Short>();
		for (Short splitNum : overalSplit) {
			if (!split.contains(splitNum)) {
				reverseSplit.add(splitNum);
			}
		}
		return reverseSplit;
	}

	private void buildDataFromHapMatrix(String url) {
		System.out.println("Reading Hap Matrix from input file");
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
		System.out.print("Adding Guide Tree Splits");
		Date startDate = new Date();
		int counter = 0;
		for (int i = 0; i < timeTrees.length; i++) {
			Tree tree = timeTrees[i];
			for (Set<Short> split : tree.getSplits()) {
				counter += addSplitTooTree(i, localTrees[i], split, false, false, true);
			}
		}
		System.out.println(" [" + (double) (new Date().getTime() - startDate.getTime()) / 1000
				+ " Seconds] ");
//		System.out.println("(" + counter + ")");
	}

	private boolean isCompatible(Set<Short> split1, Set<Short> split2) {
		if (split1.size() <= 1 || split2.size() <= 1 || split1.size() >= overalSplit.size() - 1
				|| split2.size() >= overalSplit.size() - 1) {
			return true;
		}
		boolean share = false;
		boolean coverAll = true;
		int size = 0;
		if (split1.size() <= split2.size()) {
			size = split2.size();
			for (Short integer : split1) {
				boolean contains = split2.contains(integer);
				if (!contains) {
					size++;
				}
				share = share | contains;
				coverAll = coverAll & contains;
			}
		}else{
			size = split1.size();
			for (Short integer : split2) {
				boolean contains = split1.contains(integer);
				if (!contains) {
					size++;
				}
				share = share | contains;
				coverAll = coverAll & contains;
			}
		}
		return coverAll || !share || (size == overalSplit.size() && share);
	}

	
	private void propagationRule() {
		System.out.print("Propagation Rule ");
		Date startDate = new Date();
		int counter = 0;
		for (int i = 0; i < localTrees.length; i++) {
			Tree localTree = localTrees[i];
			for (Set<Short> split : localTree.getSplits()) {
				if (split.size() > 1
						&& split.size() < localTree.getRoot().getSplit().size() - 1) {
					counter += propagateSplit(i, split);
				}
			}
		}
		System.out.println(" [" + (double) (new Date().getTime() - startDate.getTime()) / 1000
				+ " Seconds] ");
//		System.out.println("(" + counter + ")");
	}

	private void propagationRuleNew() {
		System.out.print("Propagation Rule using Guide Tree Topologies");
		Date startDate = new Date();
		for (int i = 0; i < localTrees.length; i++) {
			Tree localTree = localTrees[i];
			for (Set<Short> split : localTree.getSplits()) {
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
		System.out.print("Propagation Rule using Guide Tree TMRCAs");
		Date startDate = new Date();
		int counter = 0;
		for (int i = 0; i < localTrees.length; i++) {
			Tree localTree = localTrees[i];
			for (Set<Short> split : localTree.getSplits()) {
				if (split.size() > 1
						&& split.size() < overalSplit.size() - 1) {
					counter += propagateSplitHeight(i, split);
				}
			}
		}
		System.out.println(" [" + (double) (new Date().getTime() - startDate.getTime()) / 1000
				+ " Seconds] ");
//		System.out.println("(" + counter + ")");
	}
	
	private boolean isCompatible(Set<Short> concatenatedSplit, Tree tree) {
		boolean compatible = true;
		// checking compatibility of the split
		for (Set<Short> split : tree.getSplits()) {
			if (!isCompatible(concatenatedSplit, split)) {
				compatible = false;
				break;
			}
		}
		return compatible;
	}

	private void fullyCompatibleRegionRuleNew() {
		System.out.print("Fully Compatible Rule");
		Date startDate = new Date();
		List<Set<Set<Short>>> splitsToAdd = new ArrayList<Set<Set<Short>>>();
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
			Set<Set<Short>> splitSet = new LinkedHashSet<Set<Short>>();
			int index = 0;
			if (i - bestRegion[0] > bestRegion[1] - i) {
				index = i - bestRegion[0];
			} else {
				index = bestRegion[1] - i;
			}
			for (int j = 1; j <= index; j++) {
				Set<Set<Short>> prevSplitSet = new LinkedHashSet<Set<Short>>();
				if (bestRegion[0] <= i - j && i - j >= 0) {
					Set<Set<Short>> splits = localTrees[i - j].getSplits();
					prevSplitSet.addAll(splits);
				}
				if (i + j <= bestRegion[1] && i - j < root.length) {
					Set<Set<Short>> splits = localTrees[i + j].getSplits();
					Set<Set<Short>> splitsToRemove = new LinkedHashSet<Set<Short>>();
					Set<Set<Short>> addedSplits = new LinkedHashSet<Set<Short>>();
					for (Set<Short> split : splits) {
						if (prevSplitSet.isEmpty()) {
							splitSet.addAll(splits);
							break;
						}
						for (Set<Short> prevSplit : prevSplitSet) {
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
			Set<Set<Short>> splitSet = splitsToAdd.get(i);
			for (Set<Short> split : splitSet) {
				counter += addSplitTooTree(i, localTrees[i], split, false, true, false);
			}
		}
		System.out.println(" [" + (double) (new Date().getTime() - startDate.getTime()) / 1000
				+ " Seconds] ");
//		System.out.println("(" + counter + ")");
	}
	
	private void fullyCompatibleWithCheckNew(){
		System.out.print("Fully Compatible Rule using Guide Trees");
		Date startDate = new Date();
		initiateCompatibleMatrix();
		List<Set<Set<Short>>> splitsToAdd = new ArrayList<Set<Set<Short>>>();
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
			Set<Set<Short>> splitSet = new LinkedHashSet<Set<Short>>();
			int index = 0;
			if (i - bestRegion[0] > bestRegion[1] - i) {
				index = i - bestRegion[0];
			} else {
				index = bestRegion[1] - i;
			}
			for (int j = 1; j <= index; j++) {
				Set<Set<Short>> prevSplitSet = new LinkedHashSet<Set<Short>>();
				if (bestRegion[0] <= i - j && i - j >= 0) {
					Set<Set<Short>> splits = localTrees[i - j].getSplits();
					prevSplitSet.addAll(splits);
				}
				if (i + j <= bestRegion[1] && i - j < root.length) {
					Set<Set<Short>> splits = localTrees[i + j].getSplits();
					Set<Set<Short>> splitsToRemove = new LinkedHashSet<Set<Short>>();
					Set<Set<Short>> addedSplits = new LinkedHashSet<Set<Short>>();
					for (Set<Short> split : splits) {
						if (prevSplitSet.isEmpty()) {
							splitSet.addAll(splits);
							break;
						}
						for (Set<Short> prevSplit : prevSplitSet) {
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
			Set<Set<Short>> splitSet = splitsToAdd.get(i);
			for (Set<Short> split : splitSet) {
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
		System.out.println(" [" + (double) (new Date().getTime() - startDate.getTime()) / 1000
				+ " Seconds] ");
//		System.out.println("(" + counter + ")");
	}

	private int addSplitTooTree(int index, Tree tree, Set<Short> addedSplit,
			boolean propagate, boolean eligibilityCheck, boolean normalPropagation) {
		int counter = 0;
		if (!tree.getSplits().contains(addedSplit)
				&& !tree.getSplits().contains(reverseSplit(addedSplit))) {
			Set<Short> bestSplit = new HashSet<Short>();
			bestSplit = tree.getRoot().getSplit();
			boolean compatible = true;
			// checking compatibility of the split
			for (Set<Short> split : tree.getSplits()) {
				if (split.size() < addedSplit.size()) {
					List<Short> copyList = new ArrayList<Short>();
					copyList.addAll(split);
					copyList.retainAll(addedSplit);
					if (copyList.size() > 0 && copyList.size() < split.size()) {
						compatible = false;
						break;
					}
				} else {
					List<Short> copyList = new ArrayList<Short>();
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
					List<Node> nodesToRemove = new ArrayList<Node>();
					for (Node child : splitNode.getChildren()) {
						if (addedSplit.containsAll(child.getSplit()))
							nodesToRemove.add(child);
					}
					int capacity = 0;
					for (Node node : nodesToRemove) {
						capacity += node.getSplit().size();
					}
					Node newNode = new Node(capacity);
					newNode.setParent(splitNode);
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
	
	private int addSplitTooTreeOnly(int index, Tree tree, Set<Short> addedSplit) {
		int counter = 0;
		if (!tree.getSplits().contains(addedSplit)
				&& !tree.getSplits().contains(reverseSplit(addedSplit))) {
			Set<Short> bestSplit = new HashSet<Short>();
			bestSplit = tree.getRoot().getSplit();
			boolean compatible = true;
			// checking compatibility of the split
			for (Set<Short> split : tree.getSplits()) {
				if (split.size() < addedSplit.size()) {
					List<Short> copyList = new ArrayList<Short>();
					copyList.addAll(split);
					copyList.retainAll(addedSplit);
					if (copyList.size() > 0 && copyList.size() < split.size()) {
						compatible = false;
						break;
					}
				} else {
					List<Short> copyList = new ArrayList<Short>();
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
				List<Node> nodesToRemove = new ArrayList<Node>();
				for (Node child : splitNode.getChildren()) {
					if (addedSplit.containsAll(child.getSplit()))
						nodesToRemove.add(child);
				}
				int capacity = 0;
				for (Node node : nodesToRemove) {
					capacity += node.getSplit().size();
				}
				Node newNode = new Node(capacity);
				newNode.setParent(splitNode);
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

	private boolean isTimeEligible(int index, Tree tree, Set<Short> addedSplit) {
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
			Node taxa = new Node(0);
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

	private int propagateSplitHeight(int index, Set<Short> addedSplit) {
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
	
	private int propagateSplit(int index, Set<Short> addedSplit) {
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
	
	private int propagateSplit3(int index, Set<Short> addedSplit) {
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
		System.out.println("Initializing Trees ");
		Tree[] localTrees = new Tree[matrix[0].length];
		int counter = 0;
		originalSplits = new ArrayList<Set<Short>>();
		for (int i = 0; i < matrix[0].length; i++) {
 			Tree tree = new Tree();
			Node rootNode = new Node(NUM_TAXA);
			rootNode.setId("root");
			tree.setRoot(rootNode);
			tree.getNodes().add(rootNode);
			
			if (isInformative(i)) {
				counter ++;
				Node inputSplitNode = new Node(2);
				inputSplitNode.setParent(rootNode);
				tree.getNodes().add(inputSplitNode);
				for (int j = 0; j < matrix.length; j++) {
					Node node = new Node(1);
					node.setId(String.valueOf((short) (j+1)));
					node.getSplit().add((short) (j+1));
					rootNode.getSplit().add((short) (j+1));
					tree.getNodes().add(node);
					if (matrix[j][i] == root[i]) {
						rootNode.getChildren().add(node);
						node.setParent(rootNode);
					} else {
						inputSplitNode.getChildren().add(node);
						node.setParent(inputSplitNode);
						inputSplitNode.getSplit().add((short) (j+1));
					}
					tree.getSplits().add(node.getSplit());
				}
				rootNode.getChildren().add(inputSplitNode);
				tree.getSplits().add(inputSplitNode.getSplit());
				Set<Short> originalSplit = new HashSet<Short>();
				originalSplit.addAll(inputSplitNode.getSplit());
				originalSplits.add(originalSplit);
			} else {
				for (int j = 0; j < matrix.length; j++) {
					Node node = new Node(1);
					node.setId(String.valueOf(j + 1));
					node.getSplit().add((short) (j+1));
					rootNode.getSplit().add((short) (j+1));
					rootNode.getChildren().add(node);
					node.setParent(rootNode);
					tree.getNodes().add(node);
					tree.getSplits().add(node.getSplit());
				}
				originalSplits.add(new HashSet<Short>());
			}
			localTrees[i] = tree;
			tree.getSplits().add(rootNode.getSplit());

		}
//		System.out.println("(" + counter + ")");
		return localTrees;
	}
	
	private boolean isInformative(int i) {
		return singletons.get(i)==0;
	}
	
	private int[] buildRoot(int[][] matrix) {
		int[] root = new int[matrix[0].length];
		for (int j = 0; j < matrix[0].length; j++) {
			int zeroCounter = 0;
			int oneCounter = 0;
			int zeroValue = 0;
			int oneValue = 0;
			Set<Short> originalSplit = new HashSet<Short>();
			for (int i = 0; i < matrix.length; i++) {
				if (matrix[i][j] == 0) {
					zeroCounter++;
					zeroValue = i+1;
				} else {
					originalSplit.add((short) (i+1));
					oneCounter++;
					oneValue = i+1;
				}
			}
			if (zeroCounter == 1 || oneCounter == 1) {
				originalSplits.add(new HashSet<Short>());
				singletons.add(zeroCounter==1?zeroValue:oneValue);
			}else{
				originalSplits.add(originalSplit);
				singletons.add(0);
			}
			if (oneCounter > zeroCounter) {
				root[j] = 1;
			}
		}
		return root;
	}
	
	private void buildNewRoot() {
		System.out.println("Assigning Local Tree roots");
		for (int j = 0; j < matrix[0].length; j++) {
			Tree timeTree = timeTrees[j];
			Set<Short> zeroSplit = new HashSet<Short>();
			Set<Short> oneSplit = new HashSet<Short>();
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
			Set<Short> split1 = originalSplits.get(i);
			if (!split1.isEmpty()) {
				HashMap<Set<Short>, Integer> passedSplits = new LinkedHashMap<Set<Short>, Integer>();
				for (int j = i + 1; j < timeTrees.length; j++) {
					Set<Short> split2 = originalSplits.get(j);
					if (!split2.isEmpty()) {
						if (split1.equals(split2)
								|| split1.equals(reverseSplit(split2))) {
							break;
						} else if (!isCompatible(originalSplits.get(i), originalSplits.get(j))) {
							Set<Set<Short>> keySet = passedSplits.keySet();
							List<Set<Short>> keyList = new ArrayList<Set<Short>>();
							keyList.addAll(keySet);
							boolean takenCareOf = false;
							for (int k = keyList.size() - 1; k >= 0; k--) {
								Set<Short> split = keyList.get(k);
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
							Set<Set<Short>> keySet = passedSplits.keySet();
							List<Set<Short>> keyList = new ArrayList<Set<Short>>();
							keyList.addAll(keySet);
							for (int k = keyList.size() - 1; k >= 0; k--) {
								Set<Short> split = keyList.get(k);
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
		System.out.println(" [" + (double) (new Date().getTime() - startDate.getTime()) / 1000
				+ " Seconds] ");
//		System.out.println("(" + counter + ")");
	}

	private int resolveConflict(int i, int j) {
		Set<Short> split1 = originalSplits.get(i);
		Set<Short> split2 = originalSplits.get(j);
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
				Set<Short> split1 = upgmaNode.getChildren().get(0).getSplit();
				Set<Short> split2 = upgmaNode.getChildren().get(1).getSplit();
				Set<Short> tmp;
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
					Set<Short> childSplit = child.getSplit();
					int count = 0;
					for (Short integer : childSplit) {
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
				if (set1.size() > 1) {
					int capacity = 0;
					for (Node node_ : set1) {
						capacity += node_.getSplit().size();
					}
					Node node1 = new Node(capacity);
					node1.getChildren().addAll(set1);
					for (Node n : set1) {
						n.setParent(node1);
						node1.getSplit().addAll(n.getSplit());
					}
					node1.setParent(node);
					tree.getNodes().add(node1);
					tree.getSplits().add(node1.getSplit());
					if (node1.getChildren().size() > 2) {
						nodesToBreak.add(node1);
					}
					node.getChildren().removeAll(set1);
					node.getChildren().add(node1);
				}
				if (set2.size() > 1) {
					int capacity = 0;
					for (Node node_ : set1) {
						capacity += node_.getSplit().size();
					}
					Node node2 = new Node(capacity);
					node2.getChildren().addAll(set2);
					for (Node n : set2) {
						n.setParent(node2);
						node2.getSplit().addAll(n.getSplit());
					}
					node2.setParent(node);
					tree.getNodes().add(node2);
					tree.getSplits().add(node2.getSplit());
					if (node2.getChildren().size() > 2) {
						nodesToBreak.add(node2);
					}
					node.getChildren().removeAll(set2);
					node.getChildren().add(node2);
				}
			}
		}
	}
	
	private void inferBranchLengths() {
		Set<Set<Short>> currentSplit = localTrees[0].getSplits();
		List<int[]> regions = new ArrayList<int[]>();
		int[] region = new int[2];
		region[0] = 0;
		for (int i = 1; i < localTrees.length; i++) {
			if (!localTrees[i].getSplits().equals(currentSplit)) {
				region[1] = i-1;
				regions.add(region);
				region = new int[2];
				region[0] = i;
				currentSplit = localTrees[i].getSplits();
			}
		}
		for (int i = 0; i < timeTreeLengths.length; i++) {
			localTreeLengths[i] = timeTreeLengths[i];
		}
//		scaleTimes(timeTreeLengths);
		region[1] = localTrees.length-1;
		regions.add(region);
		Map<Set<Short>, Double> splitMap = new HashMap<Set<Short>, Double>();
		for (int[] r : regions) {
			Tree tree = localTrees[r[0]];
			splitMap.clear();
			double scaledHeight = 0.0;
			for (int i = r[0]; i <= r[1]; i++) {
				scaledHeight += timeTreeLengths[i];
			}
			scaledHeight = scaledHeight/(r[1]-r[0]+1);
			for (Set<Short> split : tree.getSplits()) {
				double height = 0.0;
				if (split.size() != 1) {
					for (int i = r[0]; i <= r[1]; i++) {
						height += Double.valueOf(timeTrees[i].getBestSplitNode(split).getInfo());
					}
				}
				height = height/(r[1]-r[0]+1);
				splitMap.put(split, height);
			}
			for (int i = r[0]; i <= r[1]; i++) {
				tree = localTrees[i];
				Tree timeTree = timeTrees[i];
				double multiplier = timeTreeLengths[i]/localTreeLengths[i];
				for (Node node : timeTree.getNodes()) {
					double branchLength = Double.valueOf(node.getBranchLength());
					double height = Double.valueOf(node.getInfo());
					node.setBranchLength(String.format( "%.2f",(branchLength)*multiplier));
					node.setInfo(String.format( "%.2f",(height)*multiplier));
				}
				double treeHeight = splitMap.get(overalSplit);
//				multiplier *= localTreeLengths[i]/treeHeight;
				multiplier = scaledHeight/treeHeight;
				for (Node node : tree.getNodes()) {
					double height = 0;
					if (node.getSplit().size() != 1) {
						height = splitMap.get(node.getSplit());
					}
					double parentHeight = 0;
					if (node.getParent() != null){
						parentHeight = splitMap.get(node.getParent().getSplit());
						node.setBranchLength(String.format( "%.2f",(parentHeight > height ? parentHeight - height:0)*multiplier));
					}else {
						node.setBranchLength("0");
						localTreeLengths[i] = height*multiplier;
					}
					if (height > parentHeight){
						height = parentHeight;
					}
					node.setInfo(String.format( "%.2f",height*multiplier));
				}
			}
		}
//		scaleTimes(localTreeLengths);
	}
	

	private void outputTMRCAs(String inputUrl) {
		try {
			BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(
			          new FileOutputStream(inputUrl + ".Tmrcas")));
			for (int i = 0; i < localTreeLengths.length; i++) {
				writer.write(String.valueOf(localTreeLengths[i]));
				writer.newLine();
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
	


	private Tree[] buildLocalTrees() {
		logMemory();
		fullyCompatibleWithCheckNew();
		System.gc();
		logMemory();
		fullyCompatibleRegionRuleNew();
		System.gc();
		logMemory();
		timeSplitOverrideRule();
		System.gc();
		logMemory();
		propagationRuleNew();
		System.gc();
		logMemory();
		propagationHeightRule();
		System.gc();
		logMemory();
		propagationRule();
		System.gc();
		logMemory();
		timeSplitRule();
		System.gc();
		logMemory();
		propagationRule();
		System.gc();
		logMemory();
		RefineNonBinaryNodes();
		System.gc();
		logMemory();
		return localTrees;
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		if (args.length == 0) {
//			String file = "15x20x720x720-ms-1";
			String file = "30x20x720x720-ms-1";
//			String file = "100x20x720x720-ms-1";
//			String file = "200x20x720x720-ms-1";
//			String file = "500x20x720x720-ms-1";
//			String file = "1000x20x720x720-ms-1";

			String url = "res/" + file + ".dat";
			String urlMs = "res/" + file + ".trace";

//			args = new String[] {url};
			args = new String[] {url, urlMs};
//			args = new String[] {"-t", url, urlMs};
//			args = new String[] {"res/test.dat"};
		}
		new Main(args);
	}
}
