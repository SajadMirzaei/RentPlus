package object;

import java.util.ArrayList;
import java.util.List;

import main.Main;

public class DistanceHolder {
//	private List<List<float[]>> distanceMatrices = new ArrayList<List<float[]>>();
	private List<float[]> distanceMatrices = new ArrayList<float[]>();
	private int totalSize = (Main.NUM_TAXA)*(Main.NUM_TAXA-1)/2;
	public DistanceHolder(int size) {
		totalSize = size*(size-1)/2;
	}
//	public void addDistance(int i, int j, int position, float value){
//		if (j < i) {
//			// swapping i and j
//			i += j;
//			j = i - j;
//			i -= j;
//		}
//		int index = totalSize - (Main.NUM_TAXA-i-1)*(Main.NUM_TAXA-i)/2;
//		index = index + j - i - 1;
//		if (distanceMatrices.size() <= index) {
//			distanceMatrices.add(new ArrayList<float[]>());
//		}
//		List<float[]> relatedList = distanceMatrices.get(index);
//		if (relatedList.isEmpty()) {
//			relatedList.add(new float[]{position, position, value});
//		}else{
//			if (relatedList.get(relatedList.size()-1)[2] == value) {
//				relatedList.get(relatedList.size()-1)[1] = position;
//			}else{
//				relatedList.add(new float[]{position, position, value});
//			}
//		}
//	}
	
	public void addDistance(int i, int j, int position, float value){
		if (j < i) {
			// swapping i and j
			i += j;
			j = i - j;
			i -= j;
		}
		int index = totalSize - (Main.NUM_TAXA-i-1)*(Main.NUM_TAXA-i)/2;
		index = index + j - i - 1;
		if (distanceMatrices.size() <= index) {
			distanceMatrices.add(new float[Main.positions.size()]);
		}
		distanceMatrices.get(index)[position] = value;
	}
	
//	public float getdistance(int i, int j, int position){
//		if (j < i) {
//			// swapping i and j
//			i += j;
//			j = i - j;
//			i -= j;
//		}
//		int index = totalSize - (Main.NUM_TAXA-i-1)*(Main.NUM_TAXA-i)/2;
//		index = index + j - i - 1;
//		List<float[]> relatedList = distanceMatrices.get(index);
//		int min = 0;
//		int max = relatedList.size()-1;
//		while(true){
//			int k = (max+min)/2;
//			if (position < relatedList.get(k)[0]) {
//				max = k-1;
//			}else if (position > relatedList.get(k)[1]) {
//				min = k+1;
//			}else{
//				return relatedList.get(k)[2];
//			}
//		}
//	}
	
	public float getdistance(int i, int j, int position){
		if (j < i) {
			// swapping i and j
			i += j;
			j = i - j;
			i -= j;
		}
		int index = totalSize - (Main.NUM_TAXA-i-1)*(Main.NUM_TAXA-i)/2;
		index = index + j - i - 1;
		return distanceMatrices.get(index)[position];
	}
	
	public float[][] getDistanceMatrixforPosition(int position){
		float[][] distanceMatrix = new float[Main.NUM_TAXA][Main.NUM_TAXA];
		for (int i = 0; i < Main.NUM_TAXA; i++) {
			for (int j = i+1; j < Main.NUM_TAXA; j++) {
				float value = getdistance(i, j, position);
				distanceMatrix[i][j] = value;
				distanceMatrix[j][i] = value;
			}
		}
		return distanceMatrix;
	}
}
