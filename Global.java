import java.util.ArrayList;
public class Global {
	public static int GAP = -2;
	public static int DIF = -1;
	public static int SAME = 1;
	
	public static ArrayList<ArrayList<Integer>> matrix = new ArrayList<ArrayList<Integer>>();
	public static ArrayList<ArrayList<ArrayList<ArrayList<String>>>> strings = new ArrayList<ArrayList<ArrayList<ArrayList<String>>>>();
	public static String seq1 = "AAG";
	public static String seq2 = "AAC";

	public static void alignment() {
		for(int i = 0; i <= seq2.length(); i++) {
			matrix.add(new ArrayList<Integer>()) ;
		//	strings.add(new ArrayList<ArrayList<ArrayList<String>>>());
			for( int j = 0; j <= seq1.length(); j++) {
			//	strings.get(i).add(new ArrayList<ArrayList<String>>());
			//	ArrayList<ArrayList<String>> possibs = strings.get(i).get(j);
				int val;
				if(i == 0 & j == 0) {
					val = 0;
				/*	possibs.add(new ArrayList<String>());
					possibs.get(0).add("");
					possibs.get(0).add("");*/
				}else if(i == 0) {
					val = matrix.get(i).get(j-1) + GAP;
				/*	possibs = strings.get(i).get(j-1);
					System.out.println(j);
					possibs.get(0).set(0, seq1.substring(0, j));
					possibs.get(0).set(1, "_");*/
				}else if(j == 0) {
					val = matrix.get(i-1).get(j) + GAP;
					/*possibs = strings.get(i-1).get(j);
					possibs.get(0).set(0, "_");
					possibs.get(0).set(1, seq2.substring(0, i));*/
				}
				else {
					int vi, vj, vd;
					val = Integer.max(Integer.max(vd = matrix.get(i-1).get(j-1) + s(i,j),
							vj = matrix.get(i).get(j-1) + GAP),
							vi = matrix.get(i-1).get(j) + GAP);
					/*if(val == vd) {
						for (ArrayList<String> s : strings.get(i-1).get(j-1)) {
							ArrayList<String> temp = new ArrayList<String>();
							temp.add(s.get(0) + seq1.substring(j-1, j));
							temp.add(s.get(1) + seq2.substring(i-1,i));
							strings.get(i).get(j).add(s);
						}
					}
					if(val == vj) {
						for (ArrayList<String> s : strings.get(i).get(j-1)) {
							ArrayList<String> temp = new ArrayList<String>();
							temp.add(s.get(0) + seq1.substring(j-1, j));
							temp.add(s.get(1) + "_");
							strings.get(i).get(j).add(s);
						}
					}
					if(val == vi) {
						for (ArrayList<String> s : strings.get(i-1).get(j)) {
							ArrayList<String> temp = new ArrayList<String>();
							temp.add(s.get(0) + "_");
							temp.add(s.get(1) + seq2.substring(i-1,i));
							strings.get(i).get(j).add(s);
						}
					}*/
				}
				matrix.get(i).add(val);
			}
		}
		
		print(matrix);
		
/*		System.out.println("Result:");
		for (ArrayList<String> s : strings.get(seq2.length()).get(seq1.length())) {
			System.out.println("-------------");
			System.out.println(s.get(0));
			System.out.println(s.get(1));
		}*/
	}

	private static void print(ArrayList<ArrayList<Integer>> matrix2) {
		for(int i = 0; i < matrix.size(); i++) {
			for(int j = 0; j < matrix.get(i).size(); j++) {
				System.out.print(matrix.get(i).get(j) + "  ");
			}
			System.out.println("");
		}
		
	}

	private static Integer s(int i, int j) {
		System.out.println(seq2.substring(i-1, i) + " " + seq1.substring(j-1, j));
		if (seq2.substring(i-1, i).equals(seq1.substring(j-1, j)))
			return SAME;
		return DIF;
	}
	
	

}
