package irreversible.util;

import java.util.*;

import beast.app.util.Application;
import beast.core.Description;
import beast.core.Input;
import beast.core.Runnable;
import beast.core.util.Log;


/*
 for 4 layers = Number of Catalan paths (nonnegative, starting and ending at 0, step +/-1) of 2*n steps with all values <= 5.
 https://oeis.org/A080937
  
1st lyr total   pattern         bits information
sites   sites   count           =log2(patterns)
 4		10		42				 5.4
 5		14		131				 7.0
 6		18		417				 8.7
 7 		22		1341			10.4
 8		26		4334			12.1
 9		30		14041			13.8
10		34		45542			15.5
11		38		147798			17.2
12		42		479779			18.9
13		46		1557649			20.6
14		50		5057369			22.3
15		54		16420730		24.0
16		58		53317085		25.7
17		62		173118414		27.4
18		66		562110290		29.1
19		70		1825158051		30.8
20		74		5926246929		32.5 verified up to here
21		78		19242396629		34.2 interpolated from here on
22		82		62479659622		35.9
23		86		202870165265	37.6
24		90		658715265222	39.3
25		94		2138834994142	41.0
26		98		6944753544643	42.7
27		102		22549473023585	44.4
 */


@Description("Counts number of states for a size of hyper cascade system")
public class HyperCascadeCounter extends Runnable {
	public Input<Integer> layersInput = new Input<>("layers", "number of layers in system", 4);
	public Input<Integer> sitesInput = new Input<>("sites", "number of sites of biggest (first) layer in system", 4);
	public Input<Boolean> verboseInput = new Input<>("verbose", "print out all states visited", false);

	int layers, sites;
	// start position of layers in memory
	int [] layerOffset;
	int [][] mutateMap;
	List<String> states;
	boolean verbose;

	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		long start = System.currentTimeMillis();
		layers = layersInput.get();
		sites = sitesInput.get();
		mutateMap = new int[layers][sites];
		states = new ArrayList<>();
		verbose = verboseInput.get();
		
		if (layers > sites) {
			throw new IllegalArgumentException("nuber of layers (" + layers + ") should not be "
					+ "larger than number of sites (" + sites + ")");
		}

		// calc total number of sites in system
		int totalSites = 0;
		layerOffset = new int[layers + 1];
		for (int i = 0; i < layers; i++) {
			layerOffset[i] = totalSites;
			totalSites += sites - i;
		}

		int [] mem = new int[totalSites]; 
		long states = countStates(mem, 0);
		
		if (verbose) {
			for (String s : this.states) {
				Log.warning(s);
			}
		}
		System.out.println("#states = " + states +
				" #bits = "  + (Math.log(states)/Math.log(2)) +
				" #sites = " + totalSites);
		long end = System.currentTimeMillis();
		Log.warning("Done in " + ((end - start)/1000) + " seconds");
	}

	private long countStates(int[] mem, int layer) {
		long count = 0;
		int layerSize = sites - layer;
		if (layer == 0) {
			// traverse all possible states
			int n = 1 << layerSize;
			for (int i = 0; i < n; i++) {
				Arrays.fill(mem, 0);
				toBinary(i, mem, layerSize);
				count++;
				if (verbose) print(mem);
				count += countStates(mem, 1);
			}
		} else {
			int [] mutateMap = this.mutateMap[layer];
			int m = 0;
			int offset = layerOffset[layer];
			int parentOffset = layerOffset[layer-1];
			for (int i = 0; i < layerSize; i++) {
				if (mem[parentOffset + i] == 1 &&  mem[parentOffset + i + 1] == 1) {
					mutateMap[m] = i;
					m++;
				}
			}
			
			if (m == 0) {
				return 0;
			}
			
			for (int i = offset; i < mem.length; i++) {
				mem[i] = 0;
			}
			
			int n = 1 << m;
			int [] buf = new int[m];
			for (int i = 1; i < n; i++) {
				toBinary(i, buf, m);
				for (int j = 0; j < m; j++) {
					mem[offset + mutateMap[j]] = buf[j];
				}
				count++;
				if (verbose) print(mem);
				if (layer +1 < layers) {
					count += countStates(mem, layer + 1);
				}
			}

			for (int i = offset; i < mem.length; i++) {
				mem[i] = 0;
			}
 		}		
		
		return count;
	}
	
	
	private void print(int[] buf) {
		int j = 1;
		StringBuilder b = new StringBuilder();
		for (int i = 0; i < buf.length; i++) {
			if (i == layerOffset[j]) {
				b.append(' ');
				j++;
			}
			b.append(buf[i] + "");
		}
		
		// System.err.println(b);
		states.add(b.toString());
	}

	/**
     * Convert the integer to an unsigned number.
     */
    private void toBinary(int val, int[] buf, int len) {
        // assert shift > 0 && shift <=5 : "Illegal shift value";
        int mag = Integer.SIZE - Integer.numberOfLeadingZeros(val);

        int charPos = 1;
        while (charPos <= mag) {
            buf[len - charPos++] = val & 1;
            val >>>= 1;
        };
    }

	public static List<String> getStates(Integer n) {
		HyperCascadeCounter h = new HyperCascadeCounter();
		h.initByName("layers", n, "sites", n, "verbose", true);
		try {
			h.run();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return h.states;
	}

	
	public static void main(String[] args) throws Exception {
		new Application(new  HyperCascadeCounter(), "Hyper Cascade Counter", args);

	}


}
