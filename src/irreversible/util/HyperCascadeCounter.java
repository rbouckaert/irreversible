package irreversible.util;

import java.util.Arrays;

import beast.app.util.Application;
import beast.core.Description;
import beast.core.Input;
import beast.core.Runnable;


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
20		74		5926246929		32.5
21		78		19242396629		34.2
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

	int layers, sites;
	// start position of layers in memory
	int [] layerOffset;
	
	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		layers = layersInput.get();
		sites = sitesInput.get();
		if (layers > sites) {
			throw new IllegalArgumentException("nuber of layers (" + layers + ") should not be "
					+ "larger than number of sites (" + sites + ")");
		}

		// calc total number of sites in system
		int totalSites = 0;
		layerOffset = new int[layers];
		for (int i = 0; i < layers; i++) {
			layerOffset[i] = totalSites;
			totalSites += sites - i;
		}

		int [] mem = new int[totalSites]; 
		long states = countStates(mem, 0);
		System.out.println("#states = " + states +
				" #bits = "  + (Math.log(states)/Math.log(2)) +
				" #sites = " + totalSites);
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
				count += countStates(mem, 1) + 1;
				// print(mem);
			}
		} else {
			boolean [] canMutate = new boolean[layerSize];
			int [] mutateMap = new int[layerSize];
			int m = 0;
			int offset = layerOffset[layer];
			int parentOffset = layerOffset[layer-1];
			for (int i = 0; i < layerSize; i++) {
				if (mem[parentOffset + i] == 1 &&  mem[parentOffset + i + 1] == 1) {
					canMutate[i] = true;
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
				if (layer +1 < layers) {
					count += countStates(mem, layer + 1);
				}
				// print(mem);
			}

 		}		
		
		return count;
	}
	
	
	private void print(int[] buf) {
		for (int i = 0; i < buf.length; i++) {
			System.err.print(buf[i]);
		}
		System.err.println();
		
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

	public static void main(String[] args) throws Exception {
		new Application(new  HyperCascadeCounter(), "Hyper Cascade Counter", args);

	}

}
