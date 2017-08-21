package com.erasmusmc.roamer.extractmethylation;
import java.io.*;
import java.util.*;
import org.broad.igv.bbfile.*;
public class em  {
	private boolean skip = false;
	private HashSet<String> contextSet = new HashSet<String>(250);
	private char space = 32;
	public em(parseParameters pp) {
		String[] coverageFiles = null;
		String[] methFiles = null;
		HashSet<String> excludeContigs = null;
		HashMap<String,HashMap<Integer,elementInfo>> ocMap = null;
		if(pp.covFiles != null)
			coverageFiles = pp.covFiles.split(",");
		else {
			System.out.println("Please specify the coverage bigWig files derived from the methylation data");
			System.exit(-1);
		}
		if(pp.methFiles != null)
			methFiles = pp.methFiles.split(",");
		else {
			System.out.println("Please specify the methylation bigWig files derived from the methylation data");
			System.exit(-1);
		}
		int i,k;
		HashMap<Integer,cpgInfo> currentList = null;
		int preceding = pp.preceding;
		int succeeding = pp.succeeding;
		if(pp.database != null)
			ocMap = readCpgDatabase(pp.database, preceding, succeeding, pp);
		else {
			System.out.println("Please specify the CG database produced with the substrateDatabase command");
			System.exit(-1);
		}
		if(pp.excludeContigs != null) {
			try {
				excludeContigs = readExcludeContigs(pp.excludeContigs);
			} catch(IOException e) {
				e.printStackTrace();
			}
		}
		if(pp.outputMethylation == null) {
			System.out.println("Please specify the output file for the 5mCG database");
			System.exit(-1);
		}
		cpgInfo cpgCurrent = null;
		BBFileReader[] covReader = new BBFileReader[coverageFiles.length];
		BBFileReader[] methReader = new BBFileReader[methFiles.length];
		BigWigIterator[] itCov = new BigWigIterator[coverageFiles.length];
		BigWigIterator[] itMeth = new BigWigIterator[methFiles.length];
		WigItem[] currentCov = new WigItem[coverageFiles.length];
		WigItem[] currentMeth = new WigItem[methFiles.length];
		BBFileHeader header = null;
		String oldString = "";
		boolean first = true;
		boolean[] boolChrom = new boolean[coverageFiles.length];
		try {
			for(i = 0; i < coverageFiles.length; i++) {
				covReader[i] = new BBFileReader(coverageFiles[i]);
				header = covReader[i].getBBFileHeader();
				if(!header.isHeaderOK()) {
					System.out.println("Header of " + coverageFiles[i]+ " is not good");
					System.exit(-1);
				}
				itCov[i] = covReader[i].getBigWigIterator();
				methReader[i] = new BBFileReader(methFiles[i]);
				header = methReader[i].getBBFileHeader();
				if(!header.isHeaderOK()) {
					System.out.println("Header of " + methFiles[i]+ " is not good");
					System.exit(-1);
				}
				itMeth[i] = methReader[i].getBigWigIterator();
			}
			while(hasNext(itCov)) {
				for(i = 0; i < covReader.length; i++) {
					if(!boolChrom[i]) {
						if(itCov[i].hasNext()) {
							currentCov[i] = itCov[i].next();
							currentMeth[i] = itMeth[i].next();
							if(!oldString.equals(currentCov[i].getChromosome())) {
								skip = (excludeContigs != null && excludeContigs.contains(currentCov[i].getChromosome())) ? true : false;
								if(!skip && pp.verbose)
									System.out.println("Currently processing chromosome for file " + i +  ": " + currentCov[i].getChromosome() + " - Skip chromosome: " + skip);
								if(skip)
									continue;
								else
									boolChrom[i] = true;
							} else {
								if(!currentList.containsKey(currentCov[i].getStartBase())) {
									cpgInfo cpg = new cpgInfo();
									cpg.start = currentCov[i].getStartBase();
									cpg.reads = currentCov[i].getWigValue();
									cpg.methreads = Math.round(cpg.reads * (currentMeth[i].getWigValue()/100.0));
									currentList.put(currentCov[i].getStartBase(),cpg);
								} else {
									cpgCurrent = currentList.get(currentCov[i].getStartBase());
									cpgCurrent.reads += currentCov[i].getWigValue();
									cpgCurrent.methreads += Math.round(currentCov[i].getWigValue() * (currentMeth[i].getWigValue()/100.0));
									currentList.put(currentCov[i].getStartBase(), cpgCurrent);
								}
							}
						}
					}
				}
				if(allTrue(boolChrom)) {
					if(!sameChromosome(currentCov)) {
						System.out.println("Error: The bigWig files are not sorted in the same order!");
						System.exit(-1);
					}
					if(currentList != null) {
						writeResults(currentList, pp.outputCoverage, pp.outputMethylation, ocMap, first, pp.verbose, oldString);
						first = false;
					}
					oldString = currentCov[0].getChromosome();
					currentList = new HashMap<Integer,cpgInfo>(5000000,0.9999f);
					for(k = 0; k < currentCov.length; k++) {
						boolChrom[k] = false;
						if(!currentList.containsKey(currentCov[k].getStartBase())) {
							cpgInfo cpg = new cpgInfo();
							cpg.start = currentCov[k].getStartBase();
							cpg.reads = currentCov[k].getWigValue();
							cpg.methreads = Math.round(cpg.reads * (currentMeth[k].getWigValue()/100.0));
							currentList.put(currentCov[k].getStartBase(),cpg);
						} else {
							cpgCurrent = currentList.get(currentCov[k].getStartBase());
							cpgCurrent.reads += currentCov[k].getWigValue();
							cpgCurrent.methreads += Math.round(currentCov[k].getWigValue() * (currentMeth[k].getWigValue()/100.0));
							currentList.put(currentCov[k].getStartBase(), cpgCurrent);
						}
					}
				}
			}
			if(currentList != null)
				writeResults(currentList, pp.outputCoverage, pp.outputMethylation, ocMap, first, pp.verbose, oldString);
		} catch(IOException e) {
			e.printStackTrace();
		}
	}
	private void writeResults(HashMap<Integer,cpgInfo> currentList, String outputCoverageFile, String outputMethFile, HashMap<String,HashMap<Integer,elementInfo>> map, boolean first, boolean verbose, String currentChromosome) throws IOException {
		cpgInfo cpgCurrent = null;
		BufferedWriter outputCoverage = null;
		StringBuffer bufferCoverage = null;
		if(outputCoverageFile != null) {
			outputCoverage = (first) ? new BufferedWriter(new FileWriter(outputCoverageFile)) : new BufferedWriter(new FileWriter(outputCoverageFile,true));
			bufferCoverage = new StringBuffer();
		}
		BufferedWriter outputMeth =  (first) ? new BufferedWriter(new FileWriter(outputMethFile)) : new BufferedWriter(new FileWriter(outputMethFile,true));
		StringBuffer bufferMeth = new StringBuffer();
		int counter = 0;
		String context = "";
		elementInfo ei = null;
		HashMap<Integer,elementInfo> currentMap = null;
		if(verbose)
			System.out.println("Currently outputting chromosome: " + currentChromosome);
		currentMap = map.get(currentChromosome);
		SortedSet<Integer> keysChrom = new TreeSet<Integer>(currentList.keySet());
		for(Integer currentLoc : keysChrom) {
			cpgCurrent = currentList.get(currentLoc);
			ei = currentMap.get(cpgCurrent.start+1);
			if(ei == null)
				continue;
			context = ei.context;
			if(outputCoverageFile != null)
				bufferCoverage.append(currentChromosome + "\t" + cpgCurrent.start + "\t" + context + "\t" + cpgCurrent.reads + "\t" + ((ei.posStrand) ? "+" : "-") + "\n");
			bufferMeth.append(currentChromosome + "\t" + cpgCurrent.start + "\t" + context + "\t" + ((cpgCurrent.methreads/cpgCurrent.reads)*100.0) + "\t" + ((ei.posStrand) ? "+" : "-") + "\n");
			counter++;
			if(counter == 1000) {
				if(outputCoverageFile != null) {
					outputCoverage.write(bufferCoverage.toString());
					outputCoverage.flush();
					bufferCoverage = new StringBuffer();
				}
				outputMeth.write(bufferMeth.toString());
				outputMeth.flush();
				bufferMeth = new StringBuffer();
				counter = 0;
			}
		}
		if(counter > 0) {
			if(outputCoverageFile != null) {
				outputCoverage.write(bufferCoverage.toString());
				outputCoverage.flush();
				bufferCoverage = new StringBuffer();
			}
			outputMeth.write(bufferMeth.toString());
			outputMeth.flush();
			bufferMeth = new StringBuffer();
			counter = 0;
		}
		if(outputCoverageFile != null)
			outputCoverage.close();
		outputMeth.close();
	}
	private boolean sameChromosome(WigItem[] cov) {
		String oldString = null;
		for(WigItem w : cov) {
			if(oldString == null)
				oldString = w.getChromosome();
			else
				if(!oldString.equals(w.getChromosome())) return false;
		}
		return true;
	}
	private boolean allTrue(boolean[] boolArray) {
		for(boolean b : boolArray) if(!b) return false;
		return true;
	}
	private boolean hasNext(BigWigIterator[] itArray) {
		for(BigWigIterator it : itArray)
			if(it.hasNext())
				return true;
		return false;
	}
	private void writeByteFile(String file, HashMap<String,HashMap<Integer,cpgInfo>> map) throws IOException {
		ByteArrayOutputStream byteBegin = new ByteArrayOutputStream(); 
		DataOutputStream dataBegin = new DataOutputStream(byteBegin);
		Collection<String> chromosomes = map.keySet();
		int maxChrom = -1;
		for(String c : chromosomes)
			if(c.length() > maxChrom) maxChrom = c.length();
		dataBegin.writeInt(chromosomes.size());
		dataBegin.writeInt(contextSet.size());
		dataBegin.writeInt(maxChrom);
		String[] chrom = new String[chromosomes.size()];
		chromosomes.toArray(chrom);
		chrom = convertList(chrom, maxChrom);
		for(String c : chrom)
			dataBegin.writeBytes(c);
		
		
	}
	private String[] convertList(String[] list, int max) {
		String[] temp = new String[list.length];
		for(int i =0; i < temp.length; i++) {
			if(list[i].length() < max) {
				char[] replace = new char[max];
				int bound = list[i].length();
				for(int j = 0; j < max; j++) {
					if(j < bound)
						replace[j] = list[i].charAt(j);
					else
						replace[j] = space;
				}
				temp[i] = new String(replace);
			}
			else 
				temp[i] = list[i];
		}
		return temp;
	}
	private HashSet<String> readExcludeContigs(String file) throws IOException {
		HashSet<String> set = new HashSet<String>(100);
		BufferedReader br = new BufferedReader(new FileReader(file));
		String line = null;
		while((line = br.readLine()) != null)
			set.add(line);
		br.close();
		return set;
	}
	private HashMap<String,HashMap<Integer,elementInfo>> readCpgDatabase(String file, int preceding, int succeeding, parseParameters pp) {
		HashMap<String,HashMap<Integer,elementInfo>> map = new HashMap<String,HashMap<Integer,elementInfo>>(50,0.9999f);
		try {
			BufferedReader br = new BufferedReader(new FileReader(file));
			String[] tokens = null;
			String line = null;
			String currentChromosome = "";
			HashMap<Integer,elementInfo> currentMap = null;
			while((line = br.readLine()) != null) {
				tokens = line.split("\t");
				tokens[0] = tokens[0].replaceFirst("chr", "");
				if(!tokens[0].equals(currentChromosome)) {
					if(pp.verbose)
						System.out.println("Currently loading chromosome " + tokens[0] + " from the CG-database");
					currentChromosome = tokens[0];
					if(!map.containsKey(currentChromosome))
						map.put(currentChromosome, new HashMap<Integer,elementInfo>(10000000));
					currentMap = map.get(tokens[0]);
				}
				elementInfo ei = new elementInfo();
				ei.context = tokens[4].substring(tokens[4].length()-preceding, tokens[4].length()) + tokens[5] + tokens[6].substring(0, succeeding);
				if(!contextSet.contains(ei.context))
					contextSet.add(ei.context);
				ei.posStrand = tokens[3].equals("+");
				currentMap.put(Integer.parseInt(tokens[1]), ei);
			}
			br.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		return map;
	}
}