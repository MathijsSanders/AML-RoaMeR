package com.erasmusmc.roamer.substratedb;
import java.io.*;
import java.util.*;

public class sdb {
	private int bufferSize = 1;
	private String substrate = null;
	private Queue<Character> buffer = null;
	private Queue<Character> substrateBuffer = null;
	private int totalLength = 0, substrateLength = 0;
	private HashSet<String> excludeContigs = null;
	private boolean skip = false;
	public sdb(parseParameters pp) {
		try {
			BufferedReader br = null;
			BufferedWriter output = null;;
			if(pp.fasta != null)
				br = new BufferedReader(new FileReader(pp.fasta));
			else {
				System.out.println("Please specify the reference FASTA file");
				System.exit(-1);
			}
			if(pp.outputDB != null)
				output = new BufferedWriter(new FileWriter(pp.outputDB));
			else {
				System.out.println("Please specify the output file for the substrate database");
				System.exit(-1);
			}
			if(pp.substrate == null)
			{
				System.out.println("Please specify the mutation substrate");
				System.exit(-1);
			}
			else
			{
				substrate = pp.substrate;
				buffer = new LinkedList<Character>();
				substrateLength = substrate.length();
			}
			if(pp.excludeFile != null)
				excludeContigs = readExcludeContigs(pp.excludeFile);
			this.bufferSize = pp.bufferSize;
			totalLength = substrate.length() + bufferSize;
			String line = null;
			String currentChromosome = "";
			int i,j;
			ArrayList<substrateInfo> pending = new ArrayList<substrateInfo>(100);
			int pos = 1;
			ArrayList<substrateInfo> currentList = null;
			ArrayList<substrateInfo> removeList = null;
			substrateInfo cpgCurrent = null;
			String substrateForward = "";
			String substrateReverse = "";
			while((line = br.readLine()) != null) {
				if(line.startsWith(">")) {
					if(currentList != null)
						writeCurrentList(currentList,output,currentChromosome);
					currentList = null;
					pending = new ArrayList<substrateInfo>(10);
					line = line.replace(">", "");
					pos = 1;
					currentList = new ArrayList<substrateInfo>(50000000);
					buffer = new LinkedList<Character>();
					substrateBuffer = new LinkedList<Character>();
					System.gc();
					currentChromosome = line;
					skip = (excludeContigs != null && excludeContigs.contains(currentChromosome)) ? true : false;
					if(pp.verbose)
						System.out.println("Chromosome: " + currentChromosome + " - Skip contig: " + skip);
				}
				else if(skip)
					continue;
				else {
					line = line.toUpperCase();
					char[] tmp = line.toCharArray();
					for(i = 0;i < tmp.length; i++) {
						if(buffer.size() == totalLength)
							buffer.remove();
						buffer.add(tmp[i]);
						if(substrateBuffer.size() == substrateLength)
							substrateBuffer.remove();
						substrateBuffer.add(tmp[i]);
						substrateForward = substrateReverse = "";
						for(Character c : substrateBuffer)
							substrateForward += c.charValue();
						substrateReverse = reverseComplement(substrateForward);
						if(pending.size() > 0) {
							removeList = new ArrayList<substrateInfo>(pending.size()+5);
							for(j = 0; j < pending.size(); j++) {
								cpgCurrent = pending.get(j);
								if(cpgCurrent.strand.equals("+")) {
									if(cpgCurrent.quadStatus == (bufferSize)) {
										cpgCurrent.succeeding += tmp[i];
										currentList.add(cpgCurrent);
										removeList.add(cpgCurrent);
									}
									else {
										cpgCurrent.succeeding += tmp[i];
										cpgCurrent.quadStatus += 1;
										pending.set(j, cpgCurrent);
									}
								}
								else {
									if(cpgCurrent.quadStatus == (bufferSize)) {
										cpgCurrent.preceding += tmp[i];
										cpgCurrent.preceding = reverseComplement(cpgCurrent.preceding);
										cpgCurrent.succeeding = reverseComplement(cpgCurrent.succeeding);
										currentList.add(cpgCurrent);
										removeList.add(cpgCurrent);
									}
									else {
										cpgCurrent.preceding += tmp[i];
										cpgCurrent.quadStatus +=1;
										pending.set(j, cpgCurrent);
									}
								}
							}
							pending.removeAll(removeList);
						}
						if(substrate.equals(substrateForward)) {
							substrateInfo forward = new substrateInfo();
							forward.pos = pos-(substrateLength-1);
							forward.strand = "+";
							if(buffer.size() == totalLength) {
								for(Character nuc : buffer)
									forward.preceding += nuc;
								forward.preceding = forward.preceding.substring(0, forward.preceding.length()-substrateLength);
								forward.context = substrate.toString();
							}
							else
								continue;
							pending.add(forward);
						}
						if(substrate.equals(substrateReverse)) {
							substrateInfo reverse = new substrateInfo();
							reverse.pos = pos;
							reverse.strand = "-";
							if(buffer.size() == totalLength) {
								for(Character nuc : buffer)
									reverse.succeeding += nuc;
								reverse.succeeding = reverse.succeeding.substring(0, reverse.succeeding.length()-substrateLength);
								reverse.context = substrate.toString();
							}
							pending.add(reverse);
						}
						pos++;
					}
				}
			}
			br.close();
			writeCurrentList(currentList,output,currentChromosome);
			output.close();
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		
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
	private void writeCurrentList(ArrayList<substrateInfo> currentList,BufferedWriter output, String chr) {
		try {
			StringBuffer buffer = new StringBuffer();
			int counter = 0;
			for(substrateInfo curr : currentList) {
				buffer.append(chr + "\t" + curr.pos + "\t"  + (curr.pos+1) + "\t" + curr.strand + "\t" + curr.preceding + "\t" + curr.context + "\t" + curr.succeeding + "\n");
				counter++;
				if(counter == 1000) {
					output.write(buffer.toString());
					output.flush();
					buffer = new StringBuffer();
					counter = 0;
				}
			}
			if(counter > 0) {
				output.write(buffer.toString());
				output.flush();
			}
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	private String reverseComplement(String seq) {
		String tmp = "";
		char a = '1';
		for(int i = (seq.length()-1); i >= 0; i--) {
			a = seq.charAt(i);
			if(a == 'A')
				tmp += "T";
			else if(a == 'T')
				tmp += "A";
			else if(a == 'G')
				tmp += "C";
			else if(a == 'C')
				tmp += "G";
			else
				tmp += "N";
		}
		return tmp;
	}
}
