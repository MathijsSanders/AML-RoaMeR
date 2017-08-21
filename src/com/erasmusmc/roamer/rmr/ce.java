package com.erasmusmc.roamer.rmr;
import java.io.*;
import java.util.*;
import com.erasmusmc.roamer.intervaltree.*;

public class ce {
	public ce(parseParameters pp) {
		try {
			HashMap<String,HashMap<String,IntervalST<segment>>> bedIntervals = null;
			int totalCases = 0;
			if(pp.vcfs == null) {
				System.out.println("Please provide VCF files to match against the database-corrected intervals!");
				System.exit(-1);
			}
			else
				totalCases = pp.vcfs.split(",").length;
			if(pp.bedFiles != null && pp.bedNameFiles != null)
				bedIntervals = readIntervals(pp.bedFiles,pp.bedNameFiles, pp.stranded, totalCases);
			else {
				System.out.println("BED files and surrogate names have to be provided for the algorithm to initiate!");
				System.exit(-1);
			}
			HashMap<String, HashMap<Integer,String>> db = null;
			if(pp.database != null)
				db = readDatabase(pp.database, bedIntervals, pp.stranded, pp.correctingFactors, pp.correctingFactorNames, pp.verbose);
			else {
				System.out.println("Please provide a database file in the appropriate format!");
				System.exit(-1);
			}
			if(pp.outputStats == null) {
				System.out.println("Please provide an output file location for saving RMR scores and statistics");
				System.exit(-1);
			}
			if(pp.sampleNames == null) {
				System.out.println("Please provide the sample names in a comma-separated format");
				System.exit(-1);
			}
			runEnrichment(pp.vcfs, pp.sampleNames, bedIntervals, db, pp.bedNameFiles, pp.outputContext, pp.outputStats, pp.stranded, pp.verbose);
		}
		catch(Exception e) {
			e.printStackTrace();
		}
	}
	private void runEnrichment(String vcfFiles, String sampleNameList, HashMap<String,HashMap<String,IntervalST<segment>>> intervals, HashMap<String, HashMap<Integer,String>> db, String bedNames, String outputContext, String outputStats, boolean stranded, boolean verbose) throws IOException {
		String[] files = vcfFiles.split(",");
		String[] names = bedNames.split(",");
		double[] stat = null;
		HashSet<String> contextSet = new HashSet<String>(50,0.99999f);
		BufferedWriter statOutput = new BufferedWriter(new FileWriter(outputStats));
		String[] sampleNames = sampleNameList.split(",");
		if(stranded) {
			statOutput.write("BEDfile");
			for(String sample : sampleNames)
				statOutput.write("\tTotal_mutation_burden_" + sample);
			statOutput.write("\tCG-abundance_same\t5mCG-abundance_same");
			for(String sample : sampleNames)
				statOutput.write("\tMutationCount_same_" + sample);
			for(String sample : sampleNames)
				statOutput.write("\tRMR_CG_same_" + sample);
			statOutput.write("\tRMR_CG_average_same\tRMR_CG_sd_same");
			for(String sample : sampleNames)
				statOutput.write("\tRMR_5mCG_same_" + sample);
			statOutput.write("\tRMR_5mCG_average_same\tRMR_5mCG_sd_same");
			statOutput.write("\tCG-abundance_diff\t5mCG-abundance_diff");
			for(String sample : sampleNames)
				statOutput.write("\tMutationCount_diff_" + sample);
			for(String sample : sampleNames)
				statOutput.write("\tRMR_CG_diff_" + sample);
			statOutput.write("\tRMR_CG_average_diff\tRMR_CG_sd_diff");
			for(String sample : sampleNames)
				statOutput.write("\tRMR_5mCG_diff_" + sample);
			statOutput.write("\tRMR_5mCG_average_diff\tRMR_5mCG_sd_diff");
		} else {
			statOutput.write("BEDfile");
			for(String sample : sampleNames)
				statOutput.write("\tTotal_mutation_burden_" + sample);
			statOutput.write("\tCG-abundance\t5mCG-abundance");
			for(String sample : sampleNames)
				statOutput.write("\tMutationCount_" + sample);
			for(String sample : sampleNames)
				statOutput.write("\tRMR_CG_" + sample);
			statOutput.write("\tRMR_CG_average\tRMR_CG_sd");
			for(String sample : sampleNames)
				statOutput.write("\tRMR_5mCG_" + sample);
			statOutput.write("\tRMR_5mCG_average\tRMR_5mCG_sd");
		}
		statOutput.write("\n");
		int i;
		String line = null;
		double[] totalMutation = new double[files.length];
		for(i = 0; i < files.length; i++) {
			BufferedReader br = new BufferedReader(new FileReader(files[i]));
			line = null;
			String[] tokens = null;
			String currentChromosome = "";
			ArrayList<IntervalST<segment>> currentTrees = null;
			IntervalST<segment> currentTree = null;
			Iterable<Interval1D> segmentSet = null;
			Iterator<Interval1D> intervalIterator = null;
			HashMap<Integer,String> currentMap = null;
			segment currentSegment = null;
			int start = -1, end = -1;
			String context = null;
			String refChar = "";
			HashMap<String,Integer> contextMutationCount = null;
			while((line = br.readLine()) != null) {
				if(line.startsWith("#"))
					continue;
				tokens = line.split("\t");
				if(tokens[3].length() > 1 || tokens[4].length() > 1)
					continue;
				tokens[0] = tokens[0].replaceFirst("chr", "");
				start = Integer.parseInt(tokens[1]);
				end = start + 1;
				refChar = tokens[3];
				if(!currentChromosome.equals(tokens[0])) {
					if(verbose)
						System.out.println("Current processing: " + tokens[0] + " - For sample: " + sampleNames[i]);
					currentChromosome = tokens[0];
					currentTrees = new ArrayList<IntervalST<segment>>(50);
					for(String bed : intervals.keySet()) {
						currentTree = intervals.get(bed).get(currentChromosome);
						if(currentTree == null)
							continue;
						else
							currentTrees.add(currentTree);
					}
					currentMap = db.get(currentChromosome);
				}
				if(currentMap == null)
					continue;
				context = currentMap.get(start - 1);
				if(context == null)
					continue;
				else {
					if(!contextSet.contains(context))
						contextSet.add(context);
					totalMutation[i] += 1.0;
				}
				for(IntervalST<segment> curr : currentTrees) {
					segmentSet = curr.searchAll(new Interval1D(start,end));
					if(size(segmentSet) > 0) {
						intervalIterator = segmentSet.iterator();
						while(intervalIterator.hasNext()) {
							currentSegment = curr.get(intervalIterator.next());
							if(currentSegment == null) {
								System.out.println("Error at line: " + line);
								System.exit(-1);
							} else if(currentSegment.mutationsSame == null) {
								System.out.println("Error at line: " + line);
								System.exit(-1);
							}
							if(stranded) {
								if((currentSegment.posStrand && refChar.equals("C")) || (!currentSegment.posStrand && refChar.equals("G"))) {
									currentSegment.mutationsSame[i] += 1.0;
									contextMutationCount = currentSegment.contextMutationCountSame.get(i);
								} else {
									currentSegment.mutationsDiff[i] += 1.0;
									contextMutationCount = currentSegment.contextMutationCountDiff.get(i);
								}
							} else {
								currentSegment.mutationsSame[i] += 1.0;
								contextMutationCount = currentSegment.contextMutationCountSame.get(i);
							}
							if(!contextMutationCount.containsKey(context))
								contextMutationCount.put(context, 1);
							else
								contextMutationCount.put(context, contextMutationCount.get(context)+1);
						}
					}
				}
			}
			br.close();
			totalMutation[i] /= 1000.0;
		}
		HashMap<String, IntervalST<segment>> chrMap = null;
		IntervalST<segment> currentTree = null;
		segment currentSegment = null;
		int j;
		List<String> contextList = new ArrayList<String>(contextSet);
		Collections.sort(contextList);
		Integer contextCount = -1;
		Double contextCountDouble = -1.0;
		BufferedWriter outputCon = null;
		if(outputContext != null) {
			outputCon = new BufferedWriter(new FileWriter(outputContext));
			outputCon.write("BEDfile");
			for(String sample : sampleNames)
				outputCon.write("\tTotal_mutation_" + sample);
			if(stranded) {
				for(String con : contextList) {
					outputCon.write("\tCG-abundance_" + con + "_same\t5mCG-abundance_" + con + "_same");
					for(String sample : sampleNames)
						outputCon.write("\tMutationCount_" + con + "_same_" + sample);
					for(String sample : sampleNames)
						outputCon.write("\tRMR_CG_" + con + "_same_" + sample);
					outputCon.write("\tRMR_CG_" + con + "same_average\tRMR_CG_" + con + "same_sd");
					for(String sample : sampleNames)
						outputCon.write("\tRMR_5mCG_" + con + "_same_" + sample);
					outputCon.write("\tRMR_5mCG_" + con + "same_average\tRMR_5mCG_" + con + "same_sd");
					outputCon.write("\tCG-abundance_" + con + "_diff\t5mCG-abundance_" + con + "_diff");
					for(String sample : sampleNames)
						outputCon.write("\tMutationCount_" + con + "_diff_" + sample);
					for(String sample : sampleNames)
						outputCon.write("\tRMR_CG_" + con + "_diff_" + sample);
					outputCon.write("\tRMR_CG_" + con + "diff_average\tRMR_CG_" + con + "diff_sd");
					for(String sample : sampleNames)
						outputCon.write("\tRMR_5mCG_" + con + "_diff_" + sample);
					outputCon.write("\tRMR_5mCG_" + con + "diff_average\tRMR_5mCG_" + con + "diff_sd");
				}
			} else {
				for(String con : contextList) {
					outputCon.write("\tCG-abundance_" + con + "\t5mCG-abundance_" + con);
					for(String sample : sampleNames)
						outputCon.write("\tMutationCount_" + con + "_" + sample);
					for(String sample : sampleNames)
						outputCon.write("\tRMR_CG_" + con + "_" + sample);
					outputCon.write("\tRMR_CG_" + con + "_average\tRMR_CG_" + con + "_sd");
					for(String sample : sampleNames)
						outputCon.write("\tRMR_5mCG_" + con + "_" + sample);
					outputCon.write("\tRMR_5mCG_" + con + "_average\tRMR_5mCG_" + con + "_sd");
				}
			}
			outputCon.write("\n");
		}
		for(String bed : names) {
			chrMap = intervals.get(bed);
			double[] caseTotalSame = new double[files.length];
			double[] caseTotalDiff = new double[files.length];
			double totalCoveredSame = 0.0;
			double correctedSame = 0.0;
			double totalCoveredDiff = 0.0;
			double correctedDiff = 0.0;
			List<HashMap<String,Double>> mutCountSame = new ArrayList<HashMap<String,Double>>();
			for(j = 0; j < files.length; j++)
				mutCountSame.add(new HashMap<String,Double>(100,0.99999f));
			HashMap<String,Double> covCountSame = new HashMap<String,Double>(50,0.99999f);
			HashMap<String,Double> corCountSame = new HashMap<String,Double>(50,0.99999f);
			List<HashMap<String,Double>> mutCountDiff = new ArrayList<HashMap<String,Double>>();
			for(j = 0; j < files.length; j++)
				mutCountDiff.add(new HashMap<String,Double>(100,0.99999f));
			HashMap<String,Double> covCountDiff = new HashMap<String,Double>(50,0.99999f);
			HashMap<String,Double> corCountDiff = new HashMap<String,Double>(50,0.99999f);
			for(String chr : chrMap.keySet()) {
				currentTree = chrMap.get(chr);
				Iterator<segment> segmentIterator = currentTree.constructIterable().iterator();
				while(segmentIterator.hasNext()) {
					currentSegment = segmentIterator.next();
					for(j = 0; j < currentSegment.mutationsSame.length; j++)
						caseTotalSame[j] += currentSegment.mutationsSame[j];
					totalCoveredSame += currentSegment.totalEnrichedSame;
					correctedSame += currentSegment.correctedEnrichedSame;
					if(stranded) {
						for(j = 0; j < currentSegment.mutationsDiff.length; j++)
							caseTotalDiff[j] += currentSegment.mutationsDiff[j];
						totalCoveredDiff += currentSegment.totalEnrichedDiff;
						correctedDiff += currentSegment.correctedEnrichedDiff;
					}
					if(outputContext != null) {
						for(String con : contextList) {
							for(j = 0; j < currentSegment.mutationsSame.length; j++) {
								contextCount = currentSegment.contextMutationCountSame.get(j).get(con);
								if(contextCount == null)
									contextCount = 0;
								if(!mutCountSame.get(j).containsKey(con))
									mutCountSame.get(j).put(con, (double)contextCount);
								else
									mutCountSame.get(j).put(con, mutCountSame.get(j).get(con) + (double)contextCount);
							}
						}
						for(String con : contextList) {
							contextCountDouble = currentSegment.contextCoveredCountSame.get(con);
							if(contextCountDouble == null)
								contextCountDouble = 0.0;
							if(!covCountSame.containsKey(con))
								covCountSame.put(con, contextCountDouble);
							else
								covCountSame.put(con, covCountSame.get(con) + contextCountDouble);
						}
						for(String con : contextList) {
							contextCountDouble = currentSegment.contextCorrectedCountSame.get(con);
							if(contextCountDouble == null)
								contextCountDouble = 0.0;
							if(!corCountSame.containsKey(con))
								corCountSame.put(con, contextCountDouble);
							else
								corCountSame.put(con, corCountSame.get(con) + contextCountDouble);
						}
						if(stranded) {
							for(String con : contextList) {
								for(j = 0; j < currentSegment.mutationsSame.length; j++) {
									contextCount = currentSegment.contextMutationCountDiff.get(j).get(con);
									if(contextCount == null)
										contextCount = 0;
									if(!mutCountDiff.get(j).containsKey(con))
										mutCountDiff.get(j).put(con, (double)contextCount);
									else
										mutCountDiff.get(j).put(con, mutCountDiff.get(j).get(con) + (double)contextCount);
								}
							}
							for(String con : contextList) {
								contextCountDouble = currentSegment.contextCoveredCountDiff.get(con);
								if(contextCountDouble == null)
									contextCountDouble = 0.0;
								if(!covCountDiff.containsKey(con))
									covCountDiff.put(con, contextCountDouble);
								else
									covCountDiff.put(con, covCountDiff.get(con) + contextCountDouble);
							}
							for(String con : contextList) {
								contextCountDouble = currentSegment.contextCorrectedCountDiff.get(con);
								if(contextCountDouble == null)
									contextCountDouble = 0.0;
								if(!corCountDiff.containsKey(con))
									corCountDiff.put(con, contextCountDouble);
								else
									corCountDiff.put(con, corCountDiff.get(con) + contextCountDouble);
							}
						}
					}
				}
			}
			int u;
			HashMap<String, Double> currentMap = null;
			if(outputContext != null) {
				outputCon.write(bed);
				for(double total : totalMutation)
					outputCon.write("\t" + (int)(total*1000));
				for(String con : contextList) {
					double[] RMR_cg_con_same = new double[caseTotalSame.length];
					double[] RMR_5mcg_con_same = new double[caseTotalSame.length];
					double[] RMR_cg_con_diff = null, RMR_5mcg_con_diff = null;
					if(stranded) {
						RMR_cg_con_diff = new double[caseTotalSame.length];
						RMR_5mcg_con_diff = new double[caseTotalSame.length];
					}
					outputCon.write("\t" + covCountSame.get(con) + "\t" + corCountSame.get(con));
					for(HashMap<String,Double> mutMap : mutCountSame)
						outputCon.write("\t" + mutMap.get(con));
					for(u = 0; u < mutCountSame.size(); u++) {
						currentMap = mutCountSame.get(u);
						RMR_cg_con_same[u] = (currentMap.get(con)/(covCountSame.get(con)/1000000.0))/totalMutation[u];
						outputCon.write("\t" + RMR_cg_con_same[u]);
					}
					stat = calcStats(RMR_cg_con_same);
					outputCon.write("\t" + stat[0] + "\t" + stat[1]);
					for(u = 0; u < mutCountSame.size(); u++) {
						currentMap = mutCountSame.get(u);
						RMR_5mcg_con_same[u] = (currentMap.get(con)/(corCountSame.get(con)/1000000.0))/totalMutation[u];
						outputCon.write("\t" + RMR_5mcg_con_same[u]);
					}
					stat = calcStats(RMR_5mcg_con_same);
					outputCon.write("\t" + stat[0] + "\t" + stat[1]);
					if(stranded) {
						outputCon.write("\t" + covCountDiff.get(con) + "\t" + corCountDiff.get(con));
						for(HashMap<String,Double> mutMap : mutCountDiff)
							outputCon.write("\t" + mutMap.get(con));
						for(u = 0; u < mutCountDiff.size(); u++) {
							currentMap = mutCountDiff.get(u);
							RMR_cg_con_diff[u] = (currentMap.get(con)/(covCountDiff.get(con)/1000000.0))/totalMutation[u];
							outputCon.write("\t" + RMR_cg_con_diff[u]);
						}
						stat = calcStats(RMR_cg_con_diff);
						outputCon.write("\t" + stat[0] + "\t" + stat[1]);
						for(u = 0; u < mutCountDiff.size(); u++) {
							currentMap = mutCountDiff.get(u);
							RMR_5mcg_con_diff[u] = (currentMap.get(con)/(corCountDiff.get(con)/1000000.0))/totalMutation[u];
							outputCon.write("\t" + RMR_5mcg_con_diff[u]);
						}
						stat = calcStats(RMR_5mcg_con_diff);
						outputCon.write("\t" + stat[0] + "\t" + stat[1]);
					}
				}
				outputCon.write("\n");
				outputCon.flush();
			}
			double[] RMR_cg_same = new double[caseTotalSame.length];
			double[] RMR_5mcg_same = new double[caseTotalSame.length];
			double[] RMR_cg_diff = null, RMR_5mcg_diff = null;
			if(stranded) {
				RMR_cg_diff = new double[caseTotalSame.length];
				RMR_5mcg_diff = new double[caseTotalSame.length];
			}
			statOutput.write(bed);
			for(u = 0; u < caseTotalSame.length; u++)
				statOutput.write("\t" + (int)(totalMutation[u]*1000));
			statOutput.write("\t" + totalCoveredSame + "\t" + correctedSame);
			for(u = 0; u < caseTotalSame.length; u++)
				statOutput.write("\t" + caseTotalSame[u]);
			for(u = 0; u < caseTotalSame.length; u++) {
				RMR_cg_same[u] = (caseTotalSame[u])/(totalCoveredSame/1000000.0)/totalMutation[u];
				RMR_5mcg_same[u] = (caseTotalSame[u]/(correctedSame/1000000.0))/totalMutation[u];
			}
			for(u = 0; u < caseTotalSame.length; u++)
				statOutput.write("\t" + RMR_cg_same[u]);
			stat = calcStats(RMR_cg_same);
			statOutput.write("\t" + stat[0] + "\t" + stat[1]);
			for(u = 0; u < caseTotalSame.length; u++)
				statOutput.write("\t" + RMR_5mcg_same[u]);
			stat = calcStats(RMR_5mcg_same);
			statOutput.write("\t" + stat[0] + "\t" + stat[1]);
			if(stranded) {
				statOutput.write("\t" + totalCoveredDiff + "\t" + correctedDiff);
				for(u = 0; u < caseTotalDiff.length; u++)
					statOutput.write("\t" + caseTotalDiff[u]);
				for(u = 0; u < caseTotalDiff.length; u++) {
					RMR_cg_diff[u] = (caseTotalDiff[u])/(totalCoveredDiff/1000000.0)/totalMutation[u];
					RMR_5mcg_diff[u] = (caseTotalDiff[u]/(correctedDiff/1000000.0))/totalMutation[u];
				}
				for(u = 0; u < caseTotalDiff.length; u++)
					statOutput.write("\t" + RMR_cg_diff[u]);
				stat = calcStats(RMR_cg_diff);
				statOutput.write("\t" + stat[0] + "\t" + stat[1]);
				for(u = 0; u < caseTotalDiff.length; u++)
					statOutput.write("\t" + RMR_5mcg_diff[u]);
				stat = calcStats(RMR_5mcg_diff);
				statOutput.write("\t" + stat[0] + "\t" + stat[1]);
			}
			statOutput.write("\n");
		}
		if(outputCon != null)
			outputCon.close();
		statOutput.close();
	}
	private double[] calcStats(double[] scores) {
		double[] stats = new double[2];
		double tot = 0.0;
		for(double score : scores)
			tot += score;
		double avrg = tot/(double)(scores.length);
		stats[0] = avrg;
		double num = 0.0;
		for(double score : scores)
			num += (score - avrg) * (score - avrg);
		stats[1] = (scores.length == 1) ? 0.0 : Math.sqrt(num / ((double)(scores.length)-1.0));
		return stats;
	}
	private HashMap<String, HashMap<String,IntervalST<segment>>> readCorrectingFiles(String files, String names,boolean stranded) throws IOException {
		HashMap<String, HashMap<String, IntervalST<segment>>> intervals = new HashMap<String, HashMap<String, IntervalST<segment>>>(25, 0.99999f);
		String[] allFiles = files.split(",");
		String[] allNames = names.split(",");
		IntervalST<segment> currentTree = null;
		int start =-1, end = -1;
		int i;
		for(i = 0; i < allFiles.length; i++) {
			HashMap<String,IntervalST<segment>> map = new HashMap<String, IntervalST<segment>>(50,0.99999f);
			String currentChromosome = "";
			BufferedReader br = new BufferedReader(new FileReader(allFiles[i]));
			String[] tokens = null;
			String line = null;
			while((line = br.readLine()) != null) {
				tokens = line.split("\t");
				tokens[0] = tokens[0].replaceFirst("chr", "");
				if(!currentChromosome.equals(tokens[0])) {
					if(!map.containsKey(tokens[0]))
						map.put(tokens[0],new IntervalST<segment>());
					currentChromosome = tokens[0];
					currentTree = map.get(currentChromosome);
				}
				start = Integer.parseInt(tokens[1]);
				end = Integer.parseInt(tokens[2]);
				segment s = new segment();
				s.start = start;
				s.end = end;
				if(stranded)
					s.posStrand = (tokens[3].equals("+"));
				if(start > end) {
					System.out.println("Error at line: " + line);
					continue;
				}
				currentTree.put(new Interval1D(start, end), s);
			}
			br.close();
			intervals.put(allNames[i], map);
		}
		return intervals;
	}
	private HashMap<String,HashMap<Integer,String>> readDatabase(String file, HashMap<String, HashMap<String,IntervalST<segment>>> intervals, boolean stranded, String correctingFiles, String correctingNames, boolean verbose) throws IOException {
		HashMap<String, HashMap<Integer,String>> map = new HashMap<String, HashMap<Integer,String>>(50,0.99999f);
		BufferedReader br = new BufferedReader(new FileReader(file));
		String[] tokens = null;
		String line = null;
		String currentChromosome = "";
		ArrayList<IntervalST<segment>> currentTrees = null;
		ArrayList<IntervalST<segment>> currentCorrectingTrees = null;
		HashMap<Integer,String> currentMap = null;
		IntervalST<segment> currentTree = null;
		IntervalST<segment> currentCorrectingTree = null;
		Iterable<Interval1D> segmentSet = null;
		Iterable<Interval1D> segmentCorrectingSet = null;
		Iterator<Interval1D> intervalIterator = null;
		Iterator<Interval1D> intervalCorrectingIterator = null;
		HashMap<String, HashMap<String, IntervalST<segment>>> correctingTrees = null;
		if(correctingFiles != null)
			correctingTrees = readCorrectingFiles(correctingFiles, correctingNames, stranded);
		segment currentSegment = null;
		int start = -1, end = -1;
		double frac = -1.0;
		HashMap<String,Double> contextCoveredCount = null;
		HashMap<String,Double> contextCorrectedCount = null;
		String segStrand = "";
		while((line = br.readLine()) != null) {
			tokens = line.split("\t");
			tokens[0] = tokens[0].replaceFirst("chr", "");
			start = Integer.parseInt(tokens[1]);
			end = start + 1;
			frac = Double.parseDouble(tokens[3])/100.0;
			if(!currentChromosome.equals(tokens[0])) {
				if(!map.containsKey(tokens[0]))
					map.put(tokens[0], new HashMap<Integer, String>(10000000,0.99999f));
				if(verbose)
					System.out.println("Currently reading chromosome: " + tokens[0]);
				currentChromosome = tokens[0];
				currentMap = map.get(currentChromosome);
				currentTrees = new ArrayList<IntervalST<segment>>(50);
				currentCorrectingTrees = new ArrayList<IntervalST<segment>>(50);
				for(String bed : intervals.keySet()) {
					currentTree = intervals.get(bed).get(currentChromosome);
					if(currentTree == null)
						continue;
					else
						currentTrees.add(currentTree);
				}
				if(correctingFiles != null) {
					for(String correctName : correctingTrees.keySet()) {
						currentCorrectingTree = correctingTrees.get(correctName).get(currentChromosome);
						if(currentCorrectingTree == null)
							continue;
						else
							currentCorrectingTrees.add(currentCorrectingTree);
					}
				}
			}
			for(IntervalST<segment> curr : currentTrees) {
				segmentSet = curr.searchAll(new Interval1D(start,end));
				if(size(segmentSet) > 0) {
					intervalIterator = segmentSet.iterator();
					while(intervalIterator.hasNext()) {
						currentSegment = curr.get(intervalIterator.next());
						segStrand = (currentSegment.posStrand) ? "+" : "-";
						if(!stranded || tokens[4].equals(segStrand)) {
							currentSegment.totalEnrichedSame += 1.0;
							currentSegment.correctedEnrichedSame += frac;
							contextCoveredCount = currentSegment.contextCoveredCountSame;
							contextCorrectedCount = currentSegment.contextCorrectedCountSame;
						} else {
							currentSegment.totalEnrichedDiff += 1.0;
							currentSegment.correctedEnrichedDiff += frac;
							contextCoveredCount = currentSegment.contextCoveredCountDiff;
							contextCorrectedCount = currentSegment.contextCorrectedCountDiff;
						}
						if(!contextCoveredCount.containsKey(tokens[2]))
							contextCoveredCount.put(tokens[2], 1.0);
						else
							contextCoveredCount.put(tokens[2], contextCoveredCount.get(tokens[2])+1.0);
						if(!contextCorrectedCount.containsKey(tokens[2]))
							contextCorrectedCount.put(tokens[2], frac);
						else
							contextCorrectedCount.put(tokens[2], contextCorrectedCount.get(tokens[2])+frac);
					}
					if(correctingFiles != null) {
						for(IntervalST<segment> correct : currentCorrectingTrees) {
							segmentCorrectingSet = correct.searchAll(new Interval1D(start,end));
							if(size(segmentCorrectingSet) > 0) {
								intervalCorrectingIterator = segmentCorrectingSet.iterator();
								while(intervalCorrectingIterator.hasNext()) {
									currentSegment = correct.get(intervalCorrectingIterator.next());
									segStrand = (currentSegment.posStrand) ? "+" : "-";
									if(!stranded || tokens[5].equals(segStrand)) {
										currentSegment.totalEnrichedSame += 1.0;
										currentSegment.correctedEnrichedSame += frac;
									} else {
										currentSegment.totalEnrichedDiff += 1.0;
										currentSegment.correctedEnrichedDiff += frac;
									}
								}
							}
						}
					}
				}
			}
			currentMap.put(start, tokens[2]);
		}
		br.close();
		if(correctingFiles != null)
			outputCorrectingStatistics(correctingTrees,stranded);
		return map;
	}
	private void outputCorrectingStatistics(HashMap<String, HashMap<String, IntervalST<segment>>> map, boolean stranded) {
		HashMap<String,IntervalST<segment>> bedMap = null;
		IntervalST<segment> currentTree = null;
		Iterator<segment> segmentIterator = null;
		segment curr = null;
		for(String bed : map.keySet()) {
			double totalEnrichedSame = 0.0;
			double correctedSame = 0.0;
			double totalEnrichedDiff = 0.0;
			double correctedDiff = 0.0;
			bedMap = map.get(bed);
			for(String chr : bedMap.keySet()) {
				currentTree = bedMap.get(chr);
				segmentIterator = currentTree.constructIterable().iterator();
				while(segmentIterator.hasNext()) {
					curr = segmentIterator.next();
					totalEnrichedSame += curr.totalEnrichedSame;
					correctedSame += curr.correctedEnrichedSame;
					if(stranded) {
						totalEnrichedDiff += curr.totalEnrichedDiff;
						correctedDiff += curr.correctedEnrichedDiff;
					}
				}
			}
			System.out.println(bed + " - " + ((!stranded) ? totalEnrichedSame + " - " + correctedSame : totalEnrichedSame + " - " + correctedSame + " - " + totalEnrichedDiff + " - " + correctedDiff));
		}
	}
	private HashMap<String,HashMap<String,IntervalST<segment>>> readIntervals(String files, String nameBedFiles, boolean stranded, int totalCases) throws IOException {
		HashMap<String,HashMap<String,IntervalST<segment>>> intervals = new HashMap<String,HashMap<String,IntervalST<segment>>>(50);
		String[] bedFiles = files.split(",");
		String[] names = nameBedFiles.split(",");
		IntervalST<segment> currentTree = null;
		int start =-1, end = -1;
		int i,j;
		for(i = 0; i < bedFiles.length; i++) {
			HashMap<String,IntervalST<segment>> map = new HashMap<String, IntervalST<segment>>(50,0.99999f);
			String currentChromosome = "";
			BufferedReader br = new BufferedReader(new FileReader(bedFiles[i]));
			String[] tokens = null;
			String line = null;
			while((line = br.readLine()) != null) {
				if(line.toLowerCase().startsWith("track") || line.toLowerCase().startsWith("browser"))
					continue;
				tokens = line.split("\t");
				tokens[0] = tokens[0].replaceFirst("chr", "");
				if(!currentChromosome.equals(tokens[0])) {
					if(!map.containsKey(tokens[0]))
						map.put(tokens[0],new IntervalST<segment>());
					currentChromosome = tokens[0];
					currentTree = map.get(currentChromosome);
				}
				start = Integer.parseInt(tokens[1]);
				end = Integer.parseInt(tokens[2]);
				segment s = new segment();
				s.start = start;
				s.end = end;
				s.contextMutationCountSame = new ArrayList<HashMap<String,Integer>>();
				for(j = 0; j < totalCases; j++)
					s.contextMutationCountSame.add(new HashMap<String,Integer>(100,0.99999f));
				if(stranded) {
					if(tokens.length != 6) {
						System.out.println("Please provide a BED6-formatted file for performing stranded analysis");
						System.exit(-1);
					}
					s.posStrand = (tokens[5].equals("+"));
					s.contextMutationCountDiff = new ArrayList<HashMap<String,Integer>>();
					for(j = 0; j < totalCases; j++)
						s.contextMutationCountDiff.add(new HashMap<String,Integer>(100,0.99999f));
					s.contextCoveredCountDiff = new HashMap<String,Double>(25, 0.9999f);
					s.contextCorrectedCountDiff = new HashMap<String,Double>(25, 0.9999f);
					s.mutationsDiff = new double[totalCases];
				}
				s.mutationsSame = new double[totalCases];
				if(start > end) {
					System.out.println("Error at line: " + line);
					continue;
				}
				currentTree.put(new Interval1D(start, end), s);
			}
			br.close();
			intervals.put(names[i], map);
		}
		return intervals;
	}
	private int size(Iterable<?> it) {
		if (it instanceof Collection)
			return ((Collection<?>)it).size();
		int i = 0;
		for (Object obj : it) i++;
		return i;
	}
}
