package com.erasmusmc.roamer.rmr;
import java.util.*;

public class segment 
{
	int start = -1;
	int end = -1;
	boolean posStrand = true;
	double totalEnrichedSame = 0;
	double correctedEnrichedSame = 0;
	double totalEnrichedDiff = 0;
	double correctedEnrichedDiff = 0;
	double[] mutationsSame = null;
	double[] mutationsDiff = null;
	List<HashMap<String,Integer>>contextMutationCountSame = null;
	HashMap<String,Double> contextCoveredCountSame = new HashMap<String,Double>(25, 0.9999f);
	HashMap<String,Double> contextCorrectedCountSame = new HashMap<String,Double>(25, 0.9999f);
	List<HashMap<String,Integer>> contextMutationCountDiff = null;
	HashMap<String,Double> contextCoveredCountDiff = null;
	HashMap<String,Double> contextCorrectedCountDiff = null;
}