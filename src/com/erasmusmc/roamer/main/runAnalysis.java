package com.erasmusmc.roamer.main;
import java.util.*;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.JCommander;
import com.erasmusmc.roamer.substratedb.substrateDatabase;
import com.erasmusmc.roamer.extractmethylation.extractMethylation;
import com.erasmusmc.roamer.rmr.calculateEnrichment;

public class runAnalysis {
	private static String versionNumber = "0.1.0";
	@Parameter
	private List<String> parameters = new ArrayList<>();
	 
	@Parameter(names = "--analysis", description = "Type of analysis (options: substrate, extractmethylation, rmr)", required = true)
	public String analysis = null;

	@Parameter(names = {"--help","-help"}, help = true, description = "Get usage information")
	private boolean help;
	
	@Parameter(names = {"--version","-version"}, description = "Get current version")
	private boolean version;
	
	public static void main(String[] args) {
		runAnalysis ra = new runAnalysis();
		JCommander jCommander = new JCommander(ra);
		jCommander.setProgramName("runAnalysis");
		JCommander.newBuilder().addObject(ra).build().parse(args);
		if(ra.version)
			System.out.println("RoaMeR " + versionNumber);
		else if(ra.help && ra.analysis == null)
			jCommander.usage();
		else {
			switch(ra.analysis.toLowerCase()) {
			case "substrate":
				new substrateDatabase().run(args);
				break;
			case "extractmethylation":
				new extractMethylation().run(args);
				break;
			case "rmr":
				new calculateEnrichment().run(args);
				break;
			default:
				 throw new IllegalArgumentException("Please specify the type of analysis: substrate, extractmethylation, rmr");
			}
		}
	}
}
