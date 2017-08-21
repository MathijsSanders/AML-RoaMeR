package com.erasmusmc.roamer.extractmethylation;
import java.util.*;
import com.beust.jcommander.Parameter;

public class parseParameters {
    @Parameter
    private List<String> parameters = new ArrayList<>();
 
    @Parameter(names = "--coverageFiles", description = "Comma-separated list of methylation coverage bigWig files", required = true)
    public String covFiles = null;
    
    @Parameter(names = "--methylationFiles", description = "Comma-separated list of methylation bigWig files", required = true)
    public String methFiles = null;
 
    @Parameter(names = "--outputCoverage", description = "Output file containing coverage information per CG (Optional)")
    public String outputCoverage = null;
    
    @Parameter(names = "--outputMethylation", description = "Output file containing methylation information per CG", required = true)
    public String outputMethylation = null;
    
    @Parameter(names = "--cgDatabase", description = "Input file containing the CG database", required = true)
    public String database = null;
    
    @Parameter(names = "--excludeContigs", description = "A file containing contigs per line to exclude (Optional)")
    public String excludeContigs = null;
    
    @Parameter(names = "--precedingBases", description = "How many preceding bases should be taking into account")
    public int preceding = 1;
    
    @Parameter(names = "--succeedingBases", description = "How many succeeding bases should be taking into account")
    public int succeeding = 0;
    
    @Parameter(names = "--verbose", description = "Provide verbose information")
    public boolean verbose  = false;
    
    @Parameter(names = {"--help","-help"}, help = true, description = "Get usage information")
	public boolean help;
}