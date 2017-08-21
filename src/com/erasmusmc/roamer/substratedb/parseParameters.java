package com.erasmusmc.roamer.substratedb;
import java.util.*;
import com.beust.jcommander.Parameter;

public class parseParameters {
    @Parameter
    private List<String> parameters = new ArrayList<>();
 
    @Parameter(names = "--inputFasta", description = "Input reference fasta file", required = true)
    public String fasta = null;
    
    @Parameter(names = "--outputDatabase", description = "Output file for the substrate database",required = true)
    public String outputDB = null;
 
    @Parameter(names = "--substrate", description = "Mutation substrate to look for (e.g., CG)", required = true)
    public String substrate = null;
    
    @Parameter(names = "--bufferSize", description = "How many bases up- and downstream of the substrate need to be included", required = true)
    public int bufferSize = 1;
    
    @Parameter(names = "--excludeContigs", description = "A file containing contigs per line to exclude (Optional)")
    public String excludeFile  = null;
    
    @Parameter(names = "--verbose", description = "Provide verbose information (Optional)")
    public boolean verbose  = false;
    
    @Parameter(names = {"--help","-help"}, help = true, description = "Get usage information")
	public boolean help;
}