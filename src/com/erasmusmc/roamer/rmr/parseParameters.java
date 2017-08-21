package com.erasmusmc.roamer.rmr;
import java.util.*;
import com.beust.jcommander.Parameter;

public class parseParameters {
    @Parameter
    private List<String> parameters = new ArrayList<>();
 
    @Parameter(names = "--database", description = "Substrate database (e.g., 5mC database)")
    public String database = null;
 
    @Parameter(names = "--vcfs", description = "Comma-separated list of VCF files")
    public String vcfs = null;
    
    @Parameter(names = "--sampleNames", description = "Comma-separated list of sample names associated with the VCF files")
    public String sampleNames = null;
 
    @Parameter(names = "--bedfiles", description = "Comma-separated list of BED files")
    public String bedFiles = null;
    
    @Parameter(names = "--namesBedFiles", description = "Comma-separated list of names associated with the BED files")
    public String bedNameFiles = null;
    
    @Parameter(names = "--outputStatistics", description = "Output RMR scores and statistics")
    public String outputStats = null;
    
    @Parameter(names = "--outputContext", description = "Output RMR context scores and statistics (Optional)")
    public String outputContext = null;
    
    @Parameter(names = "--stranded", description="Estimate bias using stranded information retrieved from the BED files (Optional)")
    public boolean stranded = false;
    
    @Parameter(names = "--correctingFactors", description="Comma-separated list of BED files used for estimating correction factors (Optional)")
    public String correctingFactors = null;
    
    @Parameter(names = "--correctingName", description="Comma-separated names of the correction factors (Optional)")
    public String correctingFactorNames = null;

    @Parameter(names = "--verbose", description = "Provide verbose information (Optional)")
    public boolean verbose  = false;
}