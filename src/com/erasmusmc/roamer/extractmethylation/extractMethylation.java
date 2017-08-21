package com.erasmusmc.roamer.extractmethylation;
import com.beust.jcommander.JCommander;

public class extractMethylation  {
	public void run(String[] args) {
		parseParameters pp = new parseParameters();
		JCommander jCommander = new JCommander(pp);
		JCommander.newBuilder().addObject(pp).build().parse(args);
		jCommander.setProgramName("extractmethylation");
		if(pp.help)
			jCommander.usage();
		else
			new em(pp);
	}
}
