package com.erasmusmc.roamer.substratedb;
import com.beust.jcommander.JCommander;

public class substrateDatabase {
	public void run(String[] args) {
		parseParameters pp = new parseParameters();
		JCommander jCommander = new JCommander(pp);
		JCommander.newBuilder().addObject(pp).build().parse(args);
		jCommander.setProgramName("substrate");
		if(pp.help)
			jCommander.usage();
		else
			new sdb(pp);
	}
}
