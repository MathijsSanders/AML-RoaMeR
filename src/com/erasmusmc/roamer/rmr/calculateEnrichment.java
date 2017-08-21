package com.erasmusmc.roamer.rmr;
import com.beust.jcommander.JCommander;

public class calculateEnrichment {
	public void run(String[] args) {
		parseParameters pp = new parseParameters();
		JCommander.newBuilder().addObject(pp).build().parse(args);
		new ce(pp);
	}
}
