package utils;

import environment.Cell;
import objects.TargetSite;
import objects.TargetSitesAndGroups;
import objects.TargetSitesGroup;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class TSfileParser {

	public TargetSitesAndGroups tsg;
	private int TSid;
	private int TSGid;

	public TSfileParser(String filename, Cell n) {
		tsg = new TargetSitesAndGroups();
		BufferedReader reader = null;
		TSid = 0;
		TSGid = 0;
		try {
			reader = new BufferedReader(new FileReader(filename));
			String text;
			while ((text = reader.readLine()) != null) {
				text = text.trim();
				if (!text.isEmpty()) {
					extractTargetSitesAndGroups(n, text);
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			try {
				if (reader != null) {
					reader.close();
				}
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}


	/**
	 * extracts the target site
	 */
	public TargetSite extractTargetSite(Cell n, int targetSiteID, String str) {
		TargetSite bufferTS = new TargetSite(n, targetSiteID, str, "", Constants.NONE, Constants.NONE);
		if (bufferTS.region.start != Constants.NONE && bufferTS.region.end != Constants.NONE) {
			n.printDebugInfo("Target site " + bufferTS + " loaded");
		} else {
			n.printDebugInfo("Error on parsing target site " + bufferTS);
			bufferTS = null;
		}
		return bufferTS;
	}


	/**
	 * parse one line from the TS file
	 */
	public void extractTargetSitesAndGroups(Cell n, String text) {
		TargetSitesGroup bufferTSG;
		bufferTSG = new TargetSitesGroup(TSGid, text);
		TargetSite bufferTS;
		String[] sitesStr;
		String bufferLogicExpression;
		int newID;

		sitesStr = text.trim().split("[ \t()]+");
		bufferLogicExpression = text;
		for (String str : sitesStr) {
			str = str.trim();
			if (!str.isEmpty() && !CellUtils.tgsl.isOperator(str)) {
				bufferTS = extractTargetSite(n, TSid, str);
				if (bufferTS != null) {
					newID = tsg.addTargetSite(bufferTS);
					bufferTSG.addTargetSite(newID);
					tsg.ts.get(newID).group.add(TSGid);
					if (newID == TSid) {
						TSid++;
					}
					bufferLogicExpression = bufferLogicExpression.replaceAll(str, newID + "");
				}
			}
		}
		bufferLogicExpression = CellUtils.tgsl.replaceOperators(bufferLogicExpression);
		bufferTSG.generateRPN(bufferLogicExpression);
		tsg.tsg.add(bufferTSG);
		TSGid++;
	}

}
