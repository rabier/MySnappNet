package beast.app.beauti;

import java.io.File;
import java.util.List;

import beast.core.BEASTInterface;

public class SpeciesNetworkAlignmentProvider extends BeautiAlignmentProvider {
 
	@Override
	public List<BEASTInterface> getAlignments(BeautiDoc doc, File[] files) {
		
		System.out.println("CE je passe dans SpeciesNetworkALignment\n");
		
		doc.autoSetClockRate = false;
		doc.beauti.autoSetClockRate.setSelected(false);
		return super.getAlignments(doc, files);
	}
	
}
