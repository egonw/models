import java.io.File;
import java.util.BitSet;
import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.io.iterator.IteratingMDLReader;
import org.openscience.cdk.nonotify.NoNotificationChemObjectBuilder;
import org.openscience.cdk.signature.AtomSignature;

fileName = "../data/CMAP_ALL_names.sdf"

bits = new HashMap<String,Integer>();
height = 1
bitCount = 0;

file = new File(fileName)
println "# " + file + ":"
iterator = new IteratingMDLReader(
  new File(fileName).newReader(),
  NoNotificationChemObjectBuilder.getInstance()
)
while (iterator.hasNext()) {
  IMolecule mol = iterator.next()
  println mol.getProperty(CDKConstants.TITLE);
}


