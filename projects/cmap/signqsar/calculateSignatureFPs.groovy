import java.io.File;
import java.util.BitSet;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.interfaces.IStereoElement;
import org.openscience.cdk.io.iterator.IteratingMDLReader;
import org.openscience.cdk.nonotify.NoNotificationChemObjectBuilder;
import org.openscience.cdk.signature.AtomSignature;

fileName = "../data/CMAP_ALL_names.sdf"

bits = new HashMap<String,Integer>();
height = 2
bitCount = 0;

file = new File(fileName)
println "# " + file + ":"
iterator = new IteratingMDLReader(
  new File(fileName).newReader(),
  NoNotificationChemObjectBuilder.getInstance()
)
while (iterator.hasNext()) {
  IMolecule mol = iterator.next()
  set = new HashMap<Integer,Integer>();
  for (int i=0; i<mol.getAtomCount(); i++) {
    signature = new AtomSignature(i, height, mol).toString();
    if (bits.containsKey(signature)) {
      bit = bits.get(signature);
      if (set.containsKey(bit)) {
        set.put(bit, set.get(bit)+1)
      } else {
        set.put(bit, 1)
      }
    } else {
      bitCount++;
      bits.put(signature, bitCount);
      set.put(bitCount, 1);
    }
  }
  println set
}


