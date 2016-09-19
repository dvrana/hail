/* Provides commands to run Unicorn and Foal */

package org.broadinstitute.hail.driver

import scala.io.Source
import org.apache.spark.RangePartitioner
import org.apache.spark.storage.StorageLevel
import org.broadinstitute.hail.methods.Unicorn
import org.broadinstitute.hail.variant.Variant
import org.apache.spark.mllib.linalg.{Vector => SVector}
import org.kohsuke.args4j.{Option => Args4jOption}
import org.broadinstitute.hail.Utils._

object Foal extends Command {
  def name = "foal"

  def description = "Generates sufficient statistics to perform the Foal association test (Unicorn without the GP)"

  class Options extends BaseOptions {
    @Args4jOption(required = true, name = "-o", aliases = Array("--output"), usage = "Output file stem")
    var output : String = _

    @Args4jOption(required = true, name = "-c", aliases = Array("--clusterfile"), usage = ".cluster file from AdaptivePCA")
    var clusterfile : String = _

  }

  def newOptions = new Options

  def supportsMultiallelic = true

  def requiresVDS = true
  
  override def hidden = true // This command is for UNICORN, and is pretty useless to anyone outside UNICORN
  def run(state: State, options: Options): State = {

    val vds = state.vds

    val clustermap : Map[String,Int] = readFile(options.clusterfile, state.hadoopConf) { reader =>
      Source.fromInputStream(reader)
        .getLines()
        .map (line => line.split("\t"))
        .filter(x => ((x.size == 2) & (x(1) != "Cluster")))
        .map(line => (line(0),line(1).toInt))
        .toMap
      }

    val k : Int = (clustermap.values.reduceLeft (Math.max(_,_))) + 1

    val clusts : Array[Set[String]] = (Array.tabulate (k) ((i : Int) => clustermap.keys filter ((id : String) => clustermap(id) == i))) map (x => x.toSet)

    val U = new Unicorn()

    val clusterPosteriorsRDD = U.foal(vds,clusts)

    val clusterPosteriors = clusterPosteriorsRDD map (x => x.collect().toMap)

    val allVariants = ((clusterPosteriors map (x => x.keys.toSet)).reduceLeft ( _ ++ _ )).toArray

    writeTextFile(options.output + ".txt", state.hadoopConf) { s =>
      val header = ((0 until k) map ((i : Int) => "\tmean_clust" + i.toString + "\tvar_clust" + i.toString)).mkString("")
      s.write("Variant" + header + "\n")
      for (v <- allVariants) {
        var row = v.toString
        for (i <- 0 until k) {
          val (mean : Double,variance : Double) = clusterPosteriors(i) getOrElse (v, (-1.0,-1.0))
          row += "\t" + mean.toString + "\t" + variance.toString
        }
        s.write(row + "\n")
      }
    }
    state
  }
}

