package org.broadinstitute.hail.io

import org.apache.hadoop.conf.Configuration
import org.apache.spark.{RangePartitioner, SparkContext}
import org.apache.spark.storage.StorageLevel
import org.broadinstitute.hail.Utils._
import org.broadinstitute.hail.variant._
import org.broadinstitute.hail.annotations._


object GenLoader {
  def apply(genFile: String, sampleFile: String, sc: SparkContext, nPartitions: Option[Int] = None): VariantSampleMatrix[Genotype] = {
    val hConf = sc.hadoopConfiguration
    val sampleIds = BgenLoader.parseSampleFile(sampleFile, hConf)
    val nSamples = sampleIds.length
    val rdd = sc.textFile(genFile, nPartitions.getOrElse(sc.defaultMinPartitions)).map{case line => readGenLine(line, nSamples)}
    val signatures = Annotations(Map("rsid" -> new SimpleSignature("String"), "varid" -> new SimpleSignature("String")))
    VariantSampleMatrix(metadata = VariantMetadata(sampleIds).addVariantAnnotationSignatures(signatures), rdd = rdd)
  }

  def convertPPsToInt(prob: Double): Int = {
    val tmp = prob * 32768
    require(tmp >= 0 && tmp < 65535.5)
    math.round(tmp).toInt
  }

  def convertPPsToInt(probArray: Array[Double]): Array[Int] = probArray.map{d => convertPPsToInt(d)}

  def readGenLine(line: String, nSamples: Int): (Variant, Annotations, Iterable[Genotype]) = {
    val arr = line.split("\\s+")
    val variant = Variant(arr(0), arr(3).toInt, arr(4), arr(5))
    val annotations = Annotations(Map[String,String]("varid" -> arr(1), "rsid" -> arr(2)))
    val dosages = arr.drop(6).map{_.toDouble}

    if (dosages.length != (3 * nSamples))
      fatal("Number of dosages does not match number of samples")

    val genotypeStream: Iterable[Genotype] = {
      for (i <- dosages.indices by 3) yield {
        val ints = convertPPsToInt(Array(dosages(i), dosages(i + 1), dosages(i + 2)))
        val pls = BgenLoader.phredScalePPs(ints(0), ints(1), ints(2))
        val gt = BgenLoader.parseGenotype(pls)
        val pls2 = if (gt == -1) null else pls //FIXME: Is this correct behavior?
        Genotype(gt = Option(gt), pl = Option(pls2))
      }
    }
    (variant, annotations, genotypeStream)
  }
}


object GenWriter {

  def apply(outputRoot: String, vds: VariantDataset, sc: SparkContext) {
    val outGenFile = outputRoot + ".gen"
    val outSampleFile = outputRoot + ".sample"
    val hConf = sc.hadoopConfiguration

    writeSampleFile(outSampleFile, hConf, vds.sampleIds)
    writeGenFile(outGenFile, vds)
  }

  def convertPlsToDosage(pl: Array[Int]): Array[Double] = {
    val tmp = pl.map{case p => math.pow(10,p / -10.0)}
    tmp.map{d => d / tmp.sum}
  }

  def appendRow(sb: StringBuilder, v: Variant, va: Annotations, gs: Iterable[Genotype]) {
    sb.append(v.contig)
    sb += ' '
    sb.append("fakeVariantID")
    sb += ' '
    sb.append("fakeRSID")
    sb += ' '
    sb.append(v.start)
    sb += ' '
    sb.append(v.ref)
    sb += ' '
    sb.append(v.alt)

    for (gt <- gs) {
      val dosages = gt.pl match {
        case Some(x) => convertPlsToDosage(x)
        case None => Array(0.0,0.0,0.0)
      }
      sb += ' '
      sb.append(dosages.mkString(" "))
    }
  }

  def writeGenFile(outFile: String, vds: VariantDataset) {
    val kvRDD = vds.rdd.map { case (v, a, gs) =>
      (v, (a, gs.toGenotypeStream(v, compress = false)))
    }
    kvRDD.persist(StorageLevel.MEMORY_AND_DISK)
    kvRDD
      .repartitionAndSortWithinPartitions(new RangePartitioner[Variant, (Annotations, Iterable[Genotype])](vds.rdd.partitions.length, kvRDD))
      .mapPartitions { it: Iterator[(Variant, (Annotations, Iterable[Genotype]))] =>
        val sb = new StringBuilder
        it.map { case (v, (va, gs)) =>
          sb.clear()
          appendRow(sb, v, va, gs)
          sb.result()
        }
      }.writeTable(outFile, None, deleteTmpFiles = true)
    kvRDD.unpersist()
  }

  def writeSampleFile(outFile: String, hConf: Configuration, sampleIds: IndexedSeq[String]){
    val header = Array("ID_1 ID_2 missing","0 0 0")
    writeTable(outFile, hConf, header ++ sampleIds.map{case s => Array(s, s, "0").mkString(" ")})
  }
}