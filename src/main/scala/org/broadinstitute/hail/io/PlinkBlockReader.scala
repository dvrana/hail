package org.broadinstitute.hail.io

import org.apache.hadoop.conf.Configuration
import org.apache.hadoop.io.LongWritable
import org.apache.hadoop.mapred.FileSplit
import org.broadinstitute.hail.variant.{GenotypeBuilder, GenotypeStreamBuilder, Variant}

import scala.collection.mutable

object PlinkBlockReader {
  final val GT_CONVERSION = Array(2, -1, 1, 0)
}

class PlinkBlockReader(job: Configuration, split: FileSplit) extends IndexedBinaryBlockReader[Int](job, split) {
  var variantIndex = 0
  val nSamples = job.getInt("nSamples", 0)
  val compressGS = job.getBoolean("compressGS", false)
  val blockLength = ((nSamples / 4.00) + .75).toInt

  val ab = new mutable.ArrayBuilder.ofByte

  seekToFirstBlock(split.getStart)

  def seekToFirstBlock(start: Long) {
    if (start > 2)
      variantIndex = ((start - 3) / blockLength).toInt

    require(variantIndex >= 0)

    pos = variantIndex * blockLength + 3
    if (pos < start) {
      variantIndex += 1
      pos = variantIndex * blockLength + 3
    }

    bfis.seek(variantIndex * blockLength + 3)
  }

  def next(key: LongWritable, value: ParsedLine[Int]): Boolean = {
    if (pos >= end)
      false
    else {
      val nullVariant = Variant("0", 0, "A", "T")
      val b = new GenotypeStreamBuilder(nullVariant, compress = compressGS)
      val genoBuilder = new GenotypeBuilder(nullVariant)


      bfis.readBytes(blockLength)
        .iterator
        .flatMap { i => Iterator(i & 3, (i >> 2) & 3, (i >> 4) & 3, (i >> 6) & 3) }
        .take(nSamples)
        .map(PlinkBlockReader.GT_CONVERSION)
        .foreach { i =>
          genoBuilder.clear()
          if (i >= 0)
            genoBuilder.setGT(i)
          b.write(genoBuilder)
        }

      require(variantIndex >= 0)
      value.setGS(b.result())
      value.setKey(variantIndex)
      variantIndex += 1
      pos += blockLength

      true
    }
  }

  def createValue(): PlinkParsedLine = {
    new PlinkParsedLine
  }
}