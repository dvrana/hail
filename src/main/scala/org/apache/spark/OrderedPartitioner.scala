package org.apache.spark

import java.io.{IOException, ObjectInputStream, ObjectOutputStream}

import org.apache.spark.rdd.{PartitionPruningRDD, RDD}
import org.apache.spark.serializer.JavaSerializer
import org.apache.spark.util.{CollectionsUtils, Utils}

import scala.collection.mutable
import scala.collection.mutable.ArrayBuffer
import scala.reflect.ClassTag
import scala.util.hashing._

/**
  * Lifted from org.apache.spark.RangePartitioner.scala in
  */
object CalculateKeyRanges {
  def apply[K : Ordering : ClassTag](rdd: RDD[K]): Array[K] = {
    // This is the sample size we need to have roughly balanced output partitions, capped at 1M.
    val sampleSize = math.min(20.0 * rdd.partitions.length, 1e6)
    // Assume the input partitions are roughly balanced and over-sample a little bit.
    val sampleSizePerPartition = math.ceil(3.0 * sampleSize / rdd.partitions.length).toInt
    val (numItems, sketched) = RangePartitioner.sketch(rdd, sampleSizePerPartition)
    if (numItems == 0L) {
      Array.empty
    } else {
      // If a partition contains much more than the average number of items, we re-sample from it
      // to ensure that enough items are collected from that partition.
      val fraction = math.min(sampleSize / math.max(numItems, 1L), 1.0)
      val candidates = ArrayBuffer.empty[(K, Float)]
      val imbalancedPartitions = mutable.Set.empty[Int]
      sketched.foreach { case (idx, n, sample) =>
        if (fraction * n > sampleSizePerPartition) {
          imbalancedPartitions += idx
        } else {
          // The weight is 1 over the sampling probability.
          val weight = (n.toDouble / sample.length).toFloat
          for (key <- sample) {
            candidates += ((key, weight))
          }
        }
      }
      if (imbalancedPartitions.nonEmpty) {
        // Re-sample imbalanced partitions with the desired sampling probability.
        val imbalanced = new PartitionPruningRDD(rdd, imbalancedPartitions.contains)
        val seed = byteswap32(-rdd.id - 1)
        val reSampled = imbalanced.sample(withReplacement = false, fraction, seed).collect()
        val weight = (1.0 / fraction).toFloat
        candidates ++= reSampled.map(x => (x, weight))
      }
      RangePartitioner.determineBounds(candidates, rdd.partitions.length)
    }
  }
}

/**
  * A [[org.apache.spark.Partitioner]] that partitions sortable records by range into roughly
  * equal ranges. The ranges are determined by sampling the content of the RDD passed in.
  *
  * Note that the actual number of partitions created by the OrderedPartitioner might not be the same
  * as the `partitions` parameter, in the case where the number of sampled records is less than
  * the value of `partitions`.
  */
case class OrderedPartitioner[K: Ordering : ClassTag](
  rangeBounds: Array[K],
  private var ascending: Boolean = true)
  extends Partitioner {

  def write(out: ObjectOutputStream) {
    out.writeBoolean(ascending)
    out.writeObject(rangeBounds)
  }

  var ordering = implicitly[Ordering[K]]

  def numPartitions: Int = rangeBounds.length + 1

  var binarySearch: ((Array[K], K) => Int) = CollectionsUtils.makeBinarySearch[K]

  def getPartition(key: Any): Int = {
    val k = key.asInstanceOf[K]
    var partition = 0
    if (rangeBounds.length <= 128) {
      // If we have less than 128 partitions naive search
      while (partition < rangeBounds.length && ordering.gt(k, rangeBounds(partition))) {
        partition += 1
      }
    } else {
      // Determine which binary search method to use only once.
      partition = binarySearch(rangeBounds, k)
      // binarySearch either returns the match location or -[insertion point]-1
      if (partition < 0) {
        partition = -partition - 1
      }
      if (partition > rangeBounds.length) {
        partition = rangeBounds.length
      }
    }
    if (ascending) {
      partition
    } else {
      rangeBounds.length - partition
    }
  }

  override def equals(other: Any): Boolean = other match {
    case r: OrderedPartitioner[_] =>
      r.rangeBounds.sameElements(rangeBounds) && r.ascending == ascending
    case _ =>
      false
  }

  override def hashCode(): Int = {
    val prime = 31
    var result = 1
    var i = 0
    while (i < rangeBounds.length) {
      result = prime * result + rangeBounds(i).hashCode
      i += 1
    }
    result = prime * result + ascending.hashCode
    result
  }
}

object OrderedPartitioner {
  def read[K: Ordering : ClassTag, T](in: ObjectInputStream)(rdd: RDD[(T, _)])(implicit ev: (T) => K): OrderedPartitioner[K] = {
    val ascending = in.readBoolean()
    val rangeBounds = in.readObject().asInstanceOf[Array[K]]
    OrderedPartitioner(rangeBounds, ascending)
  }
}
