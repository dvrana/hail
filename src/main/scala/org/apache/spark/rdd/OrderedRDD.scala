package org.apache.spark.rdd

import org.apache.spark._
import org.broadinstitute.hail.Utils._

import scala.reflect.ClassTag

object OrderedRDD {
  def apply[K, V](rdd: RDD[(K, V)])(implicit ordering: Ordering[K], tct: ClassTag[K], vct: ClassTag[V]): OrderedRDD[K, V] = {
    rdd.partitioner match {
      case Some(p: OrderedPartitioner[K, _]) => new OrderedRDD(rdd, p)
      case _ =>
        val partitioner = new OrderedPartitioner[K, V](rdd.partitions.length, rdd)
        new OrderedRDD(rdd.repartitionAndSortWithinPartitions(partitioner), partitioner)
    }
  }
}

class OrderedRDD[K, V](rdd: RDD[(K, V)], p: OrderedPartitioner[K, _])
  (implicit ordering: Ordering[K]) extends RDD[(K, V)](rdd) {

  def orderedPartitioner = p

  override val partitioner: Option[Partitioner] = Some(p)

  def getPartitions: Array[Partition] = rdd.partitions

  def compute(split: Partition, context: TaskContext): Iterator[(K, V)] = rdd.compute(split, context)

  override def getPreferredLocations(split: Partition): Seq[String] = rdd.preferredLocations(split)

  //  def leftJoin[V2](other: OrderedRDD[K, V2]): OrderedRDDLeftJoin[K, V, V2] = new OrderedRDDLeftJoin(this, other)
  def leftJoin[V2](other: OrderedRDD[K, V2]): OrderedRDDLeftJoin[K, V, V2] = new OrderedRDDLeftJoin(this, other)
}

case class OneDependency[T](rdd: RDD[T]) extends Dependency[T]

case class OrderedJoinPartition(index: Int) extends Partition

object RangeDependency {
  def getDependencies[K, V](p1: OrderedPartitioner[K, _], p2: OrderedPartitioner[K, _])(partitionId: Int): Seq[Int] = {

    val lastPartition = if (partitionId == p1.rangeBounds.length)
      p2.numPartitions - 1
    else
      p2.getPartition(p1.rangeBounds(partitionId))

    if (partitionId == 0)
      0 to lastPartition
    else {
      val startPartition = p2.getPartition(p1.rangeBounds(partitionId - 1))
      startPartition to lastPartition
    }
  }
}

class RangeDependency[K, V](p1: OrderedPartitioner[K, _], p2: OrderedPartitioner[K, _], rdd: RDD[(K, V)]) extends NarrowDependency[(K, V)](rdd) {
  override def getParents(partitionId: Int): Seq[Int] = RangeDependency.getDependencies(p1, p2)(partitionId)
}

object OrderedRDDLeftJoin {
  def apply[K, V1, V2](rdd1: RDD[(K, V1)], rdd2: RDD[(K, V2)])
    (implicit ordering: Ordering[K], tct: ClassTag[K], vct1: ClassTag[V1], vct2: ClassTag[V2]): OrderedRDDLeftJoin[K, V1, V2] = {

    new OrderedRDDLeftJoin(OrderedRDD(rdd1), OrderedRDD(rdd2))
  }
}

class OrderedRDDLeftJoin[K, V1, V2](rdd1: OrderedRDD[K, V1], rdd2: OrderedRDD[K, V2])(implicit ordering: Ordering[K])
  extends RDD[(K, (V1, Option[V2]))](rdd1.sparkContext,
    Seq(new OneToOneDependency(rdd1), new RangeDependency(rdd1.orderedPartitioner, rdd2.orderedPartitioner, rdd2)): Seq[Dependency[_]]) {

  override val partitioner: Option[Partitioner] = rdd1.partitioner

  lazy val p1 = rdd1.orderedPartitioner
  lazy val p2 = rdd2.orderedPartitioner

  val localPartitions = rdd1.getPartitions.indices.map(OrderedJoinPartition).map(_.asInstanceOf[Partition]).toArray

  def getPartitions: Array[Partition] = localPartitions

  //  override def getPartitions: Array[Partition] = {
  //    val numParts = rdd1.partitions.length
  //    Array.tabulate[Partition](numParts) { i =>
  //
  //      val prefs = rdds.map(rdd => rdd.preferredLocations(rdd.partitions(i)))
  //      // Check whether there are any hosts that match all RDDs; otherwise return the union
  //      val exactMatchLocations = prefs.reduce((x, y) => x.intersect(y))
  //      val locs = if (!exactMatchLocations.isEmpty) exactMatchLocations else prefs.flatten.distinct
  //      new ZippedPartitionsPartition(i, rdds, locs)
  //    }
  //  }

  override def getPreferredLocations(split: Partition): Seq[String] = rdd1.preferredLocations(split)

  def compute(split: Partition, context: TaskContext): Iterator[(K, (V1, Option[V2]))] = {
    val left = rdd1.compute(split, context)
    val right = RangeDependency.getDependencies(p1, p2)(split.index)
      .iterator
      .flatMap(i => rdd2.compute(rdd2.partitions(i), context))

    left.sortedLeftJoin(right)
    // FIXME there are two ways to do this, talk with Cotton
    //
    //    if (left.isEmpty)
    //      Iterator.empty
    //    else {
    //      val leftHead = left.next()
    //      val rightStart = rdd2.partitioner.get.getPartition(leftHead._1)
    //      val rightIterator = new PartitionSpanningIterator[(K, V2)](rdd2, rightStart, context)
    //      (Iterator(leftHead) ++ left).sortedLeftJoin(rightIterator)
    //    }
  }
}

class PartitionSpanningIterator[T](rdd: RDD[T], start: Int, context: TaskContext) extends Iterator[T] {
  private var i = start
  private var it = rdd.compute(rdd.partitions(i), context)

  def next(): T = {
    if (hasNext)
      it.next()
    else
      throw new java.util.NoSuchElementException("next on empty iterator")
  }

  def hasNext: Boolean = {
    if (it.hasNext)
      true
    else if (i + 1 < rdd.partitions.length) {
      i += 1
      it = rdd.compute(rdd.partitions(i), context)
      hasNext
    } else false

  }
}