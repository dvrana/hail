package org.broadinstitute.hail.rdd

import org.apache.spark.{Partition, TaskContext}
import org.apache.spark.rdd.RDD

object OrderedRDD {
  def apply[K, V](rdd: RDD[(K, V)], partitioner: OrderedPartitioner[K]): OrderedRDD[K, V] = new OrderedRDD(rdd, partitioner)

  def make[K](rdd: RDD[(K, _)]): OrderedRDD[K, _] = {
    rdd.partitioner match {
      case Some(p: OrderedPartitioner[K]) => new OrderedRDD(rdd, p)
      case _ =>

    }
  }
}

class OrderedRDD[K, V](rdd: RDD[(K, V)], partitioner: OrderedPartitioner[K]) extends RDD[(K, V)](rdd) {

  def getPartitions: Array[Partition] = rdd.partitions

  def compute(split: Partition, context: TaskContext): Iterator[(K, V)] = rdd.compute(split, context)

  //  def leftJoin[V2](other: OrderedRDD[K, V2]): OrderedRDDLeftJoin[K, V, V2] = new OrderedRDDLeftJoin(this, other)
  def leftJoin[V2](other: OrderedRDD[K, V2]): OrderedRDDLeftJoin[K, V, V2] = new OrderedRDDLeftJoin(this, other)
}

class OrderedRDDLeftJoin[K, V1, V2](rdd1: OrderedRDD[K, V1], rdd2: OrderedRDD[K, V2]) extends RDD[(K, (V1, Option[V2]))](rdd1.sparkContext, rdd1.dependencies) {

  def getPartitions: Array[Partition] = rdd1.getPartitions

  def compute(split: Partition, context: TaskContext): Iterator[(K, (V1, Option[V2]))] = {
    val res1 = rdd1.compute(split, context)

    val first = res1.head._1

    ???
  }
}