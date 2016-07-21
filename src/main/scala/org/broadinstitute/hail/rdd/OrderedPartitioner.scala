package org.broadinstitute.hail.rdd

import java.io.{ByteArrayInputStream, InputStream, ObjectInputStream}

import org.apache.spark._
import org.apache.spark.serializer.{DeserializationStream, JavaSerializer, SerializerInstance}

import scala.reflect.ClassTag
import org.broadinstitute.hail.Utils._

object OrderedPartitioner {

  def deserializeViaNestedStream[T](is: InputStream, ser: SerializerInstance)(
    f: DeserializationStream => T): T = {
    val isWrapper = ser.deserializeStream(new InputStream {
      override def read(): Int = is.read()

      override def read(b: Array[Byte], off: Int, len: Int): Int = is.read(b, off, len)
    })
    try {
      f(isWrapper)
    } finally {
      isWrapper.close()
    }
  }

  def apply[K, V](partitioner: RangePartitioner[K, V]): OrderedPartitioner[K] = {

    val bytebuffer = SparkEnv.get.serializer.newInstance().serialize(partitioner)

    val in = new ObjectInputStream(new ByteArrayInputStream(bytebuffer.array()))

    val sfactory = SparkEnv.get.serializer
    sfactory match {
      case js: JavaSerializer => fatal(s"invalid serializer: ${js.getClass.getName}")
      case _ =>
        val ascending = in.readBoolean()
        val ordering = in.readObject().asInstanceOf[Ordering[K]]
        val binarySearch = in.readObject().asInstanceOf[(Array[K], K) => Int]

        val ser = sfactory.newInstance()
        deserializeViaNestedStream(in, ser) { ds =>
          implicit val classTag = ds.readObject[ClassTag[Array[K]]]()
          val rangeBounds = ds.readObject[Array[K]]()
          new OrderedPartitioner[K](ascending, ordering, binarySearch, rangeBounds)
        }
    }
  }
}

case class OrderedPartitioner[K](ascending: Boolean, ordering: Ordering[K],
  binarySearch: (Array[K], K) => Int, rangeBounds: Array[K]) extends Partitioner {

  def numPartitions: Int = rangeBounds.length + 1

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

}
