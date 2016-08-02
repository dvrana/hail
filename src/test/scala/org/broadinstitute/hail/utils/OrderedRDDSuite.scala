package org.broadinstitute.hail.utils

import org.apache.spark.OrderedPartitioner
import org.apache.spark.rdd.OrderedRDD
import org.apache.spark.sql.Row
import org.apache.spark.sql.types._
import org.broadinstitute.hail.SparkSuite
import org.broadinstitute.hail.Utils._
import org.broadinstitute.hail.check.Arbitrary._
import org.broadinstitute.hail.check.{Gen, Prop}
import org.broadinstitute.hail.variant.Variant
import org.testng.annotations.Test

class OrderedRDDSuite extends SparkSuite {

  val g = for (uniqueInts <- Gen.buildableOf[Set[Variant], Variant](Gen.choose(1, 1000).map(i => Variant("16", i, "A", "T")))
    .map(set => set.toIndexedSeq.sorted).filter(_.nonEmpty);
               toZip <- Gen.buildableOfN[IndexedSeq[String], String](uniqueInts.size, arbitrary[String]);
               nPar <- Gen.choose(1, 10).map(i => Math.min(i, uniqueInts.size))
  ) yield {
    (nPar, uniqueInts.zip(toZip))
  }

  @Test def test() {

    val p = Prop.forAll(g, g) { case ((nPar1, it1), (nPar2, it2)) =>
      val m2 = it2.toMap

      val rdd1 = OrderedRDD(sc.parallelize(it1, nPar1))
      val rdd2 = OrderedRDD(sc.parallelize(it2, nPar2))

      println(s"""rdd1: ${rdd1.getPartitions.size}, ${rdd1.getPartitions.toIndexedSeq}""")
      println(s"""rdd2: ${rdd2.getPartitions.size}, ${rdd2.getPartitions.toIndexedSeq}""")

      val join: IndexedSeq[(Variant, (String, Option[String]))] = rdd1.orderedLeftJoin(rdd2).collect().toIndexedSeq
      //      println(join)
      //      println(join.size, join.count { case (v, (i1, i2)) => i2.isDefined })

      val check1 = it1 == join.map { case (k, (v1, _)) => (k, v1) }
      val check2 = join.forall { case (k, (_, v2)) => v2 == m2.get(k) }

      println(check1, check2)
      check1 && check2
    }

    p.check(size = 1000) // important to keep size at ~1000 to get reasonable levels of match and no match
  }

  @Test def testWriteRead() {
    val tmpPartitioner = tmpDir.createTempFile("partitioner")
    val tmpRdd = tmpDir.createTempFile("rdd", ".parquet")

    val p = Prop.forAll(g) { case (nPar, it) =>
      val rdd = OrderedRDD(sc.parallelize(it, nPar))

      val schema = StructType(Array(
        StructField("variant", Variant.schema, nullable = false),
        StructField("str", StringType, nullable = false)))
      hadoopDelete(tmpRdd, hadoopConf, recursive = true)
      val df = sqlContext.createDataFrame(rdd.map { case (v, s) => Row.fromSeq(Seq(v.toRow, s)) }, schema)
        .write.parquet(tmpRdd)

      writeObjectFile(tmpPartitioner, hadoopConf) { out =>
        rdd.partitioner.get.asInstanceOf[OrderedPartitioner[Variant, String]].write(out)
      }

      val status = hadoopFileStatus(tmpPartitioner, hadoopConf)


      val rddReadBack = sqlContext.readPartitioned.parquet(tmpRdd)
        .rdd
        .map(r => (Variant.fromRow(r.getAs[Row](0)), r.getAs[String](1)))

      val readBackPartitioner = readObjectFile(tmpPartitioner, hadoopConf) { in =>
        OrderedPartitioner.read[Variant, String](in)(rddReadBack)
      }

      val orderedRddRB = new OrderedRDD[Variant, String](rddReadBack, readBackPartitioner)


      orderedRddRB.zipPartitions(rdd) { case (it1, it2) =>
        it1.zip(it2)
      }
        .collect()
        .forall { case (v1, v2) => v1 == v2 }
    }

    p.check(size = 1000)

  }
}
