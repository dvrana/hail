package org.broadinstitute.hail.utils

import org.apache.spark.rdd.{OrderedRDD, OrderedRDDLeftJoin}
import org.broadinstitute.hail.SparkSuite
import org.broadinstitute.hail.check.{Gen, Prop}
import org.broadinstitute.hail.check.Arbitrary._
import org.broadinstitute.hail.Utils._
import org.broadinstitute.hail.variant.Variant
import org.testng.annotations.Test

class OrderedRDDSuite extends SparkSuite {

  @Test def test() {
    val g = for (uniqueInts <- Gen.buildableOf[Set[Variant], Variant](Gen.choose(1, 1000).map(i => Variant("16", i, "A", "T")))
      .map(set => set.toIndexedSeq.sorted).filter(_.nonEmpty);
                 toZip <- Gen.buildableOfN[IndexedSeq[String], String](uniqueInts.size, arbitrary[String]);
                 nPar <- Gen.choose(1, 10).map(i => Math.min(i, uniqueInts.size))
    ) yield {
      (nPar, uniqueInts.zip(toZip))
    }

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

    p.check(size = 1000, count = 100) // important to keep size at ~1000 to get reasonable levels of match and no match
  }

}
