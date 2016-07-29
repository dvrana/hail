package org.broadinstitute.hail.methods

import breeze.linalg._
import org.apache.spark.rdd.RDD
import org.broadinstitute.hail.Utils._
import org.broadinstitute.hail.annotations.Annotation
import org.broadinstitute.hail.stats.LogisticRegressionModel
import org.broadinstitute.hail.variant._

object LogisticRegression {
  def name = "LogisticRegression"

  def buildX(gs: Iterable[Genotype], C: DenseMatrix[Double]): Option[DenseMatrix[Double]] = {
    require(gs.size == C.rows)

    val (nCalled, gtSum) = gs.flatMap(_.gt).foldRight((0,0))((gt, acc) => (acc._1 + 1, acc._2 + gt))

    // allHomRef || allHomVar || allNoCall || allHet
    if (gtSum == 0 || gtSum == 2 * nCalled || nCalled == 0 || (gtSum == nCalled && gs.flatMap(_.gt).forall(_ == 1)))
      None
    else {
      val gtMean = gtSum.toDouble / nCalled
      val gtArray = gs.map(_.gt.map(_.toDouble).getOrElse(gtMean)).toArray
      Some(DenseMatrix.horzcat(new DenseMatrix(gtArray.length, 1, gtArray), C))
    }
  }

  def apply(vds: VariantDataset, y: DenseVector[Double], cov: Option[DenseMatrix[Double]]): LogisticRegression = {
    require(cov.forall(_.rows == y.size))

    // FIXME: improve message (or require Boolean and place in command)
    if (! y.forall(yi => yi == 0d || yi == 1d))
      fatal(s"For logistic regression, phenotype must be Boolean or numeric will all values equal to 0 or 1.")

    val n = y.size
    val k = if (cov.isDefined) cov.get.cols else 0
    val d = n - k - 2

    if (d < 1)
      fatal(s"$n samples and $k ${plural(k, "covariate")} with intercept implies $d degrees of freedom.")

    info(s"Running logreg on $n samples with $k sample ${plural(k, "covariate")}...")

    val C: DenseMatrix[Double] = cov match {
      case Some(dm) => DenseMatrix.horzcat(DenseMatrix.ones[Double](n, 1), dm)
      case None => DenseMatrix.ones[Double](n, 1)
    }

    val nullModel = new LogisticRegressionModel(C, y)
    val nullFit = nullModel.fit(nullModel.bInterceptOnly())
    val b0 = DenseVector.vertcat(DenseVector(0d), nullFit.b)

    val sc = vds.sparkContext
    val yBc = sc.broadcast(y)
    val CBc = sc.broadcast(C)
    val b0Bc = sc.broadcast(b0)

    new LogisticRegression(vds.rdd
      .map{ case (v, _, gs) =>
        (v, buildX(gs, CBc.value)
          .map { X =>
            val model = new LogisticRegressionModel(X, yBc.value)
            val fit = model.fit(b0Bc.value)
            fit.toAnnotation(fit.waldTest())
          }.getOrElse(Annotation.empty)
        )
      }
    )
  }
}

case class LogisticRegression(rdd: RDD[(Variant, Annotation)])