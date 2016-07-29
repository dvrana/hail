package org.broadinstitute.hail.stats

import breeze.linalg._
import breeze.numerics._
import org.apache.commons.math3.distribution.ChiSquaredDistribution

import org.broadinstitute.hail.annotations.Annotation
import org.broadinstitute.hail.expr._

object LogisticRegressionModel {

  def `waldType`: Type = TStruct(
    ("beta", TDouble),
    ("se", TDouble),
    ("zstat", TDouble),
    ("pval", TDouble),
    ("nIter", TInt),
    ("converged", TBoolean),
    ("exploded", TBoolean))

  def `likelihoodRatioType`: Type = TStruct(
    ("beta", TDouble),
    ("chi2", TDouble),
    ("pval", TDouble),
    ("nIter", TInt),
    ("converged", TBoolean),
    ("exploded", TBoolean))

  def `scoreType`: Type = TStruct(
    ("chi2", TDouble),
    ("pval", TDouble))
}

class LogisticRegressionModel(X: DenseMatrix[Double], y: DenseVector[Double]) {
  require(y.length == X.rows)
  require(y.forall(yi => yi >= 0 && yi <= 1))
  require(sum(y) > 0 && sum(y) < y.length)

  val n = X.rows
  val m = X.cols

  //possibly remove once handled by general case
  def loglkInterceptOnly(): Double = {
    val avg = sum(y) / n
    sum(log(y * avg + (1d - y) * (1d - avg)))
  }

  def bInterceptOnly(interceptCol: Int = 0): DenseVector[Double] = {
    val b = DenseVector.zeros[Double](m)
    val avg = sum(y) / n
    b(interceptCol) = math.log(avg / (1 - avg))
    b
  }

  // intended per variant, starting from fit without genotype
  def fit(b0: DenseVector[Double] = DenseVector.zeros[Double](m), maxIter: Int = 100, tol: Double = 1E-10): LogisticRegressionFit = {
    require(X.cols == b0.length)
    require(maxIter > 0)

    var b = b0.copy
    var mu = DenseVector.zeros[Double](n)
    var score = DenseVector.zeros[Double](m)
    var fisher = DenseMatrix.zeros[Double](m, m)
    var iter = 0
    var converged = false
    var exploded = false

    while (!converged && !exploded && iter < maxIter) {
      iter += 1

      mu = sigmoid(X * b)
      score = X.t * (y - mu)
      fisher = X.t * (X(::, *) :* (mu :* (1d - mu))) // would mu.map(x => x * (1 - x)) be faster?

      // for debugging:
      // println(s"b = $b")
      // println(s"mu = $mu")
      // println(s"score = $score")
      // println(s"fisher = $fisher")

      // alternative algorithm avoids both mult by X.t and direct inversion, clearly better for Firth
      // val sqrtW = sqrt(mu :* (1d - mu))
      // val qrFact = qr.reduced(X(::, *) :* sqrtW)
      // solve qrFact.R * bDiff = qrFact.Q.t * (y - mu) with R upper triangular

      // for Wald: return diagonal of inverse as well, which is diagonal of inv(R)^T * inv(R)
      // for Firth, modify score using:
      // val QQ = qrFact.Q :* qrFact.Q
      // val h = sum(QQ(*, ::))

      try {
        val bDiff = fisher \ score // could also return bDiff if last adjustment improves Wald accuracy. Conceptually better to have b, mu, and fisher correspond.

        if (norm(bDiff) < tol)
          converged = true
        else
          b += bDiff
      }
      catch {
        case e: breeze.linalg.MatrixSingularException => exploded = true
      }
    }

    LogisticRegressionFit(b, mu, fisher, iter, converged, exploded)
  }

  def scoreTest(b: DenseVector[Double], chiSqDist: ChiSquaredDistribution): Option[ScoreStat] = {
    require(X.cols == b.length)

    val mu = sigmoid(X * b)
    val score = X.t * (y - mu)

    try {
      val chi2 = score dot ((X.t * (X(::, *) :* (mu :* (1d - mu)))) \ score) // score.t * inv(X.t W X) * score

      //alternative using QR:
      //val sqrtW = sqrt(mu :* (1d - mu))
      //val Qty0 = qr.reduced.justQ(X(::, *) :* sqrtW).t * ((y - mu) :/ sqrtW)
      //val chi2 = Qty0 dot Qty0  // better to create normSq Ufunc

      val p = 1d - chiSqDist.cumulativeProbability(chi2)

      Some(ScoreStat(chi2, p))
    }
    catch {
      case e: breeze.linalg.MatrixSingularException => None
    }
  }

  def toAnnotation(scoreStat: Option[ScoreStat]): Annotation = {
    scoreStat
      .map(s => Annotation(s.chi2, s.p))
      .getOrElse(Annotation.empty)
  }
}

case class LogisticRegressionFit(
  b: DenseVector[Double],
  mu: DenseVector[Double],
  fisher: DenseMatrix[Double],
  nIter: Int,
  converged: Boolean,
  exploded: Boolean) {

  def loglk(y: DenseVector[Double]): Double = sum(log((y :* mu) + ((1d - y) :* (1d - mu))))

  def waldTest(): Option[WaldStat] = {
    if (converged) {
      try {
        val se = sqrt(diag(inv(fisher))) // breeze uses LU to invert, dgetri...for Wald, better to pass through from fit? if just gt, can solve fisher \ (1,0,...,0) or use schur complement

        val z = b :/ se
        val sqrt2 = math.sqrt(2)
        val p = z.map(zi => 1 + erf(-abs(zi) / sqrt2))

        Some(WaldStat(b, se, z, p))
      }
      catch {
        case e: breeze.linalg.MatrixSingularException => None
      }
    }
    else
      None
  }

  def likelihoodRatioTest(y: DenseVector[Double], loglk0: Double, chiSqDist: ChiSquaredDistribution): Option[LikelihoodRatioStat] = {
    if (converged) {
      val chi2 = 2 * (loglk(y) - loglk0)
      val p = 1d - chiSqDist.cumulativeProbability(chi2)

      Some(LikelihoodRatioStat(b, chi2, p)
      )
    }
    else
      None
  }

  def toAnnotation(waldStat: Option[WaldStat]): Annotation =
    waldStat
      .map(s => Annotation(s.b(0), s.se(0), s.z(0), s.p(0), nIter, converged, exploded))
      .getOrElse(
        Annotation(Annotation.empty, Annotation.empty, Annotation.empty, Annotation.empty, nIter, converged, exploded))

  def likelihoodRatioAnnotation(lrStat: Option[LikelihoodRatioStat]): Annotation =
    lrStat
      .map(s => Annotation(s.b(0), s.chi2: Double, s.p: Double, nIter, converged, exploded))
      .getOrElse(
        Annotation(Annotation.empty, Annotation.empty, Annotation.empty, nIter, converged, exploded))
}

case class WaldStat(b: DenseVector[Double], se: DenseVector[Double], z: DenseVector[Double], p: DenseVector[Double])

case class ScoreStat(chi2: Double, p: Double)

case class LikelihoodRatioStat(b: DenseVector[Double], chi2: Double, p: Double)
