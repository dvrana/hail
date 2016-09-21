package org.broadinstitute.hail.methods

import org.apache.spark.mllib.linalg.{Matrix, DenseMatrix}
import org.apache.spark.rdd.RDD
import org.broadinstitute.hail.annotations._
import org.broadinstitute.hail.variant._

class SamplePCA(k: Int, computeLoadings: Boolean, computeEigenvalues: Boolean) {
  def name = "SamplePCA"

  def apply(vds: VariantDataset): (Matrix, Option[RDD[(Variant, Array[Double])]], Option[Array[Double]])  = {

    val (variants, mat) = ToStandardizedIndexedRowMatrix(vds)
    val sc = vds.sparkContext
    val variantsBc = sc.broadcast(variants)

    val svd = mat.computeSVD(k, computeU = computeLoadings)

    val scores =
      svd.V.multiply(DenseMatrix.diag(svd.s))

    val loadings =
      if (computeLoadings)
        Some(svd.U.rows.map(ir =>
          (variantsBc.value(ir.index.toInt), ir.vector.toArray)))
      else
        None

    val eigenvalues =
      if (computeEigenvalues)
        Some(svd.s.toArray.map(x => x * x))
      else
        None

    (scores, loadings, eigenvalues)
  }

  // Projects samples into PC loadings
  def project(vds : VariantDataset, loadings : Map[Variant, Array[Double]], singular : Array[Double]) : RDD[(String,Seq[Double])] = {
    val localK = k
    val genos = (g : Genotype) => g.nNonRefAlleles match {
          case None => 0.0
          case Some(n : Int) => n.toDouble
        }

    val scores = (vds.filterVariants((v : Variant,_,_) => loadings contains v)
      .aggregateBySampleWithKeys(Vector.fill[Double](localK)(0.0))(
        (pos : Seq[Double], v : Variant, _ : String, gt : Genotype) =>
          Vector.tabulate(localK)((i : Int) => pos(i) + (loadings(v)(i) * genos(gt))),
        (pos1 : Seq[Double], pos2 : Seq[Double]) => Vector.tabulate(localK)((i : Int) => pos1(i) + pos2(i))
    ))
    scores.map( { case (s : String, a : Vector[Double]) => (s,Vector.tabulate(a.size)((i : Int) => a(i) * singular(i))) } )
  }
}
