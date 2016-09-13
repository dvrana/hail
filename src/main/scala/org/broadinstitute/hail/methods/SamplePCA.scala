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
  def project(vds : VariantDataset, loadings : Map[Variant, Array[Double]], singular : Array[Double]) : RDD[(String,Array[Double])] = {
    val scores = (vds.filterVariants((v : Variant,_,_) => loadings contains v)
      .aggregateBySampleWithKeys(Array.fill[Double](k)(0.0))(
        (pos : Array[Double], v : Variant, _ : String, gt : Genotype) =>
          Array.tabulate(k)((i : Int) => pos(i) + (loadings(v)(i) * gt.nNonRefAlleles.get)),
        (pos1 : Array[Double], pos2 : Array[Double]) => Array.tabulate(k)((i : Int) => pos1(i) + pos2(i))
    ))
    scores.map( { case (s : String, a : Array[Double]) => (s,Array.tabulate(a.size)((i : Int) => a(i) * eig(i))) } )
  }
}
