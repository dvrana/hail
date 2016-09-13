package org.broadinstitute.hail.methods

import reflect.ClassTag
import org.apache.spark.rdd.RDD
import org.broadinstitute.hail.variant._
import org.broadinstitute.hail.variant.VariantDataset

/* Implements adaptive PCA
 *
 * Adaptive PCA is a hierarchical clustering method for genotypes.  At each
 * step, PCA is used to move genotypes into a lower-dimensional space which 
 * captures most of the variance of the original.  Then, a stopcondition is
 * evaluated. If the genotypes meet the stopcondition, they are viewed as a
 * homogenous population and clustering ends.   If the stopcondition is not
 * met, a clustering algorithm is used to split the genotypes into k clust-
 * ters and adaptive PCA is run again on each of the sub-clusters.
 *
 * In this implementation, the clustering algorithm is Ward's algorithm and
 * clustering stops when the eigenvalues of the first two principal compon-
 * ents exceed p percent of the sum of eigenvalues.
 *
 * Note that this class includes a hacky definition of a tree, because real
 * binary trees aren't in Scala's standard library as far as I can tell.
 */


class AdaptivePCA(k : Int, returnTree : Boolean) {
  def name = "AdaptivePCA"
    
  abstract class Tree[A]
  case class Leaf[A](v : A) extends Tree[A]
  case class Node[A](left : Tree[A],v : A,right : Tree[A]) extends Tree[A]

  def leaves[A : ClassTag](T : Tree[A]) : Seq[A] = {
    T match {
      case Leaf(v) => Array[A](v)
      case Node(l,_,r) => leaves(l) ++ leaves(r)
    }
  }

  def apply (vds : VariantDataset, iterations : Int) : Tree[(IndexedSeq[String],Option[(Double,Double,Map[Variant, Array[Double]])])] = {
    val sampleIds = vds.sampleIds
    if (vds.nSamples <= k | iterations == 0) Leaf(sampleIds,None) else {
    val PCA = if (returnTree) new SamplePCA(k,true,true) else new SamplePCA(k,false,true)
    val (scores, loadingsOpt, Some(evalues)) = PCA(vds)
    val loadings = loadingsOpt match { case Some(l) => l.collect().toMap case None => null }
    val W = new Ward()
    val D_base = W.distMat(scores)
    val D = D_base map ((S : Array[Double]) => (S map Math.sqrt))
    val clusts = W(D,2).toSeq
    /* Note- the next line is not good Hail style, it's forced by 
     * the fact that SamplePCA returns a plain 2d array w/o IDs
     * attached.
     */
    val idclusts = clusts map ((S : Set[Int]) => S map ((i : Int) => sampleIds(i)))
    def p (i : Int) (name : String, A : Any) : Boolean = idclusts(i) contains name
    val thisNode = (sampleIds,if (returnTree) Some(0.0,0.0,loadings) else None)
    Node(apply(vds.filterSamples(p(0)),iterations-1),thisNode,apply(vds.filterSamples(p(1)),iterations-1))
  } }

  def project (vds : VariantDataset, structure : Tree[(IndexedSeq[String],Option[(Array[Double],Array[Double],Map[Variant,Array[Double]],Array[Double])])]) : Tree[(IndexedSeq[String])] = {
    val ids = vds.sampleIds
    structure match {
      case Leaf(_) => Leaf(ids)
      case Node(l,(_,Some((meanl,meanr,loadings,eigs))),r) => {
        val (_,karray) = loadings.head
        val k = karray.size
        val pca = new SamplePCA(k,false,false)
        val w = new Ward()
        val svals = eigs map (Math.sqrt _)
        val points = pca.project(vds,loadings,svals).collect()
        val assignments = w.meanJoin(Array(meanl,meanr),points map (x => x._2))
        val lids = (((assignments.zipWithIndex) filter ({ case (x,i) => x == 0})) map { case (_,i) => points(i)._1}).toSet
        val rids = (((assignments.zipWithIndex) filter ({ case (x,i) => x == 1})) map { case (_,i) => points(i)._1}).toSet
        def pl (name : String, A : Any) : Boolean = lids contains name
        def pr (name : String, A : Any) : Boolean = rids contains name
        Node(project(vds.filterSamples(pl),l),ids,project(vds.filterSamples(pr),r))
      }
    }
  }
}
