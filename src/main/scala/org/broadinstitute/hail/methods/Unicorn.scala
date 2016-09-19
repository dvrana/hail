/* Implements the core UNICORN model */

package org.broadinstitute.hail.methods

import org.broadinstitute.hail.annotations._
import org.broadinstitute.hail.expr._
import org.broadinstitute.hail.variant._

class Unicorn() {
  def name = "UNICORN"

  type Stage1Dist = Map[Variant,(Double,Double)]
  type cases = Seq[(Int,(Double,Double))]

  def alleleCountAnnotate(vds : VariantDataset, refName : String = "refCount", altName : String = "altCount") : VariantDataset = {
    val (t1,i1) = vds.vaSignature.insert(TInt,refName)
    val (t2,i2) = t1.insert(TInt,altName)
    vds.mapAnnotations((v : Variant, va : Annotation, geno) => {
      var nRef : Int = 0
      var nAlt : Int = 0
      geno foreach { (g : Genotype) =>
        g.nNonRefAlleles match {
          case None => Unit
          case Some(n) => {
            nRef += 2 - n
            nAlt += n
          }
        }
      }
      i2(i1(va,Some(nRef)),Some(nAlt))
    }).copy(vaSignature=t2)
  }

  val getP = (va : Annotation, qref : Querier, qalt : Querier) => {
    val ref : Double = qref(va) match { case Some(x : Int) => x.toDouble 
                               case _ => 0.0 }
    val alt : Double = qalt(va) match { case Some(x : Int) => x.toDouble 
                               case _ => 0.0 }
    alt / Math.max(ref + alt,1.0)
  }


  // Fst calc isn't parallelized- I found it a pretty nasty problem.  I'll swing back and touch it up when I don't have as much time pressure
  def fstAnnotate(vds : VariantDataset, subpops : Seq[VariantDataset]) : VariantDataset = {
    val (t,i) = vds.vaSignature.insert(TDouble,"Fst")
    val (_,qglobref) = vds.queryVA("va.globalRefCount")
    val (_,qglobalt) = vds.queryVA("va.globalAltCount")
    val clusterGetters = subpops map (x => {
      val (_,qref) = x.queryVA("va.refCount")
      val (_,qalt) = x.queryVA("va.altCount")
      (qref,qalt)
    })
    val getrefalt = (va : Annotation, getters : (Querier,Querier)) => {
      val (qref,qalt) = getters
      val ref : Double = qref(va) match { case Some(x : Int) => x.toDouble
                                 case _ => 0.0 }
      val alt : Double = qalt(va) match { case Some(x : Int) => x.toDouble
                                 case _ => 0.0 }
      (ref, alt)
    }
    val clustcount = subpops.size
    println("alpha")
    val subpoprefalt = (0 until clustcount) map {(i : Int) => ( (subpops(i)).variantsAndAnnotations.map( ((y : (Variant,Annotation)) => (y._1,getrefalt(y._2,clusterGetters(i)) )))).collectAsMap() }
    println("beta")
    val globalrefalt : Map[Variant,(Double,Double)] = vds.variantsAndAnnotations.map {x : (Variant,Annotation) => (x._1,getrefalt(x._2,(qglobref,qglobalt)))}.collectAsMap().toMap
    println("gamma")
    val fsts = (globalrefalt.keys map ( (v : Variant) => {
      val (gref,galt) = globalrefalt(v)
      val n = (galt + gref)
      val p = galt / n
      val totalVar = p * (1.0 - p)
      val betweenVar = ((0 until clustcount) map ((i : Int) => {
        val (ref,alt) = subpoprefalt(i)(v)
        val ni = (ref + alt)
        val pi = alt / ni
        (ni / n) * (pi - p) * (pi - p)
      })).sum

      (v,if (totalVar <= 0.0) 1.0 else betweenVar / totalVar)
      })).toMap

    println("delta")

    vds.mapAnnotations((v : Variant, va : Annotation, _) => {
      i(va,Some(fsts(v)))
    }).copy(vaSignature=t)
  }
  
  // Calculates priors for Stage 1, then updates to get posteriors (a,b)
  def calcHyperparams(vds : VariantDataset, fst : Map[Variant,Double]) : Stage1Dist = {
    val (_,qref) = vds.queryVA("va.refCount")
    val (_,qalt) = vds.queryVA("va.altCount")
    val (_,qglobref) = vds.queryVA("va.globalRefCount")
    val (_,qglobalt) = vds.queryVA("va.globalAltCount")
    vds.mapWithAll ((v, va, _, _, _) => {
      val p = getP(va,qglobref,qglobalt)
      val (a,b) = if (fst(v) == 0.0) (0.5,0.5) else (p * (1 - fst(v)) / fst(v), (1 - p) * (1 - fst(v)) / fst(v))
      val ref = qref(va) match { case Some(n : Int) => n
                                 case _ => 0 }
      val alt = qalt(va) match { case Some(n : Int) => n
                                 case _ => 0 }
      (v,(a + alt, b + ref))
    }).collect().toMap
  }

  def clusterWidePriors(data : VariantDataset, clusts : Seq[Set[String]]) : Seq[Stage1Dist] = {
    var vds = alleleCountAnnotate(data,refName = "globalRefCount",altName = "globalAltCount")
    var subvds : Array[VariantDataset] = Array.tabulate(clusts.size)((i : Int) => vds.filterSamples((name : String, A : Annotation) => clusts(i) contains name) )
    println("0")
    subvds = subvds map (g => alleleCountAnnotate(g))
    println("1")
    vds = fstAnnotate(vds,subvds)
    
    println("2")

    val (_,qfst) = vds.queryVA("va.Fst")
    val fst = ( vds.mapWithAll ((v, va, _, _, _) => qfst(va) match {
      case Some(x : Double) => (v,x)
      case _ => (v,0.0)
      } ) ).collectAsMap().toMap
    
    println("3")

    var posteriors : Seq[Stage1Dist] = subvds map (x => calcHyperparams(x,fst))
    posteriors
  }

  // Uses MCMC algorithm to compute posteriors on within-cluster GP
  def fitGP(data : VariantDataset, priors : Seq[Stage1Dist]) : Unit = {
    // Read in data & priors
    // For each variant:
    //   Generate Gaussian prior for allele frequency (Tracy 39-44)
    //   Generate priors for phi
    //   Make initial pull for phi and S
    //   For each MCMC iteration (Tracy uses 6000):
    //     d

    Unit
  }

  def betaMeanVariance(beta : (Double, Double)) : (Double,Double) = {
    val (a,b) = beta
    val mean = a / (a + b)
    val variance = (a * b) / ((a + b) * (a + b) * (a + b + 1))
    (mean,variance)
  }

  /*def nulldist(samples : cases, model : Seq[Stage1Dist]) : Map[Variant,(Double,Double)] = {
    val n = samples.size.toDouble
    val variants = model(0) mapValues (x => (0.0,0.0))
    val modelParam : Seq[Map[Variant,(Double,Double)]] = model map ((m : Stage1Dist) => m mapValues betaMeanVariance)
    val ybar = (samples foldLeft variants) ((acc,loc) => { 
      val (clust,(pc1,pc2)) = loc
      (acc.keys map ((k : Variant) => {
         val (mean,variance) : (Double,Double) = modelParam(clust)(k)
         val (aggregatemean,aggregatevariance) = acc(k)
         val meanincrement = mean * 2 / n
         val varianceincrement = ( (2 * mean * (1 - mean)) + (2 * variance)) / (n * n)
         (k,(aggregatemean + meanincrement,aggregatevariance + varianceincrement)) } )).toMap 
    } )
    ybar 
  }
  
  def chisqtop(x : Double) : Double = Gamma.regularizedGammaQ(0.5, x / 2)

  // Assumes samples and vds describe the same dataset, and ref / alt alleles
  // are the same in the cases as in the UNICORN controls 
  // Assumes null dist is accurate and factors in missingness
  def test(samples : cases, vds_samples : VariantDataset, nulldist : Map[Variant,(Double,Double)]) : Map[Variant,Option[Double]] = {
    val vds = alleleCountAnnotate(vds_samples)
    //val (_, refQuery) = vds.queryVA("va.refCount")
    val (_, altQuery) = vds.queryVA("va.altCount")
    val vam = vds.variantsAndAnnotations.collect().toMap
    (vam.keys map ((v : Variant) => {
      val va : Annotation = vam(v)
      val altCount = altQuery(va) match { case Some(n : Int) => n.toDouble
                                          case _ => 0.0 }
      val Y = altCount / samples.size
      if (! (nulldist contains v)) None else Some ( {
        val (Ybar,Ybarvar) = nulldist(v)
        val chisq = (Y - Ybar) * (Y - Ybar) / Ybarvar
        chisqtop(chisq)
        } )
    } )).toMap
  }*/

  def foal(data : VariantDataset, clusts : Seq[Set[String]]) : Seq[Stage1Dist] = {
    val stage1prior : Seq[Stage1Dist] = clusterWidePriors(data,clusts)
    println("4")
    stage1prior map ((x : Stage1Dist) => x map ((y : (Variant,(Double,Double))) => (y._1,betaMeanVariance(y._2))))
  }
}
