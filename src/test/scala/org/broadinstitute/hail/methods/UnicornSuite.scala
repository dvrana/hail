package org.broadinstitute.hail.methods

import org.broadinstitute.hail.SparkSuite
import org.broadinstitute.hail.Utils._
import org.broadinstitute.hail.variant._
import org.broadinstitute.hail.Utils._
import org.testng.annotations.Test

class UnicornSuite extends SparkSuite {

  @Test def test() {
    val U = new Unicorn()
 
    // Test alleleCountAnnote
    val vds = U.alleleCountAnnotate(LoadVCF(sc, "src/test/resources/tiny_m.vcf"))
    
    val var1 = Variant("20", 10019093, "A", "G")
    val var2 = Variant("20", 10026348, "A", "G")
    val var3 = Variant("20", 10026357, "T", "C")

    val (t1, refQuery) = vds.queryVA("va.refCount")
    val (t2, altQuery) = vds.queryVA("va.altCount")
    val variantAnnotationMap = vds.variantsAndAnnotations.collect().toMap
    
    assert(variantAnnotationMap contains var1)
    assert(variantAnnotationMap contains var2)
    assert(variantAnnotationMap contains var3)

    assert(refQuery(variantAnnotationMap(var1)) == Some(5))
    assert(refQuery(variantAnnotationMap(var2)) == Some(7))
    assert(refQuery(variantAnnotationMap(var3)) == Some(6))
    assert(altQuery(variantAnnotationMap(var1)) == Some(1))
    assert(altQuery(variantAnnotationMap(var2)) == Some(1))
    assert(altQuery(variantAnnotationMap(var3)) == Some(2))
    
    // Test nulldist
    val model1 = Array(Map((var1,(0.5,0.1)),(var2,(0.3,0.06))),Map((var1,(0.4,0.32)),(var2,(0.44,0.07))))
    val cases1 = Array((0,(0.0,0.0)),(0,(0.0,0.0)),(1,(0.0,0.0)))
    val cases2 = Array((0,(0.0,0.0)),(1,(0.0,0.0)),(0,(0.0,0.0)),(1,(0.0,0.0)))
    val nulldist1 = U.nulldist(cases1,model1)
    val nulldist2 = U.nulldist(cases2,model1)

    assert(D_==(nulldist1(var1)._1,2.8))
    assert(D_==(nulldist1(var1)._2,1.04))
    assert(D_==(nulldist1(var2)._1,2.08))
    assert(D_==(nulldist1(var2)._2,0.38))

    assert(D_==(nulldist2(var1)._1,3.6))
    assert(D_==(nulldist2(var1)._2,1.68))
    assert(D_==(nulldist2(var2)._1,2.96))
    assert(D_==(nulldist2(var2)._2,0.52))

    // Test test
    val chisq1 = U.test(cases2,LoadVCF(sc, "src/test/resources/tiny_m.vcf"),model1)
    val chisq1v1 = (3.6 - 1.0 * (8.0 / 6.0)) * (3.6 - 1.0 * (8.0 / 6.0)) / 1.68

    assert(D_==(chisq1(var1).get,chisq1v1))

    // Test fstcalc
    val clust1 = Set("C1046::HG02024","C1046::HG02025")
    val clust2 = Set("C1046::HG02026","C1047::HG00731")

    val fstvds = U.alleleCountAnnotate(LoadVCF(sc, "src/test/resources/tiny_m.vcf"),refName="globalRefCount",altName="globalAltCount")
    var fstsvds = Array(vds.filterSamples ((name : String, _ : Any) => clust1 contains name), vds.filterSamples ((name : String, _ : Any) => clust2 contains name))
    fstsvds = fstsvds map (g => U.alleleCountAnnotate(g))

    def fsts = U.fstcalc(fstvds,fstsvds)
    val var1p = 1.0 / 6.0
    val var2p = 1.0 / 8.0
    val var1fst = ( ((1.0 / 3.0) * (.5 - var1p) * (.5 - var1p)) + ((2.0 / 3.0) * (0.0 - var1p) * (0.0 - var1p)) ) / (var1p * (1.0 - var1p))
    val var2fst = ( (0.5 * (.25 - var2p) * (.25 - var2p)) + (0.5 * (0.0 - var2p) * (0.0 - var2p)) ) / (var2p * (1.0 - var2p))
    
    assert(D_==(fsts(var1),var1fst))
    assert(D_==(fsts(var2),var2fst))

    // Test clusterWidePriors

    val priors = U.clusterWidePriors(LoadVCF(sc, "src/test/resources/tiny_m.vcf"),Array(clust1,clust2)) map (x => x.collect().toMap)
    val (var1a0,var1b0) = priors(0)(var1)
    val (var1a1,var1b1) = priors(1)(var1)

    assert(D_==(var1a0,(var1p * (1.0 - fsts(var1)) / fsts(var1)) + 1.0))
    assert(D_==(var1b0,((1.0 - var1p) * (1.0 - fsts(var1)) / fsts(var1)) + 1.0))
    assert(D_==(var1a1,(var1p * (1.0 - fsts(var1)) / fsts(var1)) + 0.0))
    assert(D_==(var1b1,((1.0 - var1p) * (1.0 - fsts(var1)) / fsts(var1)) + 4.0))
  }

}
