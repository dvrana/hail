package org.broadinstitute.hail.methods


import org.apache.commons.math3.optim.univariate.BrentOptimizer
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction
import org.apache.commons.math3.optim._
import org.apache.commons.math3.analysis.UnivariateFunction
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType
import org.apache.commons.math3.distribution.ChiSquaredDistribution
import org.apache.commons.math3.special._
import org.apache.commons.math3.optim.univariate.SearchInterval


/**
  * Created by mitja on 7/28/16.
  */
object RAFTStats {

  /**
    *
    * Computes the RAFT statistic of Lim et al. AJHG 2015.
    * @param nHomCases
    * @param nCases
    * @param nHomCtrls
    * @param nCtrls
    * @param pHom probability of homozygote in a population (typically MAF squared)
    * @param disease_prevalence
    * @return tuple of (raftstatistic, genotypic relative risk, p-value)
    */
   def raftStats( nHomCases: Int, nCases:Int, nHomCtrls:Int, nCtrls: Int, pHom:Double, diseasePrevalence:Double  ): ( Option[Double], Option[Double], Option[Double] ) = {

      val rf = new RAFTFunction( nHomCases, nCases, nHomCtrls, nCtrls, pHom, diseasePrevalence)
      val of = new UnivariateObjectiveFunction(rf)

      val optim = new BrentOptimizer(1e-6, 1e-8)

      var enrichment = 1.0
      val expHomCases = pHom *( nCases )
      val expHomControls = pHom *( nCtrls)

      if (nHomCases>0 && nHomCtrls>0)
        enrichment = (nHomCases/expHomCases). / ( nHomCtrls /expHomControls)
      else if (nHomCases>0)
        enrichment = nHomCases / expHomCases

      if (enrichment < 1)
        return (None, None, None)

      println("Enrichment " + enrichment)

      val result = optim.optimize(of, GoalType.MINIMIZE, new MaxIter(1000), new InitialGuess(Array(enrichment)), new MaxEval(1000), new SearchInterval(0,100000)  )

      val minLL = result.getValue()
      val nullLL = rf.loglik(1,true)

      val raftStat = Some(-minLL + nullLL)
      val genoRR = Some(result.getPoint)

      val chiSqr = new ChiSquaredDistribution(1)

      val pVal = Some(Gamma.regularizedGammaQ( 0.5, raftStat.get.toDouble))

      return (raftStat, genoRR, pVal)
   }

}

class RAFTFunction( nCaseHom: Int, nCases:Int, nCtrlHom:Int, nCtrl:Int, pHomozygote:Double, diseasePrevalence:Double) extends UnivariateFunction {


  def value( genoRR:Double): Double = {
    return loglik(genoRR,false)
  }

  def loglik(genoRR:Double, nullModel:Boolean): Double = {
    val alleleFreq = math.sqrt(pHomozygote)

    var FStat = ((nCtrlHom / nCtrl.toFloat) - pHomozygote) / ( alleleFreq - pHomozygote ).toFloat


    if (FStat < 0) {
      FStat = 0
    }

    val corrPrHom = FStat * alleleFreq + (1- FStat) * pHomozygote

    val denom = ( 1- pHomozygote) + genoRR * corrPrHom
    val x = diseasePrevalence / denom

    var pHomCase:Double = 0.0
    var pHomControl:Double= 0.0

    if (nullModel) {
      pHomCase = corrPrHom
      pHomControl = corrPrHom
    }
    else
    {
      pHomCase = genoRR * corrPrHom / denom
      pHomControl = ((1- (genoRR *x) ) * corrPrHom ) / (1-diseasePrevalence)
    }


    val ll = ( math.log(pHomCase) * nCaseHom + math.log(pHomControl) * nCtrlHom + math.log(1-pHomControl) * (nCtrl-nCtrlHom) + math.log(1-pHomCase) * (nCases-nCaseHom) )

    return -ll
  }

}


