package org.broadinstitute.hail.driver

//import org.kohsuke.args4j.Option
import org.kohsuke.args4j.{Option => Args4jOption}
import org.broadinstitute.hail.Utils._
import org.broadinstitute.hail.variant._
import org.apache.spark.rdd.RDD

import scala.math._
import org.broadinstitute.hail.annotations.{Annotation, _}
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.PowellOptimizer
import org.apache.commons.math3.optim.univariate.BrentOptimizer
import org.apache.commons.math3.optim.univariate.UnivariateObjectiveFunction
import org.apache.commons.math3.optim._
import org.apache.commons.math3.analysis.UnivariateFunction
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType
import org.apache.commons.math3.distribution.ChiSquaredDistribution
import org.broadinstitute.hail.expr.{TBoolean, TDouble, TInt, TStruct, Type}
import org.apache.commons.math3.special._
import org.apache.commons.math3.optim.univariate.SearchInterval
import org.apache.spark.HashPartitioner
import org.apache.spark.sql.catalyst.plans.physical.HashPartitioning
import org.broadinstitute.hail.methods.{LinRegStats, LinearRegression, Phasing}
import org.broadinstitute.hail.utils.{SparseVariantSampleMatrix, SparseVariantSampleMatrixRRDBuilder}

import scala.collection.mutable
import scala.collection.mutable.ArrayBuffer



/**
  * Created by mitja on 7/22/16.
  */


object RAFT extends Command {


  def name = "RAFT"

  def description = "Performs recessive test RAFT ( Lim et al. AJHG 2014)"
  class Options extends BaseOptions {
    @Args4jOption(required = true, name = "-o", aliases = Array("--output"), usage = "Output file")
    var output: String = _

    @Args4jOption(required = true, name = "-p", aliases = Array("--phenotype"), usage = "Phenotype names (1=control, 2=case) everything else ignored")
    var pheno_name: String = _

    @Args4jOption(required = false, name = "-g", aliases = Array("--group_by"), usage = "Variant annotation to group variants by. If not specified then single variant analysis is performed.")
    var group_by: String = _

  }

  def newOptions = new Options

  override def supportsMultiallelic = false

  def requiresVDS = true

  //  seqOp: (U, Variant, Annotation, String, Annotation, T) => U,
  //  combOp: (U, U) => U)(implicit uct: ClassTag[U]): RDD[(Variant, U)] = {

  def calcRAFT(state:State, pheno:String): RDD[(Variant, RaftCounter)] = {

    val vds = state.vds

    val (base_type,pheno_query) = vds.querySA(pheno)
    val saBc = state.sc.broadcast(vds.sampleAnnotations)

    val var_rdd = vds.aggregateByVariantWithAll( new RaftCounter) ( {
      case (counter, v,va, s, sa, g ) => {
        val is_case = pheno_query(sa)

        counter.consume_geno(g, is_case)
      }
    }, { case (r1, r2) => r1.merge(r2)
    }
    ).mapValues( raftcount => raftcount.calculate_RAFT() )

    return var_rdd
  }


  def calcRAFT_by_annotation(state:State, pheno:String, annot:String, n_partitions:Int=200 ): RDD[(String,RaftAnnotationCounter) ] = {

    val vds = state.vds

    val (base_type,pheno_query) = vds.querySA(pheno)
    val (va_type,va_query) = vds.queryVA(annot)

    val partitioner = new HashPartitioner(n_partitions)

    info("Counting by anno " + annot)
    val samplesPhenos = state.sc.broadcast(for (i <- state.vds.sampleIds.indices) yield {
        pheno_query(state.vds.sampleAnnotations(i))
        }
    )

    //Array(pheno)

    SparseVariantSampleMatrixRRDBuilder.buildByVA(  vds, state.sc,   partitioner ) ( {
        case(v, va) =>  {
          info("MAPPING ANNOTATION")
          val anno = va_query(va).get.toString()
          info("printing anno")
        }
      }).map( {
      case( anno, svm) =>
        {
          (anno.toString, new RaftAnnotationCounter(anno.toString, svm, pheno))
        }
    })


    /*
    val val_rdd = vds.aggregateByAnnotation( new RaftAnnotationCounter() ) ( {

      case (rcounter, v, va, sampleid, sa, sample_indx, g) => {
        val one =1
        return Unit
      }
    }, {
      case( counter1,counter2) => counter1.merge(counter2)
    }, {
      case( v, va) =>  va_query(va)
    }
    )
    */

  }



  def run(state: State, options: Options): State = {

    //vds.rdd.map  { case (v,va,gs)  => gs.iterator.zip(bc.value.iterator).foldLeft(  ) }
    val pheno = options.pheno_name

    val group_by = options.group_by

    //val phenos = get_phenotypes(state, pheno)

    val output = options.output

  // why is this TBoolean not found even if imported?
    // state.vds.saSignature.getAsOption[ TBoolean ](pheno) match {



/*
    state.vds.saSignature.getOption(pheno) match {
      case None => {
        println("You need to have sample annotation with the given name " + pheno + " and type Boolean prior to running RAFT")
        throw new IllegalArgumentException("You need to have sample annotation with the given name " + pheno + " and type Boolean prior to running RAFT")
      }
    }
  */

    /*
    if ( group_by != null ) {
      state.vds.vaSignature.getOption(group_by) match {
        case None => {
          println("You need to have variant annotation with the given name " + group_by + " for the variants to be grouped by prior to running RAFT with group_by option")
          throw new IllegalArgumentException("You need to have variant annotation with the given name " + group_by + " for the variants to be grouped by prior to running RAFT with group_by option")
        }
      }
    }
    */

    //val (base_type,pheno_query) = state.vds.querySA(pheno)

    if ( group_by == null ) {


      info(s"Running single variant RAFT for phenotype " + pheno + " with " + state.vds.nSamples + " samples and " + state.vds.nVariants + " variants.")
      val resRDD = calcRAFT(state, pheno)

      if (output != null) {
        hadoopDelete(output, state.hadoopConf, recursive = true)
        resRDD.map { case (v, raft_stat) =>
          val sb = new StringBuilder()
          sb.append(v.contig)
          sb += '\t'
          sb.append(v.start)
          sb += '\t'
          sb.append(v.ref)
          sb += '\t'
          sb.append(v.alt)
          sb += '\t'
          raft_stat.result_string(sb)
          sb.result()
        }.writeTable(output, Some("Chrom\tPos\tRef\tAlt\t" + RaftCounter.header))
      }

      val vds = state.vds

      val (newVAS, inserter) = vds.insertVA(RaftCounter.signature, "raft")

      val newstate = state.copy(
        vds = vds.copy(
          rdd = vds.rdd.zipPartitions(resRDD) { case (it, jt) =>

            it.zip(jt).map { case ((v, va, gs), (v2, comb)) =>
              assert(v == v2)

              (v, inserter(va, comb.toAnnotation), gs)
            }
          },
          vaSignature = newVAS
        )
      )

      val sb = new StringBuilder()
      newVAS.pretty(sb, 1)
      print(sb)

      return newstate
    }
    else
    {
      info(s"Running RAFT by annotation " + group_by+ " for phenotype " + pheno + " with " + state.vds.nSamples + " samples and " + state.vds.nVariants + " variants.")
      val resRDD = calcRAFT_by_annotation( state, pheno, group_by)
      return state
    }

  }

}


class RaftAnnotationCounter( annot: String, svm:SparseVariantSampleMatrix, pheno: String ) {


  var n_cases: Int = 0
  var n_controls: Int = 0
  var n_missing_pheno: Int = 0


  info("raft counter for " + annot)

  private val phaseCache = mutable.Map[(String,String),Option[Double]]()

  // holds mapping from sample id to list of variant ids the sample carries (either homozygotes, or compound hets)
  var homvars =  Map[ String, ArrayBuffer[String] ]().withDefaultValue( new ArrayBuffer[String])

  svm.mapSamples( {

    case( sample_id, sample_idx, variant_idxs, genotypes ) => {
        val phenotype = svm.getSampleAnnotation(sample_id, pheno)

        print("sample anno of " + sample_id + " is " + phenotype)
        phenotype match {
          case Some(true) => n_cases+=1
          case Some(false) => n_controls+=1
          case None => n_missing_pheno+=1
        }

        val hetSites = new ArrayBuffer[Int]

        variant_idxs.indices.foreach({
          i =>
          val vi = variant_idxs(i)

          genotypes(i) match{
            case 1 =>
              hetSites += vi

            case 2 =>
              homvars( sample_id) += svm.variants( vi )

            case _ =>
          }
        })

        hetSites.indices.foreach( {
          i1 =>
            val v1i = hetSites(i1)
            Range(i1 + 1, hetSites.size).foreach({
              i2 =>
                val v2i = hetSites(i2)

                val var1 = svm.variants(v1i)
                val var2 = svm.variants(v2i)

                getPhase( var1, var1 ) match {
                  case Some(p_same_hap) =>
                    if( p_same_hap < 0.5 ) {
                      val comp_var_name = if (var1 < var2) var1 + "|" + var2 else var2 + "|" + var1
                      homvars( sample_id) += comp_var_name
                    }
                }
            }
            )
        }
        )



    }

  })

  private def getPhase(v1: String, v2: String) : Option[Double] = {
    phaseCache.get((v1,v2)) match{
      case Some(p) => p
      case None =>
        val pPhase = Phasing.probOnSameHaplotypeWithEM(svm.getGenotypeCounts(v1,v2))
        phaseCache.update((v1,v2),pPhase)
        pPhase
    }
  }

  def count_variants() {
    val first =1
  }

}



object RaftCounter {
  val header =
    "n_hom_ref_case\t" +
    "n_het_case\t" +
    "n_hom_alt_case\t" +
    "n_missing_case\t" +
    "n_hom_ref_control\t" +
    "n_het_control\t" +
    "n_hom_alt_control\t" +
    "n_missing_control\t" +
    "geno_RR\t" +
    "log_LR\t" +
    "p_val"

  def signature: Type = TStruct(
    ("geno_RR", TDouble),
    ("log_LR", TDouble),
    ("p_val", TDouble),
    ("n_hom_ref_case", TInt),
    ("n_het_case", TInt),
    ("n_hom_alt_case", TInt),
    ("n_missing_case", TInt),
    ("n_hom_ref_control", TInt),
    ("n_het_control", TInt),
    ("n_hom_alt_control", TInt),
    ("n_missing_control", TInt)

  )


}

class RaftCounter {

  var nHomRefCase=0
  var nHetCase=0
  var nHomAltCase=0
  var nMissingCase=0

  var nHomRefControl=0
  var nHetControl=0
  var nHomAltControl=0
  var nMissingControl=0

  var nMissingPheno=0

  var pVal: Option[Double] = None
  var raftStat: Option[Double] = None
  var genoRR: Option[Double] = None

  def toAnnotation: Some[Annotation] = Some(Annotation(genoRR, raftStat, pVal, nHomRefCase,nHetCase,nHomAltCase,nMissingCase,nHomRefControl, nHetControl,nHomAltControl, nMissingControl))

  def merge( o:RaftCounter  ):RaftCounter = {

    nHomRefCase +=  o.nHomRefCase
    nHomRefControl +=  o.nHomRefControl
    nHetCase +=  o.nHetCase
    nHetControl +=  o.nHetControl
    nHomAltCase +=  o.nHomAltCase
    nHomAltControl +=  o.nHomAltControl
    nMissingCase +=  o.nMissingCase
    nMissingControl +=  o.nMissingControl
    nMissingPheno += o.nMissingPheno
    this
  }

  def consume_geno( g:Genotype, is_case:Option[Any] ):RaftCounter = {

    if ( g.isHomRef && is_case == Some(true)) {
      nHomRefCase += 1
    }
    else if ( g.isHet && is_case == Some(true) ) {
      nHetCase += 1
    }
    else if ( g.isHomVar && is_case == Some(true) ) {
      nHomAltCase +=1
    }
    else if ( g.isHomRef && is_case == Some(false)) {
      nHomRefControl += 1
    }
    else if ( g.isHet && is_case == Some(false) ) {
      nHetControl += 1
    }
    else if ( g.isHomVar && is_case == Some(false) ) {
      nHomAltControl +=1
    }
    else if (g.isNotCalled && is_case == Some(true)  ) {
      nMissingCase +=1
    }
    else if (g.isNotCalled && is_case == Some(false)  ) {
      nMissingControl +=1
    }
    else if ( is_case.isEmpty ) {
      nMissingPheno +=1
    }

    this
  }

  def calculate_RAFT(): RaftCounter = {

    val nCtrlAlleles = (nHomRefControl+nHetControl+nHomAltControl)*2

    var afALT = (nHetControl + 2*nHomAltControl).toFloat / (( nHomRefControl+nHetControl+nHomAltControl)*2).toFloat

    if (afALT > 0.5)
    {
      // switch hom as the var of interest.
      var tmp = nHomRefControl
      nHomRefControl = nHomAltControl
      nHomAltControl = tmp

      tmp = nHomRefCase
      nHomRefCase = nHomAltCase
      nHomAltCase = tmp
      afALT = (nHetControl + 2*nHomAltControl).toFloat / (( nHomRefControl+nHetControl+nHomAltControl)*2).toFloat
    }

    if ( nCtrlAlleles==0 || nHomAltCase==0) {
      return this
    }


    if( nHomAltControl + nHetControl == 0 )
      afALT = 1 / (1 + (2*nHomRefControl+nHetControl+nHomAltControl)).toFloat

    val pHomozygote = math.pow( afALT , 2)
    val diseasePrevalence = 0.000001

    var enrichment = 1.0
    val expHomCases = pHomozygote *( nHomAltCase + nHetCase + nHomRefCase )
    val expHomControls = pHomozygote *( nHomAltControl + nHetControl + nHomRefControl )

    if (nHomAltCase>0 && nHomAltControl>0)
      enrichment = (nHomAltCase/expHomCases). / ( nHomAltControl /expHomControls)

    else if (nHomAltCase>0)
      enrichment = nHomAltCase / expHomCases

    if (enrichment<1) {
      return this
    }
    val rf = new RAFTFunction( nHomAltCase, nHomRefCase+nHetCase+nHomAltCase, nHomAltControl, nHomRefControl+nHetControl+nHomAltControl, pHomozygote, diseasePrevalence)
    val of = new UnivariateObjectiveFunction(rf)

    val optim = new BrentOptimizer(1e-6, 1e-8)

    val result = optim.optimize(of, GoalType.MINIMIZE, new MaxIter(1000), new InitialGuess(Array(enrichment)), new MaxEval(1000), new SearchInterval(1,100000)  )

    val minLL = result.getValue()
    val nullLL = rf.loglik(1,true)

    raftStat = Some(-minLL + nullLL)
    genoRR = Some(result.getPoint)

    val chiSqr = new ChiSquaredDistribution(1)

    pVal = Some(Gamma.regularizedGammaQ( 0.5, raftStat.get.toDouble))

    return this
  }

  def result_string(sb: StringBuilder) = {

    sb.append(nHomRefCase)
    sb.append("\t")
    sb.append(nHetCase)
    sb.append("\t")
    sb.append(nHomAltCase)
    sb.append("\t")
    sb.append(nMissingCase)
    sb.append("\t")

    sb.append(nHomRefControl)
    sb.append("\t")
    sb.append(nHetControl)
    sb.append("\t")
    sb.append(nHomAltControl)
    sb.append("\t")
    sb.append(nMissingControl)
    sb.append("\t")


    sb.append( genoRR.getOrElse("NA") )
    sb.append("\t")
    sb.append( raftStat.getOrElse("NA") )
    sb.append("\t")
    sb.append( pVal.getOrElse("NA") )
    sb.append("\t")

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
