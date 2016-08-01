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

import org.broadinstitute.hail.methods.RAFTStats

import scala.collection.mutable.Map
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

    /*
    val samplesPhenos = state.sc.broadcast(for (i <- state.vds.sampleIds.indices) yield {
        pheno_query(state.vds.sampleAnnotations(i))
        }
    )
    */


    val res_rdd = SparseVariantSampleMatrixRRDBuilder.buildByVAstoreAllVAandSA(  vds, state.sc,   partitioner ) ( {
        case(v, va) =>  va_query(va).get.toString
    }).map( {
      case( anno, svm) =>
        {
          (anno.toString, new RaftAnnotationCounter(anno.toString, svm, pheno))
        }
    })

    res_rdd

  }



  def run(state: State, options: Options): State = {

    //vds.rdd.map  { case (v,va,gs)  => gs.iterator.zip(bc.value.iterator).foldLeft(  ) }
    val pheno = options.pheno_name

    val group_by = options.group_by

    //val phenos = get_phenotypes(state, pheno)

    val output = options.output


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


      return newstate
    }
    else
    {
      info(s"Running RAFT by annotation " + group_by+ " for phenotype " + pheno + " with " + state.vds.nSamples + " samples and " + state.vds.nVariants + " variants.")
      val resRDD = calcRAFT_by_annotation( state, pheno, group_by)
      resRDD.mapValues( {
        case (raft_counter) =>
          print("stuff for gene. cases: " + raft_counter.n_cases + ", controls:" +   raft_counter.n_controls + ", missing_pheno" +
            + raft_counter.n_missing_pheno + " raft stat:" + raft_counter.calculateRAFT() + "cum PHom:" + raft_counter.cumulativePHom )
          raft_counter.calculateRAFT()
          raft_counter

      } )
      if (output != null) {
        hadoopDelete(output , state.hadoopConf, recursive = true)
        resRDD.map { case (g, raft_stat) =>
          val sb = new StringBuilder()
          raft_stat.result_string(sb)
          sb.result()
        }.writeTable(output, Some(RaftAnnotationCounter.header))
      }


      return state
    }

  }

}


class RaftAnnotationCounter( annot: String, svm:SparseVariantSampleMatrix, pheno: String ) {


  var n_cases: Int = 0
  var n_controls: Int = 0
  var n_missing_pheno: Int = 0


  val dis_prevalence = 0.000001

  private val phaseCache = mutable.Map[(String, String), Option[Double]]()

  var cumulativePHom = 0.0

  var pVal: Option[Double] = None
  var raftStat: Option[Double] = None
  var genoRR: Option[Double] = None

  // mapping from variant to array of (allelecount, allelenumber, missing genotypes)
  //var variant_freqs = mutable.Map[String, Array[Int]]().withDefaultValue( Array.fill[Int](3)(0))

  var variant_freqs = mutable.Map.empty[String, Array[Int]]


  var case_homvars = mutable.Map[String, ArrayBuffer[String]]()
  var ctrl_homvars = mutable.Map[String, ArrayBuffer[String]]()


  var case_hets = mutable.Set.empty[String]
  var ctrls_hets = mutable.Set.empty[String]


  var unique_case_hets: Int = 0
  var unique_ctrl_hets: Int = 0
  var unique_case_hom_ref: Int = 0
  var unique_ctrl_hom_ref: Int = 0

  var cumAF_counts = 0.0

  compute_counts()
  calculateRAFT()


  private def compute_counts() {

    svm.foreachSample({
      case (sample_id, sample_idx, variant_idxs, genotypes) => {
        val phenotype = svm.getSampleAnnotation(sample_id, pheno)

        phenotype match {
          case Some(true) => n_cases += 1
          case Some(false) => n_controls += 1
          case None => n_missing_pheno += 1
        }

        if (phenotype != None) {
          val hetSites = new ArrayBuffer[Int]

          variant_idxs.indices.foreach({
            i =>
              val vi = variant_idxs(i)
              var alleles = 0
              var chroms = 0
              var missing = 0
              genotypes(i) match {
                case 1 =>
                  hetSites += vi
                  alleles = 1
                  chroms = 2

                  if (phenotype == Some(false)) {
                    this.ctrls_hets.add(sample_id)
                  } else if (phenotype == Some(true)) {
                    this.case_hets.add(sample_id)
                  }

                case 2 =>

                  phenotype match {
                    case Some(true) =>
                      if (!case_homvars.contains(sample_id) )
                        case_homvars(sample_id) = new ArrayBuffer[String]

                      case_homvars(sample_id) += svm.variants(vi)
                    case Some(false) =>
                      if (!ctrl_homvars.contains(sample_id) )
                        ctrl_homvars(sample_id) = new ArrayBuffer[String]
                      ctrl_homvars(sample_id) += svm.variants(vi)
                  }

                  alleles = 2
                  chroms = 2
                case -1 =>
                  missing = 1

              }

              if (phenotype == Some(false)) {
                // variant freqs from controls
                val varname =svm.variants(vi)

                if(!variant_freqs.contains(varname)) {
                  variant_freqs += (varname -> Array.fill[Int](3)(0))
                }

                var v = variant_freqs(varname)
                v(0) += alleles
                v(1) += chroms
                v(2) += missing
                variant_freqs.update(varname, v)
              }

          })

          hetSites.indices.foreach({
            i1 =>
              val v1i = hetSites(i1)
              Range(i1 + 1, hetSites.size).foreach({
                i2 =>
                  val v2i = hetSites(i2)

                  val var1 = svm.variants(v1i)
                  val var2 = svm.variants(v2i)

                  getPhase(var1, var1) match {
                    case Some(p_same_hap) =>
                      if (p_same_hap < 0.5) {
                        val comp_var_name = if (var1 < var2) var1 + "|" + var2 else var2 + "|" + var1

                        //phenotype match {
                          //case Some(true) => case_homvars(sample_id) += comp_var_name
                          //case Some(false) => ctrl_homvars(sample_id) += comp_var_name
                        //}

                      }
                  }
              })
          })

        }
      }
    })
  }


  private def getPhase(v1: String, v2: String) : Option[Double] = {
    phaseCache.get((v1,v2)) match{
      case Some(p) => p
      case None =>
        val pPhase = Phasing.probOnSameHaplotypeWithEM(svm.getGenotypeCounts(v1,v2))
        phaseCache.update((v1,v2),pPhase)
        pPhase
    }
  }


  /**
    *
    * @return tuple of (raft_Stat, genotypic relative risk, p-value)
    */
  def calculateRAFT():(Option[Double], Option[Double], Option[Double]) = {

    this.cumulativePHom = this.variant_freqs.map( {
      case (v, counts) =>
        Math.pow(counts(0).toDouble / (2*( n_controls - counts(2)  )).toDouble, 2)
    }).sum

    if( this.raftStat == None ) {
      print("print calculating raft. N vars:" + this.variant_freqs.size + "case homs " + case_homvars.size + " ctrl homs" + ctrl_homvars.size)
      val res = RAFTStats.raftStats(  case_homvars.size,n_cases, n_controls, ctrl_homvars.size, cumulativePHom, dis_prevalence  )
      this.raftStat = res._1
      this.genoRR = res._2
      this.pVal = res._3


      this.unique_case_hets = this.case_hets.filter( (iid ) => !this.case_homvars.contains(iid) ).size
      this.unique_ctrl_hets = this.ctrls_hets.filter( (iid ) => !this.ctrl_homvars.contains(iid) ).size

      this.unique_case_hom_ref = this.n_cases - this.unique_case_hets - this.case_homvars.size
      this.unique_ctrl_hom_ref = this.n_controls - this.unique_ctrl_hets - this.ctrl_homvars.size

      this.cumAF_counts =  ( unique_ctrl_hets + 2* this.ctrl_homvars.size ).toDouble / (2*this.n_controls).toDouble

    }

    return ( this.raftStat ,this.genoRR, this.pVal)
  }

  def result_string(sb:StringBuilder):Unit =  {


    sb.append(this.unique_case_hom_ref )
    sb.append( "\t" )
    sb.append( this.unique_case_hets )
    sb.append( "\t" )
    sb.append( this.case_homvars.size )
    sb.append( "\t" )
    sb.append( this.unique_ctrl_hom_ref )
    sb.append( "\t" )
    sb.append( this.unique_ctrl_hets )
    sb.append( "\t" )
    sb.append( this.ctrl_homvars.size )
    sb.append( "\t" )
    sb.append( this.annot  )
    sb.append( "\t" )
    sb.append( this.variant_freqs.size.toString )
    sb.append( "\t" )
    sb.append( this.cumulativePHom.toString )
    sb.append( "\t" )
    sb.append( this.cumAF_counts.toString )
    sb.append( "\t" )
    sb.append( this.case_homvars.map{ case (iid,vars) => (iid,iid + "," + vars.result().mkString(",")) }.valuesIterator.mkString(";") )
    sb.append( "\t" )
    sb.append( this.ctrl_homvars.map{ case (iid,vars) => (iid,iid + "," + vars.result().mkString(",")) }.valuesIterator.mkString(";") )
    sb.append( "\t" )

    sb.append( this.genoRR.getOrElse("NA") )
    sb.append( "\t" )

    sb.append( this.raftStat.getOrElse("NA") )
    sb.append( "\t" )

    sb.append( this.pVal.getOrElse("NA")  )

  }

}

object  RaftAnnotationCounter {
  val header = "n_hom_ref_case\t" +
    "n_het_case\t" +
    "n_hom_alt_case\t" +
    "n_hom_ref_control\t" +
    "n_het_control\t" +
    "n_hom_alt_control\t" +
    "grouped_by\t" +
    "n_vars\t" +
    "cum_p_hom\t" +
    "cum_af_counts\t" +
    "case_carriers\t" +
    "control_carriers\t" +
    "geno_RR\t" +
    "log_LR\t" +
    "p_val"
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

    val raft = RAFTStats.raftStats( nHomAltCase, nHomRefCase+nHetCase+nHomAltCase, nHomAltControl, nHomRefControl+nHetControl+nHomAltControl, pHomozygote, diseasePrevalence)

    raftStat = raft._1
    genoRR = raft._2
    pVal = raft._3

    return this
  }




}



