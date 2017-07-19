/*
  # Pick 16S candidates

  In this first step we pick those sequences which contain a 16S gene sequence, based on their length and annotations. From the taxonomical point of view, we are only interested in sequences with at least one assignment to (a descendant of) *Bacteria* or *Archaea* which is *not* a descendant of "unclassified" taxa.

  The output of this step represents around `5%` of the sequences in RNACentral.
*/
package ohnosequences.db.cpr16s.test

import ohnosequences.db._, csvUtils._, collectionUtils._
import ohnosequences.db.rnacentral._, RNAcentral._
import ohnosequences.ncbitaxonomy._, titan._
import ohnosequences.fastarious.fasta._
import ohnosequences.statika._
import com.github.tototoshi.csv._
import better.files._

case object pick16SCandidates extends FilterData(
  RNAcentral.table,
  RNAcentral.fasta,
  ohnosequences.db.cpr16s.s3prefix
)(
  deps = ncbiTaxonomyBundle
)
{

  /* We are using the ribosomal RNA type annotation on RNACentral as a first catch-all filter. We are aware of the existence of a gene annotation corresponding to 16S, that we are **not using** due to a significant amount of 16S sequences lacking the corresponding annotation. */
  val ribosomalRNAType = "rRNA"

  /* Taxon ID for *Saccharibacteria* */
  val SaccharibacteriaTaxonID = "95818"
  /*
    ## Predicate defining a candidate

    Sequences that satisfy this predicate (on themselves together with their annotation) are included in the output of this step.
  */
  private lazy val taxonomyGraph = ohnosequences.ncbitaxonomy.ncbiTaxonomyBundle.graph

  def rowPredicate(row: Row): Boolean = {
    val taxID = row.select(tax_id)

    /* - are annotated as rRNA */
    row.select(rna_type).contains(ribosomalRNAType) &&
    /* - are *not* from the SILVA database */
    ( row.select(db).trim.toLowerCase != "silva" )  &&
    {
      taxonomyGraph.getTaxon(taxID).map(_.ancestors) match {
        case None => false // not in the DB
        case Some(ancestors) =>
    /* - is a descendant of Saccharibacteria */
          ancestors.exists { _.id  == SaccharibacteriaTaxonID }
      }
    }
  }

  /* This predicate determines whether the *sequence* value is OK, and will be kept. */
  def sequencePredicate(fastaSeq: FastaSequence): Boolean = {
    val seq = fastaSeq.value

    ( !(seq containsSlice "NNNNNNNN") )             &&
    ( (seq.count(_ == 'N') / seq.length) <= 0.01 )
  }

  // bundle to generate the DB (see the runBundles file in tests)
  def filterData(): Unit = {

    val groupedRows: Iterator[(String, Seq[Row])] =
      source.table.tsvReader.iterator.contiguousGroupBy { _.select(id) }

    val fastas: Iterator[FASTA.Value] = source.fasta.parsed.toIterator

    (groupedRows zip fastas) foreach { case ((commonID, rows), fasta) =>

      if (commonID != fasta.getV(header).id)
        sys.error(s"ID [${commonID}] is not found in the FASTA. Check RNACentral filtering.")

      val (acceptedRows, rejectedRows) =
        if ( sequencePredicate(fasta.get(sequence).value) ) {
          /* if the sequence is OK, we partition the rows based on the predicate */
          rows.partition(rowPredicate)
        } else {
          (Seq[Row](), rows)
        }

      val extendedID: String = s"gnl|${ohnosequences.db.cpr16s.dbName}|${commonID}"

      writeOutput(
        extendedID,
        acceptedRows.map{ _.select(tax_id) }.distinct,
        rejectedRows.map{ _.select(tax_id) }.distinct,
        fasta.update( header := FastaHeader(Seq(extendedID, fasta.getV(header).description).mkString(" ") ) )
      )
    }
  }
}
