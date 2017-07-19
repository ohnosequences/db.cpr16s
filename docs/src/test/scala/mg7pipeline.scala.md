
```scala
package ohnosequences.db.cpr16s.test

import ohnosequences.mg7._, loquats._
import ohnosequences.datasets._, illumina._
import ohnosequences.cosas._, types._, klists._
import ohnosequences.loquat._, utils._
import ohnosequences.statika._, aws._
import ohnosequences.blast.api._

import ohnosequences.awstools._, regions._, ec2._, s3._, autoscaling._
import com.amazonaws.services.s3.transfer._
import com.amazonaws.auth._, profile._

import better.files._
```


# BLAST reference sequences comparison

This computes the relation which serves as input for the clustering step. We require

1. close to complete query coverage
2. close to 100% identity


```scala
case object mg7 {
```

As the reference database we use the one generated from dropRedundantAssignments

```scala
  case object cpr16sRefDB extends ReferenceDB(
    ohnosequences.db.cpr16s.dbName,
    dropRedundantAssignmentsAndGenerate.s3,
    dropRedundantAssignments.output.table.s3
  )

  case object parameters extends MG7Parameters(
    splitChunkSize = 10,
    splitInputFormat = FastaInput,
    blastCommand = blastn,
    blastOutRec  = defaults.blastnOutputRecord,
    blastOptions = defaults.blastnOptions.update(
      num_threads(2)              ::
      word_size(42)               ::
      evalue(BigDecimal(1E-100))  ::
      max_target_seqs(10000)      ::
      perc_identity(99.0)         ::
      *[AnyDenotation]
    ).value,
    referenceDBs = Set(cpr16sRefDB)
  ) {
```

The only basic thing we require is at least 99% **query** coverage.

```scala
    override def blastFilter(row: csv.Row[BlastOutRecKeys]): Boolean =
      ( row.select(outputFields.qcovs).toDouble >= 99 ) &&
```

IMPORTANT: exclude the query from the results

```scala
      ( row.select(outputFields.qseqid) != row.select(outputFields.sseqid) )
  }

  case object pipeline extends MG7Pipeline(parameters) {
    override lazy val name = "db-cpr16s"

    val metadata: AnyArtifactMetadata = ohnosequences.db.generated.metadata.cpr16s
    // TODO: we should probably have a restricted role for this:
    val iamRoleName: String = "era7-projects"
    val logsS3Prefix: S3Folder = s3"era7-projects-loquats" /
```

As input we use the FASTA accepted by dropRedundantAssignments

```scala
    lazy val inputSamples: Map[ID, S3Resource] = Map(
      "refdb" -> S3Resource(ohnosequences.db.cpr16s.test.dropRedundantAssignments.output.fasta.s3)
    )

    lazy val outputS3Folder: (SampleID, StepName) => S3Folder = { (_, stepName) =>
      ohnosequences.db.cpr16s.s3prefix / "mg7" / stepName /
    }

    val splitConfig  = SplitConfig(1)
    val blastConfig  = BlastConfig(100)
    // these steps are not needed:
    val assignConfig = AssignConfig(20)
    val mergeConfig  = MergeConfig(1)
    val countConfig  = CountConfig(1)
  }
}
```

This bundle just downloads the output of the MG7 Blast step and merges the chunks

```scala
case object mg7BlastResults extends Bundle() {

  lazy val s3location: S3Folder = mg7.pipeline.outputS3Folder("", "blast") / "chunks" /

  lazy val blastChunks: File = File(s3location.key)
  lazy val blastResult: File = (blastChunks.parent / "blastResult.csv").createIfNotExists()

  def instructions: AnyInstructions = LazyTry {
    val transferManager = new TransferManager(new DefaultAWSCredentialsProviderChain())

    transferManager.downloadDirectory(
      s3location.bucket, s3location.key,
      File(".").toJava
    ).waitForCompletion

    transferManager.shutdownNow()
  } -&- LazyTry {

    loquats.mergeDataProcessing().mergeChunks(blastChunks, blastResult)
  }
}

```




[test/scala/dropRedundantAssignments.scala]: dropRedundantAssignments.scala.md
[test/scala/runBundles.scala]: runBundles.scala.md
[test/scala/mg7pipeline.scala]: mg7pipeline.scala.md
[test/scala/package.scala]: package.scala.md
[test/scala/compats.scala]: compats.scala.md
[test/scala/clusterSequences.scala]: clusterSequences.scala.md
[test/scala/dropInconsistentAssignments.scala]: dropInconsistentAssignments.scala.md
[test/scala/pick16SCandidates.scala]: pick16SCandidates.scala.md
[test/scala/releaseData.scala]: releaseData.scala.md
[main/scala/package.scala]: ../../main/scala/package.scala.md
[main/scala/data.scala]: ../../main/scala/data.scala.md