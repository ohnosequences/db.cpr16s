package ohnosequences.db.cpr16s.test

import ohnosequences.statika._, aws._
import ohnosequences.awstools._, regions._, ec2._, autoscaling._, s3._

case object compats {

  type DefaultAMI = AmazonLinuxAMI[Ireland.type, HVM.type, InstanceStore.type]
  val  defaultAMI: DefaultAMI = AmazonLinuxAMI(Ireland, HVM, InstanceStore)

  class DefaultCompatible[B <: AnyBundle](bundle: B, javaHeap: Int) extends Compatible(
    amznAMIEnv(defaultAMI, javaHeap),
    bundle,
    ohnosequences.db.generated.metadata.cpr16s
  )

  case object pick16SCandidates extends
    DefaultCompatible(ohnosequences.db.cpr16s.test.pick16SCandidates, javaHeap = 50)

  case object dropRedundantAssignmentsAndGenerate extends
    DefaultCompatible(ohnosequences.db.cpr16s.test.dropRedundantAssignmentsAndGenerate, javaHeap = 10)

  case object clusterSequences extends
    DefaultCompatible(ohnosequences.db.cpr16s.test.clusterSequences, javaHeap = 10)

  case object dropInconsistentAssignmentsAndGenerate extends
    DefaultCompatible(ohnosequences.db.cpr16s.test.dropInconsistentAssignmentsAndGenerate, javaHeap = 10)
}
