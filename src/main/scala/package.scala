package ohnosequences.db

import ohnosequences.awstools.s3._

package object cpr16s {

  val dbName = "ohnosequences.db.cpr16s"

  private val metadata = generated.metadata.cpr16s

  val s3prefix: S3Folder =
    s3"resources.ohnosequences.com" /
    metadata.organization /
    metadata.artifact /
    metadata.version /
}
