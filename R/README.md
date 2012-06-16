The *Rplinkseq* R-package will build as a self-contained package (including
*libplinkseq*). Further it will build itself in such a way that its internal
library (containing all of *libplinkseq*) can be linked against by other
downstream R packages. However, *libprotobuf*, on which Plink/Seq depends,
[does not allow itself to be dynamically loaded multiple times in the same
process](http://code.google.com/p/protobuf/issues/detail?id=128). Thus when
both *Rplinkseq* and the downstream package load *libprotobuf* the process
aborts. The solution is to statically link *libprotobuf* into *Rplinkseq* by
doing the following:

    R CMD INSTALL --configure-args='--with-static-protobuf=/usr/local/lib/libprotobuf.a' Rplinkseq

substituting the path to *libprotobuf* on your system. And to force downstream
packages to dynamically link against *libRplinkseq*. Note this is only
neccessary if you intend to link against Plink/Seq library embedded in
*Rplinkseq* in your downstream package; for most users this won't be an issue.
