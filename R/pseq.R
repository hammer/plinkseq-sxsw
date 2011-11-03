
###################################################
#                                                 #
# PLINK/SEQ R interface                           #
#                                                 #
###################################################

.First.lib <- function(libname, pkgname) 
{
  library.dynam("Rplinkseq")
}

###################################################
#                                                 #
# Open/re-open projects                           #
#                                                 #
###################################################


pseq.project <- function(s)
{
 tmp <- .Call("Rset_project",s)
}


pseq.summary <- function() {
 .Call("Rsummary");
}


###################################################
#                                                 #
# Variant database                                #
#                                                 #
###################################################


var.attach <- function(s)
{
 tmp <- .Call("Rattach",s)
}


var.detach <- function()
{	
 tmp <- .Call("Rdettach")
}


load.vcf <- function( filename , mask = "" , limit = 100000 )
{
 .Call( "Rdirect_load_vcf" , filename , mask , as.integer( limit ) )
}


###################################################
#                                                 #
# Summary functions (files, individuals, etc)     #
#                                                 #
###################################################


var.files <- function() 
{
 .Call("Rfile_list")
}

ind.fetch <- function( pheno = "" , mask = "" ) 
{
 .Call("Rind_list",mask,pheno)
}

var.headers <- function( filetag = "1" )
{
 .Call("Rhdr_list", as.character(filetag) )
}

var.meta <- function( filetag = "1" )
{

 tmpk <- .Call("Rmeta_list", as.character(filetag) )

 as.data.frame( list( NAME = unlist( lapply( tmpk , "[[" , "NAME" ) ) , 
                      NUM  = as.numeric(unlist( lapply( tmpk , "[[" , "NUM" ) ) ) , 
                      TYPE = unlist( lapply( tmpk , "[[" , "TYPE" ) ) , 
                      GRP  = unlist( lapply( tmpk , "[[" , "GRP" ) ) , 
                      DESC = unlist( lapply( tmpk , "[[" , "DESC" ) ) ) )
	  
}


###################################################
#                                                 #
# Iterate over all variants in database,          #
#  applying an R function                         #
#                                                 #
#   iterate( f_genorate )                         #
#                                                 #
###################################################

var.iterate <- function( func = NULL , mask = "" )
{
 tmp <- .Call("Riterate", as.function(func) , mask, as.integer(0), new.env() )
}

var.fetch <- function( mask = "" , limit = 1000 )
{
 .Call("Riterate", NULL , mask, as.integer(limit), new.env() )
}

meta.fetch <- function( meta , mask = "" )
{
 as.data.frame( .Call("Rgetmeta" , meta , mask , new.env() ) )
}


###################################################
#                                                 #
# Helper extraction functions (from var-lists     #
#                                                 #
###################################################


# Extract base level attribute (CHR, BP, etc)

x.core <- function(l, m) {
  as.vector( unlist(
  lapply( lapply( l , "[[" , m ) ,
                  function(x) ifelse(is.null(x),NA,x) ) ) )
}

x.meta <- function(l,m) {
 as.vector( unlist(
 lapply( lapply( lapply( l , "[[" , "META" )
         , "[[" , m ) ,
        function(x) ifelse(is.null(x),NA,x) ) ) )
}


x.consensus.core <- function(l, m) {
  as.vector( unlist(
  lapply( lapply( lapply( l , "[[" , "CON" )
                  , "[[" , m ) ,
                  function(x) ifelse(is.null(x),NA,x) ) ) )
}


x.consensus.meta <- function(l,m) {
 as.vector( unlist(
 lapply( lapply( lapply( lapply( l , "[[" , "CON" )
         , "[[" , "META" )
         , "[[" , m ) ,
        function(x) ifelse(is.null(x),NA,x) ) ) )
}


x.consensus.genotype <- function(l,m = "GT" ) {
     matrix( unlist( lapply( lapply( l , "[[" , "CON" )
                     , "[[" , c("GENO",m) ) ) ,
                nrow = length(l) , byrow=T )
}


x.sample.meta <- function(l,m) {
 fs <- sort( unique(unlist( lapply( l , function(k2) { t <- c(0) ; for( i in 8:length(k2) ) t <- c( t, k2[[i]]$FSET ) ; return(t) } ) ) ) ) 
 d <- matrix( NA , nrow = length(l) , ncol = length( fs ) )
 i <- 0  
 lapply( l , function(k2) 
   { i <<- i+1; for (j in 7:length(k2) ) 
                 { v <- as.vector(unlist( lapply( lapply( k2[j]  , "[[" , "META" ) , "[[" , m ))); 
                   if ( length(v) > 0 ) d[ i , which( fs == unlist(k2[[j]]$FSET)) ] <<- v } } ) 
  d <- as.data.frame(d)
  names(d) <- c("CON",paste("S",fs[-1],sep=""))
  return(d)
}



###################################################
#                                                 #
# Region/locus functions                          #
#                                                 #
###################################################


loc.attach <- function(s) 
{
 tmp <- .Call("Rlocdb_attach",s)
}

loc.new <- function(s) 
{
 .Call("RopenLocusDB",s);
}

loc.load.GTF <- function(x,y)
{
 tmp <- .Call("Rlocdb_load_gtf",x,y)
}

loc.load.alias <- function(x,y)
{
 tmp <- .Call("Rlocdb_load_alias",x,y)
}

loc.fetch <- function(g)
{
  .Call("Rfetch_regions", g)
}

loc.collapse <- function(x,y)
{
 tmp <- .Call("Rlocdb_collapse_subregions",x,y)
}

loc.summary <- function()
{
 tmp <- .Call("Rlocdb_summary")
}


loc.make.varset <- function(x)
{	  
  tmp <- .Call("Rvardb_make_set",x,x)
}

loc.set.fetch <- function(x,y)
{
 .Call("Rfetch_set_names", x,y )
}

loc.set.fetch.members <- function(x,y,z)
{
 .Call("Rfetch_set_members",x,y,z)
}


###################################################
#                                                 #
# Reference-database                              #
#                                                 #
###################################################

ref.attach <- function(x) {
 tmp <- .Call("Rrefdb_load",x)
}

ref.summary <- function() {
 tmp <- .Call("Rrefdb_summary")
}

ref.lookup <- function(v,g) {
 tmp <- .Call("Rrefdb_lookup",as.integer(c(v$CHR,v$BP)),g)
}


refdb_index_lookup <- function(x) {
 tmp <- .Call("Rrefdb_index_lookup",x)
}


###################################################
#                                                 #
# Sequence/annotation functions                   #
#                                                 #
###################################################


seq.load.FASTA <- function(x) {
 tmp <- .Call("Rseqdb_loadFASTA",x)
}			     

seq.attach <- function(x) {
 tmp <- .Call("Rseqdb_attach",x)
}			     

seq.lookup <- function(x,y,z) { 
 apply( cbind(x,y,z) , 1 , .seq.lookup )
}

.seq.lookup <- function(x) {
 .Call("Rseqdb_lookup", as.integer(x) );
}


annotate.var <- function(x) {
 if ( length(x) == 1 ) 
 {
  v <- x[[1]]
  tmp <- .Call("Rseqdb_annotate",as.integer(c(v$CHR,v$BP)),c(v$REF,v$ALT))
 }
}

annotate.pos <- function(x,y) {
 if ( length(x) == 2 && length(y)==2 ) 
  { 
   tmp <- .Call("Rseqdb_annotate",as.integer(x),y)
  } 
}

annotate.load <- function(x) {
 if ( length(x) == 1 )
  { 
   tmp <- .Call("Rseqdb_annotate_load",x)
  } 
}


###################################################
#                                                 #
# Visualisation                                   #
#                                                 #
###################################################

chrom2chr <-function(s)
{
  s <- sub( "chr" , "" , s )
  s[ s == "X"  ] <- 23
  s[ s == "X"  ] <- 24
  s[ s == "XY" ] <- 25
  s[ s == "M"  ] <- 26
  s[ s == "MT" ] <- 26
  as.numeric( s )
}


chrplot <- function( chr , bp , x , fltr = T , hl = F )
{
 chr <- chr[fltr]; bp <- bp[fltr]; x <- x[fltr]; hl <- hl[fltr]
 mycol <- rep( 1 , length(x) ); mycol[ hl ] <- 2

 chrlen <- matrix(
 c( 558390,247185615, 2994,242739670, 36113,199358925, 8584,191262626, 78452,180648416,
 94609,170853143, 130483,158819535, 153784,146271426, 30910,140229800, 62765,135353683,
 188510,134450734, 25919,132288869, 17929308,114125098, 18689850,106359422, 18274994,100323631,
 24170,88691089, 13043,78638558, 1523,76116152, 204938,63788762, 9795,62385675, 9887804,46924583,
 14431347,49576671 ) , byrow = 2 , nrow = 22 )

 chr <- chrom2chr(chr)
 minchr <- min(chr); maxchr <- max(chr)
 minx   <- min(x); maxx   <- max(x)
 ytot   <- 10000
 int    <- ytot / (length(unique(chr)))
 int2   <- 0.90* int

 uchr <- unique(chr); achr <- numeric()
 for (i in 1:length(uchr)) achr[ uchr[i] ] <- i

 par( mar=c(0.5,0.5,0.5,0.5), oma=c(0,0,0,0))

 plot( c(1,1) ,type="n" , 
       ylim=c(0,ytot) , 
       xlim=c(min(chrlen[uchr,]),max(chrlen[uchr,])) , 
       xlab="", ylab="" , xaxt="n" , yaxt="n" )

 for (i  in unique(chr) )
  lines( c( chrlen[i,1] , chrlen[i,2] ) , 
         c(int*(achr[i]-1),int*(achr[i]-1)) , 
         col="gray" )

 for (i in 1:length(x))
  lines( c(bp[i],bp[i]) , 
         c( int*(achr[chr[i]]-1) , int*(achr[chr[i]]-1) + (int2)*( (x[i]-minx)/(maxx-minx)) ) , 
         col= mycol[i] )
}



###################################################
#                                                 #
# Misc/legacy                                     #
#                                                 #
###################################################

variant <- function(n)
{		
 .Call("Ridx_variant",as.integer(n))
}


