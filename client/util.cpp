
#include "plinkseq/helper.h"
#include "util.h"

// so we can register commands, need to pull in other headers across PSEQ for regargs()
#include "ibs.h"


extern Pseq::Util::Commands pcomm;

enum Pseq::Util::Options::type_t types;

void Pseq::Util::Options::shortform( const std::string & sht , const std::string & lng ) 
{
  shortcuts[ sht ] = lng;
}


void Pseq::Util::populate_commands( Pseq::Util::Commands & pcomm )
{

    // command|group|description|VCF|GRP|ARG
  
    // add groups
    
    pcomm.new_group( "root"        , "All commands groups" );
    pcomm.new_group( "input"       , "Data input" );
    pcomm.new_group( "output"      , "Variant data output" );
    pcomm.new_group( "project"     , "Project functions" );
    pcomm.new_group( "stats"       , "Variant summary statistics" );
    pcomm.new_group( "tests"       , "Genotype-phenotype association tests" );
    pcomm.new_group( "qc"          , "Quality control metrics and tests" );
    pcomm.new_group( "views"       , "Viewing variant and other data" );
    pcomm.new_group( "annot"       , "Annotation functions" );
    pcomm.new_group( "varop"       , "Variant database operations" );
    pcomm.new_group( "locop"       , "Locus database operations" );
    pcomm.new_group( "refop"       , "Reference database operations" );
    pcomm.new_group( "segop"       , "Segment database operations" );
    pcomm.new_group( "seqop"       , "Sequence database operations" );
    pcomm.new_group( "indop"       , "Individual database operations" );
    pcomm.new_group( "ibd"         , "IBD analysis" );
    pcomm.new_group( "net"         , "Network-based analysis" );
    pcomm.new_group( "prot"        , "Protein domain/feature annotation" );
    pcomm.new_group( "misc"        , "Misc." );
    pcomm.new_group( "system"      , "System functions" ); // not shown in GUI
    
    pcomm.add_to_group( "root" , "input" );
    pcomm.add_to_group( "root" , "output" );
    pcomm.add_to_group( "root" , "project" );
    pcomm.add_to_group( "root" , "stats" );  
    pcomm.add_to_group( "root" , "tests" );
    pcomm.add_to_group( "root" , "qc" );
    pcomm.add_to_group( "root" , "views" );
    pcomm.add_to_group( "root" , "annot" );
    pcomm.add_to_group( "root" , "varop" );
    pcomm.add_to_group( "root" , "locop" );
//  pcomm.add_to_group( "root" , "segop" );
    pcomm.add_to_group( "root" , "seqop" );
    pcomm.add_to_group( "root" , "refop" );
    pcomm.add_to_group( "root" , "indop" );
    pcomm.add_to_group( "root" , "ibd" );
    pcomm.add_to_group( "root" , "net" );
    pcomm.add_to_group( "root" , "prot" );
    pcomm.add_to_group( "root" , "misc" );
    pcomm.add_to_group( "root" , "system" );
    
    
    // add commands (will be grouped in order of entry)
    
    pcomm << "commands|misc|list commands/groups"
	
	  << "masks|misc|list mask options"
	
	  << "new-project|project,input|set a new project|ARG:vcf,resources,locdb,refdb,vardb,seqdb,inddb,output,scratch,metameta"
	
	  << "version|project|display version information"
	
	  << "append|project,input|add a file to the project|ARG:name,file,type"
	
	  << "drop|project|drop a file from the project|ARG:name,file,type"
	
	  << "load-vcf|input|load all VCF files not already in VARDB|ARG:file,vcf,filter,include-meta,exclude-meta"
	
	  << "index-vcf|input|add index to VARDB for a BGZF-compressed VCF|ARG:vcf"
	
	  << "*index-bcf|input|add index to VARDB for a BCF|ARG:bcf"

	  << "*reload-vcf|input|clear VARDB, then reload all VCF (not implemented yet)"
	
	  << "load-plink|input|load a PLINK binary PED file (BED)|ARG:file,id,iid,fid,check-reference,fix-strand" 
      
	  << "load-dosage|input|load dosage data|ARG:file,file-list,id,check-reference,hard-call-threshold,format$space-delimited$skip-header$position-map$allele-map$dose1$dose2$prob2$prob3$as-dosage$as-posteriors,name,maf,mac"

	  << "attach-meta|input|load meta-information for existing VARDB variants|ARG:file,id,group"
	
	  << "load-pheno|input,indop|load phenotypes into INDB|ARG:file"
	
	  << "load-pedigree|input,indop|load pedigree information into INDDB|ARG:file"

	  << "swap-ids|input,indop|swap individual IDs in the VARDB and INDDB|ARG:file"

	  << "delete-meta|varop|remove meta-information|ARG:group,all" 
		
	  << "var-set|input,varop|add a variant-set to a VARDB|ARG:group,file"

	  << "var-superset|input,varop|add a super-set to a VARDB|ARG:group,file,members"

	  << "var-drop-set|varop|drop a set from a VARDB|AGR:group"

	  << "var-drop-all-sets|varop|drop all sets from a VARDB"

	  << "var-drop-superset|varop|drop a superset from a VARDB|AGR:group"

	  << "var-drop-all-supersets|varop|drop all supersets from a VARDB"


	//
	// VARDB output
	//
	
	
	  << "write-vardb|output,varop|write a new VARDB|ARG:new-vardb,new-project"
	
	  << "write-vcf|output|write a new VCF file|VCF|file,format$BGZF|OUT:vcf"

	  << "write-ped|output|write a new PLINK TPED fileset|VCF|ARG:name,use-family-id|OUT:tped,tfam" 
	
	  << "write-lik|output|write a BEALGE likelihood file|VCF|OUT:lik"

	  << "*write-haps|output|write a MaCH format haplotype file|VCF"
      
	  << "*write-bcf|output|output from VARDB to BCF|VCF|ARG:bcf|OUT:bcf"
      
      
      	
	//
	// Core view-functions for variant data
	//
	
	
	
	  << "v-view|views|view variant data|VCF|ARG:simple,vmeta,verbose,geno,gmeta,samples,hide-null,only-minor,only-alt,pheno|OUT:vars" 
	
	  << "rv-view|views|view rare alleles|VCF|ARG:pheno|OUT:vars" 
	
	  << "mv-view|views|view multiple variants|VCF|OUT:vars" 
    
	  << "mrv-view|views|view multiple rare variants|VCF|OUT:vars" 
    
	  << "g-view|views|view variants grouped by gene|GRP|ARG:vmeta,transpose,geno,gmeta,rarelist,phenotype,verbose|OUT:groups"
    
	  << "gs-view|views|view gene variants in sequence|GRP|ARG:ref,gene,protdb,domain,variant,plot"
	
	  << "i-view|views|individuals in project/file|VCF|ARG:phenotype,from-vardb|OUT:indiv"

	  << "i-matrix|views|individual phenotype matrix|ARG:phenotype|OUT:indiv.matrix"

	  << "write-phe|views|dump phenotypes in .phe file format|ARG:name|OUT:phe"

	  << "v-matrix|output|write a matrix of allele counts|VCF|OUT:matrix"

 //	  << "v-lookup|output|write a matrix of allele counts for a specific set of variants|ARG:file,region"
      
	  << "g-matrix|output|write a matix of gene-based allele counts|GRP|ARG:hide-invariant,collapse|OUT:matrix"
	
	  << "g-meta-matrix|output|matix of gene-based per-individual meta-information|GRP|ARG:name|OUT:matrix"
	
	  << "meta-matrix|output|write a matrix of variant meta-information|VCF|NOGENO|OUT:matrix"
	
	  << "v-meta-matrix|output|write a matrix of individual genotype meta-information|VCF|ARG:name|NOGENO|OUT:matrix"
	

	
	//
	// VARDB misc.
	//
	

	  << "tag-file|varop|add file-tags to VARDB|ARG:name,id"
      
	  << "var-delete|varop|remove file from VARDB|ARG:id" 
	
	  << "vacuum|varop|clean-up unused disk-space in VARDB"     

	  << "reindex|varop|drop then remake VARDB indices"
	
	
	//
	// Locus-operations
	// 
	
	  << "loc-load|input,locop|load from a .GTF or .REG file into LOCDB|ARG:file,group,keep-unmerged"
      
	  << "gtf-check|input,locop|check a .GTF (i.e. prior to loc-load'ing|ARG:file,add-frame"
      
	  << "locset-load|input,locop|load a locus-set|ARG:file,name,group,alternate-name"
	
	  << "loc-load-alias|input,locop|load a gene-alias table|ARG:file" 
	
	  << "loc-merge|locop|merge a LOCDB group|ARG:group"
	
	  << "loc-delete-alias|locop|remove gene-alias table"
	
	  << "loc-swap-names|locop|swap LOCDB names|ARG:file,group,alternate-name"

	  << "loc-update-name-table|locop|update LOCDB name table|ARG:group,alternate-name,index-name"

	  << "loc-delete|locop|remove a LOCDB group|ARG:group"
	
	  << "*loc-index|locop|index a LOCDB" 
	
	  << "*loc-drop-index|locop|remove index from LOCDB"         
	
	  << "loc-set-special|locop|set special variable in a LOCDB|ARG:key,value"
	
	  << "loc-get-special|locop|get special variable(s) from a LOCDB|ARG:key"
	
	  << "loc-delete-special|locop|remove all special variables from a LOCDB"
	

	

	  << "*intersect|views,locop|view loci from a LOCDB group with 1 or more variants|OUT:loci"
	
	  << "loc-view|views,locop|show all loci in a LOCDB group|ARG:group,alias,no-meta,show-subregions|OUT:loci"
	
	  << "loc-stats|views,locop|locus-based stats|ARG:no-subregions,show-subregions,group,ref|OUT:locstats"
	
	  << "loc-translate|locop|AA sequence of loci|OUT:loci"

	  << "loc-annotate|locop,annot|annotate loci|ARG:group,show-subregions:OUT:loci"
	
	  << "*loc-overlap|locop|show loci in groups Y, Z that overlap each locus in X|ARG:group,alias,comma,tab,row|OUT:loci"

      
	  << "seg-load|input,segop|input segment data to SEGDB" 

	  << "seg-view|views|individual segments|ARG:group|OUT:segs"

	  << "*seg-merge|segop|merge intervals in SEGDB" 
	  << "*segset-load|input,segop|load locus-sets in SEGDB" 
	  << "*seg-load-alias|input,segop|load alias-table" 
	  << "*seg-delete-alias|segop|remove alias-table"
	  << "*seg-delete|segop|remove segment group" 
	  << "*seg-index|segop|index SEGDB" 
	  << "*seg-drop-index|segop|index"            
    

    //
    // REFDB & SEQDB operations
    //
	
	
	  << "ref-load|input,refop|load data (VCF or flat-file) into REFDB|ARG:file,group,vcf"
      
	  << "*load-weights|input,refop|load weight table|ARG:name,file" 
      
	  << "*score-weights|annot|score variants for weights|VCF|ARG:name|OUT:scores"
      
	  << "seq-load|input,seqop|load FASTA into SEQDB|ARG:format$build$repeat-mode$iupac,file,name,description,"
      
	  << "lookup|misc,annot|lookup various annotatations for a list of positions|ARG:loc,alias,ref,ref-allelic,protdb,annotate,worstByAlias,worstAnnotationPriorities,titv|OUT:meta"
      
	  << "ref-view|views|view a group from a REFDB|ARG:group,vmeta,verbose|OUT:refvars"
      
	  << "seq-view|views|view regions of sequence from SEQDB|ARG:compact,region|OUT:seq"



    //
    // Database-level summaries
    //
	

	  << "summary|project|summary of all databases|ARG:ugly"

	  << "var-summary|project|summary of VARDB|ARG:ugly"

	  << "loc-summary|project|summary of LOCDB|ARG:ugly"

	  << "seg-summary|project|summary of SEQDB|ARG:ugly"

	  << "ind-summary|project|summary of INDDB|ARG:ugly"

	  << "ref-summary|project|summary of REFDB|ARG:ugly"

	  << "seq-summary|project|summary of SEQDB|ARG:ugly"

	  << "file-summary|project|summary of project files|ARG:ugly"

	  << "meta-summary|project|summary of variant meta-information|VCF|ARG:ugly"

    
    //
    // Data summary statistics
    //


	  << "v-stats|stats|variant statistics|VCF|ARG:stats$hwe$qual$dp$mac$maf$count$mean$gcount$gmean|OUT:vstats"
	
	  << "g-stats|stats|gene-based summary statistics|GRP|ARG:stats|OUT:gstats"
	
	  << "i-stats|stats|per-individual statistics|VCF|ARG:stats|OUT:istats"
	
	  << "i-pop|stats|individual posterior population probabilities|ARG:file|OUT:ipop"

	  << "v-dist|stats,tests|comparison of rare-variant group distributions|VCF|ARG:whole-sample-counts,perm|OUT:vdist"
	
	  << "v-freq|stats,qc|variant frequency data|VCF|ARG:em|OUT:vfreq"



      //
      // Genotype-phenotype association models
      //
      
      
	  << "counts|views,tests|summary/count statistics|VCF|ARG:output-vcf,name,annotate,full-annotate,show-filters,meta,phenotype|OUT:counts"
      
	  << "g-counts|views,tests|genotype summary/count statistics|VCF|ARG:output-vcf,name|OUT:gcounts"
      
	  << "means|views,tests|summary QT means per variant|VCF|ARG:annotate,full-annotate,show-filters,meta,phenotype|OUT:means"
      
	  << "assoc|tests|gene-based association tests|GRP|ARG:phenotype,tests$burden$calpha$uniq$vt$fw$sumstat$two-hit$skat$skato,info,fix-null,perm,dump-null-matrix,carriers,midpoint|OUT:assoc"
    
	  << "v-assoc|tests|single-variant association|VCF|ARG:phenotype,info,fix-null,perm,separate-chr-bp,vmeta|OUT:vassoc"
	
	  << "glm|tests|general linear models|VCF|ARG:phenotype,perm,vmeta,use-postprobs,use-dosages,covar,show-covar,show-intercept,residuals,file|OUT:glm"

	  << "unique|views,tests|view variants specific to individual groups|VCF|ARG:indiv,require,allow|OUT:uniq"

	  << "indiv-enrich|tests|test per individual for greater-than-expected burden of variants per set|ARG:phenotype,perm,locset,loc,loc.ex,reg.ex|OUT:indiv.enrich"


    //
    // Family-based operations
    //

	  << "denovo|views|filter for de-novo mutations (SNPs and indels)|VCF|ARG:minChildDP,minParDP,minChildPL,minParPL,minChild_AB_alt,minChild_AB_ref,minPar_AB_ref,minMQ,maxAAC,allowDoubleAltDeNovos,printTransmission|OUT:denovo.vars,denovo.indiv,parent_transmission.vars"
	  << "cnv-denovo|views|filter for de-novo CNV mutations|VCF|ARG:minSQ,minNQ,allowDeNovoWithoutDiscoveryInChild|OUT:denovo.cnv,denovo.indiv"

    //
    // IBD database 
    //


	  << "*ibd-load|input,ibd|load IBD segment data|ARG:ibddb,file"
	
	  << "*ibd-sharing|views,ibd|pairwise IBD sharing around rare variants|VCF|ARG:ibddb|OUT:ibd"
	
	  << "*s-assoc|tests,ibd|segment-based IBD test|ARG:perm,file|OUT:sassoc"

	  << "*mutation-screen|views,ibd|screen for new mutations given shared IBD|ARG:ibddb,indiv,region|OUT:mut"
      
      
    
    //
    // Gene-network database
    //
      
	  << "net-load|input,net|populate a NETDB|ARG:netdb,file"
      
	  << "net-view|views,net|view gene connections in a NETDB|GRP|ARG:name,group,netdb|OUT:nets"
      
	  << "net-assoc|net,tests|network-based gene-association|GRP|ARG:netdb,mask,pheno,file,outfile|OUT:net.assoc"
            
      

   //
   // Protein domain/feature annotation
   //

	  << "prot-load|input,prot|populate a PROTDB|ARG:protdb,file,group"
      
	  << "prot-view|views,prot|view entries from a PROTDB|ARG:protdb,name,group"

	  << "prot-summary|views,prot|summary of a PROTDB|ARG:protdb"
    //
    // Misc. QC etc
    //


	  << "clusters|qc|.|OUT:clusters"
	
	  << "proximity-scan|qc|VCF|ARG:distance|OUT:proximity"
	
	  << "concordance|qc|genotypic concordance checks|ARG:report-all|OUT:concord.indiv,concord.geno,concord.vars,concord"
	
	  << "group-comparison|qc,tests|OUT:group.comparison"
	
	  << "ibs-matrix|qc|IBS matrix calculation|VCF|ARG:long-format,two-counts|OUT:ibs"
	
	  << "loc-intersect|locop|intersect locus groups|ARG:file,group|OUT:loci"





    
    //
    // Misc.
    //
    
	  << "*simple-sim|misc|simple gene variant simulation|GRP";


}



std::string Pseq::Util::Options::load( int n , char ** argv )
{

    reg( "help" , NONE , "produce help message" );

    // output options

    reg( "out"    , STRING , "output root filename" );
    reg( "silent" , NONE   , "set silent output mode");
    reg( "quiet"  , NONE   , "set quiet output mode");
    reg( "noweb"  , NONE   , "skip web-check" );

    reg( "debug", NONE , "set debug mode");    


    // input sources

    reg( "vcf"       , STRING_VECTOR , "VCF file locations" );
    reg( "bcf"       , STRING_VECTOR , "BCF file locations" );
    reg( "resources" , STRING        , "central resource folder" );
    reg( "scratch"   , STRING        , "scratch folder" );
    reg( "metameta"  , STRING        , "meta-information meta-information" );

    reg( "description" , STRING        , "file description" );
    //    reg( "history"     , STRING_VECTOR , "use a .history file with GSEQ" );

    // core databases
    reg( "vardb" , STRING, "variant database location" );
    reg( "inddb" , STRING, "individual database location" );
    reg( "refdb" , STRING, "reference database location" );
    reg( "seqdb" , STRING, "sequence database location" );
    reg( "segdb" , STRING, "segent database location" );
    reg( "locdb" , STRING, "locus database location" );
    reg( "netdb" , STRING, "network database location" );  
    reg( "protdb", STRING, "protein domain/feature database" );
    reg( "ibddb" , STRING, "IBD segment database location" );  

    // misc file inputs
    reg( "file"      , STRING_VECTOR , "generic input file(s)" );
    reg( "file-list" , STRING        , "file to specify a list of files" );
    
    // dosage file inputs
    reg( "map-file"            , STRING , "map file" );  // primarily for load-dosage
    reg( "indiv-file"          , STRING , "individual ID list"); // primarily for load-dosage
    reg( "meta-file"           , STRING , "meta-information file" ); // primarily for load-dosage
    reg( "hard-call-threshold" , FLOAT , "hard-call threshold" ); // primarily for load-dosage

    reg( "variant"    , NONE          , "only show positions with variant sites in gs-view" );
    reg( "members"    , STRING_VECTOR , "super-set members" );
    reg( "ref-group"  , STRING , "REFDB group label" );
    reg( "loc-group"  , STRING , "LOCDB group label" );
    reg( "region"     , STRING_VECTOR , "region(s) ");
    reg( "alias"      , STRING_VECTOR , "locus alias group(s)" );
    reg( "name"       , STRING_VECTOR , "generic name(s) variable" );
    reg( "key"        , STRING , "key of key-value pair" );
    reg( "value"      , STRING_VECTOR , "value(s) of key-value pair" );
    reg( "type"       , STRING , "type of project entry");
    reg( "declare"    , STRING_VECTOR , "on-the-fly declaration of meta-attribute types" );

    reg( "group"      , STRING_VECTOR , "generic group label(s)" ); 
    reg( "all"        , NONE          , "generic group label(s)" );
    reg( "domain"     , STRING_VECTOR , "protein domain group from PROTDB" );
    reg( "plot"       , NONE          , "make R plotting script from gs-view" ); 
    reg( "id"         , INT_VECTOR    , "generic numeric IDs" );

    reg( "whitespace", NONE , "allow whitespace delimited input" );
    reg( "show-id" , NONE , "use chr:position:id variant format in output");

    // write-vardb
    reg( "new-project" , STRING , "new project specification filename" );
    reg( "new-vardb" , STRING , "new VARDB name, for write-vardb" );
    

    reg( "ignore-warnings" , NONE , "turn off warnings");
    reg( "early-warnings" , NONE , "display warning as soon as happens");
    reg( "limit-warnings" , INT , "print up to N warnings per topic" );
    reg( "all-warnings" , NONE , "show all warnings" );

    reg( "out-file", STRING , "set main output file");
    reg( "debug-file", STRING , "debug file name");
    reg( "prolix-file", STRING, "prolix output filename");
    reg( "long" , NONE , "set long output mode");
    reg( "long-header" , NONE , "set header/long output mode");

    // net-DB
    reg( "exclude-seed" , NONE , "exclude seed from net-assoc" ) ;

    // locus operations
    reg( "keep-unmerged" , NONE , "retain unmerged region group when loading GTF" );
    reg( "add-frame"     , NONE , "add frame to a GTF during gtf-check" );

    reg( "mask", STRING_VECTOR, "mask specification");
    
    reg( "include", STRING, "filter specification");
    reg( "eval", STRING, "eval expression");
    reg( "exclude", STRING , "filter excludes");

    reg( "v-include", STRING, "variant filter specification");
    reg( "v-eval", STRING, "variant eval expression");
    reg( "v-exclude", STRING , "variant filter excludes");

    reg( "gene", STRING_VECTOR, "gene-group gene1 gene2 ...");
    reg( "em", FLOAT , "EM calculation of P(geno|data) from GL or PL");
    
    reg( "assume-ref" , NONE , "convert null genotypes to reference homozygote");
    
    reg( "hide", STRING_VECTOR,"hide specific meta-fields");
    reg( "show", STRING_VECTOR,"show specific meta-fields" );
    reg( "force-consensus" , NONE , "set all tags to consensus" );
    
    reg( "vmeta", NONE, "show variant meta-information" );
    reg( "append-meta" , STRING_VECTOR , "on-the-fly appending of meta-information");
    reg( "samples", NONE, "show each specific sample variant");
    reg( "verbose", NONE, "verbose output");
    reg( "geno", NONE, "show genotypes");
    reg( "gmeta", NONE, "show genotype meta-information" );
    reg( "transpose", NONE, "transposed g-view output");
    reg( "simple" , NONE , "simple variant format, POS RET ALT" );
    reg( "rarelist" , NONE , "rare-list mode for v-view" );
    reg( "hide-null" , NONE , "" );
    reg( "only-minor" , NONE , "" );
    reg( "only-alt" , NONE , "" );

    reg( "variant", STRING , "show specific variant (v-view)");

    reg( "indiv", STRING_VECTOR , "specify individual(s)");
    reg( "require", INT , "require N individuals with variant" );
    reg( "allow", INT , "allow N individuals without variant" );
    
    reg( "annotate" , STRING , "transcript group for annotation" );
    reg( "titv" , NONE , "annotate locus as being transition, transversion, or neither" );
    reg( "loc" , STRING_VECTOR , "transcript group" );
    reg( "ref" , STRING_VECTOR , "reference-variant group" );
    reg( "ref-allelic" , STRING_VECTOR , "reference-variant group" );
    reg( "var-allelic" , NONE , "allele matching for var.req (ONLY)" );
    reg( "locset" , STRING_VECTOR , "locus-set group" );
    reg( "worstByAlias"      , STRING_VECTOR , "locus alias group(s) for which to split out worst annotations" );
    reg( "worstAnnotationPriorities" , STRING, "File containing an ordering of the 'worst' annotations chosen to be output by lookup" );
    
    // Genotype/phenotype inputs

    reg( "phenotype" , STRING, "phenotype specification");
    reg( "phenotype-file" , STRING_VECTOR, "phenotype specification from file");
    reg( "make-phenotype" , STRING, "dichotomise factor");

    reg( "from-vardb" , NONE , "list individuals from VARDB, not INDDB" );

    reg( "strata" , STRING,"stratifier variable");
    reg( "covar" , STRING_VECTOR , "covariate(s)");
    reg( "residuals" , STRING_VECTOR , "residualize phenotype by these covariates" );
    reg( "weights" , STRING , "name of variant weights tag");
    reg( "skat-weights" , FLOAT_VECTOR , "Beta(a,b); a=1,b=25 default");

    // De novo scan parameters:
    reg("minChildDP", INT, "Minimum child depth");
    reg("minParDP", INT, "Minimum parental depth");
    reg("minChildPL", FLOAT, "Minimum child PL (genotype likelihood) for non-called genotype");
    reg("minParPL", FLOAT, "Minimum parental PL (genotype likelihood) for non-homozygous reference genotype");
    reg("minChild_AB_alt", FLOAT, "Minimum child AB-ALT (% of reads with ALT allele)");
    reg("minChild_AB_ref", FLOAT, "Minimum child AB-REF (% of reads with REF allele)");
    reg("minPar_AB_ref", FLOAT, "Minimum parental AB-REF (% of reads with REF allele)");
    reg("minMQ", FLOAT, "Minimum MQ (read mapping quality)");
    reg("maxAAC", INT, "Maximum AAC (alternate allele count at this site)");
    reg("allowDoubleAltDeNovos", NONE, "Include de novos that consist of two alternate alleles in child");
    reg("printTransmission", NONE, "Parent variants will be printed with transmission status to the parent_transmission.vars file");

    reg("minSQ", INT, "Minimum SQ score for calling a CNV (e.g., in child)");
    reg("minNQ", INT, "Minimum NQ score for ruling out a CNV (e.g., in parent)");
    reg("allowDeNovoWithoutDiscoveryInChild", NONE, "Try to call de novo CNV in child even in cases where that CNV is not marked as having been discovered in that child");

    // Input modifiers 

    // Output modifiers
    
    reg( "output" , KEYWORD , "general output format modifiers" );
    keyword( "output" , "row" , NONE , "new row per item" );
    keyword( "output" , "tab" , NONE , "tab-delimited" );
    keyword( "output" , "comma" , NONE , "comma-delimited" );

    // generic FORMAT for 

    reg( "format" , KEYWORD , "" );
    keyword( "format" , "name" , STRING , "" );
    keyword( "format" , "build" , STRING , "" );
    keyword( "format" , "repeat-mode" , STRING , "" );
    keyword( "format" , "description" , STRING , "" ); // SEQDB, REFDB
    keyword( "format" , "iupac" , NONE , "" );
    keyword( "format" , "BGZF" , NONE , "write VCF in BGZF-compressed form" );

    keyword( "format" , "position-map" , NONE , "dosage map only contains chr/bp" );
    keyword( "format" , "allele-map" , NONE , "dosage map only contains a1/a2" );    
    keyword( "format" , "dose1" , NONE , "dosage data 0..1" );
    keyword( "format" , "dose2" , NONE , "dosage data 0..2" );
    keyword( "format" , "prob2" , NONE , "2 posterior probabilities" );
    keyword( "format" , "prob3" , NONE , "3 posterior probabilities" );
    keyword( "format" , "as-dosage" , NONE , "store as dosage" );
    keyword( "format" , "as-posteriors" , NONE , "store as 3 posterior probailities" );
    keyword( "format" , "skip-header" , NONE , "ignore header in dosage file(s)" );
    keyword( "format" , "space-delimited" , NONE , "use spaces, not tabs, as delimiters in dosag file(s)" );
    keyword( "format" , "no-ids" , NONE , "dosage file does not contain any ID column" );

    keyword( "format" , "chr" , STRING , "" ) ;
    keyword( "format" , "bp1" , STRING , "" );
    keyword( "format" , "bp2" , STRING , "" );
    keyword( "format" , "id" , STRING , "" );
    keyword( "format" , "skip" , STRING_VECTOR , "");
    
    keyword( "format" , "integer" ,STRING_VECTOR , "");
    keyword( "format" , "string" ,STRING_VECTOR , "");
    keyword( "format" , "text" ,STRING_VECTOR , "");
    keyword( "format" , "float" ,STRING_VECTOR , "");
    keyword( "format" , "flag" ,STRING_VECTOR , "");
    keyword( "format" , "bool" ,STRING_VECTOR , "");
    
    keyword( "format" , "0-based" , NONE , "genomic interval start/stop are 0-based" );
    keyword( "format" , "0-start" , NONE , "genomic interval start is 0-based" );
    keyword( "format" , "0-stop"  , NONE , "genomic interval stop is 0-based" );
    
    keyword( "format" , "include-meta" , STRING_VECTOR , "include only these tags when loading" );
    keyword( "format" , "exclude-meta" , STRING_VECTOR , "exclude these tags when loading" );


    reg( "check-reference" , NONE , "loading/indexing VCFs, report if REF != SEQDB" );    
    reg( "fix-strand" , NONE , "fix strand given a SEQDB when loading a PLINK GWAS file" ); 
    
    // Output modifiers

    reg( "family-id" , NONE , "use {FID,IID} instead of ID" );

    reg( "fid" , NONE , "use FID as ID when reading PLINK files" );
    reg( "iid" , NONE , "use IID as ID when reading PLINK files" );

    // dosage loading modifers
    reg( "maf" , FLOAT , "minimum alternate allele frequency to load a dosage entry" );
    reg( "mac" , INT , "minimum minor allele count to load a dosage entry" );

    // Association models
    
    reg( "tests"    , KEYWORD , "gene-based association tests to apply" );
    reg( "midpoint" , NONE , "report interval BP mid-point in output" );
    reg( "info"     , NONE , "report ISTAT in association test" );
    reg( "separate-chr-bp" , NONE , "report CHR-tab-BP not chr:BP" );
    reg( "fix-null" , NONE , "exclude individuals with null genotypes" );
    reg( "yates" , NONE , "use Yates' test for single-site association" );
    
    reg( "dump-null-matrix" , NONE , "output (null) gene-based statistics in matrix form" );
    reg( "carriers" , NONE , "output who carriers of ALT alleles are in .assoc.det" );
    reg( "perm" , INT, "number of permutations");
    reg( "aperm" , INT_VECTOR , "adaptive perm min, max");
    reg( "seed" , FLOAT , "seed for RNG" ); 

    reg( "use-dosages" , STRING , "genotypes are dosages in specified tag" );
    reg( "use-postprobs" , STRING , "genotype are posterior probabilities in specified tag" );

    reg( "show-covar" , NONE, "list coefficients for covariates" );
    reg( "show-intercept" , NONE, "list coefficients for b0 term in GLM");
    
    reg( "tests" , KEYWORD , "gene-based tests" );
    keyword( "tests" , "sumstat" , NONE , "sum-statistic test" );
    keyword( "tests" , "burden" , NONE , "basic burden test" );
    keyword( "tests" , "uniq" , NONE , "burden of case-unique variants" );
    keyword( "tests" , "site-burden" , NONE , "burden of variant sites (not alleles)" );
    keyword( "tests" , "mhit" , NONE , "multiple-hit association model" );
    keyword( "tests" , "vt" , NONE , "variable threshold test" );
    keyword( "tests" , "fw" , NONE , "frequency-weighted test" );
    keyword( "tests" , "calpha" , NONE , "c-alpha test" );
    keyword( "tests" , "cancor" , NONE , "canonical correlation test" );
    keyword( "tests" , "stepup" , NONE , "Hoffman-Witte step-up test" );
    keyword( "tests" , "kbac" , NONE , "KBAC test" );
    keyword( "tests" , "two-hit" , NONE, "Recessive/compound-het test" );
    keyword( "tests" , "skat" , NONE , "SKAT test" );
    keyword( "tests" , "skato" , NONE , "Optimal SKAT test" );
    
    // loading intervals/GTF 
    
    reg( "use-gene-id" , NONE , "loading GTFs" );
    reg( "index-name" , NONE , "use LOCDB primary/index name (e.g. transcript ID)" );
    reg( "alternate-name" , NONE , "use LOCDB alternate name (e.g. gene symbol)" );


    // v-stats command
    reg( "stats" , KEYWORD , "quantities calculated under (v|g|i)-stats" );
    keyword( "stats" , "hwe" , FLOAT_RANGE_VECTOR , "HWE test p-values" );
    keyword( "stats" , "ref" , STRING_VECTOR , "REFDB groups" );
    keyword( "stats" , "loc" , STRING_VECTOR , "LOCDB groups" );
    keyword( "stats" , "mac" , INT_RANGE_VECTOR , "MAC range counts" );
    keyword( "stats" , "maf" , FLOAT_RANGE_VECTOR , "MAF range counts" );
    keyword( "stats" , "qual" , FLOAT_RANGE_VECTOR , "QUAL range counts" );
    keyword( "stats" , "dp" , INT_RANGE_VECTOR , "variant read-depth counts" );
    keyword( "stats" , "mean" , STRING_VECTOR , "arbitrary mean tag values" );
    keyword( "stats" , "count" , STRING_VECTOR , "arbirary tag range counts" );
    keyword( "stats" , "filter" , STRING_VECTOR , "FILTER counts");
    keyword( "stats" , "gmean" , STRING_VECTOR , "genotype tag means" );
    keyword( "stats" , "gcount" , STRING_VECTOR , "genotype tag range counts" );

    // modifies i-stats , etc 
    reg( "alternate" , NONE , "stats for alternate, not minor, alleles" );


    // summaries
    reg( "ugly" , NONE , "do not format summary output" );
    
    
    reg( "show-subregions" , NONE , "modifies loc-view output" ); // loc-annotate
    reg( "no-subregions" , NONE , "modifies loc-view output" ); // loc-stats
    reg( "no-meta" , NONE , "do not show locus meta-information"); // loc-view, loc-
    reg( "subregions" , NONE , "show subregions from LOCDB"); //loc-view
    
    reg( "show-phase" , NONE , "" );
    reg( "acdb" , NONE , "write summary VCF in ACDB-format" );
    
    // counts and g-counts

    reg( "annotate" , STRING , "annotate variants" );
    reg( "full-annotate" , STRING , "full-annotate variants" );
    reg( "show-filters" , NONE , "show filters" );
    reg( "output-vcf" , NONE , "summary counts in VCF format" );

     reg( "filter" , STRING , "loading VCF, only load variants in specified LOCDB group" );

     

     // Misc. commands / one-off options
     
     reg( "hide-invariant" , NONE , "" );  // g-matrix
     reg( "collapse" , NONE , "" ); // g-matrix
     reg( "whole-sample-counts" , NONE , "" );  // v-dist
     reg( "report-all" , NONE , "" ); // concordance
     reg( "distance" , NONE , "" ); // ?     
     reg( "compact" , NONE , "compact seq-view output" );  // seq-view
     reg( "prev" , STRING , "disease prevalence" ); // setting prevalence for recessive/compound het tests
     reg( "func-inc" , STRING_VECTOR , "included functional annotations" ); // get list of functional annotations to include
     reg( "func-exc" , STRING_VECTOR , "excluded functional annotations" ); // get list of functional annotations to include
     reg( "singles", NONE, "include singleton-singleton compound hets" ); // if true, include compound het variants that are two singletons in the same individual
     reg( "mhit", NONE, "include multi hits on the same chromosome" ); // if true, don't perform phasing and just count multi hits



     Pseq::IBS::regargs( this );  // IBS tests (example
    

    
 
    
    


    //
    // Register some short-cuts
    //
    
    shortform( "-m" , "--mask" );
    shortform( "-f" , "--file" );
    shortform( "-g" , "--group" );
    shortform( "-p" , "--phenotype" );
    shortform( "-n" , "--name" );    
    shortform( "-o" , "--out" );

    shortform( "-h" , "--help" );  
    shortform( "help" , "--help" );

    shortform( "--string" , "--text" );
    shortform( "--pheno"  , "--phenotype" );
    shortform( "--phe"    , "--phenotype" );
    shortform( "--weight" , "--weights" );
    
    
    needs_help = n > 1 && ( std::string( argv[1] ) == "--help" 
			    || std::string( argv[1] ) == "help" 
			    || std::string( argv[1] ) == "-h" 
			    || std::string( argv[1] ) == "masks" 
			    || std::string( argv[1] ) == "mask" ) ;
    
    
    if ( needs_help ) // will be called later
    {      
	if ( n == 2 ) help_str = std::string(argv[1]).substr(0,1) == "m" ? "mask" : "." ;
	else if ( n > 2 ) help_str = argv[2];
	return "";
    }
    

    //
    //  Process actual arguments from command line
    //
    
    // create a friendly version of it for LOG
    
    std::stringstream ss;

    // arg 0 is filename
    // position 1 should be project
    // position 2 should be command
    // then commands and flags.
    
    if ( n < 2  ) Helper::halt( "no project or command specified; try 'pseq help'" );
    else if ( n < 3 ) Helper::halt( "no command specified; try 'pseq help'" );
    
    project_str = argv[1];
    command_str = argv[2];

    ss << "Project : " << project_str << "\n"
       << "Command : " << command_str << "\n";
    if ( n > 3 ) 
      ss << "Options : " ;
    
    for (int i=3 ; i < n ; i++ )
    {

      std::string s = argv[i];
	
      // Swap in long-form from short-form? e.g. -o to --options
      if ( shortcuts.find( s ) != shortcuts.end() ) s = shortcuts[s];
	
      if ( s.substr(0,2) != "--" ) Helper::halt("unknown option: " + s ); 
      s =  s.substr(2); 
      if ( ! known(s) ) Helper::halt("unknown option: " + s );
      
      // pretty-printing
      if ( i > 3 ) ss << "          ";
      ss << "--" << s ;
      const int len_first = s.size() + 2 ;
      const int pos_root = i;

      // now read values/keywords for this command, up to next arg (--x) or end
      
      std::vector<std::string> a;
      
      while ( 1 ) 
	{
	  ++i;
	  
	  if ( i == n )
	    {
	      if ( i == pos_root + 1 ) 
		ss << "\n";
	      break;
	    }

	  std::string b = argv[i];
	  
	  if ( shortcuts.find( b ) != shortcuts.end() ) b = shortcuts[b];
	  
	  if ( b.substr(0,2) == "--" ) 
	  {
	      --i; 
	      break;
	  }
	  else
	  {
	    
	    // this is an argument for 's'.  Either just add, unless 's' is of KEYWORD type, 
	    // in which case we need to check that 'b' is valid
	    
	    if ( i - pos_root > 1 ) 
	      ss << "          " << std::string( len_first , ' ' );
	    ss << " " << b << "\n";
	    
	    // --arg value value=x,y,y-2,"a,b"
	    
	    if ( arg_type[s] == KEYWORD )
	      {
		
		std::string root = b.substr( 0 , b.find( "=" ) );
		
		if ( known( s , root ) ) 
		  {
		    // comma-delimited
		    
		    
		    if ( b.find("=") == std::string::npos )
		      {
			// just set s/root pair
			data_kw[s][root].clear();
		      }
		    else
		      {
			// or add actual --s root=arg1,arg2
			
			std::string val = b.substr( b.find("=")+1 );			  
			std::vector<std::string> vals = Helper::quoted_char_split( val , ',' );
			for (int i=0;i<vals.size();i++)
			  {
			    data_kw[s][root].push_back( vals[i] );
			  }
		      }
		  }
		else Helper::halt( "unknown keyword " + root + " for argument " + s );
	      }
	    else a.push_back( b );
	    
	  }
	  
	}
      
      if ( i == pos_root ) // no options given
	ss << "\n";
      
      // store, or add in as extras
      
      if ( data.find(s) == data.end() ) data[s] = a;      
      else Helper::append( data[s] , a );      
      
    }

    return ss.str();
}


bool Pseq::Util::Options::known( const std::string & s ) const
{
    return arg_type.find(s) != arg_type.end();
}
 

bool Pseq::Util::Options::known( const std::string & s , const std::string & k ) const
{
    std::map<std::string,std::set<std::string> >::const_iterator i = keyword_map.find(s);
    if ( i == keyword_map.end() ) return false;
    return i->second.find( k ) != i->second.end();
}

void Pseq::Util::Options::reg( const std::string & s , const type_t & t , const std::string & desc )
{
  arg_type[s] = t;
  arg_desc[s] = desc;
}

void Pseq::Util::Options::keyword( const std::string & s , const std::string & k , const type_t & t , const std::string & desc )
{
    keyword_map[s].insert( k );
    keyword_type[s][k] = t;    
    keyword_desc[s][k] = desc;
}

bool Pseq::Util::Options::has( const std::string & s ) const
{
    if ( arg_type.find(s) == arg_type.end() ) Helper::halt("internal error --" + s + " not registered " );
    return data.find(s) != data.end();
}

bool Pseq::Util::Options::has( const std::string & s , const std::string & k ) const
{
    if ( arg_type.find(s) == arg_type.end() ) Helper::halt("internal error -- " + s + " not registered " );
    std::map<std::string, std::map<std::string,std::vector<std::string> > >::const_iterator i = data_kw.find(s);
    if ( i == data_kw.end() ) return false;
    return i->second.find( k ) != i->second.end();
}


std::vector<std::string> Pseq::Util::Options::as_string_vector( const std::string & a )  const
{  
  if ( ! known(a) ) Helper::halt("argument " + a + " not found" );
  return data.find(a)->second;
}

std::string Pseq::Util::Options::as_string( const std::string & a )  const 
{
  if ( ! known(a) ) Helper::halt("argument --" + a + " not found" );
  if ( data.find(a)->second.size() != 1 ) Helper::halt("expecting 1 value for --" + a );
  return data.find(a)->second[0];
}

int Pseq::Util::Options::as_int( const std::string & a )  const
{
  if ( ! known(a) ) Helper::halt("argument " + a + " not found" );
  if ( arg_type.find(a)->second != INT ) Helper::halt("incorrect return type for " + a );
  int i;
  if ( ! Helper::str2int( data.find(a)->second[0] , i ) ) 
    Helper::halt("non-integer argument for --" + a );
  return i;
}

std::vector<int> Pseq::Util::Options::as_int_vector( const std::string & a ) const
{
  if ( ! known(a) ) Helper::halt("argument " + a + " not found" );
  if ( arg_type.find(a)->second != INT_VECTOR ) Helper::halt("incorrect return type for " + a );
  std::vector<int> d;
  for (int i=0; i<data.find(a)->second.size(); i++)
    {
      int x;
      if ( ! Helper::str2int( data.find(a)->second[i] , x ) ) 
	Helper::halt("non-integer argument for --" + a );
      d.push_back(x);
    }
  return d;
}

double Pseq::Util::Options::as_float( const std::string & a ) const
{
  if ( ! known(a) ) Helper::halt("argument " + a + " not found" );
  if ( arg_type.find(a)->second != FLOAT ) Helper::halt("incorrect return type for " + a );
  double d;
  if ( ! Helper::str2dbl( data.find(a)->second[0] , d ) ) 
    Helper::halt("non-integer argument for --" + a );
  return d;
}

std::vector<double> Pseq::Util::Options::as_float_vector( const std::string & a ) const
{
  if ( ! known(a) ) Helper::halt("argument " + a + " not found" );
  if ( arg_type.find(a)->second != FLOAT_VECTOR ) Helper::halt("incorrect return type for " + a );
  std::vector<double> d;
  for (int i=0; i<data.find(a)->second.size(); i++)
    {
      double x;
      if ( ! Helper::str2dbl( data.find(a)->second[i] , x ) ) 
	Helper::halt("non-numeric argument for --" + a );
      d.push_back(x);
    }
  return d;
}



std::vector<std::string> Pseq::Util::Options::as_string_vector( const std::string & a , const std::string & b ) const
{
    std::vector<std::string> r;
    if ( ! known(a,b) ) Helper::halt("argument/keyword --" + a + " " + b + " not found" );
    
    if( data_kw.find(a) == data_kw.end() ) return r;

    std::map<std::string,std::map<std::string,std::vector<std::string> > >::const_iterator i = data_kw.find( a );
    if ( i == data_kw.end() ) return r;
    
    std::map<std::string,std::vector<std::string> >::const_iterator ii = i->second.find( b );
    if ( ii == i->second.end() ) return r;
    
    return ii->second;

}

std::string Pseq::Util::Options::as_string( const std::string & a , const std::string & b ) const
{
    std::string r;
    if ( ! known(a,b) ) Helper::halt("argument/keyword --" + a + " " + b + " not found" );
    
    if( data_kw.find(a) == data_kw.end() ) return r;

    std::map<std::string,std::map<std::string,std::vector<std::string> > >::const_iterator i = data_kw.find( a );
    if ( i == data_kw.end() ) return r;
    
    std::map<std::string,std::vector<std::string> >::const_iterator ii = i->second.find( b );
    if ( ii == i->second.end() ) return r;
    
    if ( ii->second.size() == 1 ) return ii->second[0];
    
    for (int j = 0 ; j < ii->second.size() ; j++ ) 
    {
	if ( j ) r +=",";
	r += ii->second[j];
    }
    
    return r;
}


bool Pseq::Util::Options::has( const std::string & a , const std::string & b , const std::string & k ) const
{
    std::string r;

    if ( ! known(a,b) ) Helper::halt("argument/keyword --" + a + " " + b + " not found" );
    
    if( data_kw.find(a) == data_kw.end() ) return false;

    std::map<std::string,std::map<std::string,std::vector<std::string> > >::const_iterator i = data_kw.find( a );
    if ( i == data_kw.end() ) return false;
    
    std::map<std::string,std::vector<std::string> >::const_iterator ii = i->second.find( b );
    if ( ii == i->second.end() ) return false;
    
    for (int j = 0 ; j < ii->second.size() ; j++ ) 
      if ( ii->second[j] == k ) return true;
    
    return false;
}



int Pseq::Util::Options::as_int( const std::string & a , const std::string & b ) const
{
    int r;
    Helper::halt( "as_int not implemented yet" );
    return r;
}

std::vector<int> Pseq::Util::Options::as_int_vector( const std::string & a , const std::string & b ) const
{
    std::vector<int> r;
    Helper::halt( "as_int_vector not implemented yet" );
    return r;

}

double Pseq::Util::Options::as_float( const std::string & a , const std::string & b ) const
{
    double r;
    Helper::halt( "as_float not implemented yet" );
    return r;
}



bool Pseq::Util::Options::help() const
{
  if ( ! needs_help ) return false;
  std::cout << desc( help_str ) ;
  return true;
}


std::string Pseq::Util::Options::desc( const std::string & c ) const 
{

  std::stringstream ss;
  
  // --help         { list all commands } 
  // --help v-view  { list all arguments for v-view }
  // --help masks   { list all mask options }
  
  if ( c == "." )
    {

      ss << "\nusage:\tpseq {project-file|VCF|-|.} {command} {--options}\n\n";
      
      std::vector<std::string> groups = pcomm.groups( );
      
      ss << "\tCommand groups\n";
      ss << "\t---------------------------------------------------------\n";
      for ( int g = 0 ; g < groups.size() ; g++ )
	{	  
	  std::vector<std::string> cs = pcomm.commands( groups[g] );
	  if ( cs.size() ) // do not show hidden/empty ones
	    ss << "\t" << groups[g] << "\t\t\t" << pcomm.group_description( groups[g] ) << "\n";
	}      
      
      ss << "\n\tMask groups\n";
      ss << "\t---------------------------------------------------------\n";
      ss << Mask::list_groups();
      
      ss << "\n"
	 << "pseq help all\n"
	 << "pseq help {group}\n"
	 << "pseq help {command}\n";

   }
  else if ( c == "all" )
    {
      // verbose help mode
      ss << "\nusage:\tpseq {project-file|VCF} {command} {--options}\n\n";
      
      ss << "Commands\n";
      ss << "---------------------------------------------------------\n\n";
      
      std::vector<std::string> cs = pcomm.all_commands();

      for (int c = 0 ; c < cs.size() ; c++)
	{	  
	  if ( pcomm.stable( cs[c] ) )
	    {
	      ss << "\t" << pcomm.command_description( cs[c] );
	      ss << pcomm.command_description( cs[c] , true ) << "\n";
	    }
	}
      
      
      ss << "\n\nMask groups\n";
      ss << "---------------------------------------------------------\n\n";

      ss << Mask::list_groups( true );
      
      ss << "\n"
	 << "pseq help {group}\n"
	 << "pseq help {command}\n";      

    }
  else if ( pcomm.has_group( c ) )
    {      
      ss << "\n\t" << c << " : " << pcomm.group_description( c ) << "\n";
      ss << "\t---------------------------------------------------------\n";
      std::vector<std::string> cs = pcomm.commands( c );      
      for (int c = 0 ; c < cs.size() ; c++)
	if ( pcomm.stable( cs[c]) ) ss << "\t" << pcomm.command_description( cs[c] );
      
    }
  else if ( pcomm.known( c ) )
    {
      std::vector<std::string> str = Helper::char_split( pcomm.command_description( c ) , '\t' , true );
      ss << "\n\t" << c << " : " << str[str.size()-1] ;
      ss << "\t---------------------------------------------------------\n";
      ss << pcomm.command_description( c , true );
    }
  else 
    {
      // is this a Mask group?
      
      std::string str = Mask::list_masks( c );
      if ( str == "" ) 	
	ss << "[ " << c << " ] not recognised as a command, command-group or mask-group\n";
      else
	{
	  // v. silly way to get desc... quicker for now.
	  std::string s = Mask::list_groups();	  
	  std::vector<std::string> s1 = Helper::char_split( s , '\n' );
	  std::string desc;
	  for ( int i = 0 ; i < s1.size() ; i++)
	    {
	      std::vector<std::string> s2 = Helper::char_split( s1[i] , '\t' , false );
	      if ( s2.size() == 2 && s2[0] == c ) { desc = s2[1] ; break; }
	    }
	  ss << "\n\t" << c << " : " << desc << "\n";
	  ss << "\t---------------------------------------------------------\n";
	  ss << str;
	}
    }
  ss << "\n";
  return ss.str();
}


void Pseq::Util::Options::attach( const std::string & command , const std::string & arg )
{
    // arg will be comma-delimited list. for keyword args, the acceptable list will 
    // be enclosed in [...]
    
    // ARG:arg1,arg|key1|key2|key3,arg3
    
    std::vector<std::string> opts = Helper::char_split( arg , ',' );
    for (int i=0; i<opts.size(); i++)
    {
	std::vector<std::string> x = Helper::char_split( opts[i] , '$' );
	comm2arg[ command ].insert( x[0] );
	for ( int j = 1 ; j < x.size() ; j++ ) 
	    comm2key[ command ][ x[0] ].insert( x[j] );	    
    }

}


std::string Pseq::Util::Options::attached( const std::string & command )
{
  std::map<std::string,std::set<std::string> >::iterator i = comm2arg.find( command );  
  if ( i == comm2arg.end() ) return "";

  std::string s = "";
  std::set<std::string>::iterator j = i->second.begin();
  while ( j != i->second.end() )
    {
	s += "ARG\t" + i->first + "\t" + *j + "\t" + type( *j ) + "\t" + arg_description( *j ) + "\n";
	
	std::map<std::string,std::map<std::string,std::set<std::string> > >::iterator ii = comm2key.find( command );  
	if ( ii != comm2key.end() ) 
	{
	    std::map<std::string,std::set<std::string> >::iterator iii = ii->second.find( *j );
	    if ( iii != ii->second.end() )
	    {
		std::set<std::string>::iterator si = iii->second.begin();
		while ( si != iii->second.end() )
		{
		    s += "KEY\t" + i->first + "\t" + *si + "\t" + key_type( *j , *si ) + "\t" + key_description( *j , *si ) + "\n";
		    ++si;
		}
	    }
	}	
      ++j;
    } 

  return s;

}




std::string Pseq::Util::Options::arg_description( const std::string & arg ) const
{
    std::map<std::string,std::string>::const_iterator i = arg_desc.find( arg );
    if ( i == arg_desc.end() ) return "";
    return i->second;
}

std::string Pseq::Util::Options::type( const std::string & arg ) const
{
    std::map<std::string,type_t>::const_iterator i = arg_type.find( arg );
    if ( i == arg_type.end() ) return ".";
    if ( i->second == NONE ) return "flag";
    if ( i->second == STRING ) return "str";
    if ( i->second == INT ) return "int";
    if ( i->second == KEYWORD ) return "keyword";
    if ( i->second == FLOAT ) return "float";
    if ( i->second == STRING_VECTOR ) return "str-list";
    if ( i->second == INT_VECTOR ) return "int-list";
    if ( i->second == FLOAT_VECTOR ) return "float-list";
    if ( i->second == VARIABLE_TYPE ) return "variable";
    return ".";
}


std::string Pseq::Util::Options::key_type( const std::string & a , const std::string & b ) const
{
    std::map<std::string,std::map<std::string,type_t> >::const_iterator i = keyword_type.find( a );
    if ( i == keyword_type.end() ) return ".";
    std::map<std::string,type_t>::const_iterator ii = i->second.find( b );
    if ( ii == i->second.end() ) return ".";
    if ( ii->second == NONE ) return "flag";
    if ( ii->second == STRING ) return "str";
    if ( ii->second == INT ) return "int";
    if ( ii->second == INT_RANGE ) return "int-range";
    if ( ii->second == INT_RANGE_VECTOR ) return "int-range-vector";
    if ( ii->second == KEYWORD ) return "keyword";
    if ( ii->second == FLOAT ) return "float";
    if ( ii->second == FLOAT_RANGE ) return "float-range";
    if ( ii->second == FLOAT_RANGE_VECTOR ) return "float-range-vector";

    if ( ii->second == STRING_VECTOR ) return "str-list";
    if ( ii->second == INT_VECTOR ) return "int-list";
    if ( ii->second == FLOAT_VECTOR ) return "float-list";
    if ( ii->second == VARIABLE_TYPE ) return "variable";
    return ".";
}

   
std::string Pseq::Util::Options::key_description( const std::string & a , const std::string & b ) const
{
    std::map<std::string,std::map<std::string,std::string> >::const_iterator i = keyword_desc.find( a );
    if ( i == keyword_desc.end() ) return ".";
    std::map<std::string,std::string>::const_iterator ii = i->second.find( b );
    return ii == i->second.end() ? "." : ii->second;
}

std::set<std::string> Pseq::Util::Options::get_set( const std::string & k ) const
{ 
    std::set<std::string> s; 
    std::map< std::string, std::vector< std::string > >::const_iterator i = data.find(k); 
    if ( i == data.end() ) return s; 
    for (int j=0; j < i->second.size(); j++)
	s.insert( i->second[j] );
    return s;
} 
	    
std::set<std::string> Pseq::Util::Options::get_set( const std::string & a , const std::string & b ) const
{
    std::set<std::string> s; 
    std::map<std::string, std::map<std::string,std::vector<std::string> > >::const_iterator i = data_kw.find(a); 
    if ( i == data_kw.end() ) return s; 
    
    std::map<std::string,std::vector< std::string > >::const_iterator ii = i->second.find( b );
    if ( ii == i->second.end() ) return s; 
    
    for (int j=0; j < ii->second.size(); j++)
	s.insert( ii->second[j] );
    return s;
}


std::string Pseq::Util::Options::comma_string( const std::string & k ) const 
{     
    std::map< std::string,std::vector<std::string> >::const_iterator i = data.find(k); 
    if ( i == data.end() ) return "";
    std::string s = "";		
    for (int j=0; j < i->second.size(); j++)
    {
	if ( j ) s += ",";
	s += i->second[j];
    }
    return s;
} 

std::string Pseq::Util::Commands::command_description( const std::string & c , bool show_args  ) const
{
  
  std::stringstream ss;
  
  if ( known( c ) ) 
    {
      if ( ! show_args ) 
	{
 	  ss << c << "\t";
	  if ( c.size() < 8 ) ss << "\t";
	  if ( c.size() < 16 ) ss << "\t";
	  ss << comm_desc.find( c )->second << "\n";
	}

      if ( show_args )
	{
	  std::string sarg = pargs->attached( c );
	  std::vector<std::string> sarg1 = Helper::char_split( sarg , '\n' );
	  
	  for (int i = 0 ; i < sarg1.size() ; i++ ) 
	    {
	      std::vector<std::string> v = Helper::char_split( sarg1[i] , '\t' );		    
	      if ( v.size() == 5 )
		{
		    if ( v[0] == "ARG" ) 
		    {
			ss << "\t  --" << v[2];
			if ( v[3] != "." ) ss << " { " << v[3] << " }" ;
			ss << "\t" << v[4] << "\n";
		    }
		    else if ( v[0] == "KEY" )
		    {
			ss << "\t\t";

			std::stringstream s1;
			s1 << v[2];
			if ( v[3] != "." ) s1 << " { " << v[3] << " }" ;
			std::stringstream s2;
			s2 << v[4];
			ss << std::left << std::setw(30) << s1.str() << std::setw(40) << s2.str() << "\n";
		    }
		}
	    }	  
	}      
    }
  return ss.str();
}
