
#include <iostream>
#include <iterator>
#include <limits>

#include "plinkseq.h"

#include "func.h"
#include "views.h"
#include "summaries.h"
#include "assoc.h"
#include "compare.h"
#include "ibd.h"
#include "ibs.h"
#include "extra.h"
#include "cnv.h"
#include "pops.h"

using namespace std;

GStore g;

Pseq::Util::Options args;
Pseq::Util::Commands pcomm;

std::string PSEQ_VERSION = "0.09";
std::string PSEQ_DATE    = "20-Sep-12";


int main(int argc, char ** argv)
{
  
  //
  // Get command-line options into a sensible form
  //

  // pseq {project} {command} {--options}
  
  std::string pp_args = args.load( argc , argv );


  // connect with GSEQ job tracking?
  
  // if ( args.has( "history" ) )
  //   {
  //     std::vector<std::string> tok = args.as_string_vector( "history" );
  //     if ( tok.size() != 2 ) Helper::halt( "--history {file} {job#}" );
  //     g.gseq_tracking( tok[0] , tok[1] );
  //   }



  // known PSEQ commands, and descriptions
  
  pcomm.attach( &args );

  Pseq::Util::populate_commands( pcomm );


  //
  // help message?
  //
  
  if ( args.help() ) exit(0);
  

  //
  // Misc. initialisation
  //

  if ( args.has( "seed" ) ) 
    {
      long unsigned s = static_cast<long unsigned>(args.as_float( "seed" ));      
      CRandom::srand( s );
      srand( s );
      PLINKSeq::DEFAULT_RNG_SEED() = s;      
    }

  //
  // Reporting and logging options
  //
  

  if ( args.has("ignore-warnings") )
    plog.show_warnings( false );

  if ( args.has("early-warnings") )
    plog.early_warnings( true );
  

  if ( args.has( "all-warnings" ) )
    plog.set_warning_limit( numeric_limits<int>::max() );
  else if ( args.has( "limit-warnings" ) )
    plog.set_warning_limit( args.as_int( "limit-warnings" ) );
  
             
  
  
  //
  // Default : send all output to STDOUT, warnings to STDERR, no LOG
  //
  

  if ( args.has( "silent" ) ) plog.silent( true );
  
  else if ( args.has( "quiet" ) ) plog.quiet( true );

  // Default to write to STDOUT for all output...

  Out::set_stdout( true );
  Out::set_tofile( false );

  // ... unless a file is explicitly specified, --out {file}

  std::string fileroot = args.has( "out" ) ? args.as_string("out") : "";
  
  if ( args.has( "out" ) ) 
    {      
      // inform main output class
      Out::set_stdout( false );
      Out::set_tofile( true );
      Out::set_fileroot( fileroot );

      // and LOG class
      plog.set_fileroot( fileroot );

    }
  else
    plog.close_logfile(); // ensure no output to file

  if ( args.has( "debug" ) ) 
    debug.silent( false );

  if ( args.has( "debug-file" ) )
    debug.logfile( args.as_string( "debug-file" ) ); 
  


  // ----- likely remove support for 'long-mode' from Out class
  // if ( args.has("long") )
  //   Out::longmode();
  // if ( args.has("long-header") )
  //   { Out::longmode(); Out::header(); } 


  //
  // Opening banner
  //
  
  time_t curr=time(0);
  std::string tdstamp = (std::string)ctime(&curr);

  plog >> "-------------------------------------------------------------------------------\n"    
       >> "||||||||||||||||||||||||||| PSEQ (v" 
       >> PSEQ_VERSION >> "; " >> PSEQ_DATE 
       >> ") |||||||||||||||||||||||||||\n"
       >> "-------------------------------------------------------------------------------\n\n";

  //
  // Process 'command'
  //
  
  std::string command = "";
  
      
  //
  // Always, we require a single command, and a single project to be
  // specified (by first and second arguments)
  //
  
  command = args.command(); 
  
  if ( ! pcomm.known( command ) ) 
    Helper::halt("command " + command + " not recognised" );


  //
  // Basics to the log
  //

  if ( args.has( "out" ) && ! ( args.has( "silent" ) || args.has( "quiet" ) ) )
    std::cout << "Copying this log to file [ " << fileroot << ".log ]\n";
  
  // Time-stamp 

  plog >> "Analysis started " >> tdstamp;

  //
  // Web-based version check
  //

  if ( ! args.has( "noweb" ) ) 
    Pseq::Util::webcheck( command );

  
  //
  // show actual input command lines in log
  //

  plog >> "\n"
       >> "-------------------------------------------------------------------------------\n\n"
       >> pp_args;
   
  plog >> "\n-------------------------------------------------------------------------------\n\n";

  
  //
  // Set project, attaching all relevant databases
  //
  
  std::string project_file = FileMap::tilde_expansion( args.project_file() );
  

  //
  // If a single VCF has been specified as the 'project'
  //
  
  g.single_file_mode( project_file == "-" 
		      || Helper::ends_with( project_file , ".bcf" ) 
		      || Helper::ends_with( project_file , ".bcf.gz" ) 
		      || Helper::ends_with( project_file , ".vcf"  )
		      || Helper::ends_with( project_file , ".vcf.gz" ) );
  
  if ( g.single_file_mode() )
    {

      if ( ! pcomm.single_VCF_mode( command ) )
	Helper::halt( command + " not applicable in single-VCF mode" );
  
      if ( Helper::ends_with( project_file , ".bcf" ) 
	   || Helper::ends_with( project_file , ".bcf.gz" ) )
	   g.single_file_bcf( true );      
    }


  //
  // Functions that do not depend on any databases
  //

  
  //
  // Query command table, or Mask options
  //

  if ( command == "commands" ) 
    {
      std::string n = args.as_string( "name" );      
      if ( n == "" || n == "." ) n = "root";
      pcomm.display( n );
      Pseq::finished();
    }
  
  if ( command == "masks" )
    {
      plog << Mask::describe_options();
      Pseq::finished();
    }


  //
  // Create a new project file?
  //
  
  if ( command == "new-project" ) 
    {
      Pseq::new_project( project_file , args );
      Pseq::finished();
    }
  
  
  //
  // Version information
  //
  
  if ( command == "version" ) 
    {
      Out pout( "version" , "version information" );
      pout << "PSEQ: " << PSEQ_VERSION << "\n"
	   << "PSEQ DATE: " << PSEQ_DATE << "\n";
	   
      g.show_version();
      
#ifdef ZLIB_VERNUM
      pout << "ZLIB: " << ZLIB_VERSION << "\n";
#endif

      Pseq::finished();
  
    }


  //
  // Simple routine to create dummy genic VariantGroups 
  // and feed to the genic association tests (designed 
  // for testing more than anything else)
  //

  if ( command == "simple-sim" ) 
    {
      Out output( "assoc.sim" , "output from simple-sim command" );
      Pseq::VarDB::simple_sim();
      Pseq::finished();
    }
  

  //
  // Otherwise, open existing project
  //
  
  if (  g.single_file_mode() ) 
  {
    if ( project_file != "-" && ! Helper::fileExists( project_file ) ) 
      Helper::halt( "could not open file " + project_file );      
  }
  else
    {
      if ( ! Pseq::set_project( project_file ) )
	Helper::halt("could not open project " + project_file );      
    }
  
  
  //
  // Hot-swap in alternate DBs, core folders, etc
  //
  
  if ( args.has("resources") )
    {
      std::string folder = args.as_string( "resources" );
      if ( folder.substr( folder.size()-1,1 ) != "/" )
	folder += "/";
      g.fIndex.addSpecial( RESOURCES , folder );
      g.fIndex.make_dir( g.fIndex.file( RESOURCES )->name() );
    }
  

  if ( args.has("vardb") && ! g.single_file_mode() )
    g.vardb.attach( args.as_string( "vardb" ) );

  // In single VCF mode, we still need a VARDB, but it can be in
  // memory

  if ( g.single_file_mode() ) 
    {
      g.vardb.attach( ":memory:" );
      g.inddb.attach( ":memory:" );
    }

  if ( args.has("inddb") )
    g.inddb.attach( args.as_string( "inddb" ) );
  
  if ( args.has("locdb") )
    g.locdb.attach( args.as_string( "locdb" ) );
  
  if ( args.has("segdb") )
    g.segdb.attach( args.as_string( "segdb" ) );
  
  if ( args.has("refdb") )
    g.refdb.attach( args.as_string( "refdb" ) );
  
  if ( args.has("seqdb") )
    g.seqdb.attach( args.as_string( "seqdb" ) );

  

  //
  // Add/remove files from the project file
  //

  if ( command == "drop" )
    {
      std::string s;
      if ( args.has("file") )
	s = args.as_string( "file" );
      else if ( args.has("name") )
	s = args.as_string( "name" );
      else if ( args.has( "type" ) )
	s = args.as_string( "type" );
      else
	Helper::halt("no --file, --name or --type specified");      
      
      g.fIndex.remove_from_projectfile( s );
      Pseq::finished();
    }

  if ( command == "append" ) 
    {
      std::string pname;
      if ( args.has("name") ) 
	pname = args.as_string( "name" );
      else if ( args.has("file") )
	pname = args.as_string( "file" );
      else
	Helper::halt("no --name or --file specified");
      
      if ( args.has("name") && args.has("file" ) ) 
	Helper::halt("use either --name OR --file");

      if ( ! args.has( "type" ) ) 
	Helper::halt("no --type specified");
      std::string type = args.as_string( "type" );
      
      g.fIndex.append_to_projectfile( Helper::fullpath( pname ) , type );
      Pseq::finished();
    }
  

  //
  // On-the-fly declaration of types
  //
  
  if ( args.has( "declare" ) )
    {
      // --declare DB|INFO|0|Flag|"Description a desc" DB2|FORMAT|1|GT|Float|"A value"
      
      std::vector<std::string> t = args.as_string_vector( "declare" );
      for ( int tt = 0 ; tt < t.size() ; tt++ )
	{
	  std::vector<std::string> x = Helper::char_split( t[tt] , '|' );

	  if ( x.size() != 5 ) Helper::halt("expecting ID|GROUP|N|TYPE|DESC for --declare" );

	  const std::string & name = x[0];
	  const std::string & type = x[3];
	  const std::string & desc = x[4];

	  int n = 0;
	  if ( ! Helper::str2int( x[2] , n ) ) n = -1; // variable length
	  
	  mType mt = META_UNDEFINED;	  
	  if ( Helper::is_int( type ) )  mt = META_INT;
	  else if ( Helper::is_float( type ) ) mt = META_FLOAT;
	  else if ( Helper::is_text( type ) ) mt = META_TEXT;  // text == String || Char
	  else if ( Helper::is_flag( type ) ) mt = META_FLAG;
	  
	  // Does this contain valid information?
	  
	  if ( mt == META_UNDEFINED )
	    Helper::halt( "problem defining type, " + type );

	  if ( n < -1 ) 	  
	    Helper::halt( "problem defining number, " + x[2] );
	  
	  if ( x[1] == "INFO" ) 
	    MetaInformation<VarMeta>::field( name , mt , n , desc );
	  else if ( x[1] == "FORMAT" ) 
	    MetaInformation<GenMeta>::field( name , mt , n , desc );
	  else if ( x[1] == "FILTER" ) 
	    MetaInformation<VarFilterMeta>::field( name , mt , n , desc );
	  else
	    Helper::halt( x[1] + " not one of INFO, FORMAT or FILTER" );
	}
      
    }
  
  
  //
  // Misc. formatting/display options
  //

  
  if ( args.has("whitespace") )
    {
      PLINKSeq::DELIM() = " \t\n";
    }

  if ( args.has( "show-id" ) )
    {
      g.show_id( true );
    }

  if ( args.has("hide"))
    {
      std::vector<std::string> s = args.as_string_vector( "hide" );
      for (int i=0; i<s.size(); i++) MetaMeta::hide( s[i] );
    }
  
  if ( args.has("show") )
    {
      std::vector<std::string> s = args.as_string_vector( "show" );
      for (int i=0; i<s.size(); i++) MetaMeta::show( s[i] );
    }
  
  if ( args.has( "force-consensus" ) )
    {
      MetaMeta::set_force_consensus( true );
    }


  //
  // Hint to use soft-called genotypes, if command allows
  //

  if ( args.has( "use-dosages" ) )
    {
      Genotype::using_dosage = Genotype::using_soft_calls = true;
      Genotype::soft_call_label = args.as_string( "use-dosages" );
    }
  else if  ( args.has( "use-postprobs" ) )
    {
      Genotype::using_probs = Genotype::using_soft_calls = true;
      Genotype::soft_call_label = args.as_string( "use-postprobs" );
    }


  //
  // Load reference and sequence data
  //
  
  if ( command == "ref-load" )
    {
      
      if ( ! g.refdb.attached() ) 
	Helper::halt("no refdb attached");
      
      if ( args.has("vcf") ) 
	{
	  Pseq::RefDB::load_VCF( args.as_string( "vcf" ), 
				 args.as_string( "group" ) ); 
	}
      else if ( args.has("file" ) )
	{      
	  Pseq::RefDB::load_refvar( args.as_string( "file" ) , 
				    args.as_string( "group" ) ,
				    args );
	}
      else 
	Helper::halt("no file or VCF specified");
      
      Pseq::finished();	
    }
  

  if ( command == "seq-load" )
    {	

      if ( ! args.has( "file" ) )
	Helper::halt("no file specified");

      std::vector<std::string> s = args.as_string_vector( "file" );

      if ( s.size() != 1 )
	Helper::halt("more than 1 file specified");
      
      Pseq::SeqDB::load_FASTA( s[0] );

      Pseq::finished();
    }
  
  

  //
  // PPH2 scoring
  //
  
  if ( command == "load-weights" )
    {
      Helper::halt("obsolete command score weights" );

//       std::string dbname = Pseq::Util::single_argument<std::string>( args , "name" );
//       std::string filename = Pseq::Util::single_argument<std::string>( args , "file" );
//       Pseq::PPH2DB::load( dbname , filename );
//       Pseq::finished();
    }
  


  //
  // Functions that do not depend on the mask, variant or individual database
  //
  
  if ( command == "loc-intersect" ) 
    {

      if ( ! args.has("file") )
	Helper::halt("no file specified");

      std::vector<std::string> s = args.as_string_vector( "file" );
      if ( s.size() != 1 )
	Helper::halt("more than 1 file specified");
      
      if ( ! args.has( "group" ) )
	Helper::halt("no group specified");

      std::vector<std::string> grp = args.as_string_vector( "group" );
      if ( grp.size() != 1 )
	  Helper::halt("more than 1 group specified");
      
      LocDBase * db = g.resolve_locgroup( grp[0] );
      if ( ! db ) Helper::halt("group not found");

      Out output( "loci" , "intersecting locus list" );
      
      Pseq::LocDB::intersection( s[0] , grp[0] , *db );
      
      Pseq::finished();
    }
  

  
  if ( command == "locset-load" )
    {
      if ( ! args.has("name") )
	Helper::halt("no name specified");	
      
      if ( ! args.has("file") )
	Helper::halt("no file specified");
      std::vector<std::string> s = args.as_string_vector( "file" );
      if ( s.size() != 1 )
	Helper::halt("more than 1 file specified");
      
      if ( ! args.has("group") )
	Helper::halt("no group specified");
      std::vector<std::string> grp = args.as_string_vector( "group" );
      if ( grp.size() != 1 )
	Helper::halt("more than 1 group specified");
      
      bool use_altname = args.has( "alternate-name" );
      
      if ( ! Pseq::LocDB::load_pathway( s[0],
					Pseq::Util::single_argument<std::string>( args, "name" ) ,
					grp[0] ,
					true , 
					use_altname ) )
	Helper::halt("problem loading locset");
      
      Pseq::finished();
    }
  

  if ( command == "net-load" )
    {
      if ( ! args.has( "netdb" ) ) Helper::halt( "no --netdb specified" );
      if ( ! args.has( "file" ) ) Helper::halt( "no --file specified" );
      Pseq::NetDB::loader( args.as_string( "netdb" ) , args.as_string( "file" ) );
      Pseq::finished();
    }

  if ( command == "net-view" )
    {
      if ( ! args.has( "name" ) ) Helper::halt( "no --name specified" );
      if ( ! args.has( "netdb" ) ) Helper::halt( "no --netdb specified" );
      if ( ! args.has( "group" ) ) Helper::halt( "no --group specified" );

      Out output( "partners" , "network partners of a given gene/element" );

      Pseq::NetDB::lookup( args.as_string( "netdb" ) , 
			   args.as_string( "name" ) , 
			   args.as_string( "group" ) );

      Pseq::finished();
    }


  //
  // PROTDB functions
  //

  if ( command == "prot-load" )
    {
      if ( ! args.has( "protdb" ) ) Helper::halt( "no --protdb specified" );
      if ( ! args.has( "file" ) ) Helper::halt( "no --file specified" );
      if ( ! args.has( "group" ) ) Helper::halt( "no --group specified" );
      Pseq::ProtDB::loader( args.as_string( "protdb" ) , args.as_string( "file" ) , args.as_string( "group" ) );
      Pseq::finished();
    }
  

  
  //
  // SEGDB functions
  //

  if ( command == "segset-load" )
    {

      if ( ! args.has("name") )
	Helper::halt("no name specified");	
      
      if ( ! args.has("file") )
	Helper::halt("no file specified");
      std::vector<std::string> s = args.as_string_vector( "file" );
      if ( s.size() != 1 )
	Helper::halt("more than 1 file specified");
      
      if ( ! args.has("group") )
	Helper::halt("no group specified");
      std::vector<std::string> grp = args.as_string_vector( "group" );
      if ( grp.size() != 1 )
	Helper::halt("more than 1 group specified");
      
      if ( ! Pseq::LocDB::load_pathway( s[0],
					Pseq::Util::single_argument<std::string>( args, "name" ) ,
					grp[0] ,
					false , 
					args.has("alternate-name") ) )
	Helper::halt("problem loading segset");
      
      Pseq::finished();
    }
  
  
  
  //
  // Compile mask specification
  //
  
  std::string maskspec = "";
  std::string filtspec = "";
  bool filter_T_include = true;
  
  if ( args.has( "mask" ) )
    {	
      std::vector<std::string> t = args.as_string_vector( "mask" );
      for (int i=0; i<t.size(); i++)
	{	    
	  if ( t[i].length() > 8 && t[i].substr(0,8)=="include=" ) 
	    filtspec = t[i].substr(8);
	  else if ( t[i].length() > 8 && t[i].substr(0,8)=="exclude=" ) 
	    {
	      filtspec = t[i].substr(8);
	      filter_T_include = false;
	    }
	  else
	    {
	      // to allow white-space in Mask options, put quotes around 
	      // all entries following an var=value --> var="value"
	      
	      if ( maskspec != "" ) maskspec += " ";
	      maskspec += Helper::quote_value( t[i] );		
	    }
	}
    }
  
  
  if ( args.has("include") )
    {
      filtspec = args.as_string( "include" );
    }
  else if ( args.has("eval") )
    {
      filtspec = args.as_string( "eval" ) + " ; T " ;      
    }  
  else if ( args.has("exclude") )
    {
      filtspec = args.as_string( "exclude" );
      filter_T_include = false;
    }

  if ( args.has("gene") )
    {
      
      std::vector<std::string> t = Pseq::Util::n_arguments<std::string>( args , "gene" );	
      
      if ( t.size() > 0 )
	{
	  
	  if ( pcomm.groups( command ) )  
	    maskspec += " loc.group=" + PLINKSeq::DEFAULT_LOC_GROUP();
	  else
	    maskspec += " loc=" + PLINKSeq::DEFAULT_LOC_GROUP();
	  
	  maskspec += " loc.subset=" + PLINKSeq::DEFAULT_LOC_GROUP() ;
	  
	  for (int i=0; i<t.size(); i++) 
	    {
	      // attempt to lookup up alias (symbol->refseq)
	      
	      std::set<std::string> trans_names = 
		g.locdb.targetted_lookup_alias( t[i] , 
						"symbol" , 
						PLINKSeq::DEFAULT_LOC_GROUP() ) ;
	      
	      if ( trans_names.size() == 0 ) 
		Helper::halt( "could not find " + t[i] + " in LOCDB" );

	      std::vector<std::string> tnames;
	      std::set<std::string>::iterator ii = trans_names.begin();
	      while ( ii != trans_names.end() )
		{
		  maskspec += "," + *ii;
		  ++ii;
		}		
	    }
	}   	
    }
  
    
  

  //
  // Create index and store location of project BCFs, VCFs
  //
  
  if ( command == "index-bcf" ) 
  {
    
    if ( ! args.has( "bcf" ) ) 
      Helper::halt( "no BCF files specified, use --bcf file(s)" );
    
    std::vector<std::string> t = args.as_string_vector( "bcf" );
    
    g.vardb.drop_index();
    
    for (int f=0; f<t.size(); f++)
      {
	
	if ( ! Helper::fileExists( t[f] ) )
	  {
	    plog.warn( "could not find BCF" , t[f] );
	    continue;
	  }
	
	// ensure we are using the full path
	t[f] = Helper::fullpath( t[f] );
	
	// Add to project index
	g.fIndex.append_to_projectfile( Helper::fullpath( t[f] ) , "BCF" );
	
	// Add to file-map, and create a BCF instance
	BCF * bcf = g.fIndex.add_BCF( t[f] );
	
	// Open BCF via BGZF interface	    
	bcf->reading();
	bcf->open();
	
	// Iterate through file, adding index	    	    
	g.vardb.begin();
	
	// Get header information, and add to VARDB
	bcf->read_header( &g.vardb );
	
	uint64_t inserted = 0;
	while ( bcf->index_record() )
	  {
	    if ( ++inserted % 1000 == 0 )
	      plog.counter1( "parsed " + Helper::int2str( inserted ) + " rows" );
	  }
	plog.counter1("\n");
	
	plog << "inserted " << inserted << " variants from BCF; now finishing index...\n";
	
	g.vardb.commit();
	bcf->close();
	
	// and calculate summary Ns
	int2 niv = g.vardb.make_summary( t[f] );
	
      }
    
    g.vardb.index();
      
    Pseq::finished();
  }
  
  

  
  //
  // Support PLINK binary files
  //
  
  if ( command == "load-plink" )
    {
      if ( ! args.has("file") )
	Helper::halt("no --file specified");
      
      std::string tag = "";
      if ( ! args.has("id") ) 
	tag = Pseq::Util::single_argument<std::string>( args , "id" );	
      
      Pseq::VarDB::load_PLINK( args.as_string_vector( "file" ) , args , tag );
      
      Pseq::finished();
    }
  
  
  //
  // Load .dosage file (-like) format
  //
  
  if ( command == "load-dosage" ) 
    {
      Pseq::VarDB::load_dosage();
      Pseq::finished();
    }


  //
  // Load/flush meta-information into VARDB, as independent meta-information 
  //

  // TODO  
//   if ( command == "load-meta" )
//     {
//       // either put a new meta-information, or from a file
//       if ( args.has( "file" ) ) Pseq::insert_meta_from_file( args.has( "file" ) ) ;
//       else if ( args.has( "name" ) ) Pseq::insert_meta_on_fly( args.has( "name" ) );
//       else Helper::halt( "need to specify --file or --name" );
//       Pseq::finished();
//     }
  

  if ( command == "attach-meta" )
    {
      std::string file = Pseq::Util::single_argument<std::string>( args , "file" );
      std::vector<std::string> id = Pseq::Util::n_arguments<std::string>( args , "id" );
      std::string group = Pseq::Util::single_argument<std::string>( args , "group" );
      for (int i=0;i<id.size();i++)
	g.vardb.loader_indep_meta( file , g.vardb.file_tag( id[i] ) , group );
      Pseq::finished();
    }

 
  
  if ( command == "delete-meta" )
    {
      if ( args.has( "all" ) )
	{	  
	  g.vardb.flush_indep_meta( );
	}
      else if ( args.has( "group" ) ) 
	{
	  std::string group = Pseq::Util::single_argument<std::string>( args , "group" );
	  g.vardb.flush_indep_meta( group );
	}
      else
	Helper::halt( "must specify --group {name}   OR --all " );
      
      Pseq::finished();
    }

  
  //
  // Add a descriptive 'tag' (i.e. short name) to a file in VARDB
  //
  
  if ( command == "tag-file" )
    {
      std::string newtag = Pseq::Util::single_argument<std::string>( args , "name" );
      std::string oldtag = Pseq::Util::single_argument<std::string>( args , "id" );
      int file_id = g.vardb.file_tag( oldtag );
      if ( file_id ) g.vardb.insert_file_tag( file_id , newtag );
      Pseq::finished();
    }
  

  //
  // Flush 1 or more files from VARDB, based on ID numbers
  //
  
  if ( command == "var-delete" )
    {
      if ( ! args.has("id") )
	Helper::halt("no --id specified");	      
      std::vector<std::string> f = args.as_string_vector( "id" );
      std::vector<int> fi;
      for (int i=0; i<f.size(); i++)
	{
	  int n = g.vardb.file_tag( f[i] );
	  if ( n ) fi.push_back( n );
	}
      Pseq::VarDB::flush( fi );	
      Pseq::finished();
    }
  
  
  if ( command == "vacuum" )
    {
      Pseq::VarDB::vacuum();
      Pseq::finished();
    }


  if ( command == "reindex" ) 
    {
      plog << "reindexing VARDB... this may take a while...\n";
      g.vardb.drop_index();
      g.vardb.index();
      Pseq::finished();
    }
  

  //
  // Check format and transcripts in a GTF
  //

  if ( command == "gtf-check" )
    {
      if ( ! args.has( "file" ) ) Helper::halt( "no --file specified" );
      g.locdb.check_GTF( args.as_string( "file" ) , false , false );
      Pseq::finished();
    }

  
  //
  // Load/merge GTF files into locus database
  //
  
    if ( command == "loc-load" )
      {

	if ( ! args.has("file") )
	  Helper::halt("no file specified");
	std::vector<std::string> s = args.as_string_vector( "file" );
	if ( s.size() != 1 )
	  Helper::halt("more than 1 file specified");

	if ( ! args.has("group") )
	  Helper::halt("no group specified");
	std::vector<std::string> grp = args.as_string_vector( "group" );
	if ( grp.size() != 1 )
	  Helper::halt("more than 1 group specified");

	bool remove_unmerged = ! args.has( "keep-unmerged" );

	if ( Helper::ends_with( s[0] , ".gtf" ) || Helper::ends_with( s[0] , ".gtf.gz" ) )
	  {
	      // Load GTF, merge transcripts, remove tmp group
	      const std::string realgrp = grp[0];
	      const std::string tmpgrp = grp[0] + "-unmerged";
	      Pseq::LocDB::load_GTF( s[0], tmpgrp , true );
	      Pseq::LocDB::merge( tmpgrp , realgrp , true );
	      if ( remove_unmerged ) g.locdb.flush( tmpgrp );
	      Pseq::LocDB::update_searchtable( realgrp  );
	  }
	else if ( Helper::ends_with( s[0] , ".reg" ) || Helper::ends_with( s[0] , ".reg.gz" ) )
	  {
	    Pseq::LocDB::load_generic_regions( s[0], grp[0] , args , true );
	    Pseq::LocDB::update_searchtable( grp[0]  );
	  }
	else Helper::halt("invalid file name, expecting extension: .gtf .gtf.gz .reg .reg.gz");
	
	Pseq::finished();
      }
    
    if ( command == "loc-swap-names" )
      {
	std::string filename = Pseq::Util::single_argument<std::string>( args , "file" );
	std::string group = Pseq::Util::single_argument<std::string>( args , "group" );
	Pseq::LocDB::swap_alternate_names( group , filename );
      }
    

    if ( command == "loc-update-name-table" )
      {
	if ( ! args.has("group") ) Helper::halt("no group specified");
	bool use_altname = ! args.has( "index-name" );
	Pseq::LocDB::update_searchtable( args.as_string( "group" ) , use_altname );
	Pseq::finished();
      }


    if ( command == "loc-merge" )
      {
	
	if ( ! args.has("group") )
	  Helper::halt("no group specified");
	std::vector<std::string> grp = args.as_string_vector( "group" );
	if ( grp.size() != 2 )
	  Helper::halt("need to specify two groups (old-name new-name)");

	Pseq::LocDB::merge( grp[0] , grp[1] , true );
	
	Pseq::finished();
      }


    //
    // Add locus alias table
    //

    if ( command == "loc-load-alias" )
      {	
	if ( !g.locdb.attached() ) Helper::halt("LOCDB not attached");
	if ( ! args.has("file") )
	  Helper::halt("no file specified");
	std::vector<std::string> f = args.as_string_vector( "file" );
	for (int i=0; i<f.size(); i++ )
	  g.locdb.load_alias( f[i] );
	Pseq::finished();
      }


    if ( command == "loc-delete-alias" )
      {	
	if ( ! g.locdb.attached() ) Helper::halt("LOCDB not attached");
	g.locdb.delete_aliases( );
	Pseq::finished();
      }


    //
    // Remove a group of loci from LOCDB
    //

    if ( command == "loc-delete" )
      {	
	if ( ! args.has("group") )
	  Helper::halt("no group specified");
	std::vector<std::string> grp = args.as_string_vector( "group" );
	for (int i=0; i<grp.size(); i++ )
	  g.locdb.flush( grp[i] );
	Pseq::finished();
      }

    //
    // Re/drop index LOCDB
    //

    if ( command == "loc-index" )
      {	
	g.locdb.drop_index();
	g.locdb.index();
	Pseq::finished();
      }

    if ( command == "loc-drop-index" )
      {	
	g.locdb.drop_index();
	Pseq::finished();
      }

    //
    // Set/get 'special' values from a LOCDB
    //

    if ( command == "loc-set-special" )
      {
	if ( ! args.has("key" ) )
	  Helper::halt("requires --key");
	if ( ! args.has("value" ) )
	  Helper::halt("requires --value");
	std::string key = args.as_string( "key" );

	//	Pseq::Util::n_arguments<std::string>( args , "key" );
	std::vector<std::string> value = args.as_string_vector( "value" );

	g.locdb.insert_special( key , value );
	Pseq::finished();
      }


    if ( command == "loc-get-special" )
      {
	if ( ! args.has("key" ) )
	  Helper::halt("requires --key");
	std::string key = args.as_string( "key" );
	std::vector<std::string> value = g.locdb.fetch_special( key );
	for (int i=0;i<value.size();i++)
	  plog << value[i] << "\n";
	Pseq::finished();
      }


    if ( command == "loc-delete-special" )
      {
	g.locdb.clear_special();
	Pseq::finished();
      }

    
    //
    // Remove a group of loci from SEGDB
    //

    if ( command == "seg-delete" )
      {	
	if ( ! args.has("group") )
	  Helper::halt("no group specified");
	std::vector<std::string> grp = args.as_string_vector( "group" );
	for (int i=0; i<grp.size(); i++ )
	  g.segdb.flush( grp[i] );	
	Pseq::finished();
      }

    //
    // Re/drop index SEGDB
    //

    if ( command == "seg-index" )
      {	
	g.segdb.drop_index();
	g.segdb.index();
	Pseq::finished();
      }

    if ( command == "seg-drop-index" )
      {	
	g.segdb.drop_index();
	Pseq::finished();
      }



    //
    // Load segment (locus plus individual) information into SEGDB
    // 

    if ( command == "seg-load" )
      {

	if ( ! args.has("file") )
	  Helper::halt("no file specified");
	std::vector<std::string> s = args.as_string_vector( "file" );
	if ( s.size() != 1 )
	  Helper::halt("more than 1 file specified");
	
	if ( ! args.has("group") )
	  Helper::halt("no group specified");
	std::vector<std::string> grp = args.as_string_vector( "group" );
	if ( grp.size() != 1 )
	  Helper::halt("more than 1 group specified");
	
	// Is this a GTF file, or a SEG file that contains
	// per-individual segment information? Check extension.
	
	if ( s[0].size() > 4 )
	  {	
	    std::string ext = s[0].substr( s[0].size() - 4 );
	    
// 	    if ( ext == ".gtf" )	
// 	      Pseq::LocDB::load_GTF( s[0], grp[0] , false );
	    
	    if ( ext == ".seg" )
	      Pseq::LocDB::load_segments( s[0], grp[0] , args );
	    else
	      Helper::halt( "expecting a .seg file" );
	  }
	
	Pseq::finished();
      }


    //
    // Add locus alias table to SEGDB
    //

    if ( command == "seg-load-alias" )
      {	
	if ( !g.segdb.attached() ) Helper::halt("LOCDB not attached");
	if ( ! args.has("file") )
	  Helper::halt("no file specified");
	std::vector<std::string> f = args.as_string_vector( "file" );
	for (int i=0; i<f.size(); i++ )
	  g.segdb.load_alias( f[i] );
	Pseq::finished();
      }

    if ( command == "seg-delete-alias" )
      {	
	if ( ! g.segdb.attached() ) Helper::halt("SEGDB not attached");
	g.segdb.delete_aliases( );
	Pseq::finished();
      }


    //
    // Merge in SEGDB
    //

    if ( command == "seg-merge" )
      {

	if ( ! args.has("group") )
	  Helper::halt("no group specified");
	std::vector<std::string> grp = args.as_string_vector( "group" );
	if ( grp.size() != 2 )
	  Helper::halt("need to specify two groups (old-name new-name)");
	
	Pseq::LocDB::merge( grp[0] , grp[1] , false );

	Pseq::finished();
      }



    //
    // Load individual phenotype or pedigree information into
    // individual database
    //
    
    if ( command == "load-pheno" )
      {
	if ( ! args.has("file") )
	  Helper::halt("no file specified");
	std::vector<std::string> s = args.as_string_vector( "file" );

	for ( int f=0; f<s.size(); f++)
	  {
	    Pseq::IndDB::load_phenotypes( s[f] );	
	  }
	Pseq::finished();
      }


    if ( command == "load-pedigree" )
      {
	if ( ! args.has("file") )
	  Helper::halt("no file specified");
	std::vector<std::string> s = args.as_string_vector( "file" );
	for ( int f=0; f<s.size(); f++)
	  Pseq::IndDB::load_ped_info( s[f] );
	Pseq::finished();
      }

    
    if ( command == "swap-ids" )
      {
	if ( !args.has( "file" ) )
	  Helper::halt( "no file specified" );
	Pseq::VarDB::swap_ids( args.as_string( "file" ) );	
	Pseq::finished();
      }


    
    //
    // EM calculation of P( genoytpe | data ) , from GL or PL genotype tags
    //

    if ( args.has( "em" ) )       
      maskspec += " em=" + args.as_string( "em" );
          
    
    //
    // Add indicator for single VCF-iteration mode
    //

    if ( g.single_file_mode() ) 
      {
	maskspec += " ex-vcf=" + project_file;
      }


    //
    // If a phenotype, or covariates, have been specified, by default require that we see these
    // This can be over-riden by adding phe.allow.missing to the --mask 
    //
    
    std::string phenotype_name = PLINKSeq::DEFAULT_PHENOTYPE();
    if ( ! g.single_file_mode() )
      {
	if ( args.has("phenotype") ) 
	  {
	    phenotype_name = args.as_string("phenotype");
	    if ( ! g.phmap.phenotype_exists( phenotype_name ) )
	      Helper::halt( "could not find " + phenotype_name );
	  }
	else if ( args.has("make-phenotype") )
	  {
	    std::vector<std::string> k = Helper::char_split( args.as_string("make-phenotype") , '=' );
	    if ( k.size() == 2 ) phenotype_name = k[0];
	  }
      }

    if ( phenotype_name != "" ) 
      {
	maskspec += " phe.obs=" + phenotype_name;
      }
	
    
    
    //
    // ---------------------------------------------------------------------------
    //

    // *** By this point, we must be pointing to the correct databases, etc *** 

    //
    // Process Mask(), as all other commands involve accessinng the
    // data, which can involve a mask, and possibly a phenotype
    //
    
    Mask m( maskspec , filtspec , filter_T_include , pcomm.groups( command ) );
    

    //
    // Set up individual-map and connect it to VARDB give a mask;
    // also, reg meta-types, etc
    //

    g.register_mask( m );

        
	

    //
    // In single-VCF mode, read header, set meta-types
    //

    if ( g.single_file_mode() )
      {
	if ( g.single_file_bcf() )
	  {
	    g.bcf = new BCF( project_file , 0 ); // 0=readmode
	    g.bcf->read_header( & g.vardb );
	  }
	else
	  g.vardb.vcf_iterate_read_header( m );
      }


    //
    // Set phenotype (as used by Mask phe.unobs)
    //
    
    if ( args.has( "phenotype-file" ) )
      {		

	// in project mode, this function is basically saying that we want to create a 

	if ( args.has( "phenotype" ) )
	  Helper::halt( "cannot specify both --phenotype and --phenotype-file" );

	// in single-file mode, we need to upload the phenotypes on the fly
	
	// If two cols specified, means load from file / name
	std::vector<std::string> k = args.as_string_vector( "phenotype-file" );
	
	if ( k.size() != 2 ) 
	  Helper::halt( "expecting --phenotype-file filename label" );
	
	g.phmap.direct_load( k[0] , k[1] );       	
      }
    else if ( args.has( "phenotype" ) )
      {
	if ( ! g.single_file_mode() ) 	
	  {
	    if ( ! Pseq::IndDB::set_phenotype( args.as_string( "phenotype" ) ) ) 
	      Helper::halt("no individuals selected / problem setting phenotype " + args.as_string( "phenotype" ));
	  }
	else
	  Helper::halt( "must specify --phenotype-file in single-file mode" );
      }
    else if ( args.has( "make-phenotype") )
      {
	// expect format: factor definition
	std::string k = args.as_string( "make-phenotype" );
	if ( ! Pseq::IndDB::make_phenotype( k ) )
	  Helper::halt("problem setting phenotype " + args.as_string( "make-phenotype" ) );
      }


    else if ( PLINKSeq::DEFAULT_PHENOTYPE() != "" ) 
      {
	if ( ! Pseq::IndDB::set_phenotype( PLINKSeq::DEFAULT_PHENOTYPE() ) )
	  Helper::halt("no individuals selected / problem setting (default) phenotype " + PLINKSeq::DEFAULT_PHENOTYPE() );
      }
    
    //
    // Create a residual-based phenotype (at this point, will only include people who 
    // are in the mask)
    //

    if ( args.has( "residuals" ) )
      {
	if ( ! ( g.phmap.type() == PHE_DICHOT || g.phmap.type() == PHE_QT ) )
	  Helper::halt( "no binary or continuous phenotype defeined for making residuals");
	
	Pseq::IndDB::make_residuals( args.as_string_vector( "residuals" ) );	
      }

    //
    // Set a stratifying variable?
    //

    if ( args.has( "strata" ) )
      {
	if ( ! g.phmap.set_strata( args.as_string( "strata" ) ) )
	  Helper::halt("problem setting strata");
      }




  //
  // Load VCF files into variant database
  //
  
  if ( command == "load-vcf" )
    {

      // Add any VCFs specified on the command line, via --vcf or --file
      
      if ( args.has("vcf") ) 
	{
	  std::vector<std::string> f = args.as_string_vector( "vcf" ); 
	  for (int i=0;i<f.size();i++) 
	    g.fIndex.append_to_projectfile( Helper::fullpath( f[i] ) , "VCF" );
	}
      
      if ( args.has("file") ) 
	{
	  std::vector<std::string> f = args.as_string_vector( "file" ); 
	  for (int i=0;i<f.size();i++) 
	    g.fIndex.append_to_projectfile( Helper::fullpath( f[i] ) , "VCF" );
	}
      
      Pseq::VarDB::load_VCF( m );

      Pseq::finished();
    }
 
    
  //
  // Just index a BGZF'ed VCF into the VARDB
  //

  if ( command == "index-vcf" ) 
  {
      
    if ( ! args.has( "vcf" ) ) 
      Helper::halt( "no VCF file(s) specified, use --vcf file(s)" );
    
    std::vector<std::string> t = args.as_string_vector( "vcf" );
    
    g.vardb.drop_index();
    
    for (int f=0; f<t.size(); f++)
      {
	
	if ( ! Helper::fileExists( t[f] ) )
	  {
	    plog.warn( "could not find VCF" , t[f] );
	    continue;
	  }
	
	// ensure we are using the full path
	t[f] = Helper::fullpath( t[f] );
	
	// Add to project index
	g.fIndex.append_to_projectfile( Helper::fullpath( t[f] ) , "VCFZ" );
	
	// Add to file-map, and create a BCF instance
	VCFZ * vcfz = g.fIndex.add_VCFZ( t[f] );
	
	// Open VCFZ via BGZF interface	    
	vcfz->set_vardb( &g.vardb );
	
	// Read header (uses standard plainext/GZIP reader, then closes)
	vcfz->read_header( m );

	// Re-open using BGZF interface
	vcfz->reading();
	vcfz->open();
		
	// Iterate through file, adding index	    	    
	g.vardb.begin();
		
	uint64_t inserted = 0;
	while ( vcfz->index_record() )
	  {
	    if ( ++inserted % 1000 == 0 )
	      plog.counter1( "parsed " + Helper::int2str( inserted ) + " rows" );
	  }
	plog.counter1("\n");
	
	if ( inserted == 0 ) 
	  {
	    plog.warn( "could not insert any variants: is the VCF BGZF-gzipped?" );	    
	  }
	else
	  plog << "processed " << inserted << " rows (variants & header) from compressed VCF; now finishing index...\n";
	
	g.vardb.commit();
	vcfz->close();
	
	// and calculate summary Ns
	int2 niv = g.vardb.make_summary( t[f] );
	
      }
    
    g.vardb.index();
    
    Pseq::finished();
  }
  


    //
    // Data-views
    //

    if ( command == "v-view" || 
	 command == "rv-view" || 
	 command == "mv-view" ||
	 command == "mrv-view" )
      {		

	const bool rview = command == "rv-view" || command == "mrv-view";
	const bool mview = command == "mv-view" || command == "mrv-view";
	
	OptVView opt;
	opt.simple = args.has("simple");
	opt.vmeta = args.has("vmeta") || args.has("verbose");
	opt.vexpand = args.has("verbose");
	opt.geno = args.has("geno") || args.has("gmeta") || rview || mview;
	opt.gmeta = args.has("gmeta");
	opt.show_samples = args.has("samples");
	opt.show_nonmissing_geno = ! args.has( "hide-null" );
	opt.show_only_minor      =   args.has( "only-minor" ) || rview;
	opt.show_only_alt        =   args.has( "only-alt" ) ;
	if ( opt.show_only_alt || opt.show_only_minor ) opt.show_nonmissing_geno = false;
	opt.mview = mview;
	
	// Output stream
	
	Out output( "vars" , "variant/genotype output" );
	
	IterationReport report = g.vardb.iterate( f_view , &opt , m );

	// Use g-view code to display a set of multiple variants; fix options for display here
	
	if ( mview ) 
	  {
	    output << opt.vars.dump( opt.vmeta , 
				     opt.vexpand , 
				     opt.geno , 
				     opt.gmeta , 
				     false , // transpose
				     false , // rarelist mode 
				     g.phmap.type() != PHE_NONE , 
				     rview ) << "\n";
	  }
	Pseq::finished();
      }
    

    if ( command == "g-view" )
      {
	OptGView opt;

	opt.vmeta = args.has("vmeta") || args.has("verbose");
	opt.geno = args.has("geno") || args.has("gmeta");
	opt.gmeta = args.has("gmeta");
	opt.vexpand = args.has("verbose");
	opt.transpose = args.has("transpose");
	opt.rarelist = args.has("rarelist");
	opt.show_phenotype = args.has("phenotype");

	Out output( "groups" , "variant-group output" );
	
	IterationReport report = g.vardb.iterate( g_view , &opt , m );	

	Pseq::finished();
      }
    

    if ( command == "gs-view" )
      {

	Opt_geneseq opt;
	
	if ( args.has( "ref" ) ) 
	  opt.ref = g.refdb.lookup_group_id( args.as_string( "ref" ) );
	
	if ( args.has( "variant" ) )
	  opt.only_variant_sites = true;

	if ( args.has( "plot" ) )
	  opt.R_plot = true;

	// add phenotype?
	opt.pheno = g.phmap.type() == PHE_DICHOT;
	
	// add protein annotations?
	
	if ( args.has( "protdb" ) )
	  {
	    ProtDBase * pd = new ProtDBase;
	    pd->attach( args.as_string( "protdb" ) );
	    if ( ! pd->attached() ) Helper::halt( "could not attach PROTDB" );
	    opt.protdb = pd;
	    if ( ! args.has( "domain" ) ) 
	      Helper::halt( "no --domain specified with --protdb for gs-view" );
	    opt.protdom = args.get_set( "domain" );
	  }

	if ( ! g.seqdb.attached() ) Helper::halt( "no SEQDB attached" );
	if ( ! g.locdb.attached() ) Helper::halt( "no LOCDB attached" );

	Out output( "gsview" , "variants in sequence context" );
	
	Out * outplot = args.has( "plot" )
	  ? new Out( "gsview.R" , "R script for gs-view plots" )
	  : NULL;

	IterationReport report = g.vardb.iterate( g_geneseq , &opt , m );
	
	if ( opt.protdb ) delete opt.protdb;
	
	if ( opt.R_plot ) delete outplot;

	Pseq::finished();
      }


    if ( command == "i-view" )
      {	

	Out output( "indiv" , "per-individual information" );

	if ( g.single_file_mode() )
	  Pseq::VarDB::header_VCF( false , true , m );
	else if ( args.has( "from-vardb" ) )
	  Pseq::VarDB::dump_indiv();
	else
	  Pseq::IndDB::dump_table( m );
	Pseq::finished();
      }



    if ( command == "write-phe" || command == "i-matrix" )
      {	


	std::vector<std::string> names;	

	if ( ! args.has( "name" ) ) 
	  {
	    if ( g.phmap.type() != PHE_NONE ) names.push_back( g.phmap.phenotype() );
	    else Helper::halt( "no phenotype specified: --phenotype (or --name pheno1 pheno2 ...)" );
	  }
	else names = args.as_string_vector( "name" ) ;
	
	Out output( command == "i-matrix" ? "indiv.matrix" : "phe" , "individual phenotype data" );
	
	Pseq::IndDB::dump_phenotypes( names , command == "i-matrix" ) ;
	
	Pseq::finished();
      }


    //
    // Dump individual segments from a SEGDB
    //

    if ( command == "seg-view" )
      {

	Out output( "segs" , "segments from SEGDB group" );
	
	if ( ! args.has( "group" ) ) 
	  Helper::halt("requires a locus-group to be specified");
	std::string grp = args.as_string( "group" );
	for (int i=0;i<g.indmap.size(); i++)
	  {
	    Individual * person = g.indmap(i);
	    plog << person->id() << "\t";
	    std::set<Region> s = g.segdb.get_indiv_regions( grp , person->id() );
	    plog << s.size();
	    std::set<Region>::iterator si = s.begin();
	    while ( si != s.end() ) 
	      {
		plog << "\t" << si->coordinate();
		++si;
	      }
	    plog << "\n";
	  }
	Pseq::finished();
      }
    

    //
    // Breakdown of { 0:1, 1:0 } etc
    //

    if ( command == "v-dist" )
      {	

	Out output( "vdist" , "output from v-dist command" );


	long int nrep = args.has( "perm" ) ? args.as_int( "perm" ) : 1000 ;
	Pseq::VarDB::vdist_summary( m , nrep );
	Pseq::finished();
      }


    //
    // Additional QC measures
    //

    if ( command == "v-freq" )
      {

	Out output( "vfreq" , "output from v-freq command" );

	XQCstats x;
	
	// ensure EM phasing
	
	if ( args.has( "em" ) )
	  {
	    m.EM_caller( true );
	    double threshold = args.as_float( "em" );	    
	    m.EM_threshold( threshold );
	    x.em_stats = true; 
	  }

	// header

	output.data_reset();	
	output.data_group_header( "VAR" );

	output.data_header( "CHR" );
	output.data_header( "POS" );
	output.data_header( "REF" );
	output.data_header( "ALT" );

	output.data_header( "FILTER" );
	output.data_header( "QUAL" );
	output.data_header( "TI" );  // flag to indicate Ti/Tv

	output.data_header( "GENO" );
	
	if ( Genotype::using_soft_calls )
	  {
	    output.data_header( "INFO" );
	    output.data_header( "DFRQ" );
	  }
	
	output.data_header( "MAC" );
	output.data_header( "MAF" );
	output.data_header( "REFMIN" );
	output.data_header( "HWE" );
	output.data_header( "HET" );
		
	if ( x.em_stats ) 
	  {
	    output.data_header( "MPP" );
	    output.data_header( "H" );
	    output.data_header( "HALT" );
	    output.data_header( "EMMAF" ); 
	  }

	output.data_header( "NSNP" );	

	output.data_header_done();

	// Run 

	IterationReport report 
	  = g.vardb.iterate( f_extra_qc_metrics , &x , m );
	
	// Flush out remaining items
	x.flush();
	Pseq::finished();
      }
    
    
    
    //
    // PolyPhen2 scoring of variants in database 
    //
    
    if ( command == "score-weights" )
      {
	Helper::halt("obsolete command score-weights");

	// Out output( "scores" , "output from score-weights command" );	
 	// std::string dbname = Pseq::Util::single_argument<std::string>( args , "name" );
 	// Pseq::PPH2DB::score( m , dbname );
 	// Pseq::finished();
      }


    //
    // Clusters of variants
    //

    if ( command == "clusters" )
      {
	Out output( "clusters" , "output from clusters command" );
	Pseq::VarDB::cluster_scan( m );
	Pseq::finished();
      }


    //
    // Proximity scan, for clusters of variants with LD calculation also
    //
    
    if ( command == "proximity-scan" )
      {
	Out output( "proximity" , "output from proximity-scan" );
	Pseq::VarDB::proximity_scan( m );
	Pseq::finished();
      }



    //
    // Create a summary counts file
    //

    if ( command == "counts" || command == "g-counts" )
      {
	if ( args.has( "output-vcf" ) )
	  {
	    Out output( "vcf" , command == "counts" ? "VCF from counts" : "VCF from g-counts" );
	    std::string proj = Pseq::Util::single_argument<std::string>( args , "name" );
	    Pseq::VarDB::make_counts_file( m , proj );
	  }
	else
	  {	 	    
	    Out output( command == "counts" ? "counts" : "gcounts" , 
			command == "counts" ? "VCF from counts" : "VCF from g-counts" );
	    bool genotypes = command == "g-counts" ;
	    Pseq::VarDB::simple_counts( m , genotypes , false ); // false --> not QT
	  }
	Pseq::finished();
      }


    //
    // Similar to the above, but for QTs
    //
    
    if ( command == "means" )
      {
	Out output( "means" , "QT means stratified by variant" );
	bool genotypes = command == "g-counts" ;
	Pseq::VarDB::simple_counts( m , genotypes , true ); // true --> using QT
	Pseq::finished();
      }



    //
    // Lookup DB information for a list of positions
    //

    if (command == "lookup") {
    	Out output("meta", "output from lookup command");

    	if (args.has("worstAnnotationPriorities"))
    		Annotate::setWorstAnnotationPriorities(Pseq::Util::single_argument<std::string>(args, "worstAnnotationPriorities"));

    	if (args.has("file")) {
    		std::string filename = Pseq::Util::single_argument<std::string>(args, "file");
    		Pseq::VarDB::lookup_list(filename, m);
    	}
    	if (args.has("region")) {
    		std::vector<std::string> regions = Pseq::Util::n_arguments<std::string>(args, "region");
    		std::vector<Region> regs;
    		for (int i = 0; i < regions.size(); i++) {
    			bool okay;
    			Region r(regions[i], okay);
    			if (okay)
    				regs.push_back(r);
    			else
    				plog.warn("could not parse region: " + regions[i]);
    		}
    		Pseq::VarDB::lookup_list(".", m, &regs);
    	}

    	// else directly annotate variants in VARDB
    	if (!(args.has("file") || args.has("region")))
    		Pseq::VarDB::lookup_list(".", m);
    	Pseq::finished();
    }
    

    //
    // Database summaries
    //

    if ( command == "summary" )
      {

	bool ugly = args.has( "ugly" ); 
	if ( ! ugly ) plog << "\n";
	Pseq::VarDB::summary(m , ugly );
	Pseq::IndDB::summary( ugly );
	Pseq::LocDB::summary(&g.locdb , ugly );
	Pseq::LocDB::summary(&g.segdb , ugly );
	Pseq::RefDB::summary( ugly );
	Pseq::SeqDB::summary( ugly );	
	Pseq::Util::file_summary( ugly );
	Pseq::Util::meta_summary( ugly );
	Pseq::finished();
      }


    if ( command == "var-summary" )
      {
	Pseq::VarDB::summary(m , args.has( "ugly" ) );
	Pseq::finished();
      }

    if ( command == "ind-summary" )
      {
	Pseq::IndDB::summary( args.has( "ugly" ) );
	Pseq::finished();
      }

    if ( command == "loc-summary" )
      {
	Pseq::LocDB::summary( &g.locdb , args.has( "ugly" ) );
	Pseq::finished();
      }

    if ( command == "seg-summary" )
      {
	Pseq::LocDB::summary( &g.segdb , args.has( "ugly" ) );
	Pseq::finished();
      }

    if ( command == "ref-summary" )
      {
	Pseq::RefDB::summary( args.has( "ugly" ));
	Pseq::finished();
      }

    if ( command == "seq-summary" )
      {
	Pseq::SeqDB::summary( args.has( "ugly" ) );
	Pseq::finished();
      }

    if ( command == "file-summary" )
      {
	Pseq::Util::file_summary( args.has( "ugly" ) );
	Pseq::finished();
      }

    if ( command == "meta-summary" )
      {
	Pseq::Util::meta_summary( args.has( "ugly" ) );
	Pseq::finished();
      }


    //
    // Calculated summary statistics
    //

    if ( command == "v-stats" )
      {	
	
	Out output( "vstats" , "output from v-stats command" );
	
	// Hard-code some of these values for now
	
	Pseq::VStat vstat(&g);
	
	Pseq::Util::set_default( vstat );
	
	IterationReport report = g.vardb.iterate( f_vstat , &vstat , m );	
	
	vstat.report();
	  
	Pseq::finished();
      }
    
    if ( command == "g-stats" )
      {
	
	Out output( "gstats" , "output from g-stats command" );

	if ( ! m.any_grouping() )
	  Helper::halt("no gene-grouping specified for g-stat");
	if ( ! m.group_loc() ) 
	  Helper::halt("currently g-stats only supports loc.group");

	Pseq::VStat vstat(&g);

	Pseq::Util::set_default( vstat );

	Pseq::GStat aux(&g, m.group_set(), vstat );
	
	Pseq::VarDB::gene_stats_header( vstat );

	IterationReport report = g.vardb.iterate( g_gstat , &aux , m );	
	Pseq::finished();
      }
    

    if ( command == "i-stats" )
      {
	Out output( "istats" , "output from i-stats command" );
	Pseq::IStat istat(&g);
	IterationReport report = g.vardb.iterate( f_istat , &istat , m );	
	istat.report();
	Pseq::finished();
      }



    //
    // Per-individual posterior class probabilities under simple pop. model
    //

    if ( command == "i-pop" )
      {
	if ( ! args.has( "file" ) ) Helper::halt( "no --file {pop-allele-freq-table} specified" );	
	SeqDBase * s = &g.seqdb;
	if ( ! g.seqdb.attached() ) 
	  {
	    plog.warn( "no SEQDB attached: will not be able to check REF/ALT status of A1/A2" );
	    s = NULL;
	  }
	Out output( "ipop" , "per-individual posterior class probabilities" );
	Pseq::IPop ipop( args.as_string( "file" ) , s );
	IterationReport report = g.vardb.iterate( f_ipop , &ipop , m );
	ipop.calculate();
	Pseq::finished();
      }


    //
    // Per-locus simple view, or sequence stats, e.g. GC percent 
    //

    if ( command == "loc-view" )
    {

      Out output( "loci" , "output from loc-view command" );

      if ( ! ( g.locdb.attached() | g.segdb.attached() ) ) Helper::halt("no LOCDB or SEGDB attached");
      if ( ! args.has( "group" ) ) Helper::halt("requires a locus-group to be specified");
	
      std::vector<std::string> grp = args.as_string_vector( "group" );
      if ( grp.size() != 1 ) Helper::halt("requires a single locus-group to be specified");	
      
      std::vector<std::string> alias;
      if ( args.has( "alias" ) ) alias = args.as_string_vector( "alias" );
      
      Pseq::LocDB::loc_view( grp[0] , alias , ! args.has("no-meta") , args.has("show-subregions") );
      
      Pseq::finished();
      }


    if ( command == "intersect" )
      {

	Out output( "loci" , "output of intersect command" );

	// print LOC group names that have 1+ VARDB variants, according to 
	// the mask
	
	// assumes a mask will specifiy the 'loc.group'

	g.vardb.iterate( g_loc_view , NULL , m );
	Pseq::finished();
      }


    if ( command == "loc-stats" )
      {	
	
	Out output( "locstats" , "output from loc-stats command" );

	if ( ! ( g.locdb.attached() | g.segdb.attached() ) ) Helper::halt("no LOCDB or SEGDB attached");
	
	if ( ! g.seqdb.attached() ) Helper::halt("no SEQDB attached");
	
	if ( ! args.has( "group" ) ) 
	  Helper::halt("requires a locus-group to be specified, with --group");
	
	std::vector<std::string> grp = args.as_string_vector( "group" );

	std::string refgroup = args.has( "ref" ) ? args.as_string( "ref" ) : "." ; 

	for (int i=0; i<grp.size(); i++) 
	  Pseq::SeqDB::loc_stats( grp[i] , refgroup );

	Pseq::finished();
      }
    
    //
    // View AA sequence and details
    //

    if ( command == "loc-translate" )
      {	

	Out output1( "loci" , "output from loc-translate" );
	Out output2( "bstats" , "nucleotide base-level statistics" );

	if ( ! ( g.locdb.attached() | g.segdb.attached() ) ) Helper::halt("no LOCDB or SEGDB attached");       
	if ( ! g.seqdb.attached() ) Helper::halt("no SEQDB attached");	
	std::string grp = Pseq::Util::single_argument<std::string>( args , "group" );	
	Pseq::SeqDB::loc_translate( grp );	
	Pseq::finished();
      }


    //
    // Dump RefVariants from REFDB
    //

    if ( command == "ref-view" )
      {
	Out output( "refvars" , "output from ref-view" );

	if ( ! g.refdb.attached() ) Helper::halt("no REFDB attached");
	if ( ! args.has( "group" ) ) Helper::halt("no group specified");
	std::string grp = Pseq::Util::single_argument<std::string>( args , "group" );
	bool with_meta = args.has( "vmeta" );
	bool with_verbose = args.has( "verbose" );
	g.refdb.dump( grp , with_meta , with_verbose );
	Pseq::finished();
      }


    //
    // Dump sequence from SEQDB
    //

    if ( command == "seq-view" )
      {	

	Out output( "seq" , "output from seq-view" );

	if ( ! g.seqdb.attached() ) Helper::halt("no SEQDB attached");
	
	if ( ! args.has( "region" ) )
	  Helper::halt("need to specify --region");

	bool compact = args.has( "compact" );

	std::vector<std::string> regions = Pseq::Util::n_arguments<std::string>( args, "region" );
	for (int i=0;i<regions.size(); i++)
	  {
	    bool okay;
	    Region r( regions[i] , okay );
	    if ( okay )
	      g.seqdb.dump( r , compact );
	    else	      
	      Helper::halt("could not parse region: " + regions[i] );
	  }

	Pseq::finished();
      }


    //
    // Family-based operations
    //

    if ( command == "denovo" ) 
      {
	Out output1( "denovo.vars" , "per-site output from denovo" );
	Out output2( "denovo.indiv" , "per-trio output from denovo" );
	Pseq::VarDB::denovo_scan( m );
	Pseq::finished();
      }

    //
    // Family-based operations for XHMM-like output
    //

    if ( command == "cnv-denovo" )
    {
    	Out output1( "denovo.cnv",       "per-site output from cnv-denovo" );
    	Out output2( "denovo.cnv.indiv", "per-trio output from cnv-denovo" );
    	Pseq::VarDB::cnv_denovo_scan( m );
    	Pseq::finished();
    }

    // 
    // Per indiviaul/group unique/enriched listing
    //
    
    if ( command == "unique" )
      {
	
	Out output( "uniq" , "output from unique command" );

	if ( ! args.has("indiv") ) 
	  Helper::halt("no individuals specified");

	std::vector<std::string> indiv = args.as_string_vector( "indiv" );

	OptUniq opt;

	if ( args.has( "require" ) )
	  opt.ingroup_req = args.as_int( "require" );

	if ( args.has( "allow" ) ) 
	  opt.outgroup_allow = args.as_int( "allow" );

	Pseq::VarDB::uniq_report( indiv , m , opt );
	Pseq::finished();
      }


    //
    // Concordance test
    //

    if ( command == "concordance" )
      {
	Out output1( "concord" , "summary from concordance command" );
	Out output2( "concord.geno" , "discordant genotypes" );
	Out output3( "concord.indiv" , "per-individual concordance data" );
	Out output4( "concord.vars" , "per-site concordance data" );
	
	Pseq::VarDB::check_concordance(m);
	Pseq::finished();
      }

    //
    // IBS sharing matrix
    //

    if ( command == "ibs-matrix" )
      {
	Out output( "ibs" , "output from ibs-matrix command" );
	Pseq::IBS::calculate(m);
	Pseq::finished();
      }


    // 
    // Association tests (with group or phenotype) 
    //

    if ( command == "group-comparison" )
      {
	Out output( "comp" , "output from group-comparison command" );
	Pseq::Assoc::group_comparison( m );
      }

    if ( command == "assoc" )
      {

	Out output1( "assoc" , "output from assoc command" );
	Out outputd( "assoc.det", "site-specific breakdown of included variants per gene" );
	
	Out * output2 = NULL;       
	if ( args.has( "tests", "two-hit") )
	  output2 = new Out( "twohit.vars" , "variants from two-hit test" );
	
	Out * outcarriers = args.has( "carriers" ) ? new Out( "assoc.carriers" , "ALT carriers for tests" ) : NULL ;
	
	Out * outmatrix = args.has( "dump-null-matrix" ) 
	  ? new Out( "matrix" , "permuted null statistic matrix" ) 
	  : NULL;
	
	// if no perms specified, use adaptive permutation mode
	Pseq::Assoc::set_assoc_test( m , args );
	
	if ( output2 ) delete output2;	
	if ( outmatrix ) delete outmatrix;
	if ( outcarriers ) delete outcarriers;
	
	Pseq::finished();
	
      }
    
    
    if ( command == "net-assoc" )
      {
	Out output( "net.assoc" , "output from net-assoc command" );

	if ( g.phmap.type() != PHE_DICHOT ) 
	  Helper::halt("no dichotomous phenotype specified");
	Pseq::Assoc::net_assoc_test( m , args );
	Pseq::finished();
      }


    
    //
    //
    //

    if ( command == "indiv-enrich" )
      {
	Out output( "indiv.enrich" , "output from indiv-enrich command" );
 	Pseq::Assoc::set_enrich_wrapper( m , args );
 	Pseq::finished();
      }
    

    if ( command == "s-assoc" )
      {
	Out output( "sassoc" , "output from s-assoc command" );

	if ( g.phmap.type() != PHE_DICHOT ) 
	  Helper::halt("no dichotomous phenotype specified");
	
	if ( !args.has("perm") )
	  Helper::halt("no permutations specified");

	if ( !args.has("file") ) 
	  Helper::halt("need to specify --file");
	
	std::vector<std::string> s = args.as_string_vector( "file" );
	if ( s.size() != 2 )
	  Helper::halt("need 2 files to be specified");

	Pseq::IBD::test_wrapper( s[0], s[1] , args.as_int( "perm" ) , m );

	Pseq::finished();
      }



    if ( command == "ibd-load" )
      {
	if ( !args.has("file") ) Helper::halt("need to specify --file");
	if ( !args.has("ibddb") ) Helper::halt("need to specify --ibddb");
	Pseq::IBD::load_wrapper( args.as_string( "file" ) , args.as_string( "ibddb" ) );
	Pseq::finished();
      }
    
    if ( command == "ibd-sharing" )
      {

	Out output( "ibd" , "output from ibd-sharing command" );

	if ( g.phmap.type() != PHE_DICHOT ) 
	  Helper::halt("no dichotomous phenotype specified");
	if ( !args.has("ibddb") ) Helper::halt("need to specify --ibddb");
	Pseq::IBD::sharing_wrapper( args.as_string( "ibddb" ) , m );
	Pseq::finished();
      }

    
    if ( command == "mutation-screen" )
      {
	Out output( "mut" , "output from mutation-screen command" );
	if ( ! args.has( "ibddb" ) ) Helper::halt( "need to specify --ibddb {filename}" );
	if ( ! args.has( "indiv" ) ) Helper::halt( "need to specify --indiv {ID}" );
	if ( ! args.has( "region" ) ) Helper::halt( "need to specify --region {chr1:1234567}" );
	Pseq::IBD::mutation_wrapper( args.as_string( "ibddb" ) , 
				     args.as_string( "indiv" ) , 
				     args.as_string( "region" ) , 
				     m );
	Pseq::finished();
      }
    




    //
    // Single-site association statistics
    // 
    
    if ( command == "v-assoc" )
      {

	Out output( "vassoc" , "output from v-assoc command" );

	Pseq::Assoc::Aux_vassoc_options aux;
	if ( args.has( "info" ) ) aux.show_istat = true;
	if ( args.has( "separate-chr-bp" ) ) aux.separate_chr_bp = true;
	if ( args.has( "vmeta" ) ) aux.show_meta = true;
	aux.nrep = args.has( "perm" ) ? args.as_int( "perm" ) : 0 ;

	if ( g.phmap.type() != PHE_DICHOT ) 
	  Helper::halt( "no dichotomous phenotype specified" );	

	Pseq::Assoc::variant_assoc_test( m , aux , args );

	Pseq::finished();
      }


    //
    // GLM single-variant tests
    // 
    

    if ( command == "glm" )
      {

	Out output( "glm" , "output from glm command" );	

	// allow phenotpe not to be specified yet, if reading from a file a list of tests:
	if ( g.phmap.type() != PHE_DICHOT && g.phmap.type() != PHE_QT && ! args.has("file" ) ) 
	  Helper::halt("no dichotomous or quantitative phenotype specified");	

	Pseq::Assoc::Aux_glm aux;
	
	aux.show_meta = args.has( "vmeta" );
	
	aux.dichot_pheno = g.phmap.type() == PHE_DICHOT;
	
	if ( args.has( "file" ) )
	  {
	    aux.test_list = 1;
	    aux.test_list_file = args.as_string( "file" );
	  }
	
	if ( args.has( "use-postprobs" ) )
	{
	    aux.use_postprobs = true;
	    aux.softtag = args.as_string( "use-postprobs" );
	}
	else if ( args.has( "use-dosages" ) )
	{
	    aux.use_dosage = true;
	    aux.softtag = args.as_string( "use-dosages" );
	}
	
	if ( args.has( "covar" ) )
	  {
	    aux.has_covar = true;
	    aux.covars = args.as_string_vector( "covar" );	    
	  }

	if ( args.has( "show-covar" ) )
	  aux.show_all_covar = true;

	if ( args.has( "show-intercept" ) )
	  aux.show_intercept = true;
	
	// for now, no permutations
	// aux.nrep = args.has( "perm" ) ? args.as_int( "perm" ) : 0 ;
	
	Pseq::Assoc::glm_assoc_test( m , aux );

	Pseq::finished();
      }


    
    //
    // View a ProtDB 
    //
    
    if ( command == "prot-view" )
      {
	
	if ( ! args.has( "protdb" ) ) Helper::halt( "no --protdb specified" );      
	
	Out output( "prot" , "protein domain/annotations" );
	
	if ( args.has( "group" ) && args.has("name") )
	  Pseq::ProtDB::lookup( args.as_string( "protdb" ) , 
				args.as_string( "name" ) , 
				args.as_string( "group" ) , &m );
	else if ( args.has( "name" ) )
	  Pseq::ProtDB::lookup( args.as_string( "protdb" ) , 
				args.as_string( "name" ) );
	else
	  Pseq::ProtDB::lookup( args.as_string( "protdb" ) );
	
	Pseq::finished();
      }

    
    //
    // Summary for a PROTDB
    //

    if ( command == "prot-summary" )
      {	
	if ( ! args.has( "protdb" ) ) Helper::halt( "no --protdb specified" );      	
	Pseq::ProtDB::summary( args.as_string( "protdb" ) );	
	Pseq::finished();
      }
    

    //
    // Data-dumpers
    //
    
    if ( command == "write-vardb" )
      {       
	if ( ! args.has( "new-vardb" ) ) Helper::halt("no --new-vardb for write-vardb specified");
	if ( ! args.has( "new-project" ) ) Helper::halt("no --new-project for write-vardb specified");
	std::string proj_file = Pseq::Util::single_argument<std::string>( args , "new-project" );
	std::string db_file = Pseq::Util::single_argument<std::string>( args , "new-vardb" );	
	Pseq::VarDB::write_vardb( proj_file , db_file ,  m);
	Pseq::finished();
      }
    
    
    if ( command == "write-vcf" )
      {

	Out output( "vcf" , "output write-vcf command" );
	
	bool compressed = args.has( "format" , "BGZF" );
	Pseq::VarDB::write_VCF( m , compressed );
	Pseq::finished();
      }
  
    
    if ( command == "write-bcf" )
      {	

	//	Out output( "bcf" , "BGZF-compressed BCF file from write-bcf" );
	
	if ( ! args.has( "bcf" ) ) 
	  Helper::halt( "need to specify --bcf {output.bcf}" );
	
	Pseq::VarDB::write_BCF( m , args.as_string( "bcf" ) );
	
	Pseq::finished();
      }
    
    
    if ( command == "write-ped" )
      {

	Out output1( "tped" , "transposed PED file" );
	Out output2( "tfam" , "associated FAM file" );
	
	Pseq::VarDB::write_PED( m , args.has( "family-id" ) );
	Pseq::finished();
      }


    if ( command == "write-lik" )
      {
	Pseq::VarDB::write_lik(m);
	Pseq::finished();
      }    
    

    if ( command == "write-haps" )
      {
	if ( ! args.has( "name" ) )
	  Helper::halt("requires --name to be specified");
	Pseq::VarDB::write_haps(m,args.as_string("name"));
	Pseq::finished();
      }    


//     // Create a v-matrix style output, but for a fixed set of variants
//     // (i.e. putting NA if no data)

//     if ( command == "v-lookup" )
//       {
// 	std::string filename = args.has( "file" ) ? args.as_string( "file" ) : "." ;
// 	std::vector<std::string> regs;
// 	if ( args.has( "region" ) ) regs = args.as_string_vector( "region" );
	 
// 	Pseq::VarDB::write_lookup_matrix(m,filename,regs);

// 	Pseq::finished();
//       }


    if ( command == "v-matrix" )
      {	
	Out output( "matrix" , "matrix of genotypes as allele counts" );	
	Pseq::VarDB::write_matrix(m);
	Pseq::finished();
      }


    if ( command == "meta-matrix" )
      {	
	Out output( "matrix" , "matrix of variant-level meta-data" );
	Pseq::VarDB::write_meta_matrix(m);
	Pseq::finished();
      }

    if ( command == "v-meta-matrix" )
      {	
	Out output( "matrix" , "matrix of genotype-level meta-data" );
	std::string name = Pseq::Util::single_argument<std::string>( args, "name" );
	Pseq::VarDB::write_var_meta_matrix(m,name);
	Pseq::finished();
      }

    if ( command == "loc-annotate" )
      {
	std::string grp = Pseq::Util::single_argument<std::string>( args , "group" );
	Pseq::VarDB::annotate_loc(grp,m);
	Pseq::finished();
      }

    if ( command == "loc-overlap" )
      {
	Pseq::LocDB::overlap_analysis();
	Pseq::finished();
      }

    if ( command == "g-matrix" )
      {	
	// options:
	// 1) only show genes with non-zero variance
			Out output( "matrix" , "variant/genotype matrix output" );
	OptGMatrix opt(&g);
	if ( args.has( "hide-invariant" )  )
	  opt.hide_zero_variance = true;
	if ( args.has( "collapse" ) ) 
	  opt.collapse_01 = true;
	Pseq::VarDB::write_gene_matrix(m,opt);
	Pseq::finished();
      }

    
    if ( command == "g-meta-matrix" )
      {	
				Out output( "matrix" , "variant/genotype matrix output" );
	OptGMetaMatrix opt;
	opt.name = Pseq::Util::single_argument<std::string>( args, "name" );
	opt.show_mean = true;
	Pseq::VarDB::write_gene_meta_matrix(m,opt);
	Pseq::finished();
      }


    //
    // Load variant lists into the VARDB
    //
    
    if ( command == "var-set" ) 
      {
	
	// either from a file; or all mask-passing variants
	
	if ( args.has( "file" ) )
	  Pseq::VarDB::add_to_varset( args.as_string( "file" ) ); 
	else 
	  {
	    if ( ! args.has( "group" ) ) Helper::halt( "need to specify a --group" );
	    if ( args.has( "name" ) ) 
	      Pseq::VarDB::add_to_varset( args.as_string( "group" ) , m , args.as_string( "name" ) );
	    else 
	      Pseq::VarDB::add_to_varset( args.as_string( "group" ) , m ); 
	  }
	Pseq::finished();
      }
    

    if ( command == "var-drop-set" )
      {
	if ( ! args.has( "group" ) ) Helper::halt( "need to specify a --group" );
	g.vardb.drop_set( args.as_string( "group") );
	Pseq::finished();
      }
        
    if ( command == "var-drop-all-sets" )
      {
	g.vardb.drop_set( "_ALL_" );
	Pseq::finished();
      }

    //
    // Group a bunch of variants sets
    //
    
    if ( command == "var-superset" )
      {
	std::string desc = args.has( "description" ) ? args.as_string( "description" ) : ".";
	
	// Read from file
	if ( args.has( "file" ) ) 
	  {
	    Pseq::VarDB::add_superset_from_file( args.as_string( "file" ) );
	    Pseq::finished();
	  }
	else if ( args.has( "group") )
	  {
	    if ( ! args.has( "members" ) ) Helper::halt( "need to specify --members with --group" );
	    std::vector<std::string> m = args.as_string_vector( "members" );
	    Pseq::VarDB::add_superset( args.as_string( "group" ) , m , desc );
	    Pseq::finished();
	  }
	else 
	  Helper::halt("need to specify --group and --members, or --file" );
      }


    if ( command == "var-drop-superset" )
      {
	if ( ! args.has( "group" ) ) Helper::halt( "need to specify a --group" );
	g.vardb.drop_superset( args.as_string( "group") );
	Pseq::finished();
      }
        
    if ( command == "var-drop-all-supersets" )
      {
	g.vardb.drop_superset( "_ALL_" );
	Pseq::finished();
      }
    
        
    Pseq::finished();
    return 0;
    
}

