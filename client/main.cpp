
#include <iostream>
#include <iterator>

#include "pseq.h"

#include "func.h"
#include "views.h"
#include "summaries.h"
#include "assoc.h"
#include "compare.h"
#include "ibd.h"
#include "ibs.h"
#include "extra.h"


using namespace std;

GStore g;

Pseq::Util::Options args;
Pseq::Util::Commands pcomm;

std::string PSEQ_VERSION = "0.08";
std::string PSEQ_DATE    = "25-Nov-2011";

int main(int argc, char ** argv)
{
    

  //
  // Get command-line options into a sensible form
  //

  // pseq {project} {command} {--options}
  
  args.load( argc , argv );


  // connect with GSEQ job tracking?
  
  if ( args.has( "history" ) )
    {
      std::vector<std::string> tok = args.as_string_vector( "history" );
      if ( tok.size() != 2 ) Helper::halt( "--history {file} {job#}" );
      g.gseq_tracking( tok[0] , tok[1] );
    }


  
  // known PSEQ commands, and descriptions
  
  pcomm.attach( &args );

  Pseq::Util::populate_commands( pcomm );
       

  //
  // help message?
  //
  
  if ( args.help() ) exit(0);
  

  //
  // Reporting and logging options
  //
  
  if ( args.has("debug") ) 
    debug.silent( false );

  if ( args.has("ignore-warnings") )
    plog.show_warnings( false );

  if ( args.has("early-warnings") )
    plog.early_warnings( true );

  if ( args.has("debug-file") )
    debug.logfile( args.as_string("debug-file") ); 
  
  if ( args.has("silent" ) )
    plog.silent( true );
  
  if ( args.has("out-file") ) 
    plog.logfile( args.as_string("out-file") );
  
  if ( args.has("prolix-file") ) 
    plog.prolix_logfile( args.as_string("prolix-file") );
  
  if ( args.has("long") )
    plog.longmode();

  if ( args.has("long-header") )
    { plog.longmode(); plog.header(); } 



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
  // Set project, attaching all relevant databases
  //
  
  std::string project_file = FileMap::tilde_expansion( args.project_file() );
  


  //
  // If a single VCF has been specified as the 'project'
  //
  
  g.single_file_mode( project_file == "-" 
		      || Helper::ends_with( project_file , ".vcf" ) 
		      || Helper::ends_with( project_file , ".vcf.gz" ) );
  

  if ( g.single_file_mode() && ! pcomm.single_VCF_mode( command ) )
    Helper::halt( command + " not applicable in single-VCF mode" );
  

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
      
      plog << "PSEQ" << "\t" 
	   << PSEQ_VERSION << "("
	   << PSEQ_DATE << ")\n";
      
      g.show_version();
	
      Pseq::finished();
  
    }


  //
  // Simple routine to create dummy genic VariantGroups 
  // and feed to the genic association tests (designed 
  // for testing more than anything else)
  //

  if ( command == "simple-sim" ) 
    {
      Pseq::VarDB::simple_sim();
      Pseq::finished();
    }
  
  
  //
  // Otherwise, open existing project
  //
  
  if (  g.single_file_mode() ) 
  {
      if ( project_file != "-" && ! Helper::fileExists( project_file ) ) 
	  Helper::halt( "could not open VCF file " + project_file );      
  }
  else
  {
      if ( ! Pseq::set_project( project_file ) )
	  Helper::halt("Could not open project file " + project_file );      
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
      
      if ( ! args.has( "type" ) ) 
	Helper::halt("no --type specified");
      std::string type = args.as_string( "type" );
      
      g.fIndex.append_to_projectfile( Helper::fullpath( pname ) , type );
      Pseq::finished();
    }
  


  //
  // Misc. formatting/display options
  //

  
  //
  // Misc. formatting/display options
  //
  
  if ( args.has("whitespace") )
    {
      PLINKSeq::DELIM() = " \t\n";
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
      std::string dbname = Pseq::Util::single_argument<std::string>( args , "name" );
      std::string filename = Pseq::Util::single_argument<std::string>( args , "file" );
      Pseq::PPH2DB::load( dbname , filename );
      Pseq::finished();
    }
  


  //
  // Functions that do not depend on the mask, variant or individual database
  //
  
  if ( command == "intersect" ) 
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
      Pseq::NetDB::lookup( args.as_string( "netdb" ) , args.as_string( "name" ) , args.as_string( "group" ) );
      Pseq::finished();
    }


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
  
  if ( args.has("mask") )
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
  
  if ( args.has("exclude") )
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
		  plog.counter( "parsed " + Helper::int2str( inserted ) + " rows" );
	  }
	  plog.counter("\n");
	  
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
  // Load/flush meta-information into VARDB
  //
  
  if ( command == "load-meta" )
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
      if ( ! args.has("group") )
	g.vardb.flush_indep_meta( );
      else
	{
	  std::string group = Pseq::Util::single_argument<std::string>( args , "group" );
	  g.vardb.flush_indep_meta( group );
	}
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
  
  if ( command == "delete-var" )
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
	  }
	else if ( Helper::ends_with( s[0] , ".reg" ) || Helper::ends_with( s[0] , ".reg.gz" ) )
	  Pseq::LocDB::load_generic_regions( s[0], grp[0] , args , true );
	else Helper::halt("invalid file name, expecting extension: .gtf .gtf.gz .reg .reg.gz");
	
	Pseq::finished();
      }
    
    if ( command == "loc-swap-names" )
      {
	std::string filename = Pseq::Util::single_argument<std::string>( args , "file" );
	std::string group = Pseq::Util::single_argument<std::string>( args , "group" );
	Pseq::LocDB::swap_alternate_names( group , filename );
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
    

    //
    // Any annotation specified? If so, must load transcripts 
    // prior to constructing the Mask
    //

    if ( args.has( "annotate" ) )
    {
	std::string agrp = Pseq::Util::single_argument<std::string>( args , "annotate" );	
	if ( ! Pseq::SeqDB::load_transcripts( agrp ) ) 
	    Helper::halt("problem loading annotation transcripts");
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
      maskspec += " ex-vcf=" + project_file;
    

    //
    // If a phenotype, or covariates, have been specified, by default require that we see these
    // This can be over-riden by adding phe.allow.missing to the --mask 
    //

    std::string phenotype_name = PLINKSeq::DEFAULT_PHENOTYPE();
    if ( ! g.single_file_mode() )
      {
	if ( args.has("phenotype") ) phenotype_name = args.as_string("phenotype");
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
      g.vardb.vcf_iterate_read_header( m );



    //
    // Set phenotype (as used by Mask phe.unobs)
    //
    
    if ( args.has("phenotype") )
      {		
	// in single-file mode, we need to upload the phenotypes on the fly
	
	if ( g.single_file_mode() ) 
	  {
	    // in this case, expect phenotype to be two columns	    
	    std::vector<std::string> k = args.as_string_vector( "phenotype" );
	    if ( k.size() != 2 ) Helper::halt( "expecting --phenotype filename label" );
	    g.phmap.direct_load( k[0] , k[1] );
	  }
	else
	  {
	    if ( ! Pseq::IndDB::set_phenotype( args.as_string( "phenotype" ) ) ) 
	      Helper::halt("no individuals selected / problem setting phenotype " + args.as_string( "phenotype" ));
	  }
      }

    else if ( args.has("make-phenotype") )
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
    // Set a stratifying variable?
    //

    if ( args.has("strata") )
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
	      plog.counter( "parsed " + Helper::int2str( inserted ) + " rows" );
	  }
	plog.counter("\n");
	
	plog << "inserted " << inserted << " variants from compressed VCF; now finishing index...\n";
	
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

	IterationReport report = g.vardb.iterate( f_view , &opt , m );

	// Use g-view code to display a set of multiple variants; fix options for display here
	
	if ( mview ) 
	  {
	    plog << opt.vars.dump( opt.vmeta , 
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
	
	IterationReport report = g.vardb.iterate( g_view , &opt , m );	

	Pseq::finished();
      }
    

    if ( command == "gs-view" )
      {
	Opt_geneseq opt;
	if ( args.has( "ref-variants" ) ) 
	    opt.ref = g.refdb.lookup_group_id( args.as_string( "ref-variants" ) );
	opt.pheno = g.phmap.type() == PHE_DICHOT;
	if ( ! g.seqdb.attached() ) Helper::halt( "no SEQDB attached" );
	if ( ! g.locdb.attached() ) Helper::halt( "no LOCDB attached" );
	IterationReport report = g.vardb.iterate( g_geneseq , &opt , m );
	Pseq::finished();
      }


    if ( command == "i-view" )
      {	
	if ( g.single_file_mode() )
	  Pseq::VarDB::header_VCF( false , true , m );
	else if ( args.has( "from-vardb" ) )
	  Pseq::VarDB::dump_indiv();
	else
	  Pseq::IndDB::dump_table( m );
	Pseq::finished();
      }

    //
    // Dump individual segments from a SEGDB
    //


    if ( command == "seg-view" )
      {
	if ( ! args.has( "group" ) ) 
	  Helper::halt("requires a locus-group to be specified");
	std::string grp = args.as_string( "group" );
	for (int i=0;i<g.indmap.size(); i++)
	  {
	    Individual * person = g.indmap(i);
	    std::cout << *person << "\n";
	    std::set<Region> s = g.segdb.get_indiv_regions( grp , person->id() );
	    std::set<Region>::iterator si = s.begin();
	    while ( si != s.end() ) 
	      {
		std::cout << "\t" << *si << "\n";
		++si;
	      }
	  }
	Pseq::finished();
      }
    

    //
    // Breakdown of { 0:1, 1:0 } etc
    //

    if ( command == "v-dist" )
      {	
	long int nrep = args.has( "perm" ) ? args.as_int( "perm" ) : 1000 ;
	Pseq::VarDB::vdist_summary( m , nrep );
	Pseq::finished();
      }


    //
    // Additional QC measures
    //

    if ( command == "v-freq" )
      {
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

	plog.data_reset();	
	plog.data_group_header( "VAR" );

	plog.data_header( "CHR" );
	plog.data_header( "POS" );
	plog.data_header( "REF" );
	plog.data_header( "ALT" );

	plog.data_header( "FILTER" );
	plog.data_header( "QUAL" );
	plog.data_header( "TI" );  // flag to indicate Ti/Tv

	plog.data_header( "GENO" );
	plog.data_header( "MAF" );
	plog.data_header( "REFMIN" );
	plog.data_header( "HWE" );
	plog.data_header( "HET" );
		
	if ( x.em_stats ) 
	  {
	    plog.data_header( "MPP" );
	    plog.data_header( "H" );
	    plog.data_header( "HALT" );
	    plog.data_header( "EMMAF" ); 
	  }

	plog.data_header( "NSNP" );	

	plog.data_header_done();

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
	std::string dbname = Pseq::Util::single_argument<std::string>( args , "name" );
	Pseq::PPH2DB::score( m , dbname );
	Pseq::finished();
      }


    //
    // Clusters of variants
    //

    if ( command == "clusters" )
      {
	Pseq::VarDB::cluster_scan( m );
	Pseq::finished();
      }


    //
    // Proximity scan, for clusters of variants with LD calculation also
    //
    
    if ( command == "proximity-scan" )
      {
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
	    std::string proj = Pseq::Util::single_argument<std::string>( args , "name" );
	    Pseq::VarDB::make_counts_file( m , proj );
	  }
	else
	  {	 	    
	    bool genotypes = command == "g-counts" ;
	    Pseq::VarDB::simple_counts( m , genotypes );
	  }
	Pseq::finished();
      }


    //
    // Lookup DB information for a list of positions
    //

    if ( command == "lookup" )
      {
	
	if ( args.has( "file" ) )
	  {
	    std::string filename = Pseq::Util::single_argument<std::string>( args , "file" );
	    Pseq::VarDB::lookup_list( filename , m );
	  }
	
	if ( args.has( "region" ) ) 
	  {
	    std::vector<std::string> regions = Pseq::Util::n_arguments<std::string>( args, "region" );
	    std::vector<Region> regs;
	    for (int i=0;i<regions.size(); i++)
	      {
		bool okay;
		Region r( regions[i] , okay );
		if ( okay ) regs.push_back( r );
		else plog.warn( "could not parse region: " + regions[i] );
	      }
	    Pseq::VarDB::lookup_list( "." , m , &regs );
	  }

	// else directly annotate variants in VARDB
	if ( ! ( args.has( "file" ) || args.has( "region" ) ) )
	  {
	    Pseq::VarDB::lookup_list( "." , m );
	  }
	
	Pseq::finished();
      }
    

    //
    // Database summaries
    //

    if ( command == "summary" )
      {
	Pseq::VarDB::summary(m);
	Pseq::IndDB::summary();
	Pseq::LocDB::summary(&g.locdb);
	Pseq::LocDB::summary(&g.segdb);
	Pseq::RefDB::summary();
	Pseq::SeqDB::summary();	
	Pseq::Util::file_summary();
	Pseq::Util::meta_summary();
	Pseq::finished();
      }

    
    if ( command == "var-summary" )
      {
	Pseq::VarDB::summary(m);
	Pseq::finished();
      }

    if ( command == "ind-summary" )
      {
	Pseq::IndDB::summary();
	Pseq::finished();
      }

    if ( command == "loc-summary" )
      {
	Pseq::LocDB::summary( &g.locdb );
	Pseq::finished();
      }

    if ( command == "seg-summary" )
      {
	Pseq::LocDB::summary( &g.segdb );
	Pseq::finished();
      }

    if ( command == "ref-summary" )
      {
	Pseq::RefDB::summary();
	Pseq::finished();
      }

    if ( command == "seq-summary" )
      {
	Pseq::SeqDB::summary();
	Pseq::finished();
      }

    if ( command == "file-summary" )
      {
	Pseq::Util::file_summary();
	Pseq::finished();
      }

    if ( command == "meta-summary" )
      {
	Pseq::Util::meta_summary();
	Pseq::finished();
      }


    //
    // Calculated summary statistics
    //

    if ( command == "v-stats" )
      {	

	  // Hard-code some of these values for now
	  
	  Pseq::VStat vstat(&g);
	  
	  Pseq::Util::set_default( vstat );
	  
	  IterationReport report = g.vardb.iterate( f_vstat , &vstat , m );	
	  
	  vstat.report();
	  
	  Pseq::finished();
      }
    
    if ( command == "g-stats" )
      {

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
	Pseq::IStat istat(&g);
	IterationReport report = g.vardb.iterate( f_istat , &istat , m );	
	istat.report();
	Pseq::finished();
      }



    //
    // Per-locus simple view, or sequence stats, e.g. GC percent 
    //

    if ( command == "loc-view" )
    {
	if ( ! ( g.locdb.attached() | g.segdb.attached() ) ) Helper::halt("no LOCDB or SEGDB attached");
	if ( ! args.has( "group" ) ) Helper::halt("requires a locus-group to be specified");
	
	std::vector<std::string> grp = args.as_string_vector( "group" );
	if ( grp.size() != 1 ) Helper::halt("requires a single locus-group to be specified");	
	
	std::vector<std::string> alias;
	if ( args.has( "alias" ) ) alias = args.as_string_vector( "alias" );

	Pseq::LocDB::loc_view( grp[0] , alias , ! args.has("no-meta") , args.has("show-subregions") );
	
	Pseq::finished();
      }


    if ( command == "loc-intersect" )
      {

	// print LOC group names that have 1+ VARDB variants, according to 
	// the mask
	
	// assumes a mask will specifiy the 'loc.group'

	g.vardb.iterate( g_loc_view , NULL , m );
	Pseq::finished();
      }


    if ( command == "loc-stats" )
      {	
	if ( ! ( g.locdb.attached() | g.segdb.attached() ) ) Helper::halt("no LOCDB or SEGDB attached");
	
	if ( ! g.seqdb.attached() ) Helper::halt("no SEQDB attached");
	
	if ( ! args.has( "loc-group" ) ) 
	  Helper::halt("requires a locus-group to be specified, with loc-group");
	
	std::vector<std::string> grp = args.as_string_vector( "loc-group" );
	
	for (int i=0; i<grp.size(); i++) 
	  Pseq::SeqDB::loc_stats( grp[i] , args.as_string("ref-group") );
	
	Pseq::finished();
      }
    
    //
    // View AA sequence and details
    //

    if ( command == "loc-translate" )
      {	
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
	if ( ! g.refdb.attached() ) Helper::halt("no REFDB attached");
	if ( ! args.has( "group" ) ) Helper::halt("no group specified");
	std::string grp = Pseq::Util::single_argument<std::string>( args , "group" );
	bool with_meta = args.has( "vmeta" );
	g.refdb.dump( grp , with_meta );
	Pseq::finished();
      }


    //
    // Dump sequence from SEQDB
    //

    if ( command == "seq-view" )
      {	
	
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
    // Per indiviaul/group unique/enriched listing
    //
    
    if ( command == "unique" )
      {

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
	Pseq::VarDB::check_concordance(m);
	Pseq::finished();
      }

    //
    // IBS sharing matrix
    //

    if ( command == "ibs-matrix" )
      {
	Pseq::IBS::calculate(m);
	Pseq::finished();
      }

    // 
    // Association tests (with group or phenotype) 
    //

    if ( command == "group-comparison" )
      {
	Pseq::Assoc::group_comparison( m );
      }

    if ( command == "assoc" )
      {
	if ( g.phmap.type() != PHE_DICHOT ) 
	  Helper::halt("no dichotomous phenotype specified");

	// if no perms specified, use adaptive permutation mode
	Pseq::Assoc::set_assoc_test( m , args );
	
      }
    
    if ( command == "net-assoc" )
      {
	if ( g.phmap.type() != PHE_DICHOT ) 
	  Helper::halt("no dichotomous phenotype specified");
	Pseq::Assoc::net_assoc_test( m , args );
      }

    if ( command == "s-assoc" )
      {
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
	if ( g.phmap.type() != PHE_DICHOT ) 
	  Helper::halt("no dichotomous phenotype specified");
	if ( !args.has("ibddb") ) Helper::halt("need to specify --ibddb");
	Pseq::IBD::sharing_wrapper( args.as_string( "ibddb" ) , m );
	Pseq::finished();
      }

    
    //
    // Single-site association statistics
    // 
    
    if ( command == "v-assoc" )
      {

	Pseq::Assoc::Aux_vassoc_options aux;
	if ( args.has("info") ) aux.show_istat = true;
	if ( args.has("separate-chr-bp") ) aux.separate_chr_bp = true;
	if ( args.has("vmeta") ) aux.show_meta = true;
	aux.nrep = args.has( "perm" ) ? args.as_int( "perm" ) : 0 ;

	if ( g.phmap.type() != PHE_DICHOT ) 
	  Helper::halt("no dichotomous phenotype specified");	

	Pseq::Assoc::variant_assoc_test( m , aux , args );

	Pseq::finished();
      }


    //
    // GLM single-variant tests
    // 
    

    if ( command == "glm" )
      {
	
	if ( g.phmap.type() != PHE_DICHOT && g.phmap.type() != PHE_QT ) 
	  Helper::halt("no dichotomous or quantitative phenotype specified");	

	Pseq::Assoc::Aux_glm aux;
	
	aux.show_meta = args.has( "vmeta" );
	
	aux.dichot_pheno = g.phmap.type() == PHE_DICHOT;
	
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
	    
	    if ( args.has( "show-covar" ) )
	    {
		aux.show_all_covar = true;
	      }
	}
	
	// for now, no permutations
	// aux.nrep = args.has( "perm" ) ? args.as_int( "perm" ) : 0 ;
	
	Pseq::Assoc::glm_assoc_test( m , aux );

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
	Pseq::VarDB::write_VCF(m);
	Pseq::finished();
      }

    
    if ( command == "write-bcf" )
      {	
	if ( ! args.has( "bcf" ) ) 
	  Helper::halt( "need to specify --bcf output.bcf" );
	Pseq::VarDB::write_BCF( m , args.as_string( "bcf" ) );
	Pseq::finished();
      }
    
    
    if ( command == "write-ped" )
      {
	if ( ! args.has("name") )
	  Helper::halt("no output file given, use --name");
	string filename = args.as_string( "name" );
	Pseq::VarDB::write_PED(m,filename, args.has( "family-id" ) );
	Pseq::finished();
      }


    if ( command == "write-lik" )
      {
	Pseq::VarDB::write_lik(m);
	Pseq::finished();
      }    
    

    if ( command == "v-matrix" )
      {	
	Pseq::VarDB::write_matrix(m);
	Pseq::finished();
      }


    if ( command == "meta-matrix" )
      {	
	Pseq::VarDB::write_meta_matrix(m);
	Pseq::finished();
      }

    if ( command == "v-meta-matrix" )
      {	
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
	std::vector<std::string> grps = Pseq::Util::n_arguments<std::string>( args , "group" );	
	if ( grps.size() != 2 ) Helper::halt( "exactly two groups need to be specified" );
	Pseq::LocDB::overlap_analysis( grps[0] , grps[1] );
	Pseq::finished();
      }

    if ( command == "g-matrix" )
      {	
	// options:
	// 1) only show genes with non-zero variance

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
	OptGMetaMatrix opt;
	opt.name = Pseq::Util::single_argument<std::string>( args, "name" );
	opt.show_mean = true;
	Pseq::VarDB::write_gene_meta_matrix(m,opt);
	Pseq::finished();
      }



    Pseq::finished();
    return 0;
    
}

