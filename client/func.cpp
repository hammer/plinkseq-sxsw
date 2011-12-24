
#include "func.h"
#include "views.h"
#include "summaries.h"
#include "util.h"

#include "pseq.h"

#include <vector>
#include <fstream>

using namespace std;

extern GStore g;
extern Pseq::Util::Options args;

void Pseq::finished()
{  
  if ( g.gseq_mode() )
    {
      std::ofstream O1( g.gseq_history().c_str() , std::ios::out | std::ios::app );
      O1 << "_STATUS" << "\t"
	 << g.gseq_job() << "\t"
	 << "Done" << "\n";
      O1.close();
    } 
  exit(0);
}

bool Pseq::new_project( std::string project , Pseq::Util::Options & args )
{
  
  // Ensure full paths are used
  
  project = Helper::fullpath( project );

  // write a new project file
  
  std::ofstream P( project.c_str() , std::ios::out );
  
  if ( P.bad() ) Helper::halt("problem writing new project file " + project );
  
  // Search the options for the following:
  
  // --output ( default project_out/
  // --vcf 1+ files

  // --vardb  default {output}/vardb
  // --inddb  default {output}/locdb
  // --locdb  default {output}/locdb
  // --refdb  default {output}/locdb
  // --seqdb  default {output}/locdb
  

  //
  // Core folders, for output and input
  //

  std::string output_folder = args.has("output") ?
    args.as_string( "output" ) : 
    project + "_out/";
  
  output_folder = Helper::fullpath( output_folder );
  Helper::ensure_folder( output_folder );
  
  P << output_folder << "\tOUTPUT\n";
  
  //
  // Resources folder
  //

  std::string resources_folder = args.has("resources") ?
    args.as_string( "resources" ) :
    project + "_res/";

  resources_folder = Helper::fullpath( resources_folder );
  Helper::ensure_folder( resources_folder );
  
  P << resources_folder << "\tRESOURCES\n";
  
  //
  // Scratch files
  //

  if ( args.has( "scratch" ) )
    {
      std::string scratch_folder = args.as_string( "scratch" );
      scratch_folder = Helper::fullpath( scratch_folder );
      Helper::ensure_folder( scratch_folder );
      P << scratch_folder << "\tTEMP\n";
    }


  //
  // Meta-information on meta-information
  //
  
  if ( args.has( "metameta" ) )
    {
      std::string metameta = args.as_string( "metameta" );
      metameta = Helper::fullpath( metameta );
      P << metameta << "\tMETAMETA\n";
    }


  //
  // VCF files
  //
  
  if ( args.has("vcf") )
    {
      std::vector<std::string> t = args.as_string_vector( "vcf" );
      for (int i=0; i<t.size(); i++)
	{
	  t[i] = Helper::fullpath( t[i] );
	  if ( ! Helper::fileExists( t[i] ) )
	    plog.warn( "VCF file " + t[i] + " does not exist" );	  
	  P << t[i] << "\tVCF\n";      
	}
    }
  
  //
  // Databases
  //
  
  if ( args.has( "vardb" ) ) 
    P << Helper::fullpath( args.as_string( "vardb" ) ) << "\tVARDB\n";
  else
    P << output_folder << "vardb\tVARDB\n";
  
  if ( args.has( "inddb" ) ) 
    P << Helper::fullpath( args.as_string( "inddb" ) )<< "\tINDDB\n";
  else
    P << output_folder << "inddb\tINDDB\n";
  
  if ( args.has( "locdb" ) ) 
    P << Helper::fullpath( args.as_string( "locdb" ) ) << "\tLOCDB\n";
  else
    P << resources_folder << "locdb\tLOCDB\n";
  
  if ( args.has( "refdb" ) ) 
    P << Helper::fullpath( args.as_string( "refdb" ) ) << "\tREFDB\n";
  else
    P << resources_folder << "refdb\tREFDB\n";
  
  if ( args.has("seqdb") )
    P << Helper::fullpath( args.as_string( "seqdb" ) ) << "\tSEQDB\n";
  else
    P << resources_folder << "seqdb\tSEQDB\n";
  
  P.close();
  
}

bool Pseq::set_project( std::string project )
{
  return g.set_project( project );
}
  
bool Pseq::VarDB::summary( Mask & m , bool ugly )
{
  if ( ! g.vardb.attached() ) return false;
  plog << g.vardb.summary(&m,ugly ) << "\n";
  return true;
}

bool Pseq::VarDB::attach( std::string db )
{
  return g.vardb.attach(db);
}


bool Pseq::VarDB::flush( const std::vector<int> & f )
{
  for (int i=0; i<f.size(); i++)
    g.vardb.drop(f[i]);    
  return true;
}

bool Pseq::VarDB::vacuum()
{
  g.vardb.vacuum();
  return true;
}


void f_vcf( Variant & v , void * p)
{
   plog << v.VCF();
}

void f_vcfz( Variant & v , void * p )
{
  VCFZ * vcfz = (VCFZ*)p;
  vcfz->write_record( v );
}

bool Pseq::VarDB::write_VCF(Mask & m , bool compressed )
{

  // Either write a VCF to stream, or a compressed BGZF-zipped VCF to a file

  if ( compressed ) 
    {

      std::string name = Pseq::Util::single_argument<std::string>( args, "file" );

      VCFZ vcfz( name , &g.vardb);
      vcfz.writing();
      vcfz.open();
      vcfz.write_header();

      IterationReport report = g.vardb.iterate( f_vcfz , &vcfz , m );
	
      vcfz.close();

      return true;
    }
  

  //
  // Otherwise, just return the VCF to STDOUT as a normal textual VCF file
  //
  
  // VCF headers
  
  plog << "##fileformat=" << PLINKSeq::CURRENT_VCF_VERSION() << "\n"
       << "##source=pseq\n"
       << MetaInformation<VarMeta>::headers( )
       << MetaInformation<GenMeta>::headers( META_GROUP_GEN )
       << MetaInformation<VarFilterMeta>::headers( META_GROUP_FILTER );
  
  // Header line containing Individual IDs
  
  const int n = g.indmap.size();
  plog << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";  
  for ( int i=0; i<n; i++) plog << "\t" << g.indmap(i)->id();
  plog << "\n";
  
  // Variants
  
  IterationReport report = g.vardb.iterate( f_vcf , NULL , m );
  
  

  
}


bool Pseq::VarDB::write_PED(Mask & m, string filename , bool use_family_id )
{

  // Need to write TPED and TFAM files
  
  if ( ! g.inddb.attached() ) return false;
  if ( ! g.vardb.attached() ) return false;

  const int n = g.indmap.size();
  
  ofstream TFAM( ( filename + ".tfam").c_str(), std::ios::out );
  
  for (int i = 0 ; i < n ; i++) 
    {      

      Individual * pmapped = g.indmap.ind( i );      

      if ( use_family_id ) 
	TFAM << pmapped->fid() << "\t"
	     << pmapped->iid() << "\t";
      else
	TFAM << pmapped->id() << "\t"
	     << "." << "\t";

      if ( pmapped->father() == "." )
	TFAM << "0\t";
      else
	TFAM << pmapped->father() << "\t";
      
      if ( pmapped->mother() == "." )
	TFAM << "0\t";
      else
	TFAM << pmapped->mother() << "\t";

      // 1=male, 2=female, 0=unknown coding for sex
      TFAM << pmapped->sex() << "\t";
      
      // Phenotype
      
      if ( g.phmap.type() == PHE_DICHOT )
	{
	  if ( pmapped->affected() == CASE ) TFAM << "2\n";
	  else if ( pmapped->affected() == CONTROL ) TFAM << "1\n";
	  else TFAM << "0\n";
	}
      else if ( g.phmap.type() == PHE_QT )
	TFAM << pmapped->qt() << "\n";
      else if ( g.phmap.type() == PHE_FACTOR )
	TFAM << pmapped->group_label() << "\n";
      else
	TFAM << "0\n";
    }
 
  TFAM.close();
 
  // Now MAP and genotype data
  
  std::ofstream TPED( ( filename + ".tped" ).c_str() , std::ios::out );
  
  IterationReport report = g.vardb.iterate( f_view_tped , &TPED , m );
  
  TPED.close();

  return true;
}

bool Pseq::VarDB::dump_indiv()
{
  int cnt = 0;
  if ( ! g.vardb.attached() ) return false;
  std::map<int,std::string> files = g.vardb.fetch_files();
  std::map<int,std::string>::iterator i = files.begin();
  while ( i != files.end() )
    {
      std::vector<std::string> inds = g.vardb.fetch_individuals( i->first );
      for (int j = 0 ; j < inds.size(); j++)
	{
	  plog << ++cnt << "\t" 
	       << j+1 << "\t"
	       << inds[j] << "\t"
	       << g.vardb.file_tag( i->first ) << "\t"
	       << i->second << "\n";
	}
      ++i;
    }
  return true;
}

bool Pseq::VarDB::write_lik(Mask & m)
{
  std::vector<std::string> ids = g.indmap.ind_id();
  plog << "VAR\tREF\tALT";
  for (int i=0; i<ids.size(); i++) 
    plog << "\t" << ids[i] 
	 << "\t" << ids[i] 
	 << "\t" << ids[i] ;
  plog << "\n";
  
  // Variants

  IterationReport report = g.vardb.iterate( f_view_lik , &g , m );

  return true;

}

bool Pseq::VarDB::write_matrix(Mask & m)
{
  std::vector<std::string> ids = g.indmap.ind_id();
  plog << "VAR\tREF\tALT";
  for (int i=0; i<ids.size(); i++) 
    plog << "\t" << ids[i] ; 
  plog << "\n";
  
  IterationReport report = g.vardb.iterate( f_view_matrix , &g , m );
  return true;
}

bool Pseq::VarDB::write_meta_matrix( Mask & m )
{
  
  std::map<int,string> files = g.vardb.fetch_files( &m );    
  
  // Display consensus meta-info, if >1 file
  
  plog << "CHR\tPOS\tID";

  // Variant level information

  int nelem_static = MetaInformation<VarMeta>::n_visible_keys_static();
  int nelem_nonstatic = MetaInformation<VarMeta>::n_visible_keys_nonstatic();

  if ( nelem_static ) plog << "\t" << MetaInformation<VarMeta>::display_header_static();
  
  // Consensus sample-variant
  if ( nelem_nonstatic ) plog << "\t" << MetaInformation<VarMeta>::display_header_nonstatic();
  
  if ( files.size() > 1 ) 
    {      
      std::map<int,std::string>::iterator i = files.begin();
      while ( i != files.end() ) 
	{ 
	  std::string m = MetaInformation<VarMeta>::display_header_nonstatic( "F"+Helper::int2str(i->first)+"_" ) ;
	  std::string f = MetaInformation<VarFilterMeta>::display_header( "F"+Helper::int2str(i->first)+"_" );
	  if ( m != "" ) plog << "\t" << m;
	  if ( f != "" ) plog << "\t" << f;
	  ++i;
	}  
    }
  else
    {      
      std::string f = MetaInformation<VarFilterMeta>::display_header( );
      if ( f != "" ) plog << "\t" << f;
    }
  plog << "\n";  
  
  IterationReport report = g.vardb.iterate( f_view_meta_matrix , &files , m );
  
  return true;
}


bool Pseq::VarDB::write_var_meta_matrix(Mask & m, std::string & name)
{
  std::vector<std::string> ids = g.indmap.ind_id();
  plog << "MARKER";
  for (int i=0; i<ids.size(); i++) 
    plog << "\t" << ids[i] ; 
  plog << "\n";
  IterationReport report = g.vardb.iterate( f_view_var_meta_matrix , &name , m );
  return true;
}


bool Pseq::VarDB::write_gene_matrix(Mask & m, OptGMatrix & opt)
{
  std::vector<std::string> ids = g.indmap.ind_id();
  plog << "GENE\tNV";
  for (int i=0; i<ids.size(); i++) 
    plog << "\t" << ids[i] ; 
  plog << "\n";
  IterationReport report = g.vardb.iterate( f_view_gene_matrix , &opt , m );
  return true;
}


bool Pseq::VarDB::write_gene_meta_matrix(Mask & m, OptGMetaMatrix & opt)
{
  std::vector<std::string> ids = g.indmap.ind_id();
  plog << "GENE\tNV";
  for (int i=0; i<ids.size(); i++) 
    plog << "\t" << ids[i] ; 
  plog << "\n";
  IterationReport report = g.vardb.iterate( f_view_gene_meta_matrix , &opt , m );
  return true;
}


bool Pseq::VarDB::consolidate(Mask & m, std::string label)
{
  
}


void f_uniq_report( Variant & v , void * p )
{
  
  OptUniq * opt = (OptUniq*)p;
  
  int a = 0;
  int c = 0;
  bool altmin = v.n_minor_allele( &a , &c );
  
  //  std::cout << opt->ingroup_req << " " << opt->outgroup_allow << "\n";
 
  if ( a == 0 || a == c ) return;
  
  int obs_group = 0;
  int obs_outgroup = 0;

  // Is this variant minor allele uniq to this/these people? 

  for (int i=0; i<v.size(); i++)
    {
      if ( ! v(i).null() )
	{
	  // Seen in non-group?
	  if ( opt->indiv.find( v.ind(i) ) == opt->indiv.end() )
	    {
	      if ( v(i).minor_allele( altmin ) > 0 )
		{
		  if ( ++obs_outgroup > opt->outgroup_allow ) return;
		}
	    }
	  else
	    {
	      if ( v(i).minor_allele( altmin ) ) obs_group++;	      
	    }
	}
      
    }


  // Did we see in all group members?
  
  if ( opt->ingroup_req == -1 )
    {
      if ( obs_group == opt->indiv.size() )
	plog  << obs_group << "\t"
	      << obs_outgroup << "\t"
	      << v << "\t" << v.meta << "\n";
    }
  else
    {     
      if ( obs_group >= opt->ingroup_req ) 
 	plog << obs_group << "\t"
	     << obs_outgroup << "\t"
	     << v << "\t" << v.meta << "\n";
    } 


//   if ( obs_group > 0 )
//     plog << "DET\t"
// 	      << obs_group << "\t"
// 	      << obs_outgroup << "\t"
// 	      << v << "\t"
// 	      << v.meta << "\n";

}


bool Pseq::VarDB::uniq_report( std::vector<std::string> & indiv , 
			       Mask & m , 
			       OptUniq & opt )
{
  
  for (int i = 0 ; i < indiv.size(); i++) 
    {
      std::vector< std::string> indiv2 = Helper::parse( indiv[i] , "," );
      for (int j = 0 ; j < indiv2.size(); j++)
	{
	  Individual * person = g.indmap.ind( indiv2[j] );
	  if ( person ) 
	    opt.indiv.insert( person ) ;
	  else 
	    plog.warn( "coud not find individual " + indiv2[j] );
	}
    }
  
  IterationReport rep = g.vardb.iterate( f_uniq_report , &opt , m );
  
  return true;
  
}
  
    
bool Pseq::VarDB::variant_stats(Mask & m)
{

}
	    


void f_write_into_vardb( Variant & var , void * p )
{
  
  VarDBase * db = (VarDBase*)p;
  
  //
  // Insert this into a single file in the new database
  //

  db->insert_consensus( 1 , var );  
  
}


bool Pseq::VarDB::write_vardb( const std::string & new_project , 
			       const std::string & new_vardb , 
			       Mask & mask )
{

  //
  // Create the new VARDB. Using a Mask, place all data into a single
  // new "file", but track the original files in the new header of
  // this dataset
  //

  if ( Helper::fileExists( new_vardb ) )
    Helper::halt( new_vardb + " already exists" );
  
  VarDBase newdb( g.indmap );
  newdb.attach( new_vardb );
  newdb.drop_index();
  newdb.begin();
  
  //
  // Add header and meta-information 
  //
  
  newdb.insert( "newvardb" , "1" );
  newdb.insert_header( 1, "format" , "VARDB" ); 
  newdb.insert_header( 1, "source" , "pseq" ); 
  
  //
  // Add individuals
  //
  
  std::vector<std::string> ids = g.indmap.ind_id();  
  for (int i = 0 ; i < ids.size() ; i++ )
    newdb.insert( 1 , Individual( ids[i] ) );

  //
  // Add meta-information descriptors
  //

  std::map<meta_name_t,meta_index_t> midx_var = MetaInformation<VarMeta>::dump_types();
  std::map<meta_name_t,meta_index_t> midx_flt = MetaInformation<VarFilterMeta>::dump_types();
  std::map<meta_name_t,meta_index_t> midx_gen = MetaInformation<GenMeta>::dump_types();

  std::map<meta_name_t,meta_index_t>::iterator mi = midx_var.begin();
  while ( mi != midx_var.end() )
  {
      meta_index_t & midx = mi->second;
      newdb.insert_metatype( 1 , midx.name , midx.mt , midx.len , META_GROUP_VAR , midx.description );
      ++mi;
  }

  mi = midx_flt.begin();
  while ( mi != midx_flt.end() )
  {
      meta_index_t & midx = mi->second;
      newdb.insert_metatype( 1 , midx.name , midx.mt , midx.len , META_GROUP_FILTER , midx.description );
      ++mi;
  }

  mi = midx_gen.begin();
  while ( mi != midx_gen.end() )
  {
      meta_index_t & midx = mi->second;
      newdb.insert_metatype( 1 , midx.name , midx.mt , midx.len , META_GROUP_GEN , midx.description );
      ++mi;
  }



  
  //
  // Add selected variants (and meta-fields)
  //
  
  g.vardb.iterate( f_write_into_vardb , &newdb , mask );
  
  //
  // All done, finish up and write new project file
  //
  
  newdb.commit();  
  newdb.index();

  // Insert summary Ni, Nv into database for this file  
  int2 niv = newdb.make_summary( "newvardb" ) ;  
  plog << new_vardb << " : inserted " << niv.p2 << " variants\n";
  
  // Update / create new project specification file  
  g.fIndex.addSpecial( VARDB , Helper::fullpath( new_vardb )  );
  g.fIndex.write_new_projectfile( new_project );
  
  return true;

}    
   
bool Pseq::LocDB::summary( LocDBase * db , bool ugly )
{
  if ( ! db ) return false;
  if ( ! db->attached() ) return false;
  plog << db->summary( ugly ) << "\n";
  return true;
}


bool Pseq::LocDB::attach( std::string db )
{
  return g.locdb.attach(db);
}

bool Pseq::LocDB::load_GTF( std::string file , std::string label , bool locdb )
{

  // use transcript_id field by default to index genes, unless 'gene-id' option
  // is given

  if ( locdb ) 
    return g.locdb.load_GTF( file , label , ! args.has("use-gene-id") ) != 0;
  else
    return g.segdb.load_GTF( file , label , ! args.has("use-gene-id") ) != 0;
}



bool Pseq::LocDB::merge( std::string label1 , std::string label2 , bool locdb )
{
  if (locdb ) 
    return g.locdb.merge( label1, label2 ) != 0;
  else 
    return g.segdb.merge( label1, label2 ) != 0;
}


bool Pseq::LocDB::index()
{  
  if ( ! g.locdb.attached() ) return false;
  g.locdb.index();
  return true;
}
      

bool Pseq::LocDB::swap_alternate_names( const std::string & group , const std::string & filename )
{
  
  int gid = g.locdb.lookup_group_id( group );
  
  if ( gid == 0 ) 
    {
      plog.warn("group not found in LOCDB");
      return false;
    }

  Helper::checkFileExists( filename );
  
  InFile infile( filename );

  int cnt = 0;

  while ( ! infile.eof() )
    {
      // expect two fields per line, tab-delimited
      std::vector<std::string> tok = infile.tokenizeLine( "\t" );
      
      if ( tok.size() == 0 ) continue;
      
      if ( tok.size() != 2 ) 
	{
	  plog.warn("bad format of gene-name file (not 2 tab-delimited fields) ");
	  continue;
	}

      g.locdb.replace_real_names( gid , tok[0] , tok[1] , args.has( "alternate-name" ) );

      if ( ++cnt % 1000 == 0 ) 
	plog.counter( "replaced " + Helper::int2str( cnt ) + " locus names" );
    }
  
  plog.counter("\n");
  
  infile.close();
 
  return true;
}

bool Pseq::LocDB::overlap_analysis( )
{

  std::vector<std::string> grps = Pseq::Util::n_arguments<std::string>( args , "group" );	
  if ( grps.size() != 2 ) Helper::halt( "exactly two groups need to be specified" );
  if ( grps[0] == grps[1] ) Helper::halt( "not going to calculate overlap of group with self" );

  std::string a = "";
  if ( args.has("alias") ) a = args.as_string( "alias" ); 
  
  std::string ol = "";
  if ( args.has("output","tab")) ol="tab";
  else if ( args.has("output","comma")) ol="comma";
  else if ( args.has("output","row")) ol="row";

  if ( ! g.locdb.attached() ) Helper::halt( "no LOCDB attached" );

  g.locdb_overlap_analysis( grps[0] , grps[1] , a , ol );
  return true;
}


bool Pseq::LocDB::intersection( std::string filename , std::string group , LocDBase & db )
{

    if ( ! db.attached() ) return false;

    if ( ! Helper::fileExists( filename ) ) 
	Helper::halt("could not find " + filename );
    
    if ( db.lookup_group_id( group ) == 0 )
      Helper::halt("could not find group " + group );

  //
  // Read a vector of regions
  //

  std::vector<Region> regions;
  
  std::ifstream IN1( filename.c_str() , std::ios::in );

  int inv = 0, blank = 0;

  while ( ! IN1.eof() )
    {

      std::vector<std::string> tok = Helper::tokenizeLine( IN1 );

      std::string region ;

      const int sz = tok.size();

      if ( sz == 1 ) 
	region = tok[0];
      else if ( sz == 2 ) 
	region = tok[0] + ":" + tok[1] + ".." + tok[1] ;
      else if ( sz == 3 )       
	region = tok[0] + ":" + tok[1] + ".." + tok[2] ;
      else
	continue;
      
      bool valid = true;
      
      Region r(region,valid);

      // Get overlap
      
      if ( valid )
	{
	  
	  set<Region> olaps = db.get_regions( group ,  r ) ;	  

	  plog << r.coordinate() << "\t";
	  
	  if ( olaps.size() == 0 ) ++blank;
	  set<Region>::iterator i = olaps.begin();
	  set<string> genes;
	  while ( i != olaps.end() )
	    {
	      genes.insert( i->altname );
	      ++i;
	    }	  
	  
	  
	  if ( genes.size() == 0 ) 
	    plog << ".";
	  else
	    {
	      set<string>::iterator j = genes.begin();
	      while ( j != genes.end() )
		{
		  if ( j != genes.begin() ) plog << "|";
		  plog << *j;
		  ++j;
		}
	    }
	  
	  plog << "\n";
	  
	  
	}
      else
	++inv;
    }
  
  plog << "read " << regions.size() << " regions\n";
  
  if( inv > 0 ) 
    plog.warn( "found " + Helper::int2str( inv ) +  " invalid regions" );
  
  if( blank > 0 ) 
    plog.warn( "found " + Helper::int2str( blank ) + " empty regions" );

  IN1.close();
  
  return true;
}


bool Pseq::LocDB::load_pathway( std::string file , std::string label , std::string group , bool locdb , bool use_altname )
{
  if ( locdb )
    return g.locdb.load_set( file , label, group , use_altname );
  else
    return g.segdb.load_set( file , label, group , use_altname );
}


bool Pseq::RefDB::summary( bool ugly )
{
  if ( ! g.refdb.attached() ) return false;
  plog << g.refdb.summary( ugly ) << "\n";
  return true;
}

bool Pseq::RefDB::attach( std::string db )
{
  g.refdb.attach( db );
}
      

bool Pseq::SeqDB::load_transcripts( std::string label )
{
  return Annotate::load_transcripts( LOCDB , label );
}

  
bool Pseq::IndDB::attach( std::string db )
{
  return g.inddb.attach(db);
}

bool Pseq::IndDB::summary( bool ugly )
{
  if ( ! g.inddb.attached() ) return false;
  plog << g.inddb.summary( ugly ) << "\n";
  return true;
}

bool Pseq::IndDB::set_phenotype( std::string label )
{
  return g.phmap.set_phenotype( label );  
}

bool Pseq::IndDB::make_phenotype( std::string def )
{
  return g.phmap.make_phenotype( def );
}


bool Pseq::IndDB::load_ped_info( std::string file )
{
  return g.inddb.load_ped_info( file );
}

bool Pseq::IndDB::load_phenotypes( std::string file )
{
  return g.inddb.load_phenotypes( file );  
}


bool Pseq::IndDB::dump_table(Mask & m)
{


  if ( ! g.inddb.attached() ) 
    {
      plog.warn("no INDDB attached");
      return false;
    }
  

  // List all information in IndDB, one person per line
  
  std::map<std::string,std::vector<std::string> > pheno = g.inddb.fetch_phenotype_info();
  std::map<std::string,std::vector<std::string> >::iterator p = pheno.begin();
  
  while ( p != pheno.end() )
    {
      plog << "#" << p->first << " "
	   << "(" << p->second[0]  << ") "
	   << p->second[1] << "\n";
      ++p;
    }

  // Explicitly note phenotype/strata
  
  plog << "#PHE\t" << ( g.phmap.phenotype_set() ? g.phmap.phenotype() : "." ) << "\n";
  plog << "#STRATA\t" << ( g.phmap.strata_set() ? g.phmap.strata() : "." ) << "\n";


  //
  // Get phenotype information
  //
  

  //
  // Select who to look at. The idea is this will create a dummy
  // person on the fly, if they do not exist; they will have a missing
  // phenotype, etc, so they will be appropriately handled downstream,
  // but this way we never encounter a NULL pointer Have to decide: if
  // we create them, then make sure all get deleted at end...  
  //
  
  const int n = g.indmap.size();
  
  // track any people in the phmap who are not in the mask 

  plog << "#ID\tFID\tIID\tMISS\tSEX\tPAT\tMAT\tMETA";
  if ( g.phmap.phenotype_set() && g.phmap.type() == PHE_DICHOT ) plog << "\tPHE";
  if ( g.phmap.strata_set() ) plog << "\tSTRATA";
  plog << "\n";

  std::set<Individual*> observed;
  
  for (int i = 0 ; i < n ; i++) 
    {

      // note: should never return a null, if based on 0..(n-1) where n is 
      // returned from create_from_mask()
      
      Individual * pmapped = g.indmap.ind( i );
            
      plog << *pmapped;
      
      if ( g.phmap.type() == PHE_DICHOT )
	{
	  if ( pmapped->affected() == CASE ) plog << "\tCASE";
	  else if ( pmapped->affected() == CONTROL ) plog << "\tCONTROL";
	  else plog << "\t.";
	}
      
      if ( g.phmap.strata_set() )
	{
	  plog << "\t" << pmapped->group_label();
	}
      
      
      plog << "\n";      

      observed.insert( pmapped );
      
    }


  //
  // List the people we are not looking at
  //

  // note -- currently, this function only lists the people who are in VARDB, and 
  // specified by the Mask.  If they also exist in INDDB, then the phenotypic info. 
  // is attacged.  Currently, there is no option to get the list of people in INDDB
  // who are not in VARDB

//   string ex = "";
//   int x = 0;
//   set<Individual*> all = g.phmap.all_in_map();
//   set<Individual*>::iterator i = all.begin();
//   while ( i != all.end() )
//     {
//       if ( observed.find( *i ) == observed.end() )
// 	{
// 	  ex += " " + (*i)->id();
// 	  if ( ++x % 8 ) ex += "\n";
// 	}
//       ++i;
//     }
  
//   if ( x != 0 ) 
//     plog << x << " individuals in INDB not included by this mask:" << ex << "\n";
  
  return true;
}


bool Pseq::SeqDB::summary( bool ugly )
{
  if ( ! g.seqdb.attached() ) return false;
  plog << g.seqdb.summary( ugly ) << "\n";
  return true;
}




Mask Pseq::Util::construct_mask( std::string desc )
{
  return Mask( desc );
}


void Pseq::Util::set_default( VStat & vstat )
{
    
    if ( args.has( "stats" , "hwe" ) )
	vstat.add_hwe_p( args.as_string( "stats" , "hwe" ) ); // implicitly a comma-string
    
    if ( args.has( "stats" , "ref" ) )
    {
	std::vector<std::string> s = args.as_string_vector( "stats" , "ref" );
	std::vector<std::string>::iterator i = s.begin();
	while ( i != s.end() )
	{
	  vstat.add_refgroup( *i );
	  ++i;
	}
    }

    if ( args.has( "stats" , "loc" ) )
    {
	std::vector<std::string> s = args.as_string_vector( "stats" , "loc" );
	std::vector<std::string>::iterator i = s.begin();
	while ( i != s.end() )
	{
	    vstat.add_locgroup( *i );
	    ++i;
	}
    }
    
    if ( args.has( "stats" , "mac" ) )
    {
	std::vector<std::string> s = args.as_string_vector( "stats" , "mac" );
	std::vector<std::string>::iterator i = s.begin();
	while ( i != s.end() )
	{
	    vstat.add_mac( *i );
	    ++i;
	}
    }
    
    if ( args.has( "stats" , "maf" ) )
    {
	std::vector<std::string> s = args.as_string_vector( "stats" , "maf" );
	std::vector<std::string>::iterator i = s.begin();
	while ( i != s.end() )
	{
	    vstat.add_maf( *i );
	    ++i;
	}
    }
    
    if ( args.has( "stats" , "qual" ) )
    {
	std::vector<std::string> s = args.as_string_vector( "stats" , "qual" );
	std::vector<std::string>::iterator i = s.begin();
	while ( i != s.end() )
	{
	    vstat.add_qual( *i );
	    ++i;
	}
    }
    
    if ( args.has( "stats", "dp" ) )
    {
	std::vector<std::string> s = args.as_string_vector( "stats" , "dp" );
	std::vector<std::string>::iterator i = s.begin();
	while ( i != s.end() )
	{
	    vstat.add_depth( *i );
	    ++i;
	}
    }
    
    if ( args.has( "stats" , "mean" ) )
    {
	std::vector<std::string> s = args.as_string_vector( "stats" , "mean" );
	std::vector<std::string>::iterator i = s.begin();
	while ( i != s.end() )
	{
	    vstat.add_mean_tag( *i );
	    ++i;
	}
    }
    
  // --options count=AN:2:3,DP
  //           count=AN,2:3  
  //           count=AN      //table

    if ( args.has( "stats" , "count" ) )
    {
	std::vector<std::string> s = args.as_string_vector( "stats" , "count" );
	std::vector<std::string>::iterator i = s.begin();
	while ( i != s.end() )
	{
	    std::vector<std::string> tok = Helper::char_split( *i , ':' );
	    if ( tok.size() == 3 )  
		vstat.add_count_tag( tok[0] , tok[1] + ":" + tok[2] );
	    if ( tok.size() == 2 )  
		vstat.add_count_tag( tok[0] , tok[1] );
	    else if ( tok.size() == 1 )  // i.e. list all values
		vstat.add_table_tag( tok[0] );
	    ++i;
	}
    }
    
    bool all_istat_filter = false;
    
    if ( args.has( "stats" , "filter" ) )
    {
	std::vector<std::string> s = args.as_string_vector( "stats" , "filter" );
	std::vector<std::string>::iterator i = s.begin();
	while ( i != s.end() )
	{
	    if ( *i != PLINKSeq::PASS_FILTER() )
	    {
		if ( *i == "." ) all_istat_filter = true;
		else vstat.add_istat_filter( *i );
	    }
	    ++i;
	}
	
    }
    
  // Add all non-PASS filters to the report fields
  std::vector<std::string> f = MetaInformation<VarFilterMeta>::field_names();
  for ( int i = 0 ; i < f.size() ; i++ ) 
    {
      if ( f[i] != PLINKSeq::PASS_FILTER() )
	{
	  vstat.add_filter( f[i] );
	  if ( all_istat_filter ) 
	    vstat.add_istat_filter( f[i] );
	}
    }


  // G-TAGs
  
  if ( args.has( "stats" , "gmean" ) )
  {
      std::vector<std::string> s = args.as_string_vector( "stats" , "gmean" );
      std::vector<std::string>::iterator i = s.begin();
      while ( i != s.end() )
      {
	  vstat.add_mean_gtag( *i );
	  ++i;
      }
    }
  
  
  if ( args.has( "stats" , "gcount" ) )
  {
      std::vector<std::string> s = args.as_string_vector( "stats" , "gcount" );
      std::vector<std::string>::iterator i = s.begin();
      while ( i != s.end() )
      {
	  std::vector<std::string> tok = Helper::char_split( *i , ':' );
	  if ( tok.size() == 3 )  
	      vstat.add_count_gtag( tok[0] , tok[1] + ":" + tok[2] );
	  if ( tok.size() == 2 )  
	      vstat.add_count_gtag( tok[0] , tok[1] );
	  ++i;
	}
  }
  
}


bool Pseq::Util::file_summary( bool ugly )
{
  plog << g.fIndex.summary( ugly ) << "\n";
  return true;
}

bool Pseq::Util::meta_summary( bool ugly )
{
  plog << Helper::metatype_summary( ugly ) << "\n";
  return true;
}


bool f_next_refvar( Region & reg , void * p )
{
  RefVariant rv;
  if ( ! g.refdb.iterate(&rv) ) return false;
  reg = Region( rv );
  return true;  
}


void f_loc_stats_range_iterator_report( const Region & a , const Region & b , void * p ) 
{
  //  plog << "intersect " << a << "\t" << b << "\n";  
  ++(*(int*)p);
} 


bool Pseq::SeqDB::loc_stats( const std::string & grp , const std::string & refgroup )
{

  
  //
  // Pull the whole locus group of interest into memory
  //

  int gid = g.locdb.lookup_group_id( grp );

  if ( gid == 0 ) return false;

  std::set<Region> regions = g.locdb.get_regions( gid ); 
  
  if ( regions.size() == 0 ) return false;


  //
  // Main report. Note -- currently we ignore overlapping subregions... 
  //

  //   number of subregions
  //   % of sites in subregions
  //   GC pct
  //   N (repeat-masked) pct
  //   (optionally) # of ref-variants in a given ref-group (REFDB)
  

  //
  // Optional to display subregions 
  //

  bool sub = args.has("show-subregions");
  bool ignore_sub = args.has("no-subregions");

  //
  // Do we wish to count the instances of intersecting reference-group variants?
  //

  bool refvars = refgroup != "";
  
  if ( refvars )
    {
      if ( ! g.refdb.attached() ) 
	Helper::halt("REFDB not attached");
      if ( ! g.refdb.init_iterate( refgroup ) )
	Helper::halt("could not initiate specified reference group " + refgroup );
    }
  
  RangeIntersector range_intersector( f_next_refvar , f_loc_stats_range_iterator_report );
  

  //
  // Header
  //
  
  plog << "GROUP" << "\t"
       << "ALIAS" << "\t"
       << "NSUB" << "\t"
       << ( sub ? "SUB\t" : "" ) 
       << "POS" << "\t"
       << "L" << "\t"
       << "L_ALL " << "\t"
       << "GC" << "\t"
       << "GC_ALL" << "\t"
       << "N" << "\t"
       << "N_ALL"  
       << ( refvars ? "\tNREF" : "" ) 
       << ( refvars ? "\tNREF_ALL" : "" ) 
       << "\n";
    


  //
  // Start iterating over locus-group selected
  //
  
  std::set<Region>::const_iterator i = regions.begin();

  while ( i != regions.end() )
    {
            
      std::string alias = g.locdb.alias( i->name , true );
      if ( alias == "" ) alias = ".";

      //
      // Also take whole region (must do this first, beause of how
      // RangeIntersector works, i.e. whole transcript < last exon )
      //
      
      
      int gc0 = 0 , n0 = 0 , tot0 = 0;
      int refcnt_whole = 0;
      
      bool okay = g.seqdb.GC( *i , gc0 , tot0 );      
      g.seqdb.N( *i , n0 , tot0 ); 
      
      if ( refvars ) 
	range_intersector.intersect( *i , &refcnt_whole );
      

      //
      // Ignore subregions in calculations (i.e. either if none exist, or told to do so)
      //
      
      bool nosubs = i->subregion.size() == 0 || ignore_sub;

      // # site GC, # repeat-masked (N) and total number
      int gc = 0 , n = 0 , tot = 0;

      // number of refvariants
      int refcnt = 0 ;
      
      // number of valid subregions
      int nreg = 0;
      

      // Consider each subregion 

      for ( int s = 0 ; s < ( nosubs ? 0 : i->subregion.size() ) ; s++ )
	{

	  //	  std::cout << "S = " << s << "\n";

	  int sgc = 0 , sn = 0 , stot = 0;
	  
	  Region sr( i->subregion[s] );
	  
	  bool okay = g.seqdb.GC( sr , sgc, stot );
	  
	  g.seqdb.N( sr , sn, stot) ;
	  
	  if ( okay )
	    {
	      ++nreg;
	      gc += sgc;
	      n += sn;
	      tot += stot;
	    }
	  
	  
	  if ( sub ) 
	    {	      
	      plog << i->name << "\t" 
		   << alias << "\t"
		   << "NA\t"
		   << s+1 << "\t"
		   << sr.coordinate() << "\t";
	      
	      plog << stot << "\t"
		   << "NA" << "\t";
	      
	      if ( okay ) 
		plog 
		  << (double)sgc/(double)(stot-sn) << "\t"
		  << "NA" << "\t"
		  << (double)sn/(double)stot << "\t"
		  << "NA";
	      else
		plog << "NA\tNA\tNA\tNA";
	    }
	  

	  //
	  // Count intersector reference variants
	  //

	  if ( refvars ) 
	    {	      
	      int c = 0; 
	      range_intersector.intersect( sr , &c );
	      refcnt += c;
	      // display?
	      if ( sub ) plog << "\t" << c << "\tNA";
	    }
	  
	  if ( sub ) plog << "\n";
	  
	}

      
	  
      //
      // Print summary of all subregions (S=0 code in output)
      //
      

      plog << i->name << "\t"
	   << alias << "\t";
      plog << ( nosubs ? 1 : nreg ) << "\t";
      if (sub) plog << "0\t";	  
      plog << i->coordinate() << "\t";
      
      if ( nosubs ) plog << ".\t";
      else plog << tot << "\t";
      plog << tot0 << "\t";	  
      
      bool only0 = tot == 0;
      
      if ( okay && ( nosubs || nreg > 0 ) )
	{
	  if ( only0 ) plog << ".\t"; 
	  else
	    {
	      if ( tot-n != 0 )
		plog << (double)gc/(double)(tot-n) << "\t";
	      else
		plog << "NA\t";
	    }
	  
	  if ( tot0-n0 ) 
	    plog << (double)gc0/(double)(tot0-n0) << "\t";
	  else
	    plog << "NA\t";
	  
	  if ( only0 ) plog << ".\t"; else plog << (double)n/(double)tot << "\t";
	  plog << (double)n0/(double)tot0 ;
	}
      else
	plog << "NA\tNA\tNA\tNA";
      
      if ( refvars ) 
	plog << "\t" << refcnt 
	     << "\t" << refcnt_whole ;
      
      plog << "\n";

      //
      // next region
      //

      ++i;

    }

  return true;
}



bool Pseq::PPH2DB::load( const std::string & dbname , const std::string & filename )
{
  PPH2DBase ppdb;
  ppdb.attach( dbname );
  ppdb.set_locdb( &g.locdb );
  ppdb.load( filename );
  return true;
}


void f_pph2_scoring( Variant & v , void * p )
{
  PPH2DBase * pph2 = (PPH2DBase*)p;
  double score = 0;   
  int prediction = 0;
  // only output scores for missense variants found in PPH2
  if ( pph2->score(v,score,prediction) )
    plog << v << "\t" << score << "\t" << prediction << "\n";
  else
    plog << v << "\tNA\t" << prediction << "\n";
}

bool Pseq::PPH2DB::score( Mask & m , const std::string & dbname )
{
  PPH2DBase ppdb;
  ppdb.attach( dbname );
  ppdb.set_locdb( &g.locdb );
  Annotate::load_transcripts( LOCDB , PLINKSeq::DEFAULT_LOC_GROUP() );
  IterationReport report = g.vardb.iterate( f_pph2_scoring , &ppdb , m );
}



struct AuxCountReport 
{
  std::string label;
  bool by_case_control;
  bool unphased;
  bool acdb_format;
};
  


void f_counts_report( Variant & v , void * p )
{  
  
  // one line per sample, per phenotype class, per allele
  
  AuxCountReport * aux = (AuxCountReport*)p;
  
  const int ns = v.n_samples();
  
  for (int s = 0 ; s < ns ; s++ )
    {

      // Attach sample, along with storage, etc
      
      SampleVariant &             svar          =  v.sample( s );
      SampleVariant &             svar_meta     =  v.sample_metainformation( svar );      
      std::string                 sample_label  =  g.vardb.file_tag( svar.fileset() );
      

      // phenotype
      
      int pi = aux->by_case_control ? 1 : 0;
      
      // collect genotypes, counts and groups

      std::stringstream ginfo;
      std::stringstream cinfo;
      std::stringstream pinfo;      

      // ACDB format      

      if ( aux->acdb_format ) 
	{

	  ginfo << "_GENO=";
	  cinfo << "_GCNT=";
	  pinfo << "_PGRP=";
	  
	  bool first = true;
	  
	  while ( 1 ) 
	    {
	      
	      // all/case/control  (0/2/1)
	      std::string phe_label = pi == 0 ? "0" : pi == 1 ? "2" : "1" ;
	      affType aff = pi == 0 ? UNKNOWN_PHE : pi == 1 ? CASE : CONTROL ;
	      
	      std::map<std::string,int> c = v.genotype_counts( svar , aff , aux->unphased );
	      
	      std::map<std::string,int>::iterator i = c.begin();
	      while ( i != c.end() )
		{
		  if ( ! first ) 
		    {
		      ginfo << ",";
		      cinfo << ",";
		      pinfo << ",";
		    }
		  first = false;
		  
		  ginfo << i->first;
		  cinfo << i->second;
		  pinfo << phe_label;
		  ++i;
		}
	      
	      if ( pi == 0 ) break;
	      ++pi;
	      if ( pi == 3 ) break;
	      
	    } // next phenotype group
	}
      else // pretty-print VCF format 
	{
	  ginfo << "GENO=";
	  
	  std::map<std::string,int> gtyp;
	  std::map<std::string,std::map<std::string,int> > cnts;	  
	  
	  while ( 1 ) 
	    {
	      
	      // all/case/control  (0/2/1)
	      std::string phe_label = pi == 0 ? "0" : pi == 1 ? "2" : "1" ;
	      affType aff = pi == 0 ? UNKNOWN_PHE : pi == 1 ? CASE : CONTROL ;
	      
	      std::map<std::string,int> c = v.genotype_counts( svar , aff , aux->unphased );
	      
	      
	      std::map<std::string,int>::iterator i = c.begin();
	      while ( i != c.end() )
		{
		  gtyp[ i->first ] += i->second;
		  cnts[ i->first ][phe_label] = i->second;
		  ++i;
		}
	      
	      // next 'phenotype'
	      if ( pi == 0 ) break;
	      ++pi;
	      if ( pi == 3 ) break;

	    }
	  
	  // display
	  
	  if ( pi == 0 ) 
	    {
	      bool first = true;    
	      std::map<std::string,int>::iterator ii = gtyp.begin();
	      while ( ii != gtyp.end() )
		{		  
		  if ( ! first ) ginfo << ",";
		  else first = false;		  
		  ginfo << ii->second << "(" << ii->first << ")";
		  ++ii;
		}
	    }
	  else // implies a split by a dichotomous phenotype
	    {
	      bool first = true;    
	      std::map<std::string,std::map<std::string,int> >::iterator ii = cnts.begin();	  
	      while ( ii != cnts.end() )
		{		  
		  if ( ! first ) ginfo << ",";
		  else first = false;		  
		  ginfo << ii->second["2"] << ":" << ii->second["1"] << "(" << ii->first << ")";
		  ++ii;
		}
	  
	    }
	}

      // write to VCF
      
      plog << v.chromosome() << "\t"
	   << v.position() << "\t"
	   << v.name() << "\t"
	   << svar_meta.reference() << "\t"
	   << svar_meta.alternate() << "\t";
      
      if ( svar_meta.quality() < 0 ) 
	plog << ".\t";
      else
	plog << svar_meta.quality() << "\t";
      
      plog << svar_meta.filter() << "\t";
      
      if ( aux->acdb_format ) plog << "_S=" << g.vardb.file_tag( svar.fileset() ) << ";" ;
      else plog << "SAMPLE=" << g.vardb.file_tag( svar.fileset() ) << ";" ;
 
      plog << ginfo.str() << ";";

      if ( aux->acdb_format ) plog << cinfo.str() << ";"
				   << pinfo.str() << ";";

      plog << svar_meta.meta << "\n";
      
      // next sample
    }

}



bool Pseq::VarDB::make_counts_file( Mask &m, const std::string & name  )
{
  
  AuxCountReport aux;
  aux.label = name;
  aux.by_case_control = g.phmap.type() == PHE_DICHOT;
  aux.unphased = ! args.has( "show-phase" );
  aux.acdb_format = args.has( "acdb" );

  // VCF header
  
  plog << "##fileformat=VCFv4.0\n"
       << "##source=pseq\n"
       << "##_PROJ=" << name << "\n";
  
  std::set<int> w = g.indmap.samples();
  std::set<int>::iterator i = w.begin();
  while ( i != w.end() )
    {
      plog << "##_N=" << *i << "," << g.indmap.size( *i ) << "\n";
      ++i;
    }
  
  
  if ( aux.acdb_format ) 
    plog << "##INFO=<ID=_S,Number=1,Type=String,Description=\"Sample tag\">\n"
	 << "##INFO=<ID=_GENO,Number=.,Type=String,Description=\"Genotypes\">\n"
	 << "##INFO=<ID=_GCNT,Number=.,Type=Integer,Description=\"Genotype counts\">\n"
	 << "##INFO=<ID=_PGRP,Number=.,Type=Integer,Description=\"Phenotypic groups (0=all;1=control;2=case)\">\n";
  else
    plog << "##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample tag\">\n"
	 << "##INFO=<ID=GENO,Number=.,Type=String,Description=\"Genotype counts\">\n";

      plog   << MetaInformation<VarMeta>::headers( )
	 << MetaInformation<VarFilterMeta>::headers( META_GROUP_FILTER );
  
  plog << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

  g.vardb.iterate( f_counts_report , &aux , m );

  return true;
}


bool Pseq::VarDB::simple_counts( Mask & m , bool genotypes )
{
  
  OptSimpleCounts opt;
  opt.apply_annot = args.has( "annotate" );
  opt.apply_full_annot = args.has( "full-annotate" );
  if ( opt.apply_full_annot ) opt.apply_annot = true;
  opt.dichot_pheno = g.phmap.type() == PHE_DICHOT;
  opt.show_filter = args.has( "show-filters" );
  opt.genotypes = genotypes;
  
  plog << "VAR";
  plog << "\tREF/ALT";
  
  if ( genotypes )
    {
      plog << "\tGENOTYPES";
    }
  else
    {      
      plog << "\tMINOR";

      // binary phenotype?
      if ( g.phmap.type() == PHE_DICHOT )
	plog << "\tCNTA"
	     << "\tCNTU"
	     << "\tTOTA\tTOTU";
      else
	plog << "\tCNT\tTOT";
    }

  if ( opt.apply_annot )
    plog << "\tFUNC"
	 << "\tGENE";
  
  if ( opt.apply_full_annot )
    plog << "\tFULL"
	 << "\tPROTEIN"	 
	 << "\tALIAS";
  
  
  //
  // Load transcripts for annotation
  //

  if ( opt.apply_annot )
    Annotate::load_transcripts( LOCDB, PLINKSeq::DEFAULT_LOC_GROUP() );
  
  
  //
  // Any optional variant meta-fields to be displayed?
  //

  opt.meta = args.get_set( "meta" );
  std::set<std::string>::iterator i = opt.meta.begin();
  while ( i != opt.meta.end() )
    {
      plog << "\t" << *i ;
      ++i;
    }

  if ( opt.show_filter )
    plog << "\tFILTER";
  
  plog << "\n";


  //
  // Process counts
  //

  g.vardb.iterate( f_simple_counts , &opt, m ); 
  
  return true;
}



void f_write_to_bcf( Variant & var , void * p )
{
  if ( ! var.valid() ) return; 
  BCF * bcf = (BCF*)p;  
  bcf->write_record( var );  
}


bool Pseq::VarDB::write_BCF( Mask & mask , const std::string & bcffile )
{
  
  // We will either be in single VCF mode, or we will be dumping 
  // a VARDB into a single BCF file. Either way, we will have read
  // the header (from VARDB, or from the VCF) by now.  Thus, first 
  // create that.
  
  BCF bcf( bcffile );  
  bcf.create_header();		    
  IterationReport report = g.vardb.iterate( f_write_to_bcf , &bcf , mask );  
  bcf.close();  
  return true;
}




