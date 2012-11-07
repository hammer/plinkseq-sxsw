
#include "func.h"
#include "views.h"
#include "summaries.h"
#include "util.h"

#include "plinkseq.h"

#include <vector>
#include <fstream>

using namespace std;

extern GStore g;
extern Pseq::Util::Options args;

void Pseq::finished()
{  
  
  time_t curr=time(0);
  std::string tdstamp = (std::string)ctime(&curr);
  plog >> "\nAnalysis finished " >> tdstamp    
       >> ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::\n";

  // if ( g.gseq_mode() )
  //   {
  //     std::ofstream O1( g.gseq_history().c_str() , std::ios::out | std::ios::app );
  //     O1 << "_STATUS" << "\t"
  // 	 << g.gseq_job() << "\t"
  // 	 << "Done" << "\n";
  //     O1.close();
  //   } 

  exit(0);
}

bool Pseq::new_project( std::string project , Pseq::Util::Options & args )
{

  //
  // Make end in .pseq for actual filename
  //
  
  std::string project_filename = Helper::ends_with( project , ".pseq" ) ? project : project + ".pseq" ; 

  
  //
  // Ensure full paths are used
  //
  
  project = Helper::fullpath( project );
 

  //
  // Core folders, for output and input
  //

  std::string output_folder = args.has("output") ?
    args.as_string( "output" ) : 
    project + "_out/";
  
  output_folder = Helper::fullpath( output_folder );


  // Fail if VARDB or INDDB already exist, unless they have been explicitly 
  // mentioned on the command line
    
  if ( ! args.has( "vardb" ) ) 
    {
      if ( Helper::fileExists( output_folder + "vardb" ) )
	Helper::halt( "VARDB [ " + output_folder + "vardb ] already exists, delete file first" );
    }

  if ( ! args.has( "inddb" ) ) 
    {
      if ( Helper::fileExists( output_folder + "inddb" ) )
	Helper::halt( "INDDB [ " + output_folder + "inddb ] already exists, delete file first" );
    }
  
  

  // write a new project file

  plog << "Creating new project specification file [ " << project_filename << " ]\n";

  std::ofstream P( project_filename.c_str() , std::ios::out );
  
  if ( P.bad() ) Helper::halt("problem writing new project file " + project );
  
  // Search the options for the following:
  
  // --output ( default project_out/
  // --vcf 1+ files

  // --vardb  default {output}/vardb
  // --inddb  default {output}/locdb
  // --locdb  default {output}/locdb
  // --refdb  default {output}/locdb
  // --seqdb  default {output}/locdb
  
  
  // As of 0.09, use KEY  VALUE  ordering 
  // and specify PROJN
  
  P << "PROJN\t" << PLINKSeq::PROJECT_VERSION_NUMBER() << "\n";

  //
  // Main output folder
  //
  
  Helper::ensure_folder( output_folder );
  
  P << "OUTPUT\t" << output_folder << "\n";
  
  //
  // Resources folder
  //

  std::string resources_folder = args.has("resources") ?
    args.as_string( "resources" ) :
    project + "_res/";

  resources_folder = Helper::fullpath( resources_folder );
  Helper::ensure_folder( resources_folder );
  
  P << "RESOURCES\t" << resources_folder << "\n";
  
  //
  // Scratch files
  //

  if ( args.has( "scratch" ) )
    {
      std::string scratch_folder = args.as_string( "scratch" );
      scratch_folder = Helper::fullpath( scratch_folder );
      Helper::ensure_folder( scratch_folder );
      P << "TEMP\t" << scratch_folder << "\n";
    }


  //
  // Meta-information on meta-information
  //
  
  if ( args.has( "metameta" ) )
    {
      std::string metameta = args.as_string( "metameta" );
      metameta = Helper::fullpath( metameta );
      P << "METAMETA\t" << metameta << "\n";
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
	  P << "VCF\t" << t[i] << "\n";      
	}
    }
  
  //
  // Databases
  //

  if ( args.has( "vardb" ) ) 
    P << "VARDB\t" << Helper::fullpath( args.as_string( "vardb" ) ) << "\n";
  else
    P << "VARDB\t" << output_folder << "vardb\n";
  
  if ( args.has( "inddb" ) ) 
    P << "INDDB\t" << Helper::fullpath( args.as_string( "inddb" ) ) << "\n";
  else
    P << "INDDB\t"<< output_folder << "inddb\n";
  
  if ( args.has( "locdb" ) ) 
    P << "LOCDB\t" << Helper::fullpath( args.as_string( "locdb" ) ) << "\n";
  else
    P << "LOCDB\t" << resources_folder << "locdb\n";
  
  if ( args.has( "refdb" ) ) 
    P << "REFDB\t" << Helper::fullpath( args.as_string( "refdb" ) ) << "\n";
  else
    P << "REFDB\t" << resources_folder << "refdb\n";
  
  if ( args.has("seqdb") )
    P << "SEQDB\t" << Helper::fullpath( args.as_string( "seqdb" ) ) << "\n";
  else
    P << "SEQDB\t" << resources_folder << "seqdb\n";
  
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
  Out & pout = Out::stream( "vcf" );
  pout << v.VCF();
}

void f_vcfz( Variant & v , void * p )
{
  VCFZ * vcfz = (VCFZ*)p;
  vcfz->write_record( v );
}

bool Pseq::VarDB::write_VCF(Mask & m , bool compressed )
{

  Out & pout = Out::stream( "vcf" );

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
  
  pout << "##fileformat=" << PLINKSeq::CURRENT_VCF_VERSION() << "\n"
       << "##source=pseq\n"
       << MetaInformation<VarMeta>::headers( )
       << MetaInformation<GenMeta>::headers( META_GROUP_GEN )
       << MetaInformation<VarFilterMeta>::headers( META_GROUP_FILTER );
  
  // Header line containing Individual IDs
  
  const int n = g.indmap.size();
  pout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";  
  for ( int i=0; i<n; i++) pout << "\t" << g.indmap(i)->id();
  pout << "\n";
  
  // Variants
  
  IterationReport report = g.vardb.iterate( f_vcf , NULL , m );
  
  

  
}


bool Pseq::VarDB::write_PED(Mask & m, bool use_family_id )
{

  // Need to write TPED and TFAM files
  
  if ( ! g.inddb.attached() ) return false;
  if ( ! g.vardb.attached() ) return false;

  const int n = g.indmap.size();

  Out & tfam = Out::stream( "tfam" );
  
  for (int i = 0 ; i < n ; i++) 
    {      

      Individual * pmapped = g.indmap.ind( i );      

      if ( use_family_id ) 
	tfam << pmapped->fid() << "\t"
	     << pmapped->iid() << "\t";
      else
	tfam << pmapped->id() << "\t"
	     << "." << "\t";
      
      if ( pmapped->father() == "." )
	tfam << "0\t";
      else
	tfam << pmapped->father() << "\t";
      
      if ( pmapped->mother() == "." )
	tfam << "0\t";
      else
	tfam << pmapped->mother() << "\t";
      
      // 1=male, 2=female, 0=unknown coding for sex
      tfam << pmapped->sex() << "\t";
      
      // Phenotype
      
      if ( g.phmap.type() == PHE_DICHOT )
	{
	  if ( pmapped->affected() == CASE ) tfam << "2\n";
	  else if ( pmapped->affected() == CONTROL ) tfam << "1\n";
	  else tfam << "0\n";
	}
      else if ( g.phmap.type() == PHE_QT )
	tfam << pmapped->qt() << "\n";
      else if ( g.phmap.type() == PHE_FACTOR )
	tfam << pmapped->group_label() << "\n";
      else
	tfam << "0\n";
    }
 
  tfam.close();
  
  // Now MAP and genotype data
  
  IterationReport report = g.vardb.iterate( f_view_tped , NULL , m );
  
  return true;
}

bool Pseq::VarDB::dump_indiv()
{
  Out & pout = Out::stream( "indiv" );

  int cnt = 0;
  if ( ! g.vardb.attached() ) return false;
  std::map<int,std::string> files = g.vardb.fetch_files();
  std::map<int,std::string>::iterator i = files.begin();
  while ( i != files.end() )
    {
      std::vector<std::string> inds = g.vardb.fetch_individuals( i->first );
      for (int j = 0 ; j < inds.size(); j++)
	{
	  pout << ++cnt << "\t" 
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
  Out & pout = Out::stream( "lik" );

  std::vector<std::string> ids = g.indmap.ind_id();
  pout << "VAR\tREF\tALT";
  for (int i=0; i<ids.size(); i++) 
    pout << "\t" << ids[i] 
	 << "\t" << ids[i] 
	 << "\t" << ids[i] ;
  pout << "\n";
  
  // Variants

  IterationReport report = g.vardb.iterate( f_view_lik , &g , m );

  return true;

}



//
//  -------- WRITE HAPS -------------
//

struct aux_haps_filemap 
{ 

  aux_haps_filemap( const std::string & r ) { filename_root = r; } 
  
  std::ofstream * handle( const std::string & id, const std::string & chr , int hap )
  {
    const std::string corename = id + "-" + ( hap==1 ? "1" : "2" ) + "-" + chr;
    std::map<std::string,std::ofstream*>::iterator ii = m.find( corename );
    if ( ii != m.end() ) return ii->second;
    std::ofstream * nf = new std::ofstream;
    nf->open( ( filename_root + "-" + corename + ".haps" ).c_str() );
    if ( ! nf->is_open() ) Helper::halt( "trouble opening file to write to : " + corename );
    *nf << id << "\t" << "HAPLO" << hap << "\t";
    m[corename] = nf;
    return nf;
  }
  
  ~aux_haps_filemap()
  {
    std::map<std::string,std::ofstream*>::iterator ii = m.begin();
    while ( ii != m.end() )
      {
	if ( ii->second ) 
	  {
	    ii->second->close();
	    delete ii->second;
	  }
	++ii;
      }
  }
  
  std::string filename_root;
  std::map<std::string,std::ofstream*> m;

};


void f_write_haps( Variant & v , void * p )
{
  aux_haps_filemap * files = (aux_haps_filemap*)p;
  const int n = v.size();
  const std::string chr = Helper::chrCode( v.chromosome() );
  for (int i=0;i<n;i++)
    {
      std::ofstream * f1 = files->handle( v.ind(i)->id() , chr , 1 ) ;
      std::ofstream * f2 = files->handle( v.ind(i)->id() , chr , 2 );
      const Genotype & g = v(i);
      *f1 << v.allele1_label( v(i) );
      *f2 << v.allele1_label( v(i) ); 
      //      std::cout << "want to write " << v.allele1_label( v(i) ) << " and " << v.allele2_label( v(i) ) << "\n";
    }
}


bool Pseq::VarDB::write_haps(Mask & m , const std::string & rt )
{
  // helper struct to organise filemap
  aux_haps_filemap filemap( rt );
  
  // write
  IterationReport report = g.vardb.iterate( f_write_haps , &filemap , m );
  
  std::cerr << "done... now, for each chromosome, concatenate, e.g:\n  cat " << rt << "-*-chr1 > all-chr1.haps\n";
}


bool Pseq::VarDB::write_matrix(Mask & m)
{
  Out & pout = Out::stream( "matrix" );
  std::vector<std::string> ids = g.indmap.ind_id();
  pout << "VAR\tREF\tALT";
  for (int i=0; i<ids.size(); i++) 
    pout << "\t" << ids[i] ; 
  pout << "\n";
  
  IterationReport report = g.vardb.iterate( f_view_matrix , &g , m );
  return true;
}

bool Pseq::VarDB::write_meta_matrix( Mask & m )
{

  Out & pout = Out::stream( "matrix" );

  std::map<int,string> files = g.vardb.fetch_files( &m );    
  
  // Display consensus meta-info, if >1 file
  
  pout << "CHR\tPOS\tID";

  // Variant level information

  int nelem_static = MetaInformation<VarMeta>::n_visible_keys_static();
  int nelem_nonstatic = MetaInformation<VarMeta>::n_visible_keys_nonstatic();

  if ( nelem_static ) pout << "\t" << MetaInformation<VarMeta>::display_header_static();
  
  // Consensus sample-variant
  if ( nelem_nonstatic ) pout << "\t" << MetaInformation<VarMeta>::display_header_nonstatic();
  
  if ( files.size() > 1 ) 
    {      
      std::map<int,std::string>::iterator i = files.begin();
      while ( i != files.end() ) 
	{ 
	  std::string m = MetaInformation<VarMeta>::display_header_nonstatic( "F"+Helper::int2str(i->first)+"_" ) ;
	  std::string f = MetaInformation<VarFilterMeta>::display_header( "F"+Helper::int2str(i->first)+"_" );
	  if ( m != "" ) pout << "\t" << m;
	  if ( f != "" ) pout << "\t" << f;
	  ++i;
	}  
    }
  else
    {      
      std::string f = MetaInformation<VarFilterMeta>::display_header( );
      if ( f != "" ) pout << "\t" << f;
    }
  pout << "\n";  
  
  IterationReport report = g.vardb.iterate( f_view_meta_matrix , &files , m );
  
  return true;
}


bool Pseq::VarDB::write_var_meta_matrix(Mask & m, std::string & name)
{
  Out & pout = Out::stream( "matrix" );

  std::vector<std::string> ids = g.indmap.ind_id();
  pout << "MARKER";
  for (int i=0; i<ids.size(); i++) 
    pout << "\t" << ids[i] ; 
  pout << "\n";

  IterationReport report = g.vardb.iterate( f_view_var_meta_matrix , &name , m );
  return true;
}


bool Pseq::VarDB::write_gene_matrix(Mask & m, OptGMatrix & opt)
{
  Out & pout = Out::stream( "matrix" );
  std::vector<std::string> ids = g.indmap.ind_id();
  pout << "GENE\tNV";
  for (int i=0; i<ids.size(); i++) 
    pout << "\t" << ids[i] ; 
  pout << "\n";
  IterationReport report = g.vardb.iterate( f_view_gene_matrix , &opt , m );
  return true;
}


bool Pseq::VarDB::write_gene_meta_matrix(Mask & m, OptGMetaMatrix & opt)
{
  Out & pout = Out::stream( "matrix" );
  std::vector<std::string> ids = g.indmap.ind_id();
  pout << "GENE\tNV";
  for (int i=0; i<ids.size(); i++) 
    pout << "\t" << ids[i] ; 
  pout << "\n";
  IterationReport report = g.vardb.iterate( f_view_gene_meta_matrix , &opt , m );
  return true;
}


bool Pseq::VarDB::consolidate(Mask & m, std::string label)
{
  return false;
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

  Out & pout = Out::stream( "uniq" );

  // Did we see in all group members?
  
  if ( opt->ingroup_req == -1 )
    {
      if ( obs_group == opt->indiv.size() )
	pout  << obs_group << "\t"
	      << obs_outgroup << "\t"
	      << v << "\t" << v.meta << "\n";
    }
  else
    {     
      if ( obs_group >= opt->ingroup_req ) 
 	pout << obs_group << "\t"
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
  return false;
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
	plog.counter1( "replaced " + Helper::int2str( cnt ) + " locus names" );
    }
  
  plog.counter1("\n");
  
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

    Out & pout = Out::stream( "loci" );
    
    //
    // Read a vector of regions
    //
    //std::vector<Region> regions;
    
    std::ifstream IN1( filename.c_str() , std::ios::in );
    
    int readRegion = 0, inv = 0, blank = 0;
    
    while ( ! IN1.eof() )
      {

	std::vector<std::string> tok = Helper::tokenizeLine( IN1 );
	if (IN1.fail())
		break;

	const int sz = tok.size();
	if (sz > 0) {
		const std::string& firstChar = tok[0].substr(0,1);
		if (firstChar == "#" || firstChar == "@")
		continue; // ignore comment lines
	}

    ++readRegion;
	std::string region;
	
	if ( sz == 1 ) 
	  region = tok[0];
	else if ( sz == 2 ) 
	  region = tok[0] + ":" + tok[1] + ".." + tok[1] ;
	else if ( sz >= 3 )
	  region = tok[0] + ":" + tok[1] + ".." + tok[2] ;
	else {
	  ++inv;
	  continue;
	}
	
	bool valid = true;
	
	Region r(region,valid);

	// Get overlap
	
	if ( valid )
	  {
	    
	    std::set<Region> olaps = db.get_regions( group ,  r ) ;	  
	    
	    pout << r.coordinate() << "\t";
	    
	    if ( olaps.size() == 0 ) ++blank;
	    std::set<Region>::iterator i = olaps.begin();
	    std::set<string> genes;
	    while ( i != olaps.end() )
	      {
		genes.insert( i->altname );
		++i;
	      }	  
	    
	    // report # of loci overlapping
	    pout << genes.size() << "\t";

	    if ( genes.size() == 0 ) 
	      pout << ".";
	    else
	      {
		std::set<std::string>::iterator j = genes.begin();
		while ( j != genes.end() )
		  {
		    if ( j != genes.begin() ) pout << "|";
		    pout << *j;
		    ++j;
		  }
	      }
	    
	    pout << "\n";
	  
	}
	else
	  ++inv;
      }
    
    plog << "read " << readRegion << " regions\n";
    
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
      

// bool Pseq::SeqDB::load_transcripts( std::string label )
// {
//   return Annotate::load_transcripts( LOCDB , label );
// }

  
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
  
  Out & pout = Out::stream( "indiv" );
  
  if ( ! g.inddb.attached() ) 
    {
      plog.warn("no INDDB attached");
      return false;
    }
  
  bool verbose = args.has( "verbose" );
  
  // List all information in IndDB, one person per line
  
  std::map<std::string,std::vector<std::string> > pheno = g.inddb.fetch_phenotype_info();
  std::map<std::string,std::vector<std::string> >::iterator p = pheno.begin();
  
  while ( p != pheno.end() )
    {
      pout << "#" << p->first << " "
	   << "(" << p->second[0]  << ") "
	   << p->second[1] << "\n";
      ++p;
    }

  // Explicitly note phenotype/strata
  
  pout << "#PHE\t" << ( g.phmap.phenotype_set() ? g.phmap.phenotype() : "." ) << "\n";
  pout << "#STRATA\t" << ( g.phmap.strata_set() ? g.phmap.strata() : "." ) << "\n";


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
  
  if ( ! verbose ) 
    {
      pout << "#ID\tFID\tIID\tMISS\tSEX\tPAT\tMAT\tMETA";
      if ( g.phmap.phenotype_set() && g.phmap.type() == PHE_DICHOT ) pout << "\tPHE";
      if ( g.phmap.strata_set() ) pout << "\tSTRATA";
      pout << "\n";
    }


  std::set<Individual*> observed;
  
  for (int i = 0 ; i < n ; i++) 
    {
      
      // note: should never return a null, if based on 0..(n-1) where n is 
      // returned from create_from_mask()
      
      Individual * pmapped = g.indmap.ind( i );
      
      if ( ! verbose ) 
	{
	  pout << *pmapped;
	  
	  if ( g.phmap.type() == PHE_DICHOT )
	    {
	      if ( pmapped->affected() == CASE ) pout << "\tCASE";
	      else if ( pmapped->affected() == CONTROL ) pout << "\tCONTROL";
	      else pout << "\t.";
	    }
	  
	  if ( g.phmap.strata_set() )
	    {
	      pout << "\t" << pmapped->group_label();
	    }
	}
      else // VERBOSE output mode
	{
	  
	  pout << "ID:\t" << pmapped->id() << "\n";
	  
	  if ( g.phmap.type() == PHE_DICHOT )
	    {
	      if ( pmapped->affected() == CASE ) pout << "\tPHENO:\tCASE\n";
	      else if ( pmapped->affected() == CONTROL ) pout << "\tPHENO:\tCONTROL\n";
	      else pout << "\tPHENO:\t.\n";
	    }
	  
	  if ( g.phmap.strata_set() )
	    {
	      pout << "\tSTRATA\t" << pmapped->group_label() << "\n";
	    }

	  // Now all phenotypes, in long form:
	  pout << pmapped->meta.display();
	  
	}
      
      pout << "\n";      
      
      observed.insert( pmapped );
      
    }

  
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

  bool refvars = ! ( refgroup == "" || refgroup == "." ) ;
  
  if ( refvars )
    {
      if ( ! g.refdb.attached() ) 
	Helper::halt("REFDB not attached");
      if ( ! g.refdb.init_iterate( refgroup ) )
	Helper::halt("could not initiate specified reference group " + refgroup );
    }
  
  RangeIntersector range_intersector( f_next_refvar , f_loc_stats_range_iterator_report );
  
  Out & pout = Out::stream( "locstats" );
  
  //
  // Header
  //
  
  pout << "LOCUS" << "\t"
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
  // di- and tri-nucleotide counts per transcript
  //

  std::map<std::string,int> dic;
  std::map<std::string,int> tri;
  

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
      
      // count di- and tri-nuc contexts
      g.seqdb.dinucleotide( *i , dic );
      g.seqdb.trinucleotide( *i , tri );      
      
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
	  
	  g.seqdb.N( sr , sn, stot ) ;
	  
	  if ( okay )
	    {
	      ++nreg;
	      gc += sgc;
	      n += sn;
	      tot += stot;
	    }
	  
	  
	  if ( sub ) 
	    {	      
	      pout << i->name << "\t" 
		   << alias << "\t"
		   << "NA\t"
		   << s+1 << "\t"
		   << sr.coordinate() << "\t";
	      
	      pout << stot << "\t"
		   << "NA" << "\t";
	      
	      if ( okay ) 		
		{
		  if ( stot-sn > 0 ) 
		    pout << (double)sgc/(double)(stot-sn) << "\t";
		  else 
		    pout << "NA" << "\t";
		  
		  pout << "NA" << "\t";

		  if ( stot > 0 ) 
		    pout << (double)sn/(double)stot << "\t";
		  else
		    pout << "NA" << "\t";

		  pout << "NA";
		}
	      else
		pout << "NA\tNA\tNA\tNA";
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
	      if ( sub ) pout << "\t" << c << "\tNA";
	    }
	  
	  if ( sub ) pout << "\n";
	  
	}

      
	  
      //
      // Print summary of all subregions (S=0 code in output)
      //
      

      pout << i->name << "\t"
	   << alias << "\t";
      pout << ( nosubs ? 1 : nreg ) << "\t";
      if (sub) pout << "0\t";	  
      pout << i->coordinate() << "\t";
      
      if ( nosubs ) pout << ".\t";
      else pout << tot << "\t";
      pout << tot0 << "\t";	  
      
      bool only0 = tot == 0;
      
      if ( okay && ( nosubs || nreg > 0 ) )
	{
	  if ( only0 ) pout << ".\t"; 
	  else
	    {
	      if ( tot-n != 0 )
		pout << (double)gc/(double)(tot-n) << "\t";
	      else
		pout << "NA\t";
	    }
	  
	  if ( tot0-n0 ) 
	    pout << (double)gc0/(double)(tot0-n0) << "\t";
	  else
	    pout << "NA\t";
	  
	  if ( only0 ) pout << ".\t"; else pout << (double)n/(double)tot << "\t";
	  pout << (double)n0/(double)tot0 ;
	}
      else
	pout << "NA\tNA\tNA\tNA";
      
      if ( refvars ) 
	pout << "\t" << refcnt 
	     << "\t" << refcnt_whole ;
      
      pout << "\n";

      //
      // next region
      //

      ++i;

    }


  // DINUCLEOTIDE CALCS

//   std::map<std::string,int>::iterator ii = dic.begin();
//   while ( ii != dic.end() ) 
//     {
//       std::cout << ii->first << "\t" << ii->second << "\n";
//       ++ii;
//     }


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
  Out & pout = Out::stream( "scores" );

  PPH2DBase * pph2 = (PPH2DBase*)p;
  double score = 0;   
  int prediction = 0;
  // only output scores for missense variants found in PPH2
  if ( pph2->score(v,score,prediction) )
    pout << v << "\t" << score << "\t" << prediction << "\n";
  else
    pout << v << "\tNA\t" << prediction << "\n";
}

bool Pseq::PPH2DB::score( Mask & m , const std::string & dbname )
{
  PPH2DBase ppdb;
  ppdb.attach( dbname );
  ppdb.set_locdb( &g.locdb );

  Annotate::setDB( LOCDB );
  Annotate::set_transcript_group( PLINKSeq::DEFAULT_LOC_GROUP() );

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
  
  Out & pout = Out::stream( "vcf" );

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
      
      pout << v.chromosome() << "\t"
	   << v.position() << "\t"
	   << v.name() << "\t"
	   << svar_meta.reference() << "\t"
	   << svar_meta.alternate() << "\t";
      
      if ( svar_meta.quality() < 0 ) 
	pout << ".\t";
      else
	pout << svar_meta.quality() << "\t";
      
      pout << svar_meta.filter() << "\t";
      
      if ( aux->acdb_format ) pout << "_S=" << g.vardb.file_tag( svar.fileset() ) << ";" ;
      else pout << "SAMPLE=" << g.vardb.file_tag( svar.fileset() ) << ";" ;
 
      pout << ginfo.str() << ";";

      if ( aux->acdb_format ) pout << cinfo.str() << ";"
				   << pinfo.str() << ";";

      pout << svar_meta.meta << "\n";
      
      // next sample
    }

}



bool Pseq::VarDB::make_counts_file( Mask &m, const std::string & name  )
{
  
  Out & pout = Out::stream( "vcf" );

  AuxCountReport aux;
  aux.label = name;
  aux.by_case_control = g.phmap.type() == PHE_DICHOT;
  aux.unphased = ! args.has( "show-phase" );
  aux.acdb_format = args.has( "acdb" );

  // VCF header
  
  pout << "##fileformat=VCFv4.0\n"
       << "##source=pseq\n"
       << "##_PROJ=" << name << "\n";
  
  std::set<int> w = g.indmap.samples();
  std::set<int>::iterator i = w.begin();
  while ( i != w.end() )
    {
      pout << "##_N=" << *i << "," << g.indmap.size( *i ) << "\n";
      ++i;
    }
  
  
  if ( aux.acdb_format ) 
    pout << "##INFO=<ID=_S,Number=1,Type=String,Description=\"Sample tag\">\n"
	 << "##INFO=<ID=_GENO,Number=.,Type=String,Description=\"Genotypes\">\n"
	 << "##INFO=<ID=_GCNT,Number=.,Type=Integer,Description=\"Genotype counts\">\n"
	 << "##INFO=<ID=_PGRP,Number=.,Type=Integer,Description=\"Phenotypic groups (0=all;1=control;2=case)\">\n";
  else
    pout << "##INFO=<ID=SAMPLE,Number=1,Type=String,Description=\"Sample tag\">\n"
	 << "##INFO=<ID=GENO,Number=.,Type=String,Description=\"Genotype counts\">\n";

      pout   << MetaInformation<VarMeta>::headers( )
	 << MetaInformation<VarFilterMeta>::headers( META_GROUP_FILTER );
  
  pout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

  g.vardb.iterate( f_counts_report , &aux , m );

  return true;
}


bool Pseq::VarDB::simple_counts( Mask & m , bool genotypes , bool qt )
{

  // Wrapper for either QT means, genotypic counts, or allelic counts
  Out & pout = Out::stream( qt ? "means" : ( genotypes ? "gcounts" : "counts" ) );
  
  OptSimpleCounts opt;
  opt.apply_annot = args.has( "annotate" );
  opt.apply_full_annot = args.has( "full-annotate" );
  if ( opt.apply_full_annot ) opt.apply_annot = true;
  opt.dichot_pheno = g.phmap.type() == PHE_DICHOT;
  opt.qt_pheno = g.phmap.type() == PHE_QT;
  opt.show_filter = args.has( "show-filters" );
  opt.genotypes = genotypes;
  
  if ( qt && ! opt.qt_pheno ) Helper::halt( "no QT specified" );  
  
  pout << "VAR";
  pout << "\tREF/ALT";
  
  if ( genotypes )
    {
      pout << "\tGENOTYPES";
    }
  else
    {      
      pout << "\tMINOR";

      // binary phenotype?
      if ( opt.dichot_pheno )
	pout << "\tCNTA"
	     << "\tCNTU"
	     << "\tTOTA\tTOTU";
      else if ( opt.qt_pheno ) // collapse HET and HOM alternate genotypes for now
	pout << "\tREFMEAN\tREFSD\tREFN"
	     << "\tALTMEAN\tALTSD\tALTN";
      else
	pout << "\tCNT\tTOT";
    }

  if ( opt.apply_annot )
    pout << "\tFUNC"
	 << "\tGENE";
  
  if ( opt.apply_full_annot )
    pout << "\tFULL"
	 << "\tPROTEIN"	 
	 << "\tALIAS";
  
  
  //
  // Load transcripts for annotation
  //

  if ( opt.apply_annot )
    {
      Annotate::setDB( LOCDB );
      Annotate::set_transcript_group( PLINKSeq::DEFAULT_LOC_GROUP() );
    }
  
  //
  // Any optional variant meta-fields to be displayed?
  //

  opt.meta = args.get_set( "meta" );
  std::set<std::string>::iterator i = opt.meta.begin();
  while ( i != opt.meta.end() )
    {
      pout << "\t" << *i ;
      ++i;
    }

  if ( opt.show_filter )
    pout << "\tFILTER";
  
  pout << "\n";


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

  plog << "Writing to BCF [ " << bcffile << " ]\n";

  BCF bcf( bcffile );  
  bcf.writing();
  bcf.open();
  bcf.write_header();		    
  IterationReport report = g.vardb.iterate( f_write_to_bcf , &bcf , mask );  
  bcf.close();  
  return true;
}





bool Pseq::IndDB::dump_phenotypes( const std::vector<std::string> & pheno , bool as_matrix )
{

  Out & pout = Out::stream( as_matrix ? "indiv.matrix" : "phe" );

  std::vector<bool> m( pheno.size() , true );
  
  if ( ! as_matrix ) 
    {
      
      for (int p=0;p<pheno.size();p++)
	{
	  pType ptype = g.phmap.type( pheno[p] );
	  if      ( ptype == PHE_DICHOT ) pout << "##" << pheno[p] << ",Integer,.,\"Integer phenotype\"\n";
	  else if ( ptype == PHE_QT )     pout << "##" << pheno[p] << ",Float,.,\"Float phenotype\"\n";
	  else if ( ptype == PHE_FACTOR ) pout << "##" << pheno[p] << ",String,.,\"Factor phenotype\"\n";
	  else { m[p] = false; continue; } 
	}
      pout << "#";
    }
  else
    {
      for (int p=0;p<pheno.size();p++)
	{
	  pType ptype = g.phmap.type( pheno[p] );
	  if ( ptype != PHE_DICHOT &&
	       ptype != PHE_QT && 
	       ptype != PHE_FACTOR ) m[p] = false;
	}
    }

  pout << "ID";
  for (int p=0;p<pheno.size();p++)
    if ( m[p] ) pout << "\t" << pheno[p];
  pout << "\n";
  
  const int n = g.indmap.size();

  for (int i=0;i<n;i++)
    {
      pout << g.indmap(i)->id() ;
      for (int p=0;p<pheno.size();p++)
	if ( m[p] ) 
	  {
	    Individual * person = g.indmap(i);
	    if ( person->meta.has_field( pheno[p] ) )
	      {
		pType ptype = g.phmap.type( pheno[p] );
		if      ( ptype == PHE_DICHOT ) pout << "\t" << person->meta.get1_int( pheno[p] ) ;
		else if ( ptype == PHE_QT )     pout << "\t" << person->meta.get1_double( pheno[p] ) ;
		else if ( ptype == PHE_FACTOR ) pout << "\t" << person->meta.get1_string( pheno[p] ) ;
	      }
	    else 
	      {
		pType ptype = g.phmap.type( pheno[p] );
		if      ( as_matrix && ( ptype == PHE_DICHOT || ptype == PHE_QT ) ) pout << "\tNA";
		pout << "\t.";
	      }
	  }
      pout << "\n";
    }
}
