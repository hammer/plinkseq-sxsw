#include "loaders.h"
#include "func.h"
#include "netdb.h"

extern GStore g;
extern Pseq::Util::Options options;
extern Log plog;

bool Pseq::VarDB::load_VCF()
{
  // are any optinal locus groups specified (this will act as a mask, such that
  // we only read in variants that fall in this region
  
  bool region_mask = options.key("filter");
  
  std::string mask;
  if ( region_mask ) mask = options.as<std::string>("filter");

  // filters on meta-fields?
  std::set<std::string> includes, excludes;
  if ( options.key( "meta" ) ) 
    includes = options.get_set( "meta" );
  if ( options.key( "meta.ex" ) ) 
    excludes = options.get_set( "meta.ex" );

  if ( ! options.key("check-reference") ) g.seqdb.dettach();
       
  return g.vardb_load_vcf( includes , excludes , region_mask ? & mask : NULL );
}

bool Pseq::RefDB::load_VCF( const std::string & filename , const std::string & grp ) 
{

  bool region_mask = options.key("filter");
  
  std::set<Region> filter;
  if ( region_mask ) filter = GP->locdb.get_regions( options.as<std::string>("filter") );  
  std::set<Region> * pfilter = region_mask ? &filter : NULL  ;
  
  // filters on meta-fields?
  std::set<std::string> includes, excludes;
  if ( options.key( "meta" ) ) 
    includes = options.get_set( "meta" );
  if ( options.key( "meta.ex" ) ) 
    excludes = options.get_set( "meta.ex" );
  
  std::string comment;
  if ( options.key("description") ) comment = options.as<std::string>( "description" );
  
  if ( ! options.key("check-reference") ) g.seqdb.dettach();

  if ( ! Helper::valid_name( grp ) ) 
    Helper::halt( grp + " is not a valid name" );

  return g.refdb.load_VCF( filename , 
			   grp , 
			   comment , 
			   includes , 
			   excludes , 
			   pfilter );
  
}


bool Pseq::VarDB::load_PLINK( const std::vector<std::string> & name , const Pseq::Util::Options & opt , const std::string & tag )
{
  if ( name.size() == 0 ) return false; 
  std::string fileroot = name[0];
  
  BEDReader b( &g.vardb );
  b.set_root( fileroot );
  if ( ! options.key("check-reference") ) g.seqdb.dettach();
  if ( g.seqdb.attached() ) b.use_seqdb( &g.seqdb );
  b.store_phenotypes( &g.inddb, options.key("phenotype") ? options.as<std::string>("phenotype") : "phe1" );
  if ( tag != "" ) b.set_tag( tag );
  bool okay = b.read_bed();
  if ( ! okay ) plog.warn( "problems detected loading BED file, " + fileroot );
  g.fIndex.append_to_projectfile( Helper::fullpath( fileroot ) , "PLINK" );
  return okay;
}


bool Pseq::SeqDB::load_FASTA( const std::string & filename )
{
  // process necessary options that describe the data loaded in

  std::map<std::string,std::string> meta;
	
  if ( options.key("name") ) meta[ PLINKSeq::SEQDB_NAME_KEY() ] = options.as<std::string>("name");
  else Helper::halt("no name specified in options");

  if ( options.key("build") ) meta[ PLINKSeq::SEQDB_GENOME_BUILD_KEY() ] = options.as<std::string>("build");
  else Helper::halt("no build specified in options");
  
  if ( options.key("repeat-mode") ) 
    {
      std::string m = options.as<std::string>("repeat-mode");
      
      // should be none  = none
      //           N     = N-masked
      //           lower = lower-case
      
      if ( m != "none" && m != "N" && m != "lower" )
	Helper::halt("repeat-mode should be 'none', 'N' or 'lower'");
      
      meta[ PLINKSeq::SEQDB_REPEAT_MODE_KEY() ] = m;
    }
  else Helper::halt("no repeat-mode specified in options");
  
  if ( options.key("description") ) meta[ PLINKSeq::SEQDB_DESCRIPTION_KEY() ] = options.as<std::string>("description");
  else Helper::halt("no description specified in options");
  
  if ( options.key("iupac") || options.key("IUPAC") ) 
    meta[ PLINKSeq::SEQDB_IUPAC_KEY() ] = "1";
  else
    meta[ PLINKSeq::SEQDB_IUPAC_KEY() ] = "0";
  
  g.seqdb.create( filename );
  g.seqdb.loadFASTA( filename , meta );
  
  return true;
}


bool Pseq::RefDB::load_refvar( const std::string & filename , 
			       const std::string & label , 
			       Pseq::Util::Options & options )
{
  
  //
  // Load a generic rectangular flat file into the REFDB
  //  
  

  // We require a header field, that starts #

  InFile ftmp( filename );
  std::vector<std::string> h = ftmp.tokenizeLine();
  ftmp.close();
  ftmp.clear();

  if ( h.size() == 0 ) return false;
  
  if ( h[0].substr(0,1) != "#" ) 
    Helper::halt("reference variant file requires #header as first row");
  
  // skip leading #
  h[0] = h[0].substr(1);

  int col_chr = -1;
  int col_bp1 = -1;
  int col_bp2 = -1;	
  int col_name = -1;
  
  std::map<std::string,int> mf;  // meta-field information
  std::set<int> flags;
  std::set<int> skip;

  //
  // Handle fixed field names (and allow alternatives)
  //

  std::string f_chr = "CHR";
  std::string f_bp1 = "BP1";
  std::string f_bp2 = "BP2";		       
  std::string f_id  = "ID";
  
  // allow command-line swap-ins

  if ( options.key("CHR") ) f_chr   = options.as<std::string>("CHR");
  if ( options.key("BP1") )   f_bp1   = options.as<std::string>("BP1");
  if ( options.key("BP2") )   f_bp2   = options.as<std::string>("BP2");
  if ( options.key("ID") )    f_id    = options.as<std::string>("ID");
  
  for (int i=0; i<h.size(); i++)
    {
      // fixed fields?
      if ( h[i] == f_chr  )        col_chr   = i; 
      else if ( h[i] == f_bp1   )  col_bp1   = i;
      else if ( h[i] == f_bp2   )  col_bp2   = i;
      else if ( h[i] == f_id    )  col_name  = i;
      else // a user-defined type
	{
	  
	  // Should we skip this field?
	  
	  if ( options.value( h[i] , "Skip" ) || 
	       options.value( "Skip" , h[i] ) ) 
	    {
	      skip.insert(i);
	      continue;
	    }

	  
	  // If this is a "Flag" field, it means that values
	  // themselves will be the key, so we do not register the
	  // column name with the meta-information class.
	  
	  if ( options.value( h[i] , "Flag" ) )
	    flags.insert(i);
	  else
	    {
	      mType mt = META_TEXT;	  
	      if ( options.value( h[i] , "Integer" ) ) mt = META_INT;
	      else if ( options.value( h[i] , "Float" ) ) mt = META_FLOAT;
	      else if ( options.value( h[i] , "Bool" ) ) mt = META_BOOL;
	      
	      // ignore description for now
	      std::string desc = "";

	      mf[ h[i] ] = i;
	      registerMetatype( h[i] , mt , 1 , META_GROUP_REF , desc );
	    }  

	}

    }


  //
  // Misc. options
  //

  bool zero1 = options.key("0-start") || options.key("0-based");
  bool zero2 = options.key("0-stop") || options.key("0-based");
  

  //
  // Load data
  //

  uint64_t id = g.refdb.loadRefVariants( filename , 
					 label,
					 col_name, 
					 col_chr , 
					 col_bp1, 
					 mf.size() == 0 ? NULL : &mf ,
					 flags.size() == 0 ? NULL : &flags , 
					 skip.size() == 0 ? NULL : &skip ,
					 col_bp2 , 
					 zero1 , 
					 zero2 );   
  

   return true;
 }


//
// Load a generic file into the SEGDB
//


bool Pseq::LocDB::load_segments( std::string filename , std::string label , Pseq::Util::Options & options )
{
  
  // Format for ".seg" file as follows
  
  // HEADER 
  //  CHR, BP1 BP2 (or POS1, POS2) 
  //   Optional  SEGID
  //   ID  Individual ID
  //   No subregions allowed (for now)
  //   Any other field is treated as meta-information, with type specifiers in Options


  // We require a header field, that starts #

  InFile ftmp( filename );
  std::vector<std::string> h = ftmp.tokenizeLine();
  ftmp.close();
  ftmp.clear();
  
  if ( h.size() == 0 ) return false;
  
  if ( h[0].substr(0,1) != "#" ) 
    Helper::halt("reference variant file requires #header as first row");
  
  // skip leading #
  h[0] = h[0].substr(1);

  int col_chr = -1;
  int col_bp1 = -1;
  int col_bp2 = -1;	
  int col_name = -1;
  int col_indiv = -1;
  int col_sub = -1;
    



  std::map<std::string,int> mf;  // meta-field information
  
  for (int i=0; i<h.size(); i++)
    {
      
      // fixed fields?
      if ( h[i] == "CHR" || h[i] == "CHROM" ) col_chr = i; 
      else if ( h[i] == "BP1" || h[i] == "POS1" ) col_bp1 = i;
      else if ( h[i] == "BP2" || h[i] == "POS2" ) col_bp2 = i;
      else if ( h[i] == "SEGID" ) col_name = i;
      else if ( h[i] == "ID" ) col_indiv = i;
      else // a user-defined type
	{
	  mType mt = META_TEXT;
	  if ( options.value( h[i] , "Integer" ) ) mt = META_INT;
	  else if ( options.value( h[i] , "Float" ) ) mt = META_FLOAT;
	  else if ( options.value( h[i] , "Flag" ) ) mt = META_FLAG;

	  // ignore description for now
	  std::string desc = "";
	  
	  mf[ h[i] ] = i;
	  
	  registerMetatype( h[i] , mt , 1 , META_GROUP_REF , desc );
	}
    }
  

  //
  // Did we observe the required fields?
  //

  if ( col_chr == -1 || col_bp1 == -1 || col_bp2 == -1 || col_indiv == -1 ) 
    Helper::halt("required headers not observed: CHR BP1 BP2 ID)");
  
  return g.segdb.load_regions( filename , 
			       label , 
			       -1 ,  // no col_pos here yet
			       col_chr , col_bp1 , col_bp2 , 
			       col_name , col_sub , 
			       col_indiv , 
			       mf.size() == 0 ? NULL : &mf ) != 0;
  
}


//
// Load generic regions into LOCDB or SEGDB
//

bool Pseq::LocDB::load_generic_regions( std::string & filename , const std::string & label , Pseq::Util::Options & , bool locdb )
{

  LocDBase * db = locdb ? &g.locdb : &g.segdb ;

  // Format for ".reg" file as follows
  
  // HEADER 
  //  CHR, BP1 BP2 (or POS1, POS2) 
  //   or a single POS
  //   ID  Individual ID   ( --> NAME of regions)
  //   Any other field is treated as meta-information, with type specifiers in Options


  // We require a header field, that starts #

  InFile ftmp( filename );
  std::vector<std::string> h = ftmp.tokenizeLine();
  ftmp.close();
  ftmp.clear();
  
  if ( h.size() == 0 ) return false;
  
  if ( h[0].substr(0,1) != "#" ) 
    Helper::halt("reference variant file requires #header as first row");
  
  // skip leading #
  h[0] = h[0].substr(1);

  int col_chr = -1;
  int col_bp1 = -1;
  int col_bp2 = -1;	
  int col_name = -1;
  int col_indiv = -1;
  int col_sub = -1;
  int col_pos = -1;

  std::map<std::string,int> mf;  // meta-field information
  
  for (int i=0; i<h.size(); i++)
    {

      // fixed fields?
      if ( h[i] == "CHR" || h[i] == "CHROM" ) col_chr = i; 
      else if ( h[i] == "BP1" || h[i] == "POS1" ) col_bp1 = i;
      else if ( h[i] == "BP2" || h[i] == "POS2" ) col_bp2 = i;
      else if ( h[i] == "POS" ) col_pos = i;
      else if ( h[i] == "ID" ) col_name = i;
      else // a user-defined type
	{
	  mType mt = META_TEXT;
	  if ( options.value( h[i] , "Integer" ) ) mt = META_INT;
	  else if ( options.value( h[i] , "Float" ) ) mt = META_FLOAT;
	  else if ( options.value( h[i] , "Flag" ) ) mt = META_FLAG;

	  // ignore description for now
	  std::string desc = "";
	  
	  mf[ h[i] ] = i;
	  
	  registerMetatype( h[i] , mt , 1 , META_GROUP_REF , desc );
	}
    }
  
  
  // Did we observe the required fields?
  
  if ( col_pos != -1 ) 
    {
      if ( col_name == -1 ) 
	Helper::halt("required headers not observed: POS ID");
      if ( col_chr != -1 ) 
	Helper::halt("cannot specific POS and CHR/BP1/BP2");
    }
  else 
    {
      if ( col_chr == -1 || col_bp1 == -1 || col_bp2 == -1 || col_name == -1 ) 
	Helper::halt("required headers not observed: CHR BP1 BP2 ID, or POS ID");
    }

    

  //
  // All good to go...
  //

  return db->load_regions( filename , 
			   label , 
			   col_pos,
			   col_chr , col_bp1 , col_bp2 , 
			   col_name , col_sub , 
			   col_indiv , 
			   mf.size() == 0 ? NULL : &mf ) != 0;
  

}




bool Pseq::NetDB::loader( const std::string & db , const std::string & file )
{
  NetDBase netdb;
  netdb.attach( db );
  if ( ! netdb.attached() ) Helper::halt( "problem creating/attaching database" );
  netdb.load( file );
  return true;
}


