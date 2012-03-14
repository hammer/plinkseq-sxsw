#include "loaders.h"
#include "util.h"
#include "netdb.h"

extern GStore g;
extern Pseq::Util::Options args;
extern Log plog;

bool Pseq::VarDB::load_VCF( Mask & mask )
{

    // are any optinal locus groups specified (this will act as a mask, such that
    // we only read in variants that fall in this region

    bool region_mask = args.has("filter");
    
    std::string rmask;
    if ( region_mask ) rmask = args.as_string( "filter" );

    // filters on meta-fields?

    std::set<std::string> includes, excludes;
    if ( args.has( "format" , "include-meta" ) ) includes = args.get_set( "format" , "include-meta" );
    if ( args.has( "format" , "exclude-meta" ) ) excludes = args.get_set( "format" , "exclude-meta" );

    if ( ! args.has("check-reference") ) g.seqdb.dettach();

    return g.vardb_load_vcf( mask , includes , excludes , region_mask ? & rmask : NULL );

}

bool Pseq::RefDB::load_VCF( const std::string & filename , const std::string & grp ) 
{

  bool region_mask = args.has( "filter" );
  
  std::set<Region> filter;
  if ( region_mask ) filter = g.locdb.get_regions( args.as_string( "filter" ) );  
  std::set<Region> * pfilter = region_mask ? &filter : NULL  ;
  
  // filters on meta-fields?
  std::set<std::string> includes, excludes;
  if ( args.has( "format" , "include-meta" ) ) includes = args.get_set( "format" , "include-meta" );
  if ( args.has( "format" , "exclude-meta" ) ) excludes = args.get_set( "format" , "exclude-meta" );
  
  std::string comment;
  if ( args.has( "description" ) ) comment = args.as_string( "description" );
  
  if ( ! args.has( "check-reference" ) ) g.seqdb.dettach();

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
  if ( args.has("iid") ) b.use_iid();
  else if ( args.has("fid") ) b.use_fid();
  b.set_root( fileroot );

  if ( args.has("fix-strand") && g.seqdb.attached() ) b.fix_strand();  
  else if ( ! args.has("check-reference") ) g.seqdb.dettach();
  
  if ( g.seqdb.attached() ) b.use_seqdb( &g.seqdb );
  b.store_phenotypes( &g.inddb, args.has("phenotype") ? args.as_string("phenotype") : "phe1" );
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
	
  if ( args.has("name") ) meta[ PLINKSeq::SEQDB_NAME_KEY() ] = Pseq::Util::single_argument<std::string>( args , "name" ); 
  else Helper::halt("no name specified");

  if ( args.has("format", "build") ) meta[ PLINKSeq::SEQDB_GENOME_BUILD_KEY() ] = args.as_string( "format" , "build" );
  else Helper::halt("need to specify --format build=[build]");
  
  if ( args.has("format", "repeat-mode") )
    {

	std::string m = args.as_string( "format", "repeat-mode" );
      
	// should be none  = none
	//           N     = N-masked
	//           lower = lower-case
	
	if ( m != "none" && m != "N" && m != "lower" )
	    Helper::halt("repeat-mode should be 'none', 'N' or 'lower'");
      
	meta[ PLINKSeq::SEQDB_REPEAT_MODE_KEY() ] = m;
    }
  else Helper::halt("no repeat-mode specified in options");
  
  if ( args.has( "description" ) ) meta[ PLINKSeq::SEQDB_DESCRIPTION_KEY() ] = args.as_string( "description" );
  else Helper::halt("no description specified in options");
  
  if ( args.has( "format", "iupac" ))
    meta[ PLINKSeq::SEQDB_IUPAC_KEY() ] = "1";
  else
    meta[ PLINKSeq::SEQDB_IUPAC_KEY() ] = "0";
  
  if (!args.has( "seqdb" ))
     Helper::halt("Need to specify SEQDB path with --seqdb");
  std::string seqdb_filename = args.as_string( "seqdb" );
  g.seqdb.create( seqdb_filename );

  plog << "loading from FASTA..\n";
  
  g.seqdb.loadFASTA( filename , meta );

  return true;
}


bool Pseq::RefDB::load_refvar( const std::string & filename , 
			       const std::string & label , 
			       Pseq::Util::Options & args )
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
  std::string f2_chr = "CHROM";

  std::string f_bp1 = "BP1";
  std::string f2_bp1 = "POS";
  std::string f3_bp1 = "POS1";   

  std::string f_bp2 = "BP2";		       
  std::string f2_bp2 = "POS2";		       

  std::string f_id  = "ID";
  
  // allow command-line swap-ins

  if ( args.has( "format" , "chr" ) ) f_chr = args.as_string( "format" , "chr" );
  if ( args.has( "format" , "bp1" ) ) f_bp1 = args.as_string( "format" , "bp1" );
  if ( args.has( "format" , "bp2" ) ) f_bp2 = args.as_string( "format" , "bp2" );
  if ( args.has( "format" , "id" ) )  f_id  = args.as_string( "format" , "id" );

  for (int i=0; i<h.size(); i++)
    {

      // fixed fields?
      if ( h[i] == f_chr || h[i] == f2_chr )                          col_chr   = i; 
      else if ( h[i] == f_bp1 || h[i] == f2_bp1 || h[i] == f3_bp1  )  col_bp1   = i;
      else if ( h[i] == f_bp2 || h[i] == f2_bp2 )                     col_bp2   = i;
      else if ( h[i] == f_id    )                                     col_name  = i;
      else // a user-defined type
	{
	    
	  // Should we skip this field?
	  // --skip col
	  
	  std::set<std::string> s = args.get_set( "format" , "skip" );
	  if ( s.find( h[i] ) != s.end() )
	    {
	      skip.insert(i);
	      continue;
	    }
	  
	  
	  // If this is a "Flag" field, it means that values
	  // themselves will be the key, so we do not register the
	  // column name with the meta-information class.

	    // --flags col2
	  
	  std::set<std::string> fs = args.get_set( "format" , "flag" );
	  if ( fs.find( h[i] ) != fs.end() ) flags.insert(i);
	  else
	    {
	      mType mt = META_TEXT;	  
	      if      ( args.has( "format" , "integer" , h[i] ) ) mt = META_INT;
	      else if ( args.has( "format" , "float"   , h[i] ) ) mt = META_FLOAT;
	      else if ( args.has( "format" , "bool"    , h[i] ) ) mt = META_BOOL;
	      else if ( args.has( "format" , "text"    , h[i] ) ) mt = META_TEXT;
	      else if ( args.has( "format" , "string"  , h[i] ) ) mt = META_TEXT;
	      
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

  bool zero1 = args.has( "format" , "0-start") || args.has( "format" , "0-based");
  bool zero2 = args.has( "format" , "0-stop")  || args.has( "format" , "0-based");
  

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


bool Pseq::LocDB::load_segments( std::string filename , std::string label , Pseq::Util::Options & args )
{

  if ( ! g.segdb.attached() ) 
    Helper::halt( "no SEGDB attached" );

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
  int col_meta = -1;



  std::map<std::string,int> mf;  // meta-field information
  
  for (int i=0; i<h.size(); i++)
    {
      
      // fixed fields?
      if ( h[i] == "CHR" || h[i] == "CHROM" ) col_chr = i; 
      else if ( h[i] == "BP1" || h[i] == "POS1" ) col_bp1 = i;
      else if ( h[i] == "BP2" || h[i] == "POS2" ) col_bp2 = i;
      else if ( h[i] == "SEGID" ) col_name = i;
      else if ( h[i] == "ID" ) col_indiv = i;
      else if ( h[i] == "META" ) col_meta = i;
      else // a user-defined type
	{

	  mType mt = META_TEXT;

	  if      ( args.has( "integer" , h[i] ) ) mt = META_INT;
	  else if ( args.has( "float"   , h[i] ) ) mt = META_FLOAT;
	  else if ( args.has( "flag"    , h[i] ) ) mt = META_FLAG;
	  else if ( args.has( "bool"    , h[i] ) ) mt = META_BOOL;
	  
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
			       col_meta ,
			       col_indiv , 
			       mf.size() == 0 ? NULL : &mf ) != 0;
  
}


bool Pseq::LocDB::update_searchtable( const std::string & group , bool locdb )
{
  LocDBase * db = locdb ? &g.locdb : &g.segdb ;
  if (!db) return false;
  db->populate_searchname_table( group );
  return true;
}


//
// Load generic regions into LOCDB or SEGDB
//

bool Pseq::LocDB::load_generic_regions( std::string & filename , const std::string & label , Pseq::Util::Options & args , bool locdb )
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
  int col_meta = -1;

  std::map<std::string,int> mf;  // meta-field information
  
  for (int i=0; i<h.size(); i++)
    {

      // fixed fields?
      if ( h[i] == "CHR" || h[i] == "CHROM" ) col_chr = i; 
      else if ( h[i] == "BP1" || h[i] == "POS1" ) col_bp1 = i;
      else if ( h[i] == "BP2" || h[i] == "POS2" ) col_bp2 = i;
      else if ( h[i] == "POS" ) col_pos = i;
      else if ( h[i] == "ID" ) col_name = i;
      else if ( h[i] == "META" ) col_meta = i;
      else // a user-defined type
	{
	  mType mt = META_TEXT;
	  if      ( args.has( "integer" , h[i] ) )  mt = META_INT;
	  else if ( args.has( "float"   , h[i] ) )  mt = META_FLOAT;
	  else if ( args.has( "flag"    , h[i] ) )  mt = META_FLAG;
	  else if ( args.has( "bool"    , h[i] ) )  mt = META_BOOL;

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

  bool state = db->load_regions( filename , 
				 label , 
				 col_pos,
				 col_chr , col_bp1 , col_bp2 , 
				 col_name , col_sub , 
				 col_meta , 
				 col_indiv , 
				 mf.size() == 0 ? NULL : &mf ) != 0 ;
  
  return state;
}




bool Pseq::NetDB::loader( const std::string & db , const std::string & file )
{
  NetDBase netdb;
  netdb.attach( db );
  if ( ! netdb.attached() ) Helper::halt( "problem creating/attaching database" );
  netdb.load( file );
  return true;
}



struct aux_addvar  {
  std::string group;
  std::string mtag;
  bool tagvalue;  
  std::set<std::string> sets;
};


void f_add_to_varset( Variant & var , void * p )
{
  aux_addvar * aux = (aux_addvar*)p;
  
  if ( aux->tagvalue )
    {
      // to to group-name defined by tag-value
      std::string val = "";
      
      // only allow this for text and integer tag values
      meta_index_t midx = MetaInformation<VarMeta>::field( aux->mtag );
      if ( !(  midx.mt == META_TEXT || midx.mt == META_INT ) ) return;

      MetaInformation<VarMeta> * m = NULL;
      if ( var.meta.has_field( aux->mtag ) ) m = &var.meta;
      else if ( var.consensus.meta.has_field( aux->mtag ) ) m = &var.consensus.meta; 
      if ( m == NULL ) return;
      
      if ( midx.mt == META_TEXT )
	{
	  std::vector<std::string> t = m->get_string( aux->mtag );
	  for (int i=0; i<t.size(); i++ ) 
	    {
	      const std::string n = aux->mtag + "[" + t[i] + "]";
	      g.vardb.add_var_to_set( n , var  );
	      if ( aux->sets.find( n ) == aux->sets.end() )
		{
		  g.vardb.add_set_to_superset( aux->group , n );
		  aux->sets.insert( n );
		}
	    }
	}
      
      if ( midx.mt == META_INT )
	{
	  std::vector<int> t = m->get_int( aux->mtag );
	  for (int i=0; i<t.size(); i++ ) 
	    {
	      const std::string n = aux->mtag + "[" + Helper::int2str( t[i] ) + "]";
	      g.vardb.add_var_to_set( n , var  );
	      if ( aux->sets.find( n ) == aux->sets.end() )
		{
		  g.vardb.add_set_to_superset( aux->group , n );
		  aux->sets.insert( n );
		}
	    }
	}
    }
  else
    g.vardb.add_var_to_set( aux->group , var );
}


bool Pseq::VarDB::add_to_varset( const std::string & group , Mask & mask , const std::string & mtag , const std::string & desc )
{

  // Create set if it does not exist

  // if setting sets based on the value of a meta-field, then create an umbrella super-set 
  // based on the --group given
  
  if ( mtag == "" ) g.vardb.add_set( group );
  else g.vardb.add_superset( group );

  if ( desc != "." ) g.vardb.add_set_description( group , desc );

  // Add all mask-passing variants

  aux_addvar av;
  av.group = group;
  if ( mtag != "" ) { av.tagvalue = true; av.mtag = mtag; } 
  g.vardb.begin();
  g.vardb.iterate( f_add_to_varset , &av , mask );
  g.vardb.commit();

  return true;
}


bool Pseq::VarDB::add_to_varset( const std::string & filename  )
{

  // format ( tab-delim, 2/3/4 cols ) : 

  //  #name        description
  //  snp/region  var-pos   name   { alternate-allele } 
       

  Helper::checkFileExists( filename );

  InFile f( filename );

  bool vardb_exists = g.vardb.attached();
  
  while ( ! f.eof() )
    {

      std::vector<std::string> h = f.tokenizeLine("\t");
      
      // a description / add a var-set name 
      if ( h.size() == 2 ) 
	{
	  g.vardb.add_set_description( h[0] , h[1] ); 
	}
      else 
	{
	  if ( h.size() != 3 && h.size() != 4 ) continue;

	  bool add_region = false;
	  Helper::str2upper( h[0] );
	  
	  if ( h[0].substr(0,3) == "REG" ) add_region = true;
	  else if ( h[0].substr(0,3) == "SNP" ) add_region = false;
	  else 
	    {
	      plog.warn( "invalid region/snp code " , h[0] );
	      continue; 
	    }

	  // add a set of SNPs in a region, or as a single variant ? 
	  
	  if ( add_region ) 
	    {

	      // this only works when a VARDB is attached
	      if ( ! vardb_exists ) continue;
	      
	      // REG name var
	      
	      bool okay = true;
	      Region reg( h[2] , okay );
	      if ( ! okay ) continue;
	      
	      std::string group = h[1];

	      std::set<Variant> vs = g.vardb.key_fetch( reg );
	      std::set<Variant>::iterator vi = vs.begin();
	      while ( vi != vs.end() ) 
		{		  
		  g.vardb.add_var_to_set( group , (Variant&)*vi );
		  ++vi;
		}
	      
	    }
	  else
	    {
	      
	      // SNP  name  var {alt}
	      
	      // add variant to 
	      Variant v;
	      bool okay = true;
	      Region reg( h[2] , okay );
	      if ( ! okay ) continue;

	      // if a VARDB attached, only add SNP if it exists in the VARDB
	      if ( vardb_exists ) 
		{
		  std::set<Variant> vs = g.vardb.key_fetch( reg );
		  if ( vs.size() == 0 ) continue;
		  v.consensus.index( vs.begin()->consensus.index() );
		}
	      
	      v.chromosome( reg.chromosome() );
	      v.position( reg.start.position() );
	      v.stop( reg.stop.position() );

	      // also specified an alternate allele?
	      if ( h.size() == 4 ) 
		v.consensus.alternate( h[3] );
	      
	      // if a VARDB is attached, only add this variant to the
	      // set if it actually exists in 1 or more files
	      
	      // add SNP (group-name, variant)
	      g.vardb.add_var_to_set( h[1] , v );
	      
	    } // end add-region
	} // end add-snp
    } // next row
  
  f.close();    
  return true;
}

bool Pseq::VarDB::add_superset_from_file( const std::string & filename )
{
  Helper::checkFileExists( filename );
  InFile f( filename );
  while ( ! f.eof() )
    {
      std::string superset, set;
      f >> superset >> set ;
      if ( superset == "" ) continue;
      g.vardb.add_set_to_superset( superset , set );
    }
  f.close();
  return true;
}

bool Pseq::VarDB::add_superset( const std::string & group , const std::vector<std::string> & members , const std::string & desc )
{
  g.vardb.add_superset( group );
  
  //  g.vardb.add_superset_description( group , desc );

  for (int m=0;m<members.size(); m++) g.vardb.add_set_to_superset( group , members[m] );
  return true;
}


// Insert meta-information, either from a file or from a new meta-tag created on-the-fly

bool Pseq::VarDB::insert_meta_from_file( const std::string & filename )
{

}

bool Pseq::VarDB::insert_meta_on_fly( const std::string & name )
{

}

