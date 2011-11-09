
#include "vardb.h"
#include "genotype.h"
#include "regions.h"
#include "gstore.h"
#include "em.h"
#include "regions.h"
#include "sqlz.h"

extern GStore * GP;

bool VarDBase::newDB( std::string n )
{
  
  sql.open(n);
  
  // register compression functions 
  
  sqlite3_create_function( sql.pointer(), "mycompress", 1, SQLITE_UTF8, 0, &compressFunc, 0, 0);
  sqlite3_create_function( sql.pointer(), "myuncompress", 1, SQLITE_UTF8, 0, &uncompressFunc, 0, 0);
  
  sql.synchronous(false);

  sql.query("PRAGMA encoding='UTF-8'");
  
  using_compression = true;

  //
  // DB version, and a place for various other meta-information in future
  //
  
  sql.query(" CREATE TABLE IF NOT EXISTS dbmeta("
	    "   varname      VARCHAR(20) NOT NULL , "
	    "   varvalue    VARCHAR(20) NOT NULL , "
	    " CONSTRAINT uMeta UNIQUE (varname ) ); " );  
  

  //
  // Headers 
  //
  
  sql.query(" CREATE TABLE IF NOT EXISTS headers("
	    "   header_id INTEGER PRIMARY KEY , "
	    "   file_id  INTEGER NOT NULL , "
	    "   name     VARCHAR(20) , "
	    "   value    VARCHAR(20) ) ; " );
  

  //
  // Meta-value type information
  //

  sql.query(" CREATE TABLE IF NOT EXISTS metatypes("
	    "   file_id      INTEGER NOT NULL , "
	    "   name         VARCHAR(20) , "
	    "   type         VARCHAR(8) , "
	    "   number       INTEGER , "
	    "   grp          INTEGER , "
	    "   description  VARCHAR(20) , "
	    " CONSTRAINT mtCon UNIQUE (file_id,name,grp) ); " );
  

  //
  // Files/summary
  //
  
  sql.query(" CREATE TABLE IF NOT EXISTS files("
	    "   file_id  INTEGER PRIMARY KEY , "
	    "   name     VARCHAR(20) NOT NULL , "
	    "   tag      VARCHAR(20) , "
	    "   ni       INTEGER , "
	    "   nv       INTEGER ); " );

  
  //
  // Chromosome codes, and populate with basics
  //
  
  bool populate_chr_codes = ! sql.table_exists( "chrcodes" );
  
  sql.query( " CREATE TABLE IF NOT EXISTS chrcodes("
             "   chr_id  INTEGER PRIMARY KEY , "
             "   name    VARCHAR(20) NOT NULL , "
	     "   ploidy  INTEGER NOT NULL , "
             " CONSTRAINT cChr UNIQUE ( name ) ) ; " );
  
  
  //
  // We only want to do this once (otherwise write to db will cause lock)
  //
  
  if ( populate_chr_codes ) 
    {
      
      stmt_insert_chr_code = 
	sql.prepare( " INSERT OR REPLACE INTO chrcodes ( name , ploidy ) values( :name , :ploidy ); " );      
      stmt_fix_chr_code = 
	sql.prepare( " INSERT OR REPLACE INTO chrcodes ( chr_id, name , ploidy ) values( :chr_id , :name , :ploidy ); " );      
      stmt_fetch_chr_code = 
	sql.prepare( " SELECT chr_id , ploidy FROM chrcodes WHERE name == :name ; " );      
      stmt_fetch_chr_name = 
	sql.prepare( " SELECT name , ploidy FROM chrcodes WHERE chr_id == :chr_id ; " );

            
      chr_code( 1, "chr1"  , PLOIDY_AUTOSOMAL );  
      chr_code( 2, "chr2"  , PLOIDY_AUTOSOMAL );  
      chr_code( 3, "chr3"  , PLOIDY_AUTOSOMAL );  
      chr_code( 4, "chr4"  , PLOIDY_AUTOSOMAL );
      chr_code( 5, "chr5"  , PLOIDY_AUTOSOMAL );  
      chr_code( 6, "chr6"  , PLOIDY_AUTOSOMAL );  
      chr_code( 7, "chr7"  , PLOIDY_AUTOSOMAL );  
      chr_code( 8, "chr8"  , PLOIDY_AUTOSOMAL );
      chr_code( 9, "chr9"  , PLOIDY_AUTOSOMAL );  
      chr_code(10, "chr10" , PLOIDY_AUTOSOMAL );  
      chr_code(11, "chr11" , PLOIDY_AUTOSOMAL );  
      chr_code(12, "chr12" , PLOIDY_AUTOSOMAL );
      chr_code(13, "chr13" , PLOIDY_AUTOSOMAL );  
      chr_code(14, "chr14" , PLOIDY_AUTOSOMAL );  
      chr_code(15, "chr15" , PLOIDY_AUTOSOMAL );  
      chr_code(16, "chr16" , PLOIDY_AUTOSOMAL );
      chr_code(17, "chr17" , PLOIDY_AUTOSOMAL );  
      chr_code(18, "chr18" , PLOIDY_AUTOSOMAL );  
      chr_code(19, "chr19" , PLOIDY_AUTOSOMAL );  
      chr_code(20, "chr20" , PLOIDY_AUTOSOMAL );
      chr_code(21, "chr21" , PLOIDY_AUTOSOMAL );  
      chr_code(22, "chr22" , PLOIDY_AUTOSOMAL );       
      chr_code(23, "chrX"  , PLOIDY_X ); 
      chr_code(24, "chrY"  , PLOIDY_Y ); 
      chr_code(25, "chrM"  , PLOIDY_HAPLOID );      
      
      sql.finalise( stmt_insert_chr_code );
      sql.finalise( stmt_fix_chr_code );
      sql.finalise( stmt_fetch_chr_code );
      sql.finalise( stmt_fetch_chr_name );
    }


  // Individuals
  
  sql.query( " CREATE TABLE IF NOT EXISTS individuals("
	     "   indiv_id INTEGER PRIMARY KEY , "
	     "   file_id  INTEGER NOT NULL , "
	     "   name     VARCHAR(20) , "
	     " CONSTRAINT cInd UNIQUE ( file_id, indiv_id ) ); ");
  

  // Variants
  
  sql.query(" CREATE TABLE IF NOT EXISTS variants("
	    "   var_id   INTEGER PRIMARY KEY , "
	    "   file_id  INTEGER NOT NULL , "
	    "   name     VARCHAR(40) , "
	    "   chr      INTEGER NOT NULL , "
	    "   bp1      INTEGER NOT NULL , "
	    "   bp2      INTEGER  , "
	    "   offset   INTEGER ) ; " );

  // Variant data-store
  
  sql.query(" CREATE TABLE IF NOT EXISTS vdat("
	    "   meta_id INTEGER PRIMARY KEY , "
	    "   var_id  INTEGER NOT NULL , "
	    "   data    BLOB , "      // basic variant data (alleles, etc)
	    "   vdata   BLOB , "      // variant meta-information
	    "   gdata   BLOB , "      // genotype calls
	    "   gmdata  BLOB ); " );  // genotype meta-information

  
  // file-id <--> BCF and BGZF-compressed VCFs
  //  type, 1 = VCF w/ BGZF compression
  //        2 = BCF

  sql.query(" CREATE TABLE IF NOT EXISTS bcfs("
	    "   file_id    INTEGER NOT NULL , "
	    "   type       INTEGER NOT NULL , "
	    "   filepath   VARCHAR(20) NOT NULL , "
            "   nind       INTEGER  ) ; " );
    
  // Independent meta-information (appended outside of main VCF)
  // Not automatically attached to the variant
    
  sql.query(" CREATE TABLE IF NOT EXISTS indep_meta_groups("
	    "   group_id    INTEGER PRIMARY KEY , "
	    "   file_id     INTEGER NOT NULL , "
	    "   name        VARCHAR(20) NOT NULL , "
	    " CONSTRAINT c1 UNIQUE ( file_id , name ) ) ; ");
  
  
     sql.query(" CREATE TABLE IF NOT EXISTS indep_meta_types("
 	      "   meta_id     INTEGER PRIMARY KEY , "
 	      "   name        VARCHAR(20) NOT NULL , "
 	      "   length      INTEGER NOT NULL , "
 	      "   type        INTEGER NOT NULL , "
 	      "   desc        VARCHAR(20) , "
               " CONSTRAINT c2 UNIQUE ( name ) ) ; ");
  
  
     sql.query(" CREATE TABLE IF NOT EXISTS indep_meta_data("
	       " group_id     INTEGER NOT NULL , "
	       "   var_id     INTEGER NOT NULL , "
	       "   meta_id    INTEGER NOT NULL , "
               "   value      VARCHAR(20) , "
               " CONSTRAINT c3 UNIQUE ( var_id, meta_id ) ); " );
     

    
    // Sets
    
    sql.query(" CREATE TABLE IF NOT EXISTS set_groups("
	      "   group_id    INTEGER PRIMARY KEY , "
	      "   name        VARHCAR(20) , "
	      "   temp        CHAR(1) , "
	      "   description VARCHAR(20) ) ; " );
    
    sql.query(" CREATE TABLE IF NOT EXISTS set_members("
	      "   set_id      INTEGER PRIMARY KEY , "
	      "   group_id    INTEGER NOT NULL , "
	      "   name        VARCHAR(20) ) ; " );
    
    sql.query(" CREATE TABLE IF NOT EXISTS set_data("
	      "   var_id      INTEGER NOT NULL , "
	      "   set_id      INTEGER NOT NULL ) ; " );
    

    // Super-sets 

    sql.query(" CREATE TABLE IF NOT EXISTS superset_groups("
	      "   group_id    INTEGER PRIMARY KEY , "
	      "   set_group_id INTEGER NOT NULL , "
	      "   name        VARHCAR(20) , "
	      "   temp        CHAR(1) , "
	      "   description VARCHAR(20) ) ; " );

    sql.query(" CREATE TABLE IF NOT EXISTS superset_members("
	      "   superset_id  INTEGER PRIMARY KEY , "
	      "   group_id     INTEGER NOT NULL , "
	      "   name         VARCHAR(20) ) ; " );

    sql.query(" CREATE TABLE IF NOT EXISTS superset_data("
	      "   superset_id  INTEGER NOT NULL , "
	      "   set_id       INTEGER NOT NULL ) ; " );


    //
    // i.e. relationship between sets and super-sets
    
    // Set = CCDS
    //  set members -> gene names
    //  set data    -> which variants are in which genes

    // Super-set = KEGG 
    //   based on CCDS ( linked by set_group_id )
    //    set members --> pathway names
    //    set data    --> which genes (set_id) are in which pathway (superset_id)


  return true;

}


bool VarDBase::attach( std::string n )
{

  if ( n == "." ) { dettach(); return false; } 
  
  //
  // Close any existing database
  //
  
  if ( attached() ) release();


  //
  // If DB already exists, this function will do nothing; 
  // otherwise, create the DB
  //

  newDB(n);
     
  
  //
  // Prepared statements
  //

  init();


  //
  // Check version number 
  //

  check_version();


  //
  // Load new meta-information table
  //

  set_metatypes( ); 
  
  populate_indep_metadata_map();
  
  populate_bcf_map();

  return true;
}


bool VarDBase::dettach()
{
  if ( attached() )
    {
      release();
      sql.close();
    }
  return true;
}


bool VarDBase::init()
{
  
    // 
    // Insertions
    // 
    
    stmt_insert_header = 
	sql.prepare(" INSERT OR IGNORE INTO headers ( file_id, name, value ) "
		    " values( :file_id, :name, :value ); " );

    // Chromosome codes

    stmt_insert_chr_code = 
      sql.prepare( " INSERT OR REPLACE INTO chrcodes ( name , ploidy ) values( :name , :ploidy ); " );

    stmt_fix_chr_code = 
      sql.prepare( " INSERT OR REPLACE INTO chrcodes ( chr_id, name , ploidy ) values( :chr_id , :name , :ploidy ); " );
    
    stmt_fetch_chr_code = 
      sql.prepare( " SELECT chr_id , ploidy FROM chrcodes WHERE name == :name ; " );
    
    stmt_fetch_chr_name = 
      sql.prepare( " SELECT name , ploidy FROM chrcodes WHERE chr_id == :chr_id ; " );


    // Meta-value type information
    
    stmt_insert_metatype = 
	sql.prepare(" INSERT OR IGNORE INTO metatypes "
		    " ( file_id, name, type, number, grp, description ) "
		    " values ( :file_id , :name , :type, :number, :group, :description ) ; " );
    
    stmt_insert_file = 
	sql.prepare(" INSERT OR IGNORE INTO files ( name , tag ) values ( :name , :tag ) ; " );
    
    stmt_insert_file_summary =
	sql.prepare(" UPDATE files "
		    " SET ni = :ni, nv = :nv "
		    " WHERE name == :name ; ");

    // BCFs

    stmt_insert_bcf_n = 
      sql.prepare( " INSERT OR IGNORE INTO bcfs ( file_id , type , filepath , nind ) "
                  " values ( :file_id , :type , :filepath, :nind ) ; " );
    
    stmt_fetch_bcf = 
      sql.prepare( " SELECT filepath FROM bcfs WHERE file_id == :file_id ; " );

    stmt_fetch_bcfs = 
      sql.prepare( " SELECT * FROM bcfs; " );
    
    stmt_insert_bcf_idx = 
      sql.prepare(" INSERT OR IGNORE INTO variants "
		  "          ( file_id, name , chr, bp1 , bp2 , offset ) "
		  "   values ( :file_id, :name , :chr, :bp1 , :bp2, :offset ) ; " );
    

    // File tags

    stmt_fetch_file_from_tag =
      sql.prepare( " SELECT file_id FROM files WHERE tag == :tag; " );

    stmt_fetch_tag_from_file =
      sql.prepare( " SELECT tag FROM files WHERE file_id == :file_id; " );

    stmt_insert_file_tag =
      sql.prepare( " UPDATE files SET tag = :tag WHERE file_id == :file_id ; " );

    stmt_insert_variant_key = 
	sql.prepare(" INSERT OR IGNORE INTO variants "
		    "          ( file_id, name, chr, bp1 , bp2 ) "
		    "   values ( :file_id, :name, :chr, :bp1 , :bp2 ) ; " );
    
    if ( using_compression ) 
      {

	stmt_insert_variant_data = 
	  sql.prepare(" INSERT OR IGNORE INTO vdat "
		      "          ( var_id, data , vdata , gdata , gmdata ) "
		      "   values ( :var_id , :data , :vdata , mycompress( :gdata ) , mycompress( :gmdata ) ) ; " );	    
      }
    else
      stmt_insert_variant_data = 
	sql.prepare(" INSERT OR IGNORE INTO vdat "
		    "          ( var_id, data , vdata , gdata , gmdata ) "
		    "   values ( :var_id, :data , :vdata , :gdata , :gmdata ) ; " );
    
    stmt_insert_individual = 
      sql.prepare(" INSERT OR IGNORE INTO individuals "
		  "          ( file_id, indiv_id, name ) "
		  "   values ( :file_id, :indiv_id, :name ); " );
    
    
    //
    // Queries
    //

    // Var-ID lookups

    stmt_fetch_var_from_position = 
      sql.prepare(" SELECT var_id,file_id FROM variants WHERE chr == :chr AND bp1 == :bp1; " );
    
    stmt_fetch_var_from_position2 = 
      sql.prepare(" SELECT var_id,file_id FROM variants WHERE chr == :chr AND bp1 == :bp1 AND bp2 == :bp2; " );

    stmt_fetch_var_from_name = 
      sql.prepare(" SELECT var_id,file_id FROM variants WHERE name == :name ; " );


    //
    // Main data queries
    //

    stmt_fetch_files = 
	sql.prepare(" SELECT file_id, name FROM files "
		    " ORDER BY file_id; " );

    stmt_fetch_file_id = 
	sql.prepare(" SELECT file_id FROM files "
		    " WHERE name == :name; " );

    stmt_fetch_file_summary = 
	sql.prepare(" SELECT ni,nv FROM files "
		    " WHERE file_id == :file_id; " );

    stmt_fetch_headers = 
	sql.prepare(" SELECT name, value FROM headers "
		    " WHERE file_id == :file_id ORDER BY header_id; " );

    stmt_fetch_metatypes = 
	sql.prepare(" SELECT name , type , number, grp, description "
		    " FROM metatypes WHERE file_id == :file_id ; " );
        
    stmt_fetch_variant_key = 
	sql.prepare(" SELECT * FROM variants WHERE var_id == :var_id ; " );

    stmt_fetch_variant_pos = 
      sql.prepare(" SELECT * FROM variants WHERE chr == :chr AND bp1 == :bp1 ;" );

    stmt_fetch_variant_range = 
      sql.prepare(" SELECT * FROM variants WHERE chr == :chr AND bp1 >= :rstart AND bp1 <= :rend;" );

    if ( using_compression )
      {
	stmt_fetch_variant_data_all = 
	  sql.prepare(" SELECT data , vdata , myuncompress(gdata) , myuncompress(gmdata) FROM vdat "
		      " WHERE var_id == :var_id ; " );
	
	stmt_fetch_variant_data_vmeta_geno = 
	  sql.prepare(" SELECT data , vdata , myuncompress(gdata) FROM vdat "
		      " WHERE var_id == :var_id ; " );
	
	stmt_fetch_variant_data_vmeta = 
	  sql.prepare(" SELECT data , vdata FROM vdat "
		      " WHERE var_id == :var_id ; " );

	stmt_fetch_variant_data_geno = 
	  sql.prepare(" SELECT data , gdata FROM vdat "
		      " WHERE var_id == :var_id ; " );

      }
    else
      {
	stmt_fetch_variant_data_all = 
	  sql.prepare(" SELECT data , vdata , gdata , gmdata FROM vdat WHERE var_id == :var_id ; " );

	stmt_fetch_variant_data_vmeta_geno = 
	  sql.prepare(" SELECT data , vdata , gdata FROM vdat WHERE var_id == :var_id ; " );

	stmt_fetch_variant_data_vmeta = 
	  sql.prepare(" SELECT data , vdata FROM vdat WHERE var_id == :var_id ; " );

	stmt_fetch_variant_data_geno = 
	  sql.prepare(" SELECT data , gdata FROM vdat WHERE var_id == :var_id ; " );

      }


    stmt_fetch_individual = 
	sql.prepare(" SELECT * FROM individuals WHERE indiv_id == :indiv_id  ; " );
    
    stmt_fetch_individuals = 
	sql.prepare(" SELECT * FROM individuals WHERE file_id == :file_id "
		    " ORDER BY indiv_id;");

    stmt_iterate_variants = 
	sql.prepare(" SELECT * FROM variants ORDER BY chr,bp1,bp2 ; " );

    stmt_fetch_sets = 
	sql.prepare(" SELECT group_id, name FROM set_groups "
		    " ORDER BY group_id; " );
    

    //
    // Aux. meta-data
    //


    stmt_insert_indep_meta_group = sql.prepare( " INSERT OR IGNORE INTO indep_meta_groups ( file_id, name ) "
						"  values( :file_id, :name ) ; " ); 
    
    stmt_fetch_indep_meta_group = sql.prepare( " SELECT group_id FROM indep_meta_groups "
					       " WHERE name == :name AND file_id == :file_id; " );
    
    stmt_dump_indep_meta_group = sql.prepare( " SELECT group_id,name,file_id FROM indep_meta_groups; " );


    stmt_insert_indep_meta_type = sql.prepare( " INSERT OR IGNORE INTO indep_meta_types ( name, length, type, desc ) "
					       "  values( :name , :length , :type , :desc ) ; " );
    
    stmt_fetch_indep_meta_type = sql.prepare( " SELECT * FROM indep_meta_types ; " );
    
    stmt_insert_indep_meta_value = sql.prepare( " INSERT OR IGNORE INTO indep_meta_data ( group_id, var_id, meta_id, value ) "
						"  values( :group_id, :var_id, :meta_id, :value ) ; " );
    
    stmt_fetch_indep_meta_value = sql.prepare( " SELECT meta_id,value FROM indep_meta_data "
					       " WHERE var_id == :var_id ; " );


    //
    // Sets
    //

    stmt_insert_group =
	sql.prepare( " INSERT OR IGNORE INTO set_groups( name, temp, description ) "
		     " values( :name, :temp, :description) ; " );
    stmt_lookup_group =
	sql.prepare( " SELECT group_id FROM set_groups WHERE name == :name ; " );

    stmt_lookup_group_name =
	sql.prepare( " SELECT name FROM set_groups WHERE group_id == :group_id ; " );

    stmt_insert_group_member =
	sql.prepare( " INSERT OR IGNORE INTO set_members ( group_id , name ) "
		     " values( :group_id , :name ) ; " );

    stmt_insert_group_variant =
	sql.prepare( " INSERT OR IGNORE INTO set_data ( set_id , var_id ) "
		     " values( :set_id , :var_id ) ; " );
    
    stmt_iterate_group = 
	sql.prepare( " SELECT * FROM variants "
		     " WHERE var_id IN "
		     "  ( SELECT var_id FROM set_data WHERE set_id IN "
		     "  ( SELECT set_id FROM set_members WHERE group_id == :group_id ) ) "
		     " ORDER BY chr,bp1,bp2 ; " );
	
    stmt_fetch_set_names1 = 
	sql.prepare( " SELECT name FROM set_members WHERE set_id IN "
		     "  ( SELECT set_id FROM set_data WHERE var_id == :var_id ) ; " ); 

    stmt_fetch_set_names2 = 
	sql.prepare( " SELECT name FROM set_members "
		     " WHERE group_id == :group_id AND set_id IN "
		     "  ( SELECT set_id FROM set_data WHERE var_id == :var_id ) ; " ); 
    
    
    //
    // Misc
    //
    
    stmt_vcount = 
	sql.prepare("SELECT count(*) FROM variants WHERE file_id == :file_id; ");

    stmt_indcount = 
	sql.prepare("SELECT count(*) FROM individuals WHERE file_id == :file_id; ");

    stmt_totvcount = 
	sql.prepare("SELECT count(*) FROM ( SELECT DISTINCT chr,bp1 FROM variants ) ; ");

    stmt_setcount = 
	sql.prepare("SELECT count(*) FROM set_members WHERE group_id == :group_id; ");

  return true;

} 

void VarDBase::attachMemoryDB()
{
    sql.query(" ATTACH \":memory:\" AS tmp ; " );
    sql.query(" CREATE TABLE tmp.tbl ( name VARCHAR(20) ) ; " );
    stmt_tmp_insert = 
	sql.prepare( " INSERT INTO tmp.tbl ( name ) values ( :name ) ; " );
}


void VarDBase::insertMemoryDB(const std::string & name)
{
  sql.bind_text( stmt_tmp_insert , ":name" , name );
  sql.step( stmt_tmp_insert );
  sql.reset( stmt_tmp_insert );
}

void VarDBase::detachMemoryDB()
{
  if ( attached() )
    {
      sql.query(" DETACH DATABASE tmp; ");    
      sql.finalise( stmt_tmp_insert );
    }
}

bool VarDBase::release()
{
  sql.finalise( stmt_dump_indep_meta_group );
  sql.finalise( stmt_fetch_chr_code );
  sql.finalise( stmt_fetch_chr_name );
  sql.finalise( stmt_fetch_file_from_tag );
  sql.finalise( stmt_fetch_file_summary  );
  sql.finalise( stmt_fetch_tag_from_file );
  sql.finalise( stmt_fix_chr_code );
  sql.finalise( stmt_insert_chr_code );
  sql.finalise( stmt_insert_file_summary );
  sql.finalise( stmt_insert_file_tag );
  sql.finalise( stmt_lookup_group_name );    
  
  sql.finalise( stmt_insert_header );
  sql.finalise( stmt_insert_metatype ); 
  sql.finalise( stmt_insert_file ); 
  sql.finalise( stmt_insert_variant_key ); 
  sql.finalise( stmt_insert_variant_data ); 
  sql.finalise( stmt_insert_individual ); 
  
  sql.finalise( stmt_fetch_var_from_position );
  sql.finalise( stmt_fetch_var_from_position2 );
  sql.finalise( stmt_fetch_var_from_name );
  
  sql.finalise( stmt_fetch_headers ); 
  sql.finalise( stmt_fetch_metatypes ); 
  sql.finalise( stmt_fetch_variant_key ); 
  sql.finalise( stmt_fetch_variant_pos ); 
  sql.finalise( stmt_fetch_variant_range ); 

  sql.finalise( stmt_fetch_variant_data_all ); 
  sql.finalise( stmt_fetch_variant_data_vmeta_geno ); 
  sql.finalise( stmt_fetch_variant_data_vmeta ); 
  sql.finalise( stmt_fetch_variant_data_geno ); 

  sql.finalise( stmt_fetch_files );
  sql.finalise( stmt_fetch_file_id );
  
  sql.finalise( stmt_insert_bcf_n );
  sql.finalise( stmt_fetch_bcf );
  sql.finalise( stmt_fetch_bcfs );
  sql.finalise( stmt_insert_bcf_idx ); 
  
  sql.finalise( stmt_insert_indep_meta_group );
  sql.finalise( stmt_fetch_indep_meta_group );
  sql.finalise( stmt_insert_indep_meta_type );
  sql.finalise( stmt_fetch_indep_meta_type );
  sql.finalise( stmt_insert_indep_meta_value );
  sql.finalise( stmt_fetch_indep_meta_value );

  sql.finalise( stmt_fetch_individual ); 
  sql.finalise( stmt_fetch_individuals ); 
  sql.finalise( stmt_iterate_variants ); 
  
  sql.finalise( stmt_insert_group );
  sql.finalise( stmt_lookup_group );
  sql.finalise( stmt_insert_group_member );
  sql.finalise( stmt_insert_group_variant );
  sql.finalise( stmt_iterate_group ); 
  sql.finalise( stmt_fetch_set_names1 ); 
  sql.finalise( stmt_fetch_set_names2 ); 
  sql.finalise( stmt_fetch_sets );
  
  sql.finalise( stmt_vcount );
  sql.finalise( stmt_totvcount );
  sql.finalise( stmt_indcount );
  sql.finalise( stmt_setcount );
  
  return true;
}


bool VarDBase::index()
{
    sql.query( "CREATE INDEX IF NOT EXISTS pos_var ON variants(chr,bp1,bp2);" );
    //sql.query( "CREATE INDEX IF NOT EXISTS file_idx ON variants(file_id);");
    sql.query( "CREATE INDEX IF NOT EXISTS vIndx1 ON vdat( var_id ) ; ");
    
    sql.query( "CREATE INDEX IF NOT EXISTS bcfIdx ON bcfs( file_id ); ");
    
    sql.query( "CREATE INDEX IF NOT EXISTS set1 ON set_data( set_id ) ; ");    
    sql.query( "CREATE INDEX IF NOT EXISTS set2 ON set_members( group_id ); ");
    sql.query( "CREATE INDEX IF NOT EXISTS set3 ON set_data( var_id ) ; ");
    sql.query( "CREATE INDEX IF NOT EXISTS set4 ON set_members( name ); ");

    sql.query( "CREATE INDEX IF NOT EXISTS meta1 ON indep_meta_data( var_id ) ; ");
    sql.query( "CREATE INDEX IF NOT EXISTS filetags ON files( tag ) ; " );

    release();
    init();
}

bool VarDBase::drop_index()
{
  sql.query( "DROP INDEX IF EXISTS pos_var;");
  sql.query( "DROP INDEX IF EXISTS vIndx1; ");
  sql.query( "DROP INDEX IF EXISTS set1; ");    
  sql.query( "DROP INDEX IF EXISTS set2; ");
  sql.query( "DROP INDEX IF EXISTS set3; ");
  sql.query( "DROP INDEX IF EXISTS set4; ");
  sql.query( "DROP INDEX IF EXISTS meta1; ");
  sql.query( "DROP INDEX IF EXISTS filetags; " );
}

int VarDBase::variant_count( uint64_t file_id )
{
    int n = 0;
    sql.bind_int64( stmt_vcount , ":file_id" , file_id );
    if ( sql.step( stmt_vcount ) )
	n = sql.get_int( stmt_vcount , 0 );
    sql.reset( stmt_vcount );
    return n;
}

int VarDBase::variant_count()
{
    int n = 0;
    if ( sql.step( stmt_totvcount ) )
	n = sql.get_int( stmt_totvcount , 0 );
    sql.reset( stmt_totvcount );
    return n;
}

int VarDBase::indiv_count( uint64_t file_id )
{
    int n = 0;
    sql.bind_int64( stmt_indcount , ":file_id" , file_id );
    if ( sql.step( stmt_indcount ) )
	n = sql.get_int( stmt_indcount , 0 );
    sql.reset( stmt_indcount );
    return n;
}

int VarDBase::set_count( uint64_t group_id )
{
    int n = 0;
    sql.bind_int64( stmt_setcount , ":group_id" , group_id );
    if ( sql.step( stmt_setcount ) )
	n = sql.get_int( stmt_setcount , 0 );
    sql.reset( stmt_setcount );
    return n;
}


uint64_t VarDBase::insert( const std::string & filename , 
			   const std::string & tag )
{
  sql.bind_text( stmt_insert_file , ":name" , filename );
  sql.bind_text( stmt_insert_file , ":tag" , tag );
  sql.step( stmt_insert_file );
  sql.reset( stmt_insert_file );
  return sql.last_insert_rowid();
}

// Headers

void VarDBase::insert_header( uint64_t file_id , const std::string & name , std::string value )
{
  sql.bind_int64( stmt_insert_header , ":file_id" , file_id );
  sql.bind_text( stmt_insert_header , ":name" , name );
  sql.bind_text( stmt_insert_header , ":value" , value );
  sql.step( stmt_insert_header );
  sql.reset( stmt_insert_header );
}

std::string VarDBase::print_headers( uint64_t file_id )
{
  std::string hd;
  sql.bind_int64( stmt_fetch_headers , ":file_id" , file_id );
  while ( sql.step( stmt_fetch_headers ) )
    {
      std::string key = sql.get_text( stmt_fetch_headers , 0 );
      std::string value = sql.get_text( stmt_fetch_headers , 1 );
      hd += key + "=" + value + "\n";
    }    
  sql.reset( stmt_fetch_headers );
  return hd;
}

void VarDBase::insert_metatype( uint64_t file_id , 
				const std::string & name , 
				mType mt , 
				int number, 
				int group,
				const std::string & description )
{
  sql.bind_int64( stmt_insert_metatype , ":file_id" , file_id );
  sql.bind_text( stmt_insert_metatype , ":name" , name );
  sql.bind_int( stmt_insert_metatype , ":type" , (int)mt );
  sql.bind_int( stmt_insert_metatype , ":number" , number );
  sql.bind_int( stmt_insert_metatype , ":group" , group );
  sql.bind_text( stmt_insert_metatype , ":description" , description );
  sql.step( stmt_insert_metatype );
  sql.reset( stmt_insert_metatype );
}

void VarDBase::set_metatypes( bool clear )
{
  if ( ! attached() ) return;  
  std::map<int,std::string> files = fetch_files();  
  std::map<int,std::string>::iterator i = files.begin();  
  while ( i != files.end() ) { set_file_metatypes( i->first , clear ); ++i; }
}


//
// and some special fields
//

void VarDBase::set_mask_metatypes( const Mask & mask )
{
  
  MetaInformation<VarFilterMeta>::field( PLINKSeq::PASS_FILTER() , META_FLAG , 1 , "Passed filters" );
  
  if ( mask.var() || mask.var_append() )
    {
      MetaInformation<VarMeta>::field( PLINKSeq::META_VSET() , META_TEXT , -1 , "Variant set name");
      MetaInformation<VarMeta>::field( PLINKSeq::META_VGRP() , META_INT  , -1 , "Variant set group");
      MetaMeta::set_static( PLINKSeq::META_VSET() );
      MetaMeta::set_static( PLINKSeq::META_VGRP() );
    }
  
  if ( mask.loc() || mask.loc_append() )
    {
      MetaInformation<VarMeta>::field( PLINKSeq::META_LSET() , META_TEXT , -1 , "Locus name");
      MetaInformation<VarMeta>::field( PLINKSeq::META_LGRP() , META_INT  , -1 , "Locus group");
      MetaMeta::set_static( PLINKSeq::META_LSET() );
      MetaMeta::set_static( PLINKSeq::META_LGRP() );
    }
  
  if ( mask.loc_set() || mask.loc_set_append() )
    {
      MetaInformation<VarMeta>::field( PLINKSeq::META_LSGRP() , META_TEXT , -1 , "Locus set name");
      MetaInformation<VarMeta>::field( PLINKSeq::META_LSGRP() , META_INT  , -1 , "Locus set group");
      MetaMeta::set_static( PLINKSeq::META_LSSET() );
      MetaMeta::set_static( PLINKSeq::META_LSGRP() );
    }

}

void VarDBase::set_file_metatypes( uint64_t file_id, bool clear )
{

  if ( clear ) 
    {            
      MetaInformation<VarMeta>::reset();
      MetaInformation<VarFilterMeta>::reset();
      MetaInformation<GenMeta>::reset();
    }

    sql.bind_int64( stmt_fetch_metatypes , ":file_id" , file_id );

    while ( sql.step( stmt_fetch_metatypes ) )
    {      
      std::string name = sql.get_text( stmt_fetch_metatypes , 0 );
      mType mt = (mType)sql.get_int( stmt_fetch_metatypes , 1 );
      int num = sql.get_int( stmt_fetch_metatypes , 2 );
      int grp = sql.get_int( stmt_fetch_metatypes , 3 );
      std::string desc = sql.get_text( stmt_fetch_metatypes , 4 );
      registerMetatype( name, mt , num, grp, desc );
    }
    sql.reset( stmt_fetch_metatypes );
}



uint64_t VarDBase::insert( uint64_t file_id , const Individual & person )
{  
  sql.bind_int64( stmt_insert_individual , ":file_id" , file_id );
  sql.bind_text( stmt_insert_individual , ":name" , person.id() );
  sql.step( stmt_insert_individual );
  sql.reset( stmt_insert_individual );

  uint64_t indiv_id = sql.last_insert_rowid();
  indiv.push_back( indiv_id );		      
  return indiv_id ;
}



uint64_t VarDBase::insert_consensus( uint64_t file_id , const Variant & var )
{
  
  // Insert a single, consensus SampleVariant, i.e. as read from a 
  // VCF file.  We will need another variant-insert function to insert
  // the multiple SampleVariants (i.e. when writing variants to the VARDB
  // from sources other than the VCF.  

  //
  // For now, allow the insert, but just put the consensus in  
  //

//   if ( ! var.single_sample() ) 
//     Helper::halt("Trying to save a MultiSample Variant in insert_consensus()");
  
  // 1) Insert variant information in the variant table
    
  sql.bind_int64(  stmt_insert_variant_key , ":file_id" , file_id );
  sql.bind_text(   stmt_insert_variant_key , ":name" , var.name() );
  sql.bind_int(    stmt_insert_variant_key , ":chr" , var.chromosome() );
  sql.bind_int(    stmt_insert_variant_key , ":bp1" , var.position() );
  sql.bind_int(    stmt_insert_variant_key , ":bp2" , var.stop() );
  sql.step( stmt_insert_variant_key );
  sql.reset( stmt_insert_variant_key );
    
  uint64_t var_id = sql.last_insert_rowid();
  
  // Obtain binary representation of variant and insert
  
  blob data = var.consensus.encode_var_BLOB();
  blob vdata = var.consensus.encode_vmeta_BLOB();
  blob gdata = var.consensus.encode_geno_BLOB();
  blob gmdata = var.consensus.encode_gmeta_BLOB();

  sql.bind_int64( stmt_insert_variant_data , ":var_id" , var_id );

  sql.bind_blob( stmt_insert_variant_data , ":data" , data );
  sql.bind_blob( stmt_insert_variant_data , ":vdata" , vdata );
  sql.bind_blob( stmt_insert_variant_data , ":gdata" , gdata );
  sql.bind_blob( stmt_insert_variant_data , ":gmdata" , gmdata );

  sql.step( stmt_insert_variant_data );
  sql.reset( stmt_insert_variant_data );
  return var_id;

}


SampleVariant & VarDBase::construct( Variant & var , sqlite3_stmt * s ,  IndividualMap * align )
{
  
  var.valid( true );  
  
  // fileset is : sql.get_int(  s , 1 ) 
  
  SampleVariant & sample = var.add( sql.get_int(  s , 1 ) );  
  
  // Sample attributes
  sample.index( sql.get_int64( s , 0 ) ) ; 
  
  // Core variant attributes
  
  var.name( sql.get_text( s , 2 ) );
  var.chromosome( sql.get_int( s , 3 ) );
  var.position( sql.get_int( s , 4 ) );
  var.stop( sql.get_int( s , 5 ) );
  
  // Are data stored in VARDB? 
  // Or look-up this variant from a compressed VCF or BCF? 
  
  int64_t offset_idx = sql.get_int64( s , 6 );
  
  if ( offset_idx )
  {
    
      int file_id = sql.get_int( s , 1 );
      
      // Is this a VCF or BCF?
      
      VCFZ * vcfz = vcfzmap[ file_id ];
      
//       std::cout << vcfzmap.size() << " is vcfzmap.size()\n";
//       std::map<int,VCFZ*>::iterator ii = vcfzmap.begin();
//       while ( ii != vcfzmap.end() )
// 	{
// 	  std::cout << ii->first << "\t" << ii->second->name() << "\n";
// 	  ++ii;
// 	}
//       if ( vcfz == NULL ) std::cout << "null VCFZ for " << file_id << "\n";

      if ( vcfz ) 
	{
	  SampleVariant & target = ( ! align->multi_sample() ) ? var.consensus : sample ;      
	  SampleVariant & genotype_target = align->flat() ? var.consensus : sample ;
	  vcfz->read_record( var , target , genotype_target, offset_idx ); 
      }
      else 
      {
	  BCF * bcf = bcfmap[ file_id ];	  

	  if ( bcf ) 
	  {
	      SampleVariant & target = ( ! align->multi_sample() ) ? var.consensus : sample ;      
	      SampleVariant & genotype_target = align->flat() ? var.consensus : sample ;
	      bcf->read_record( var , target , genotype_target, offset_idx ); 
	  }
	  else
	      Helper::halt( "requested compressed-VCF or BCF not attached" );
      }
      
  }
  else
    {

      // Attach BLOB from VARDB, which will be later expanded (if needed)
      
      sqlite3_stmt * s = stmt_fetch_variant_data_all;
      if ( fetch_mode == NO_GMETA ) s = stmt_fetch_variant_data_vmeta_geno;
      else if ( fetch_mode == ONLY_VMETA ) s = stmt_fetch_variant_data_vmeta;
      else if ( fetch_mode == ONLY_GENO ) s = stmt_fetch_variant_data_geno;
      
      sql.bind_int64( s , ":var_id" , sample.index() );
      sql.step( s );
      
      // Only store BLOB in raw form: do not parse until we 
      // know we are definitely interested in this variant      
      
      blob var_blob = sql.get_blob( s , 0 );
      
      blob * vmeta_blob  = NULL;
      blob * geno_blob   = NULL;
      blob * gmeta_blob  = NULL;
      

      if ( fetch_mode == ALL || fetch_mode == NO_GMETA )
	{	       
	  vmeta_blob  = new blob( sql.get_blob( s , 1 ) );
	  geno_blob   = new blob( sql.get_blob( s , 2 ) );
	  if ( fetch_mode == ALL ) gmeta_blob  = new blob( sql.get_blob( s , 3 ) ); 
	}
      else if ( fetch_mode == ONLY_VMETA )
	{
	  vmeta_blob  = new blob( sql.get_blob( s , 1 ) );
	}
      else if ( fetch_mode == ONLY_GENO )
	{
	  geno_blob  = new blob( sql.get_blob( s , 1 ) );
	}

      // this extracts out a PB string, so we can remove BLOBs afterwards
      
      sample.store_BLOBs( &var_blob , vmeta_blob , geno_blob , gmeta_blob );      
      
      if ( vmeta_blob ) delete vmeta_blob;
      if ( geno_blob ) delete geno_blob;
      if ( gmeta_blob ) delete gmeta_blob;
      
      sql.reset( s );
    }

  return sample;
}





//
// Sets
//


uint64_t VarDBase::set_group_id( const std::string & name , bool temp , const std::string & desc)
{
  
  uint64_t group_id = 0;
  
  sql.bind_text(stmt_lookup_group, ":name" , name ); 
  
  if ( sql.step( stmt_lookup_group ) ) 
    {
      group_id = sql.get_int64( stmt_lookup_group , 0 ) ;
      sql.reset( stmt_lookup_group );
    }
  else
    {
      sql.reset( stmt_lookup_group );
      
      sql.bind_text( stmt_insert_group , ":name" , name );
      sql.bind_int( stmt_insert_group , ":temp" , 1 );
      sql.bind_text( stmt_insert_group , ":description" , name + " (default name)" );
      sql.step( stmt_insert_group );      
      group_id = sql.last_insert_rowid();
      sql.reset( stmt_insert_group );
    }
  
  return group_id;
  
}

uint64_t VarDBase::lookup_group_id( const std::string & name )
{
  uint64_t group_id = 0;
  sql.bind_text(stmt_lookup_group, ":name" , name ); 
  if ( sql.step( stmt_lookup_group ) ) 
    group_id = sql.get_int64( stmt_lookup_group , 0 ) ;
  sql.reset( stmt_lookup_group );
  return group_id;
}

uint64_t VarDBase::set_member_id( const uint64_t grp_id , const std::string & name )
{
  sql.bind_int64( stmt_insert_group_member , ":group_id" , grp_id );
  sql.bind_text( stmt_insert_group_member , ":name" , name );
  sql.step( stmt_insert_group_member );      
  uint64_t new_id = sql.last_insert_rowid();
  sql.reset( stmt_insert_group_member );
  return new_id;
}


void VarDBase::set_add_variant( const uint64_t set_id , const uint64_t var_id )
{
  sql.bind_int64( stmt_insert_group_variant , ":set_id" , set_id );
  sql.bind_int64( stmt_insert_group_variant , ":var_id" , var_id );
  sql.step( stmt_insert_group_variant );      
  sql.reset( stmt_insert_group_variant );
}


std::vector<std::string> VarDBase::fetch_individuals(uint64_t file_id)
{
  std::vector<std::string> res;
  sql.bind_int64( stmt_fetch_individuals, ":file_id" , file_id );
  while ( sql.step( stmt_fetch_individuals ) )
    {	
      std::string i = sql.get_text( stmt_fetch_individuals , 2 );	
      res.push_back(i);
    }
  sql.reset( stmt_fetch_individuals );
  return res;        
}


std::vector<std::map<std::string,std::string> > VarDBase::fetch_headers(uint64_t file_id)
{
  std::vector<std::map<std::string,std::string> > res;
  sql.bind_int64( stmt_fetch_headers, ":file_id" , file_id );
  while ( sql.step( stmt_fetch_headers ) )
    {	
      std::string k = sql.get_text( stmt_fetch_headers , 0 );
      std::string v = sql.get_text( stmt_fetch_headers , 1 );
      std::map<std::string,std::string> m;	
      m["KEY"] = k;
      m["VALUE"] = v;
      res.push_back(m);
    }
  sql.reset( stmt_fetch_headers );
  return res;    
}


std::vector<std::map<std::string,std::string> > VarDBase::fetch_metatypes(uint64_t file_id)
{
  std::vector< std::map<std::string,std::string> > res;
  sql.bind_int64( stmt_fetch_metatypes, ":file_id" , file_id );
  
  while ( sql.step( stmt_fetch_metatypes ) )
    {
      std::string name = sql.get_text( stmt_fetch_metatypes , 0 );
      int t    = sql.get_int(  stmt_fetch_metatypes , 1 );
      std::string n = sql.get_text( stmt_fetch_metatypes , 2 );
      int    g = sql.get_int( stmt_fetch_metatypes , 3 );
      std::string d = sql.get_text( stmt_fetch_metatypes , 4 );
      std::map<std::string,std::string> m;
      m[ "NAME" ] = name;
	
      switch ( t ) 
	{
	case META_FLAG :
	  m[ "TYPE" ] = "Flag";
	  break;
	case META_UNDEFINED :
	  m[ "TYPE" ] = "Undefined";
	  break;
	case META_TEXT :
	  m[ "TYPE" ] = "String";
	  break;
	case META_INT :
	  m[ "TYPE" ] = "Integer";
	  break;
	case META_FLOAT :
	  m[ "TYPE" ] = "Float";
	  break;
	case META_BOOL :
	  m[ "TYPE" ] = "Bool";
	  break;
	case META_CHAR :
	  m[ "TYPE" ] = "Char";
	}
            
      if      ( g == 1 ) m[ "GRP" ] = "Variant";
      else if ( g == 2 ) m[ "GRP" ] = "Genotype";
      else if ( g == 8 ) m[ "GRP" ] = "Variant Filter";
      else               m[ "GRP" ] = "?";

      m[ "NUM" ] = n;	
      m[ "DESC" ] = d;
      
      res.push_back(m);
    }

  sql.reset( stmt_fetch_metatypes );
  
  return res;

}


int VarDBase::fileID( const std::string & filename )
{
  int r = 0;
  sql.bind_text( stmt_fetch_file_id , ":name" , filename );
  if ( sql.step( stmt_fetch_file_id ) )
    r = sql.get_int( stmt_fetch_file_id , 0 );
  sql.reset( stmt_fetch_file_id );
  return r;		
}


std::map<int,std::string> VarDBase::fetch_files( Mask * mask )
{
  std::map<int,std::string> res;
  while( sql.step( stmt_fetch_files ) )
    {
      int f = sql.get_int( stmt_fetch_files , 0 );
      std::string n = sql.get_text( stmt_fetch_files , 1 );      
      if ( mask && ! mask->use_file( f ) ) { continue; }
      res[f] = n;
    }
  sql.reset( stmt_fetch_files );
  return res;
}


int VarDBase::n_files( Mask * mask )
{
  std::map<int,std::string> f = fetch_files( mask );
  return f.size();
}


std::map<int,std::string> VarDBase::fetch_sets()
{
  std::map<int,std::string> res;
  while( sql.step( stmt_fetch_sets ) )
    {
      int f = sql.get_int( stmt_fetch_sets , 0 );
      std::string n = sql.get_text( stmt_fetch_sets , 1 );
      res[f] = n;
    }
  sql.reset( stmt_fetch_sets );
  return res;
}

std::string VarDBase::group_name( const uint64_t group_id )
{
  sql.bind_int64( stmt_lookup_group_name , ":group_id" , group_id );
  if ( sql.step( stmt_lookup_group_name ) )
    return sql.get_text( stmt_lookup_group_name , 0 );
  else
    return "";
}



void VarDBase::addMetaFields( Variant & var, sqlite3_stmt * s, Mask & mask)
{
  
  // Determine the specified ordering of meta-information
  
  // 0   var_id
  // 1   file_id
  // 2   name
  // 3   chr
  // 4   bp1  
  // 5   bp2
  // 6   offset

  // 7   type 1=var, 2=loc, 3=locset, 0=ignore
  // 8   name
  // 9   grp

  
  // types

  int type = sql.get_int( s , 7 ) ;
  if ( type == 0 ) return;

  std::string name = sql.get_text( s, 8 );
  if ( name == "" ) return;
  
  int grp = sql.get_int( s, 9 );
  
  if ( type == 1 ) // variants
    {
      if ( var.meta.add_if_unique( PLINKSeq::META_VSET() , name ) )
	var.meta.add( PLINKSeq::META_VGRP() , grp );	          
    }
  else if ( type == 2 ) // locdb    
    {
      if ( var.meta.add_if_unique( PLINKSeq::META_LSET() , name ) )
	var.meta.add( PLINKSeq::META_LGRP() , grp );	          
    }
  else if ( type == 3 ) // pathways "locus-set"
    {
      if ( var.meta.add_if_unique( PLINKSeq::META_LSSET() , name ) )
	var.meta.add( PLINKSeq::META_LSSET() , grp );	    
    }
  
}


Variant VarDBase::fetch( uint64_t var_id )
{

  // Return a single SampleVariant (wrapped in a Variant)

  Variant var;
  if ( ! attached() ) { var.valid( false ); return var; }

  sql.bind_int64( stmt_fetch_variant_key , ":var_id" , var_id );
  if ( sql.step( stmt_fetch_variant_key ) )
    {      
      SampleVariant & sample = construct( var , stmt_fetch_variant_key , &indmap );

      sample.decode_BLOB( &var ,     // parent
			  &indmap ,  // indiv-map
			  NULL );    // mask

      var.make_consensus( &indmap );  
    }
  sql.reset( stmt_fetch_variant_key );
  return var;
}

Variant VarDBase::fetch( int chr , int bp1 )
{
  // Return 1+ SampleVariants, based on physical position, 
  //  in the form of a single Variant

  Variant var;
  if ( ! attached() ) { var.valid( false ); return var; }

  sql.bind_int( stmt_fetch_variant_pos , ":chr" , chr );
  sql.bind_int( stmt_fetch_variant_pos , ":bp1" , bp1 );
  
  while ( sql.step( stmt_fetch_variant_pos ) )
    {      
      SampleVariant & sample = construct( var , stmt_fetch_variant_pos , &indmap );
      sample.decode_BLOB( &var , &indmap , NULL );      
    }
  
  var.make_consensus( &indmap );      

  sql.reset( stmt_fetch_variant_pos );
  return var;
}


std::set<Variant> VarDBase::fetch( const Region & region )
{
  
  std::set<Variant> s;
 
  if ( ! attached() ) return s; 

  sql.bind_int( stmt_fetch_variant_range , ":chr" , region.chromosome() );
  sql.bind_int( stmt_fetch_variant_range , ":rstart" , region.start.position() );
  sql.bind_int( stmt_fetch_variant_range , ":rend" , region.stop.position() );
    
  std::map<int,Variant> vmap;

  while ( sql.step( stmt_fetch_variant_range ) )
    {
      // extract BP position on this chromosome
      int pos = sql.get_int( stmt_fetch_variant_range , 4 );      
      SampleVariant & sample = construct( vmap[pos] , stmt_fetch_variant_range , &indmap );
      sample.decode_BLOB( &vmap[pos] , &indmap , NULL );
    } 
  sql.reset( stmt_fetch_variant_range ) ;  

  std::map<int,Variant>::iterator i = vmap.begin();
  while ( i != vmap.end() )
    {
      i->second.make_consensus( &indmap );
      s.insert(i->second);
      ++i;
    }

  return s;
}



int2 VarDBase::make_summary( std::string filename )
{

  sql.bind_text( stmt_fetch_file_id, ":name" , filename );
  
  int file_id;
  
  if ( sql.step( stmt_fetch_file_id ) )
    {
      file_id = sql.get_int( stmt_fetch_file_id , 0 ); 
      sql.reset( stmt_fetch_file_id );
    }
  else 
    {
      sql.reset( stmt_fetch_file_id );
      return int2(0,0);
    }
  
  return make_summary( file_id );
}

int2 VarDBase::make_summary( int file_id )  
{
  std::map<int,std::string> f = fetch_files();
  int2 niv( indiv_count(file_id ) , variant_count(file_id ) );
  sql.bind_text( stmt_insert_file_summary, ":name" , f[file_id] );
  sql.bind_int( stmt_insert_file_summary , ":ni" ,  niv.p1 );
  sql.bind_int( stmt_insert_file_summary , ":nv" ,  niv.p2 );
  sql.step( stmt_insert_file_summary );
  sql.reset( stmt_insert_file_summary );
  return niv;
}


std::string VarDBase::summary( Mask * mask )
{
  
  std::stringstream ss;
  
  std::map<int,std::string> f = fetch_files( mask );
  
  ss << "VARDB\t"
     << "N_UNIQ_VAR=" << variant_count() << "\n";
  
  
  std::map<int,std::string>::iterator i = f.begin();
  while ( i != f.end() )
    {
      
      int ni = 0, nv = 0;
      sql.bind_int64( stmt_fetch_file_summary , ":file_id" , i->first );
      if ( sql.step( stmt_fetch_file_summary ) )
	{
	  ni = sql.get_int( stmt_fetch_file_summary , 0 );
	  nv = sql.get_int( stmt_fetch_file_summary , 1 );
	}
      sql.reset( stmt_fetch_file_summary );
      
      ss << "VARDB\t"
	 << "FILE_N=" << i->first << "\t" 
	 << "TAG=" << file_tag( i->first ) << "\t"
	 << "N_INDIV=" << ni << "\t"
	 << "N_VAR=" << nv << "\t"
	 << "FILE_NAME=" << i->second << "\n";
      ++i;
    }
    
    // Sets
    
  std::map<int,std::string> s = fetch_sets();
  std::map<int,std::string>::iterator j = s.begin();
  while ( j != s.end() )
    {
      ss << "VARDB\t"
	 << "SET_N=" << j->first << "\t" 
	 << "N_VAR=" << set_count( j->first ) << "\t"
	 << "SET_NAME=" << j->second << "\n";
      ++j;
    }
  
  // Indep meta-information
  
  std::map<std::string,int>::iterator k = indep_metamap.begin();
  while ( k != indep_metamap.end() )
    {
      ss << "VARDB\t"
	 << "ADD_META\t"
	 << "NAME=" << k->first 
	 << "\n";
      ++k;
    }
  
  return ss.str();
}


void VarDBase::append_metainformation( Variant & v , const std::set<int> & grp )
{
  Helper::halt( "append_metainformation not implemented yet..." ); 
}


void VarDBase::drop(int g)
{
  
  // Obtain list of variants that need to be deleted from 'vdat' and 'set_data' 
  
  sql.query(" ATTACH \":memory:\" AS tmpdel ; " );

  sql.query("CREATE TABLE tmpdel.tmp AS SELECT var_id FROM variants WHERE file_id == " 
	    + Helper::int2str( g ) + " ; " );
   
  sql.query( "DELETE FROM headers     WHERE file_id == " + Helper::int2str( g ) + ";" );
  sql.query( "DELETE FROM metatypes   WHERE file_id == " + Helper::int2str( g ) + ";" );
  sql.query( "DELETE FROM files       WHERE file_id == " + Helper::int2str( g ) + ";" );
  sql.query( "DELETE FROM individuals WHERE file_id == " + Helper::int2str( g ) + ";" );
  sql.query( "DELETE FROM vdat        WHERE var_id IN ( SELECT var_id FROM tmpdel.tmp ) ; " );
  sql.query( "DELETE FROM set_data    WHERE var_id IN ( SELECT var_id FROM tmpdel.tmp ) ; " );
  sql.query( "DELETE FROM variants    WHERE file_id == " + Helper::int2str( g ) + ";" );
  
}


void VarDBase::vacuum()
{
  sql.query("VACUUM;");
}

uint64_t VarDBase::lookup_file_id( const std::string & tag ) 
{
  uint64_t id = 0;
  sql.bind_text( stmt_fetch_file_from_tag , ":tag" , tag );
  if ( sql.step( stmt_fetch_file_from_tag ) )
    id = sql.get_int64( stmt_fetch_file_from_tag , 0 );
  else
    {
      // otherwise, is tag a valid file number?
      int t = 0;
      if ( Helper::str2int(tag , t ) ) id = t;
    }
  sql.reset( stmt_fetch_file_from_tag );
  return id;
}

void VarDBase::insert_file_tag( uint64_t id , const std::string & tag )
{
  sql.bind_int64( stmt_insert_file_tag , ":file_id" , id );
  sql.bind_text( stmt_insert_file_tag , ":tag" , tag );
  sql.step( stmt_insert_file_tag );
  sql.reset( stmt_insert_file_tag );
}



std::string VarDBase::file_tag( uint64_t id )
{

  if ( id == 0 ) return ".";

  // use cache 
  std::map<int,std::string>::iterator i = file_tag_map.find( id );
  if ( i != file_tag_map.end() ) return i->second;

  std::string ftag = "";
  
  sql.bind_int64( stmt_fetch_tag_from_file, ":file_id" , id );
  if ( sql.step( stmt_fetch_tag_from_file ) )
    ftag = sql.get_text( stmt_fetch_tag_from_file , 0 );
  sql.reset( stmt_fetch_tag_from_file );  

  // default to actual number if nothing else found
  if ( ftag == "" ) 
    ftag = Helper::int2str( id );

  file_tag_map[ id ] = ftag;
  return ftag;
}

uint64_t VarDBase::file_tag( const std::string & filetag )
{
  // use cache 
  std::map<std::string,int>::iterator i = reverse_file_tag_map.find( filetag );
  if ( i != reverse_file_tag_map.end() ) return i->second;

  uint64_t id = lookup_file_id( filetag ); 
  if ( id != 0 ) reverse_file_tag_map[ filetag ] = id;
  return id;    
}

bool VarDBase::chr_code( const int c , const std::string & n , const ploidy_t ploidy )
{
  sql.bind_text( stmt_fix_chr_code , ":name" , n );
  sql.bind_int( stmt_fix_chr_code , ":chr_id" , c );
  sql.bind_int( stmt_fix_chr_code , ":ploidy" , ploidy );
  bool okay = sql.step( stmt_fix_chr_code );
  sql.reset( stmt_fix_chr_code );
  chr_name_map[c]=n;
  chr_code_map[n]=c;
  chr_ploidy_map[c]=ploidy;
  return okay;  
}

int VarDBase::chr_code( const std::string & n , ploidy_t * ploidy ) 
{
  std::map<std::string,int>::iterator i = chr_code_map.find(n);
  if ( i != chr_code_map.end() ) 
    {
      if ( ploidy ) *ploidy = chr_ploidy_map[ i->second ];
      return i->second;
    }

  sql.bind_text( stmt_fetch_chr_code , ":name" , n );
  if ( sql.step( stmt_fetch_chr_code ) )
    {
      int c = sql.get_int( stmt_fetch_chr_code , 0 );
      chr_code_map[n]=c;
      chr_ploidy_map[c] = (ploidy_t)sql.get_int( stmt_fetch_chr_code , 1 );
      if ( ploidy ) *ploidy = chr_ploidy_map[c];
      sql.reset( stmt_fetch_chr_code );
      return c;
    }

  //otherwise, we need to insert
  sql.bind_text( stmt_insert_chr_code , ":name" , n );
  if ( ploidy ) sql.bind_int( stmt_insert_chr_code , ":ploidy" , *ploidy );
  else sql.bind_int( stmt_insert_chr_code , ":ploidy" , PLOIDY_UNKNOWN );
  sql.step( stmt_insert_chr_code );
  int c = sql.last_insert_rowid();
  sql.reset( stmt_insert_chr_code );
  return c;
}

ploidy_t VarDBase::ploidy( const int c )
{
  std::map<int,ploidy_t>::iterator i = chr_ploidy_map.find( c );
  return i == chr_ploidy_map.end() ? PLOIDY_UNKNOWN : i->second;
}


std::string VarDBase::chr_name( const int c ) 
{
  std::map<int,std::string>::iterator i = chr_name_map.find(c);
  if ( i != chr_name_map.end() ) return i->second;
  sql.bind_int( stmt_fetch_chr_name , ":chr_id" , c );
  std::string n = ".";
  if ( sql.step( stmt_fetch_chr_name ) )
    {
      n = sql.get_text( stmt_fetch_chr_name , 0 );
      chr_name_map[c]=n;
      chr_ploidy_map[c]= (ploidy_t)sql.get_int( stmt_fetch_chr_name , 1 );
    }
  sql.reset( stmt_fetch_chr_name );
  return n;  
}

bool VarDBase::check_version()
{
  
  if ( ! sql.table_exists( "dbmeta" ) )
    Helper::halt( "old database format, expecting VARDB v" 
		  + Helper::int2str( PLINKSeq::VARDB_VERSION_NUMBER() ) 
		  + " : to fix, remake your project" );

  // expected version # is given by  PLINKSeq::VARDB_VERSION_NUMBER()
  int v = 0;
  
  sqlite3_stmt * s = sql.prepare( "SELECT varvalue FROM dbmeta WHERE varname == 'VERSION'; " );
  
  if ( sql.step(s) ) 
    {
      v = sql.get_int( s , 0 );
      sql.finalise(s);
    }
  else // implies a new database, as version note yet set -- so add one
    {
      sql.finalise(s);
      sqlite3_stmt * si = sql.prepare("INSERT OR REPLACE INTO dbmeta(varname, varvalue ) values( :x , :y ) ; " );
      std::string vn = "VERSION";
      v = PLINKSeq::VARDB_VERSION_NUMBER();
      sql.bind_text( si , ":x" , vn );
      sql.bind_int( si , ":y" , v );
      sql.step(si);
      sql.finalise(si);
    }

  if ( v != PLINKSeq::VARDB_VERSION_NUMBER() ) 
    Helper::halt("VARDB version " 
		 + Helper::int2str( v ) + " but expected " 
		 + Helper::int2str( PLINKSeq::VARDB_VERSION_NUMBER() ) 
		 + " : to fix, remake your project" );
  
  return true;
}


void VarDBase::populate_bcf_map()
{
    
  bcfmap.clear();
  vcfzmap.clear();
  
  while ( sql.step( stmt_fetch_bcfs ) )
    {
      
	int fid              = sql.get_int( stmt_fetch_bcfs , 0 );
	int type             = sql.get_int( stmt_fetch_bcfs , 1 );
	std::string filename = sql.get_text( stmt_fetch_bcfs , 2 );
	int nind             = sql.get_int( stmt_fetch_bcfs , 3 );
	
	
	if ( type == 1 ) // compressed VCF
	  {
	
	    VCFZ * vcfz = GP->fIndex.vcfz( filename );
	    
	    if ( vcfz ) 
	      {
		vcfzmap[ fid ] = vcfz;	  
		vcfz->reading();
		vcfz->open();
	      }
	    else
	      plog.warn( "could not find compressed VCF " , filename );	
	  }
	else if ( type == 2 ) // BCF
	  {
	    
	    BCF * bcf = GP->fIndex.bcf( filename );
	    
	    if ( bcf ) 
	    {
	      bcfmap[ fid ] = bcf;	  
	      bcf->set_nind( nind );
	      bcf->reading();
	      bcf->open();
	      bcf->set_types();
	    }
	    else
	      plog.warn( "could not find BCF " , filename );	
	  }
	
    } // next file
    sql.reset( stmt_fetch_bcfs );
}


void VarDBase::store_bcf_n( uint64_t f , const std::string & filename , int type , int nind )
{
  // type, 1 = VCF(compressed), 2 = BCF
  sql.bind_int64( stmt_insert_bcf_n , ":file_id" , f ); 
  sql.bind_int( stmt_insert_bcf_n , ":type" , type ); 
  sql.bind_int( stmt_insert_bcf_n , ":nind" , nind ); 
  sql.bind_text( stmt_insert_bcf_n , ":filepath" , filename );
  sql.step( stmt_insert_bcf_n );
  sql.reset( stmt_insert_bcf_n );
}

void VarDBase::insert_bcf_index( uint64_t file_id , const Variant & var , int64_t offset )
{
  // add row to variants table
  sql.bind_int64( stmt_insert_bcf_idx , ":file_id" , file_id );
  sql.bind_text( stmt_insert_bcf_idx , ":name" , var.name() );
  sql.bind_int( stmt_insert_bcf_idx , ":chr" , var.chromosome() );
  sql.bind_int( stmt_insert_bcf_idx , ":bp1" , var.position() );
  sql.bind_int( stmt_insert_bcf_idx , ":bp2" , var.stop() == var.position() ? 0 : var.stop() );
  sql.bind_int64( stmt_insert_bcf_idx , ":offset" , offset );
  sql.step( stmt_insert_bcf_idx );
  sql.reset( stmt_insert_bcf_idx );  
}


