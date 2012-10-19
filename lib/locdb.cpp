#include "plinkseq/locdb.h"
#include "plinkseq/defs.h"
#include "plinkseq/regions.h"
#include "plinkseq/gstore.h"
#include "plinkseq/annot.h"

#include <iostream>
#include <set>
#include <utility>

extern GStore * GP;

bool LocDBase::wipe( const std::string & n )
{
  if ( Helper::fileExists(n) ) 
    {
      Helper::remove_file( n );
      return true;
    }
  else
    return false;
}



bool LocDBase::attach( const std::string & n )
{
  
  if ( attached() ) dettach();

  if ( n == "-" || n == "." ) { dettach(); return false; } 
  
  sql.open(n); 
  
  sql.synchronous(false);
  
  fname = n;
  
  //
  // DB version, and a place for various other meta-information in future
  //
  
  sql.query(" CREATE TABLE IF NOT EXISTS dbmeta("
            "   varname      VARCHAR(20) NOT NULL , "
            "   varvalue    VARCHAR(20) NOT NULL , "
            " CONSTRAINT uMeta UNIQUE (varname ) ); " );


  // Main locus table
  
  sql.query(" CREATE TABLE IF NOT EXISTS loci("
	    "   loc_id   INTEGER PRIMARY KEY , "
	    "   name     VARCHAR(20) , "	      
	    "   group_id INTEGER NOT NULL , "
	    "   chr      INTEGER  , "
	    "   bp1      INTEGER  , "
	    "   bp2      INTEGER  , "
	    "   altname  VARCHAR(20)  ); " );
  
  // Simple name-search table (for dumping gene-lists)

  sql.query(" CREATE TABLE IF NOT EXISTS searchnames("
	    "   group_id INTEGER NOT NULL , "
	    "   name     VARCHAR(20) ); " );
  
  // Sub-region/exon table (with strand/frame as fixed fields
  
  sql.query(" CREATE TABLE IF NOT EXISTS subloci("
	    "   sub_id   INTEGER PRIMARY KEY , "
	    "   loc_id   INTEGER NOT NULL , "
	    "   name     VARCHAR(20) , "
	    "   chr      INTEGER  , "
	    "   bp1      INTEGER  , "
	    "   bp2      INTEGER  , "
	    "   strand   INTEGER  , "
	    "   frame    INTEGER  ); " );    
  

  //
  // Region-set information
  //

  sql.query(" CREATE TABLE IF NOT EXISTS set_groups("
	    "   group_id       INTEGER PRIMARY KEY , "
	    "   loc_group_id   INTEGER NOT NULL , "
	    "   name           VARHCAR(20) , "
	    "   temp           CHAR(1) , "
	    "   description    VARCHAR(20) ) ; " );
  
  sql.query(" CREATE TABLE IF NOT EXISTS set_members("
	    "   set_id       INTEGER PRIMARY KEY , "
	    "   group_id     INTEGER NOT NULL , "
	    "   name         VARCHAR(20) ); " );    
  
  sql.query(" CREATE TABLE IF NOT EXISTS set_data("
	    "   loc_id   INTEGER NOT NULL , "
	    "   set_id   INTEGER NOT NULL ); " );
  

  //
  // Meta-information tables
  //

  sql.query( " CREATE TABLE IF NOT EXISTS metatypes("
             "   field_id     INTEGER PRIMARY KEY , "
             "   name         VARCHAR(8) , "
	     "   type         VARCHAR(8) , "
	     "   number       INTEGER , "
	     "   description  VARCHAR(20) ); ");

  sql.query( " CREATE TABLE IF NOT EXISTS loc_meta("
	     "   loc_id    INTEGER NOT NULL , "
	     "   value     VARCHAR(20) ); " );
  
  sql.query( " CREATE TABLE IF NOT EXISTS subloc_meta("
	     "   sub_id    INTEGER NOT NULL , "
	     "   value     VARCHAR(20) ); ");


  //
  // Name alias tables
  //

  sql.query(" CREATE TABLE IF NOT EXISTS aliases("
	    "   group_id1 INTEGER      NOT NULL , "
	    "   name1     VARCHAR(20)  NOT NULL , "
	    "   group_id2 INTEGER      NOT NULL , "
	    "   name2     VARCHAR(20)  , "
	    " CONSTRAINT uniqNames UNIQUE ( group_id1, name1, group_id2, name2 ) ); " );
  
  sql.query(" CREATE TABLE IF NOT EXISTS alias_groups("
	    "   group_id    INTEGER PRIMARY KEY , "
	    "   group_name  VARCHAR(20) NOT NULL , "
            " CONSTRAINT uniqID UNIQUE ( group_name ) ); " );
  

  // Allow a locus to be assigned to 1 or more individuals
  // This will be empty for the reference LOCDB, but for the 
  // user LOCDB could be used to represent CNVs, shared segments, 
  // ROH, etc
  
  sql.query(" CREATE TABLE IF NOT EXISTS segments("
            "   indiv_id INTEGER NOT NULL , "
	    "   loc_id   INTEGER NOT NULL  ); " );
  
  sql.query(" CREATE TABLE IF NOT EXISTS individuals("
            "   indiv_id INTEGER PRIMARY KEY , "
	    "   name     VARCHAR(20) NOT NULL , "
            " CONSTRAINT uniqID UNIQUE ( name ) ); " );


  //
  // Special region table, e.g. X chromosome
  //

  sql.query( "CREATE TABLE IF NOT EXISTS special("
	     "  loc_id   INTEGER PRIMARY KEY , "
	     "  name     VARCHAR(20) , "
	     "  value    VARCHAR(20) ); ");


  //
  // Group information table
  //

  sql.query( " CREATE TABLE IF NOT EXISTS groups("
	     "   group_id     INTEGER PRIMARY KEY , "
	     "   name         VARCHAR(20) NOT NULL , "
             "   temp         CHAR(1) , "
	     "   description  TEXT ); " );


  //
  // Region overlap table (i.e. all pairs of loci that overlap)
  //

  sql.query( " CREATE TABLE IF NOT EXISTS overlaps("
	     "   loc_id1    INTEGER NOT NULL , "
	     "   loc_id2    INTEGER NOT NULL , "
	     "   val_inter  INTEGER NOT NULL , "
	     "   val_union  INTEGER NOT NULL ); ");


  // Prepate some key queries
  
  init();

  index();
  
  set_metatypes();

  read_alias_groups();

  
  return true;

}

uint64_t LocDBase::insert_alias_group( const std::string & alias )
{
  sql.bind_text( stmt_loc_alias_group_insert , ":group_name" , alias );
  sql.step( stmt_loc_alias_group_insert );
  sql.reset( stmt_loc_alias_group_insert );

  // add to internal map, and return the code
  read_alias_groups();
  return alias_group_table[alias];
}

void LocDBase::clear_alias_groups()
{
  alias_group_table.clear();
  alias_group_table[""] = 0;
  alias_group_reverse_table[0] = "";
}

void LocDBase::read_alias_groups()
{
  while ( sql.step( stmt_loc_alias_group_dump ) )
    {      
      int n = sql.get_int( stmt_loc_alias_group_dump , 0 ) ;
      std::string g = sql.get_text( stmt_loc_alias_group_dump , 1 );      
      alias_group_table[ g ] = n;
      alias_group_reverse_table[ n ] = g;
    }      
  sql.reset( stmt_loc_alias_group_dump );
}

void LocDBase::set_metatypes( bool clear )
{
    if ( clear ) 
    {
	MetaInformation<LocMeta>::reset();    
    }
    
    while ( sql.step( stmt_fetch_metatypes ) )
    {	
      std::string name = sql.get_text( stmt_fetch_metatypes , 0 );
      mType mt = (mType)sql.get_int( stmt_fetch_metatypes , 1 );
      int num = sql.get_int( stmt_fetch_metatypes , 2 );
      std::string desc = sql.get_text( stmt_fetch_metatypes , 3 );
      registerMetatype( name, mt , 1, META_GROUP_LOC, desc );
    }
    sql.reset( stmt_fetch_metatypes );
}


bool LocDBase::dettach()
{
  release();
  sql.close();
}


bool LocDBase::init()
{

    populate_meta_field_map();

    stmt_insert_special = 
      sql.prepare(" INSERT INTO special ( name , value ) values( :name , :value ) ; ");
    
    stmt_fetch_special = 
      sql.prepare(" SELECT value FROM special WHERE name == :name ; ");

    stmt_loc_insert_group_name = 
	sql.prepare("INSERT OR REPLACE INTO groups ( name, temp, description ) "
		" values( :name, :temp, :description ) ; ");
  

    stmt_loc_lookup_group_name = 
      sql.prepare("SELECT group_id FROM groups WHERE name == :name ; " );

    stmt_loc_lookup_group_id = 
      sql.prepare("SELECT name FROM groups WHERE group_id == :group_id ; " );

    stmt_loc_lookup_real_name = 
      sql.prepare("SELECT * FROM loci WHERE group_id == :group_id AND altname == UPPER(:altname) ; " );
    
    stmt_loc_lookup_real_name_only = 
      sql.prepare("SELECT name FROM loci WHERE group_id == :group_id AND altname == UPPER(:altname) ; " );

    stmt_loc_replace_real_name = 
	sql.prepare("UPDATE OR IGNORE loci SET altname = :altname WHERE group_id == :group_id AND name == :name ; " );

     stmt_loc_replace_real_name_alternate = 
 	sql.prepare("UPDATE OR IGNORE loci SET altname = :newname WHERE group_id == :group_id AND altname == :oldname  ; " );

//     stmt_loc_replace_real_name_alternate = 
// 	sql.prepare("SELECT * FROM loci WHERE group_id == :group_id AND altname == :oldname ; " );

    stmt_loc_update_temp_status = 
      sql.prepare("INSERT OR REPLACE INTO groups ( group_id, temp ) "
		  " values( :group_id, :temp ) ; " );
    
    stmt_loc_lookup_temp_status = 
      sql.prepare("SELECT temp FROM groups WHERE group_id == :group_id ; " );
    
    stmt_loc_remove_group1 = 
      sql.prepare("DELETE FROM loci WHERE group_id == :group_id ; " );
    
    stmt_loc_remove_group2 = 
      sql.prepare("DELETE FROM groups WHERE group_id == :group_id ; " );
    
    stmt_loc_lookup_name = 
      sql.prepare(" SELECT * FROM loci WHERE name == :name ; " );    
    
    stmt_loc_lookup_group = 
      sql.prepare(" SELECT * FROM loci WHERE group_id == :group_id ORDER BY chr,bp1; " );    
    
    stmt_loc_name_dump = 
      sql.prepare(" SELECT name FROM searchnames WHERE group_id == :group_id;");

    stmt_loc_lookup_group_and_name = 
      sql.prepare(" SELECT * FROM loci WHERE group_id == :group_id AND name == :name ; " );
    
    stmt_loc_lookup_id = 
      sql.prepare(" SELECT * FROM loci WHERE loc_id == :loc_id ; " );    
        
    stmt_loc_lookup_range = 
      sql.prepare(" SELECT * FROM loci WHERE chr == :chr AND bp1 <= :end AND bp2 >= :start ; " );
    
    stmt_loc_lookup_group_and_range = 
      sql.prepare(" SELECT * FROM loci WHERE group_id == :group_id "
		  " AND chr == :chr AND bp1 <= :end AND bp2 >= :start ORDER BY altname; " );
    
    stmt_loc_lookup_id_group_and_range = 
      sql.prepare(" SELECT * FROM loci WHERE group_id == :group_id "
		  " AND chr == :chr AND bp1 <= :end AND bp2 >= :start ; " );
    
    stmt_loc_insert = 
      sql.prepare("INSERT OR REPLACE INTO loci ( name, group_id, chr, bp1, bp2 , altname ) "
		  " values( :name , :group_id, :chr, :bp1, :bp2, :altname ) ; " );
    
    stmt_loc_intersect = 
      sql.prepare(" SELECT * FROM loci WHERE group_id == :group1_id "
		  " AND name IN ( SELECT name FROM loci WHERE group_id == :group2_id ) ; " );    
    
    stmt_loc_iterate =
      sql.prepare("SELECT * FROM loci ORDER BY chr,bp1;");
    
    stmt_loc_iterate_group =
      sql.prepare("SELECT * FROM loci WHERE group_id == :group_id ORDER BY chr,bp1 ;");
    
    stmt_loc_iterate_two_groups =
      sql.prepare("SELECT * FROM loci WHERE group_id == :group_id1 OR group_id == :group_id2 ORDER BY chr,bp1 ;");


    stmt_loc_fetch_altnames = 
      sql.prepare(" SELECT altname FROM loci WHERE group_id == :group_id AND chr == :chr AND bp1 <= :bp AND bp2 >= :bp ; " );

    stmt_loc_fetch_altnames_indel = 
      sql.prepare(" SELECT altname FROM loci WHERE group_id == :group_id "
		  " AND chr == :chr AND bp2 >= :start AND bp1 <= :stop ; " );

    stmt_loc_group_list = 
      sql.prepare("SELECT * FROM groups;");
    
    stmt_loc_iterate_overlap =
      sql.prepare("SELECT * FROM overlaps;");
    
    stmt_loc_overlap_insert = 
      sql.prepare("INSERT OR REPLACE INTO overlaps ( loc_id1, loc_id2, val_inter, val_union ) "
		  " values( :loc1, :loc2, :vint, :vun ) ; " );
    
    stmt_fetch_metatypes = 
      sql.prepare(" SELECT name , type , number, description "
		  " FROM metatypes ; " );
    
    stmt_loc_overlap_lookup = 
      sql.prepare("SELECT * FROM overlaps WHERE loc_id1 == :loc_id OR loc_id2 == :loc_id ;" );
    
    stmt_loc_subregion_insert =
      sql.prepare("INSERT OR REPLACE INTO subloci ( loc_id, name, chr, bp1, bp2, strand, frame ) "
		  " values ( :loc_id, :name, :chr, :bp1, :bp2 , :strand , :frame ) ; " );
    
    
    stmt_loc_subregion_lookup = 
      sql.prepare("SELECT * FROM subloci WHERE loc_id == :loc_id ; ");
    
    stmt_loc_lookup_group_with_overlap = 
      sql.prepare("SELECT * FROM loci a , loci b, overlaps o  "
		  "WHERE a.group_id == :group_id AND ( ( a.loc_id == o.loc_id1 AND b.loc_id == o.loc_id2 )"
		  "OR ( a.loc_id == o.loc_id2 AND b.loc_id == o.loc_id1 ) ) ; ");
    
    stmt_loc_lookup_group_with_overlap_p1 = 
      sql.prepare("SELECT * FROM loci a , overlaps o  "
		  "WHERE a.loc_id == o.loc_id1  "
		  "ORDER BY o.ROWID ;" );
    
    stmt_loc_lookup_group_with_overlap_p2 = 
      sql.prepare("SELECT * FROM loci a , overlaps o  "
		  "WHERE a.loc_id == o.loc_id2  "
		  "ORDER BY o.ROWID ;" );
    
    stmt_loc_meta_insert_type = 
      sql.prepare(" INSERT INTO metatypes (name,type,number,description) values( :name, :type, :number, :description ); " );
  

    stmt_loc_meta_insert = 
      sql.prepare(" INSERT OR REPLACE INTO loc_meta (loc_id,value)"
		  "  values( :reg_id, :value ); ");
    
    stmt_loc_submeta_insert = 
      sql.prepare(" INSERT OR REPLACE INTO subloc_meta (sub_id,value)"
		  "  values( :reg_id, :value ); ");
    
    stmt_loc_get_meta = 
      sql.prepare(" SELECT value FROM loc_meta WHERE loc_id == :loc_id " );
    
    stmt_loc_get_submeta = 
      sql.prepare(" SELECT value FROM subloc_meta WHERE sub_id == :sub_id ; " );
    
    stmt_loc_alias_insert = 
      sql.prepare(" INSERT OR REPLACE INTO aliases (group_id1,name1,group_id2,name2)"
		  "  values( :group_id1, :name1, :group_id2, :name2 ); ");
    
    stmt_loc_alias_group_insert = 
      sql.prepare(" INSERT OR IGNORE INTO alias_groups (group_name)"
		  "  values( :group_name ); ");
  
  stmt_loc_alias_lookup = 
    sql.prepare(" SELECT group_id2, name2 FROM aliases WHERE name1 == :name ;" );
  
  stmt_loc_group_alias_lookup = 
    sql.prepare(" SELECT group_id2, name2 FROM aliases WHERE group_id1 == :group_id AND name1 == :name ;" );

  stmt_loc_targetted_group_alias_lookup = 
    sql.prepare(" SELECT name2 FROM aliases WHERE group_id1 == :group_id1 AND name1 == :name AND group_id2 == :group_id2 ;" );

  stmt_loc_alias_group_dump = 
    sql.prepare(" SELECT * FROM alias_groups; " );



  // Set functions  
  
  stmt_set_group_lookup = 
    sql.prepare( "SELECT group_id FROM set_groups WHERE name == :name AND loc_group_id == :loc_group_id ; ");
  
  stmt_set_group_insert =
      sql.prepare("INSERT OR REPLACE INTO set_groups ( name, loc_group_id, temp, description ) "
		  " values( :name, :loc_group_id, :temp, :description ) ; ");
  
  stmt_set_member_lookup = 
    sql.prepare( "SELECT set_id FROM set_members WHERE name == :name AND group_id == :group_id ; ");
  
  stmt_set_names_fetch = 
    sql.prepare( "SELECT name FROM set_members WHERE group_id == :group_id ; ");
  
  stmt_set_names_and_id_fetch = 
    sql.prepare( "SELECT name,set_id FROM set_members WHERE group_id == :group_id ; ");

  stmt_set_data_dumper = 
    sql.prepare( "SELECT loc_id , set_id FROM set_data; ");

  stmt_set_members_fetch = 
    sql.prepare( "SELECT l.name FROM loci AS l , set_data AS sd "
		 " WHERE l.loc_id == sd.loc_id "
		 "   AND sd.set_id == :set_id ; " );

  stmt_set_members_fetch_regions = 
    sql.prepare( "SELECT l.loc_id,l.name,l.group_id,l.chr,l.bp1,l.bp2,l.altname FROM loci AS l , set_data AS sd "
		 " WHERE l.loc_id == sd.loc_id "
		 "   AND sd.set_id == :set_id ; " );

  stmt_set_member_insert =
    sql.prepare( " INSERT OR REPLACE INTO set_members ( group_id , name )  "
		 " values ( :group_id , :name ) ; " );
  
  stmt_set_data_insert =
    sql.prepare( "INSERT OR IGNORE INTO set_data ( set_id , loc_id ) "
		 "  values( :set_id , :loc_id ) ; " );

  stmt_loc_lookup_set = 
    sql.prepare(" SELECT * FROM loci WHERE loc_id IN ( SELECT loc_id FROM set_data WHERE set_id == :set_id ) ; " ); 

  stmt_dump_all_sets = 
    sql.prepare( "SELECT loc_group_id , name , description FROM set_groups ;");

  
  // Individual/segment functions

  stmt_insert_indiv =
    sql.prepare( "INSERT OR REPLACE INTO individuals ( name ) values ( :name ) ; " );
  
  stmt_lookup_indiv_id = 
    sql.prepare(" SELECT indiv_id FROM individuals WHERE name == :name ; " ) ;

  stmt_insert_segment =
    sql.prepare( "INSERT OR REPLACE INTO segments ( loc_id , indiv_id ) values ( :loc_id , :indiv_id ) ; " );

  stmt_fetch_segment =
    sql.prepare( "SELECT * FROM loci WHERE group_id == :group_id AND loc_id IN "
		 " ( SELECT loc_id FROM segments "
		 "      WHERE indiv_id == :indiv_id ) ; " );
  

  stmt_loc_fetch_id_given_name = 
    sql.prepare( "SELECT loc_id FROM loci WHERE group_id == :group_id AND name == :name ; ") ;


}

bool LocDBase::release()
{

  sql.finalise(stmt_loc_insert_group_name );
  sql.finalise(stmt_loc_lookup_group_name ); 
  sql.finalise(stmt_loc_lookup_group_id ); 
  sql.finalise(stmt_loc_update_temp_status ); 
  sql.finalise(stmt_loc_lookup_temp_status ); 
  sql.finalise(stmt_loc_remove_group1 );
  sql.finalise(stmt_loc_remove_group2 ); 
  sql.finalise(stmt_loc_lookup_name );
  sql.finalise(stmt_loc_lookup_group ); 
  sql.finalise(stmt_loc_lookup_group_and_name ); 
  sql.finalise(stmt_loc_group_list );

  sql.finalise(stmt_loc_fetch_id_given_name);

  sql.finalise(stmt_loc_lookup_real_name);
  sql.finalise(stmt_loc_lookup_real_name_only);

  sql.finalise( stmt_loc_replace_real_name );
  sql.finalise( stmt_loc_replace_real_name_alternate );

  sql.finalise( stmt_loc_alias_group_insert );
  sql.finalise( stmt_loc_alias_group_dump );

  sql.finalise(stmt_dump_all_sets);
  sql.finalise(stmt_loc_intersect );

  sql.finalise(stmt_loc_lookup_id );
  sql.finalise(stmt_loc_lookup_range ); 
  sql.finalise(stmt_loc_insert );
  sql.finalise(stmt_loc_iterate ); 
  sql.finalise(stmt_loc_iterate_group ); 
  sql.finalise(stmt_loc_iterate_overlap);
  sql.finalise(stmt_loc_overlap_insert);
  sql.finalise(stmt_loc_overlap_lookup);
  sql.finalise(stmt_loc_subregion_insert);
  sql.finalise(stmt_loc_subregion_lookup); 
  sql.finalise(stmt_loc_meta_insert_type );
  sql.finalise(stmt_loc_meta_insert );
  sql.finalise(stmt_loc_get_meta );

  sql.finalise(stmt_loc_alias_insert);
  sql.finalise(stmt_loc_alias_lookup);
  sql.finalise(stmt_loc_group_alias_lookup);

  sql.finalise(stmt_set_names_fetch);
  sql.finalise(stmt_set_members_fetch);
  sql.finalise(stmt_set_members_fetch_regions);

  sql.finalise(stmt_fetch_metatypes);
  sql.finalise(stmt_fetch_segment);
  sql.finalise(stmt_insert_indiv);
  sql.finalise(stmt_insert_segment);
  sql.finalise(stmt_loc_get_submeta);
  sql.finalise(stmt_loc_lookup_group_and_range);
  sql.finalise(stmt_loc_lookup_group_with_overlap_p1);
  sql.finalise(stmt_loc_lookup_group_with_overlap_p2);
  sql.finalise(stmt_loc_lookup_group_with_overlap);

  sql.finalise(stmt_loc_lookup_set);
  sql.finalise(stmt_loc_name_dump);
  sql.finalise(stmt_loc_submeta_insert);
  sql.finalise(stmt_loc_targetted_group_alias_lookup);
  sql.finalise(stmt_lookup_indiv_id);
  sql.finalise(stmt_set_data_insert);
  sql.finalise(stmt_set_group_insert);
  sql.finalise(stmt_set_group_lookup);
  sql.finalise(stmt_set_member_insert);
  sql.finalise(stmt_set_member_lookup);

  sql.finalise(stmt_loc_lookup_id_group_and_range ); 

}


void LocDBase::temporary(const uint64_t id, const bool b)
{ 
  sql.bind_int( stmt_loc_update_temp_status , ":temp" , (int)b );
  sql.step( stmt_loc_update_temp_status );
  sql.reset( stmt_loc_update_temp_status );
}  


bool LocDBase::temporary(const uint64_t id)
{
  bool r = false;
  sql.bind_int64( stmt_loc_lookup_temp_status , ":group_id" , id );
  if ( sql.step( stmt_loc_lookup_temp_status ) )
    r = sql.get_int( stmt_loc_lookup_temp_status , 0 );
  sql.reset( stmt_loc_lookup_temp_status );  
  return r;
}


void LocDBase::flush( const std::string & grp )
{
  uint64_t loc_group_id = lookup_group_id( grp );
  if ( loc_group_id == 0 ) return;
  flush( loc_group_id );
}

void LocDBase::flush()
{
  // Remove all groups that are temporary
  sql.query("DELETE FROM loci WHERE group_id IN ( SELECT group_id FROM groups WHERE temp == 1 ); ");
  sql.query("DELETE FROM groups WHERE temp == 1 ;");
}


void LocDBase::flush(const uint64_t group_id)
{
  sql.bind_int64( stmt_loc_remove_group1 , ":group_id", group_id );
  sql.bind_int64( stmt_loc_remove_group2 , ":group_id", group_id );
  sql.step( stmt_loc_remove_group1 );
  sql.step( stmt_loc_remove_group2 );
  sql.reset( stmt_loc_remove_group1 );
  sql.reset( stmt_loc_remove_group2 );
}
   

bool LocDBase::index()
{
  if ( ! attached() ) return false;
  
  sql.query( "CREATE INDEX IF NOT EXISTS groupPositionIndex ON loci(group_id,chr, bp1); " );  
  sql.query( "CREATE INDEX IF NOT EXISTS nameIndex ON loci(group_id,name);" );
  sql.query( "CREATE INDEX IF NOT EXISTS altNameIndex ON loci(group_id,altname);" );
  
  sql.query( "CREATE INDEX IF NOT EXISTS searchNameIdx ON searchnames(group_id); " );

  sql.query( "CREATE INDEX IF NOT EXISTS setmem ON set_members(group_id,name);");

  sql.query( "CREATE INDEX IF NOT EXISTS indivIndex ON segments(indiv_id); ");
  sql.query( "CREATE INDEX IF NOT EXISTS indivIndex2 ON individuals(name);");

  sql.query( "CREATE INDEX IF NOT EXISTS aliasIndex ON aliases(name1);");
  sql.query( "CREATE INDEX IF NOT EXISTS groupAliasIndex ON aliases(group_id1,name1);");

  sql.query( "CREATE INDEX IF NOT EXISTS subRegIndex ON subloci(loc_id); " );  
  sql.query( "CREATE INDEX IF NOT EXISTS overlapIndex ON overlaps(loc_id1, loc_id2);" );
  sql.query( "CREATE INDEX IF NOT EXISTS metaIndex ON loc_meta(loc_id);" );
  sql.query( "CREATE INDEX IF NOT EXISTS submetaIndex ON subloc_meta(sub_id);" );
 
  // Schema changed, so update prepared queries
  release();
  init();

  return true;
}


bool LocDBase::drop_index()
{
  if ( ! attached() ) return false;

  sql.query( "DROP INDEX IF EXISTS groupPositionIndex;");
  sql.query( "DROP INDEX IF EXISTS nameIndex;");
  sql.query( "DROP INDEX IF EXISTS altNameIndex;");

  sql.query( "DROP INDEX IF EXISTS searchNameIdx;" );

  sql.query( "DROP INDEX IF EXISTS indivIndex;");
  sql.query( "DROP INDEX IF EXISTS indiv2Index;");

  sql.query( "DROP INDEX IF EXISTS alias1Index;");
  sql.query( "DROP INDEX IF EXISTS alias2Index;");
  sql.query( "DROP INDEX IF EXISTS groupAlias1Index;");
  sql.query( "DROP INDEX IF EXISTS groupAlias2Index;");

  sql.query( "DROP INDEX IF EXISTS subRegIndex;");
  sql.query( "DROP INDEX IF EXISTS overlapIndex; ");
  sql.query( "DROP INDEX IF EXISTS metaIndex;");
  sql.query( "DROP INDEX IF EXISTS submetaIndex;");
  
  // Schema changed, so update prepared queries
  release();
  init();

  return true;
}


bool LocDBase::clear_overlaps()
{
  if ( ! attached() ) return false;
  sql.query( "DELETE FROM overlaps; ");
  return true;
}

bool LocDBase::clear_overlaps(uint64_t id1)
{
    sql.query( "DELETE FROM overlaps "
	       " WHERE loc_id1 IN ( SELECT loc_id FROM loci WHERE group_id == " + Helper::int2str(id1) + " ) "
	       "    OR loc_id2 IN ( SELECT loc_id FROM loci WHERE group_id == " + Helper::int2str(id1) + " ) ; " );

}

bool LocDBase::clear_overlaps(uint64_t id1, uint64_t id2)
{
    sql.query( "DELETE FROM overlaps "
	       " WHERE loc_id1 IN ( SELECT loc_id FROM loci "
	       "                    WHERE group_id == " + Helper::int2str(id1) + " "
	       "                       OR group_id == " + Helper::int2str(id2) + " ) "
	       "   AND loc_id2 IN ( SELECT loc_id FROM loci "
	       "                    WHERE group_id == " + Helper::int2str(id1) + " "
	       "                       OR group_id == " + Helper::int2str(id2) + " ) " );

}


Region LocDBase::construct_region( sqlite3_stmt * s  )
{

  // From loci table

  // 0  loc_id
  // 1  name
  // 2  group_id
  // 3  chr
  // 4  bp1
  // 5  bp2

  uint64_t id = sql.get_int64( s , 0 );
  int chr = sql.get_int( s , 3 ) ;
  int bp1 = sql.get_int( s , 4 ) ;
  int bp2 = sql.get_int( s , 5 ) ;

  std::string name = sql.get_text( s , 1 ) ; 
  std::string altname = sql.get_text( s , 6 ) ; 
  uint64_t gid = sql.get_int64( s , 2 ); 

  Region r(id,chr,bp1,bp2,name,altname,gid);
  
 
  /////////////////////////////////
  //                             //
  // Optionally, add subregions? //
  //                             //
  /////////////////////////////////


  if ( vget_subregions )
    {

      // From subloci table:
      
      // 0  sub_id
      // 1  loc_id
      // 2  name
      // 3  chr
      // 4  bp1
      // 5  bp2
      // 6  strand
      // 7  frame

	sql.bind_int64( stmt_loc_subregion_lookup, ":loc_id" , id );

	while ( sql.step( stmt_loc_subregion_lookup ) )
	{

	    uint64_t id = sql.get_int64( stmt_loc_subregion_lookup , 0 );
	    std::string name = sql.get_text( stmt_loc_subregion_lookup , 2 ) ;
	    int bp1 = sql.get_int( stmt_loc_subregion_lookup , 4 ) ;
	    int bp2 = sql.get_int( stmt_loc_subregion_lookup , 5 ) ;
	    int strand = sql.get_int( stmt_loc_subregion_lookup , 6 ) ;
	    int frame = sql.get_int( stmt_loc_subregion_lookup , 7 ) ;
	    
	    r.addSubRegion( id, name, chr, bp1, bp2 , strand , frame );
	    
  	    if ( vget_meta ) 
		r.subregion.back().meta = submeta( id );
	}
      sql.reset( stmt_loc_subregion_lookup );
    }



  /////////////////////////////////
  //                             //
  // Optionally, add meta-info?  //
  //                             //
  /////////////////////////////////
  
  if ( vget_meta )
       r.meta = meta( id );



  /////////////////////////////////
  //                             //
  // Return the final region     //
  //                             //
  /////////////////////////////////
  
  return r;

}


void LocDBase::add_overlap_table( uint64_t group1_id, uint64_t group2_id )
{
  
  // For each overlapping pair of regions, give intersection and union
  // distances From these, any subsequent definition of overlap can be
  // computed 
  

  // Iterate over all regions, sorted by position
  
  std::set<Region> current;

  int cnt = 0;

  sql.bind_int64( stmt_loc_iterate_two_groups , ":group_id1" , group1_id );
  sql.bind_int64( stmt_loc_iterate_two_groups , ":group_id2" , group2_id );

  sql.begin();

  while ( sql.step( stmt_loc_iterate_two_groups ) )
    {
      
      Region r = construct_region( stmt_loc_iterate_two_groups );
      
      std::set<Region>::iterator i = current.begin();
      
      // Perform all inserts within a single transaction

      while ( i != current.end() )
	{

	    if ( i->overlaps( r ) )
	    {
		// add overlap entry to list	      
		

		sql.bind_int64( stmt_loc_overlap_insert , ":loc1" , r.id );
		sql.bind_int64( stmt_loc_overlap_insert , ":loc2" , i->id );
		
		int a1 = r.start.position();
		int a2 = r.stop.position();
		
		int b1 = i->start.position();
		int b2 = i->stop.position();
		
		int consensusStart = a1 > b1 ? a1 : b1;
		int consensusStop = a2 < b2 ? a2 : b2;
		int intersection = consensusStop - consensusStart + 1;
		
		int unionStart = a1 < b1 ? a1 : b1;
		int unionStop = a2 > b2 ? a2 : b2;
		int vunion = unionStop - unionStart + 1;
		
		sql.bind_int( stmt_loc_overlap_insert , ":vint" , intersection );
		sql.bind_int( stmt_loc_overlap_insert , ":vun" , vunion );
		
		sql.step( stmt_loc_overlap_insert );
		sql.reset( stmt_loc_overlap_insert );
		
		// Advance to next possible segment
		
		++i;
		
	    }
	    else
	    {
	      // otherwise, we no longer need to consider this 
	      // region
		
	      std::set<Region>::iterator t = i;
	      ++i;
	      current.erase(t);
	    }
	}
      
      
      // Add this region to the pool
      
      current.insert(r);
      
    }
  
  
  sql.reset( stmt_loc_iterate_two_groups );
  
  sql.commit();
  
}



bool LocDBase::range_insertion(const Region & region , uint64_t indiv_id )
{
  
  sql.bind_text(  stmt_loc_insert , ":name"     , region.name               );
  sql.bind_int64( stmt_loc_insert , ":group_id" , (uint64_t)region.group    );
  sql.bind_int(   stmt_loc_insert , ":chr"      , region.start.chromosome() );
  sql.bind_int(   stmt_loc_insert , ":bp1"      , region.start.position()   );
  sql.bind_int(   stmt_loc_insert , ":bp2"      , region.stop.position()    );
  sql.bind_text(  stmt_loc_insert , ":altname"  , region.altname            );
  
  sql.step( stmt_loc_insert );
  sql.reset( stmt_loc_insert );
  
  // Subregions?

  uint64_t loc_id = sql.last_insert_rowid();
  int chr = region.start.chromosome();
  
  for (int s = 0 ; s < region.subregion.size(); s++)
    {
      
      sql.bind_int64( stmt_loc_subregion_insert , ":loc_id" , loc_id );
      sql.bind_text( stmt_loc_subregion_insert  , ":name"   , region.subregion[s].name );
      sql.bind_int( stmt_loc_subregion_insert   , ":chr"    , chr );
      sql.bind_int( stmt_loc_subregion_insert   , ":bp1"    , region.subregion[s].start.position() );
      sql.bind_int( stmt_loc_subregion_insert   , ":bp2"    , region.subregion[s].stop.position() );      
      sql.bind_int( stmt_loc_subregion_insert   , ":strand" , region.subregion[s].strand );
      sql.bind_int( stmt_loc_subregion_insert   , ":frame"  , region.subregion[s].frame );
      
      sql.step( stmt_loc_subregion_insert );
      sql.reset( stmt_loc_subregion_insert );

      uint64_t sub_id = sql.last_insert_rowid();
      
      // Subregion meta-information?
      insertMeta( stmt_loc_submeta_insert , region.subregion[s].meta , sub_id );
    }
  
  
  // Meta-information?
  
  insertMeta( stmt_loc_meta_insert , region.meta, loc_id );
  
  
  // Individual information ?

  if ( indiv_id )
    {
      sql.bind_int64( stmt_insert_segment , ":loc_id" , loc_id );
      sql.bind_int64( stmt_insert_segment , ":indiv_id" , indiv_id );
      sql.step( stmt_insert_segment );
      sql.reset( stmt_insert_segment );
    }
  
  
  return true;
}


void LocDBase::insertMeta( sqlite3_stmt * s , const MetaInformation<LocMeta> & meta, const int id )
{
  
  // Just insert as a string/INFO format, e.g.  "ID=name1;VALIDATED;SC=0.99,0.1"
  // But track meta-types in LocDB as before
  
  // 1) Insert any unknown types in LOCDB register 
  
  std::vector<std::string> keys = meta.keys();  

  for (int j=0; j<keys.size(); j++)
    {
      meta_index_t midx = MetaInformation<LocMeta>::field( keys[j] );
      
      if ( midx.mt == META_UNDEFINED ) 
	{
	
	  MetaInformation<LocMeta>::field( keys[j] , META_TEXT , -1 , "undeclared tag");
	  
	  sql.bind_text( stmt_loc_meta_insert_type , ":name"        , keys[j] );
	  sql.bind_int(  stmt_loc_meta_insert_type , ":type"        , midx.mt );
	  sql.bind_int(  stmt_loc_meta_insert_type , ":number"      , midx.len );
	  sql.bind_text( stmt_loc_meta_insert_type , ":description" , midx.description );

	  sql.step( stmt_loc_meta_insert_type );		
	  sql.reset( stmt_loc_meta_insert_type );

	}
    }
  
  
  // 2) Insert actual meta-information as a single string
  
  std::stringstream ss;
  ss << meta;
  
  // NOTE: Need to allocate ssStr here on the stack so that it stays in scope until the step/reset that actually use its contents (since bind_text takes it as a by-reference parameter):
  std::string ssStr = ss.str();

  sql.bind_int(  s , ":reg_id" , id );
  sql.bind_text( s , ":value"  , ssStr );
  sql.step( s );
  sql.reset( s );	  
    
}



uint64_t LocDBase::set_group_id(const std::string & grp, const bool temp , const std::string & desc )
{

  uint64_t group_id = 0;

  sql.bind_text(stmt_loc_lookup_group_name, ":name" , grp ); 

  if ( sql.step( stmt_loc_lookup_group_name ) ) 
    {
      group_id = sql.get_int64( stmt_loc_lookup_group_name , 0 ) ;
      sql.reset( stmt_loc_lookup_group_name );
    }
  else
    {
      sql.reset( stmt_loc_lookup_group_name );

      sql.bind_text( stmt_loc_insert_group_name , ":name" , grp );
      sql.bind_int( stmt_loc_insert_group_name , ":temp" , 1 );
      sql.bind_text( stmt_loc_insert_group_name , ":description" , desc );
      sql.step( stmt_loc_insert_group_name );      
      group_id = sql.last_insert_rowid();
      sql.reset( stmt_loc_insert_group_name );
    }
  
  return group_id;
}

uint64_t LocDBase::lookup_group_id(const std::string & grp)
{
  uint64_t group_id = 0;
  if ( !attached() ) return 0;
  sql.bind_text(stmt_loc_lookup_group_name, ":name" , grp );     
  if ( sql.step( stmt_loc_lookup_group_name ) ) 
    group_id = sql.get_int64( stmt_loc_lookup_group_name , 0 ) ;    
  sql.reset( stmt_loc_lookup_group_name );
  return group_id;
}

std::string LocDBase::lookup_group_id(const int group_id)
{
  std::string grp = "";
  sql.bind_int(stmt_loc_lookup_group_id, ":group_id" , group_id );     
  if ( sql.step( stmt_loc_lookup_group_id ) ) 
    grp = sql.get_text( stmt_loc_lookup_group_id , 0 ) ;    
  sql.reset( stmt_loc_lookup_group_id );
  return grp;
}


uint64_t LocDBase::set_set_id(const std::string & name, 
			      const int loc_group_id , 
			      const bool temp , 
			      const std::string & desc )
{
  
  uint64_t group_id = 0;
  
  sql.bind_text( stmt_set_group_lookup, ":name" , name );
  sql.bind_int( stmt_set_group_lookup, ":loc_group_id" , loc_group_id );

  if ( sql.step( stmt_set_group_lookup ) ) 
    {
      group_id = sql.get_int64( stmt_set_group_lookup , 0 ) ;
      sql.reset( stmt_set_group_lookup );
    }
  else
    {
      sql.reset( stmt_set_group_lookup );

      sql.bind_text( stmt_set_group_insert , ":name" , name );
      sql.bind_int( stmt_set_group_insert , ":loc_group_id" , loc_group_id );
      sql.bind_int( stmt_set_group_insert , ":temp" , 1 );
      sql.bind_text( stmt_set_group_insert , ":description" , desc );
      sql.step( stmt_set_group_insert );
      group_id = sql.last_insert_rowid();
      sql.reset( stmt_set_group_insert );
    }
  
  return group_id;
}


uint64_t LocDBase::lookup_set_id(const std::string & grp, const std::string & name )
{
  uint64_t loc_group_id = lookup_group_id( grp );
  if ( loc_group_id == 0 ) return 0;

  sql.bind_text( stmt_set_group_lookup, ":name" , name );
  sql.bind_int( stmt_set_group_lookup, ":loc_group_id" , loc_group_id );

  uint64_t group_id = 0;
  
  if ( sql.step( stmt_set_group_lookup ) ) 
    group_id = sql.get_int64( stmt_set_group_lookup , 0 ) ;      

  sql.reset( stmt_set_group_lookup );

  return group_id;
}


uint64_t LocDBase::load_GTF( const std::string & filename, const std::string & grp, bool use_transcript_id)
{

  if ( ! attached() ) Helper::halt( "no LOCDB attached" );

  if ( ! Helper::fileExists( filename ) ) return 0;
  
  InFile f( filename );
  
  // Expect GTF2.2 format, as documented 
  //  http://mblab.wustl.edu/GTF22.html
  
  // Register meta-types  

  uint64_t group_id = set_group_id( grp );
  

  /////////////////////////////////////////
  //                                     //
  // Begin SQL transaction               //
  //                                     //
  /////////////////////////////////////////
  
  sql.begin();

  int inserted = 0;
  
  std::set<std::string> names;
  

  /////////////////////////////////////////
  //                                     //
  // Process each row of input           //
  //                                     //
  /////////////////////////////////////////

  while ( ! f.eof() )
    {

      // tab-delimited line

      std::vector<std::string> tok = Helper::char_split( f.readLine() , '\t' );
      
      if ( tok.size() == 0 ) continue;
      
      // Should contain exactly 9 tab-delimited elements      

      if ( tok.size() != 9 ) 
	{
	  plog.warn("line in GTF column count != 9");
	  continue;
	}
      
      
      // <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]

      // AB000381 Twinscan  CDS          380   401   .   +   0  gene_id "001"; transcript_id "001.1";

            
      // Start/stop positions

      int p1,p2;
      if ( ! Helper::str2int( tok[ 3 ] , p1 ) ) continue;
      if ( ! Helper::str2int( tok[ 4 ] , p2 ) ) continue;
      
      // Name (from gene_id

      std::vector<std::string> tok2 = Helper::char_split( tok[8] , ' ' , false );
      
      // requires atleast 4 manadatory fields: gene_id XXX transcript_id XXX
      
      if ( tok2.size() < 4 ) 
	{
	  plog.warn("badly formed GTF, col 9");
	  continue;
	}
      
      // remove quotes and last semi-colon

      if ( tok2[0] != "gene_id" ) plog.warn("expecting gene_id in GTF :" + tok2[0] );
      if ( tok2[2] != "transcript_id" ) plog.warn("expecting transcript_id in GTF :" + tok2[2]);

      if ( tok2[1].substr( tok2[1].size()-1 ) == ";" ) tok2[1] = tok2[1].substr( 0, tok2[1].size()-1 );
      else plog.warn("no ; after gene_id in GTF");

      if ( tok2[3].substr( tok2[3].size()-1 ) == ";" ) tok2[3] = tok2[3].substr( 0, tok2[3].size()-1 );
      else plog.warn("no ; after transcript_id in GTF");
      
      tok2[1] = Helper::unquote( tok2[1] );
      tok2[3] = Helper::unquote( tok2[3] );

      // tok2[0] should equal  gene_id
      // use either gene-name or transcript name as 'name' 
      
      std::string name = use_transcript_id ? tok2[3] : tok2[1];
      
      int chromosome = Helper::chrCode( tok[0] ) ;
      
      if ( chromosome == 0 ) continue;
	  
      Region r( chromosome ,
		p1 , p2 , 
		name , 
		(int)group_id ); 

      
      // Track gene-name, if unique name is a transcript
      
      if ( use_transcript_id ) r.altname = tok2[1];
      else r.altname = tok2[3]; 
      

      // Always expect gene_id and transcript_id 
      // Use gene_id

      std::string & feature = tok[2];

      r.meta.set( "feature" , feature );
      

      // Track strand for CDS (-1,1,  0=missing)
      
      // Use 'strand' to encode two distinct things: a) whether
      //  transcript is on +ve or -ve strand, b) but also some
      //  meta-data about the type of feature. Not particularly
      //  elegant, but will work for now...

      //         <other>     0
      //         CDS         -1,+1
      //         exon        -2,+2
      //         start_codon -3,+3
      //         stop_codon  -4,+4

      int strand = 0;
      if      ( tok[6] == "-" ) strand = -1;
      else if ( tok[6] == "+" ) strand =  1;

      if      ( feature == "exon" ) strand *= 2;
      else if ( feature == "start_codon" ) strand *= 3;
      else if ( feature == "stop_codon" ) strand *= 4;

      r.meta.set( PLINKSeq::TRANSCRIPT_STRAND() , strand );
      
      
      //
      // Frame: for CDS 0, 1 or 2 
      //   for start_codon
      //

      int frame = 0;
      
      if ( Helper::str2int( tok[7] , frame ) ) 
	{	
	  // GTF encodes frame 021,021 (as number of bases to next start of codon (inc. this one)
	  if ( frame == 2 ) frame = 1;
	  else if ( frame == 1 ) frame = 2;
	  r.meta.set( PLINKSeq::TRANSCRIPT_FRAME() , frame );
	}

      ///////////////////////////
      // Add region to database
      
      range_insertion(r);

      ++inserted;

    }




  /////////////////////////////////////////
  //                                     //
  // Finish transaction                  //
  //                                     //
  /////////////////////////////////////////


  sql.commit();
  
  f.close();

  plog << "inserted " << inserted << " rows\n";

  populate_meta_field_map();


  /////////////////////////////////////////
  //                                     //
  // Return group ID                     //
  //                                     //
  /////////////////////////////////////////
  
  return group_id;


}


void LocDBase::populate_searchname_table( const std::string & g , bool altname )
{
  
  if ( ! attached() ) return;
  int grp_id = lookup_group_id( g );
  if ( grp_id == 0 ) return;

  // clear out any existing entries
  
  sql.query( "DELETE FROM searchnames WHERE group_id == " + Helper::int2str( grp_id ) + " ; " ) ;

  // get new list of genes to add in

  sqlite3_stmt * s = altname ? 
    sql.prepare( "SELECT altname FROM loci WHERE group_id == " + Helper::int2str( grp_id ) + " ; " ) :
    sql.prepare( "SELECT    name FROM loci WHERE group_id == " + Helper::int2str( grp_id ) + " ; " ) ;

  std::set<std::string> names;
  while ( sql.step(s) )
    {
      std::string n = sql.get_text( s , 0 ); 
      names.insert( n );
    }
  sql.reset(s);
  sql.finalise(s);


  // and now insert them

  sql.begin();
  s = sql.prepare( "INSERT OR IGNORE INTO searchnames ( group_id , name ) values ( :group_id , :name ) ; " );
  sql.bind_int64( s , ":group_id" , grp_id );  
  
  std::set<std::string>::iterator ii = names.begin();
  while ( ii != names.end() )
    {
      sql.bind_text( s , ":name" , *ii );
      sql.step(s);
      sql.reset(s);
      ++ii;
    }
  sql.finalise(s);
  sql.commit();

  plog << "inserted " << names.size() << " into the name-table\n";

}



uint64_t LocDBase::load_GFF( const std::string & filename, const std::string & grp, const std::string & name  )
{
  
  Helper::halt("GFF support not yet implemented");

  // if name == "" || ".", then use 'source' as the set name
  // otherwise, look for a meta-field named 'name' and use that as the name
  
  if ( ! attached() ) Helper::halt( "no LOCDB attached" );
  
  if ( ! Helper::fileExists( filename ) ) return 0;
  
  InFile f( filename );
  
  // Expect GFF v2 format, as documented 
  //  http://www.sanger.ac.uk/resources/software/gff/spec.html

  // Tab-delimited
  // <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
  
  // expect 8 fields always; no whitespace

  // missing data = '.'
  // sequence numbering starts at 1
  // if score missing, '.'
  // strand: + - .
  // frame: same as GTF
  
  // comments -- everything after # is ignored
  // note -- # could come at start, which means whole line is ignored;
  //         OR after first mandatory fields
  //  typically may have special ## lines as the header
   // comments, semi-colon delimited
  //   can be key 
  //    or    key=value
  //    or    key val1 val2 
  // Free test must be within double-quotes

    
  // Register meta-types
    
  uint64_t group_id = set_group_id( grp );
  
  /////////////////////////////////////////
  //                                     //
  // Begin SQL transaction               //
  //                                     //
  /////////////////////////////////////////
  
  sql.begin();

  int inserted = 0;
  
  
  /////////////////////////////////////////
  //                                     //
  // Process each row of input           //
  //                                     //
  /////////////////////////////////////////

  while ( ! f.eof() )
    {

      // tab-delimited line

      std::vector<std::string> tok = Helper::char_split( f.readLine() , '\t' );

      // empty line, or comment line?
      if ( tok.size() == 0 || tok[0].substr(0,1) == "#" ) continue;
      
      // Should contain at least 8 tab-delimited elements      
      if ( tok.size() < 8 ) 
	{
	  plog.warn( "invalid GFF entry", tok );
	  continue;
	}

      // GFF2 format:
      // <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]

       
      // chr5    nhgri_gwas_catalog      sequence_variant        150240076       150240076       3E-7    +       .       ID=ngc00001;Note=Crohn's disease rs1000113*T 3E-7 17554300;rsid=rs1000113;pmid=17554300;link=http://www.ncbi.nlm.nih.gov/pubmed/17554300
      
      // Start/stop positions

      int p1,p2;
      if ( ! Helper::str2int( tok[ 3 ] , p1 ) ) continue;
      if ( ! Helper::str2int( tok[ 4 ] , p2 ) ) continue;
      
      
      // Name (from gene_id)

      std::vector<std::string> tok2 = Helper::char_split( tok[8] , ' ' , false );
            
//       std::string name = use_transcript_id ? tok2[3] : tok2[1];
      
//       // Chromosome
//       int chromosome = Helper::chrCode( tok[0] ) ;      
//       if ( chromosome == 0 ) continue;
	  
//       Region r( chromosome ,
// 		p1 , p2 , 
// 		name , 
// 		(int)group_id ); 

//       // Track gene-name, if unique name is a transcript
      
//       if ( use_transcript_id ) r.altname = tok2[1];
//       else r.altname = tok2[3]; 


//       // Always expect gene_id and transcript_id 
//       // Use gene_id

//       r.meta.set( "source" , tok[1] );
//       r.meta.set( "feature" , tok[2] );
      

//       // Track strand and frame

//       r.meta.set( PLINKSeq::TRANSCRIPT_STRAND() , tok[6] );
      
//       int frame = 0;

//       // implies "." --> 0  (i.e. for stop_codon )

//       if ( Helper::str2int( tok[7] , frame ) ) 
// 	r.meta.set( PLINKSeq::TRANSCRIPT_FRAME() , frame );
      

//       ///////////////////////////
//       // Add region to database
      
//       range_insertion(r);

//       ++inserted;

    }


  /////////////////////////////////////////
  //                                     //
  // Finish transaction                  //
  //                                     //
  /////////////////////////////////////////


  sql.commit();
  
  f.close();

  plog << "inserted " << inserted << " rows\n";

  populate_meta_field_map();


  /////////////////////////////////////////
  //                                     //
  // Return group ID                     //
  //                                     //
  /////////////////////////////////////////
  
  return group_id;


}



uint64_t LocDBase::load_regions( const std::string & filename,
				 const std::string & grp, 
				 int col_pos,
				 int col_chr, 
				 int col_bp1, 
				 int col_bp2,
				 int col_name, 				 
				 int col_sub,
				 int col_meta,
				 int col_indiv,
				 std::map<std::string,int> * meta )
{

  if ( ! attached() ) Helper::halt( "no LOCDB attached" );

  bool subregions = col_sub >= 0;
  bool names = col_name >= 0;
  bool chr = col_chr >= 0;
  bool bp1 = col_bp1 >= 0;
  bool bp2 = col_bp2 >= 0;
  bool individuals = col_indiv >= 0;
  bool single_locus_pos = col_pos >= 0;

  int maxcol = col_name;
  if ( names && col_name > maxcol ) maxcol = col_name;
  if ( single_locus_pos )
    {
      if ( col_pos > maxcol ) maxcol = col_pos;
    }
  else
    {
      if ( col_chr > maxcol ) maxcol = col_chr;
      if ( col_bp1 > maxcol ) maxcol = col_bp1;
      if ( col_bp2 > maxcol ) maxcol = col_bp2;
    }
  if ( subregions && col_sub > maxcol ) maxcol = col_sub;
  if ( col_meta > maxcol ) maxcol = col_meta;


  // user defined meta-data (i.e. col name == "META")
  // and format is assumed to be String (way to specify in advance?)
  
  bool free_meta = col_meta >= 0 ;

  // fixed field meta columns?
  if ( meta )
    {
      std::map<std::string,int>::iterator i = meta->begin();
      while ( i != meta->end() )
	{
	  if ( i->second > maxcol ) 
	    maxcol = i->second;
	  ++i;
	}
    }
  ++maxcol;
  

  // No data specified? 
  if ( maxcol == 0 ) 
      return 0;


  InFile f( filename );
  
  /////////////////////////////////////////
  //                                     //
  // Begin SQL transaction               //
  //                                     //
  /////////////////////////////////////////

  uint64_t group_id = set_group_id( grp );

  sql.begin();

  int inserted = 0;
  int inserted2 = 0;


  /////////////////////////////////////////
  //                                     //
  // Process each row of input           //
  //                                     //
  /////////////////////////////////////////

  while ( ! f.eof() )
    {
      
      strList buffer = f.tokenizeLine();

      if( buffer.size() < maxcol ) 
	continue;

      // Comment line?  
      if ( buffer[0].substr(0,1) == "#" ) 
	continue;
      
      int chromosome = -1;
      int p1 = -1;
      int p2 = -1;

      if ( single_locus_pos )
	{
	  bool okay = true;
	  Region t( buffer[ col_pos ] , okay );
	  if ( ! okay ) continue;
	  chromosome = t.chromosome() ;
	  p1 = t.start.position();
	  p2 = t.stop.position();
	}
      else
	{
	  // CHR, BP1 and BP2 specified separately
	  if ( chr )
	    chromosome = Helper::chrCode( buffer[ col_chr ] ) ;	    	  
	  if ( bp1 && ! Helper::str2int( buffer[ col_bp1 ] , p1 ) ) continue;
	  if ( bp2 && ! Helper::str2int( buffer[ col_bp2 ] , p2 ) ) continue;
	}

      std::string name = names ? buffer[ col_name ] : "-";      
      
      
      Region r( chromosome ,
		p1 , p2 , 
		name , 
		(int)group_id ); 



      ////////////////////
      // Sub-regions
      
      if ( subregions )
	{
	  /// Assume 100-110;200-210;300-305
	  // For now, no range checking wrt to parent region
	  // Also, automatically assign to same chromosome
	  
	  
	  std::vector<std::string> tok = Helper::char_split( buffer[col_sub] , '-' , ';' );

	  // Should contain even number of elements
	  
	  if ( tok.size() % 2 != 0 ) 
	    {
	      plog.warn("badly formed subregion info");
	      continue;
	    }

	  for (int s=0; s<tok.size(); s+=2)
	    {
	      int p1, p2;
	      if ( ! Helper::str2int( tok[s] , p1 ) ) continue;
	      if ( ! Helper::str2int( tok[s+1] , p2 ) ) continue;	      
	      r.addSubRegion( chromosome, p1, p2 );
	    }
	  
	  inserted2 += (int)(tok.size() / 2 );
	}



      ////////////////////
      // Meta-information
      
      // Fixed fields

      if ( meta )
	{
	  std::map<std::string,int>::iterator m = meta->begin();
	  while ( m != meta->end() )
	    {
	      r.meta.set( m->first , buffer[ m->second ] );
	      ++m;
	    }

	}
      
      
      // free fields
      
      if ( free_meta )
	{

	  // expecting key=value,value;key=value;key  format
	  
	  const std::string & s = buffer[ col_meta ];
	  
	  if ( s == "." || s == "" ) continue;
	  
	  r.meta.parse( s , ';' , true ); 
	  
	}
 
    
      ///////////////////////////
      // Add region to database
      
      if ( individuals ) 
	range_insertion( r , insert_indiv( buffer[ col_indiv ] ) );
      else 
	range_insertion(r);

      ++inserted;

    }


  /////////////////////////////////////////
  //                                     //
  // Finish transaction                  //
  //                                     //
  /////////////////////////////////////////

  
  sql.commit();
  
  plog << "inserted " << inserted << " rows, " << inserted2 << " subregions\n";
  
  f.close();

  populate_meta_field_map();

  /////////////////////////////////////////
  //                                     //
  // Return group ID                     //
  //                                     //
  /////////////////////////////////////////
  
  return group_id;
  
}


std::set<Region> LocDBase::get_overlaps(uint64_t loc_id)
{
  OverlapDefinition o;
  return get_overlaps(loc_id, o);
}


std::set<Region> LocDBase::get_overlaps(uint64_t loc_id, OverlapDefinition & od)
{

  // OLD CODE...

  std::set<Region> regions;

  sql.bind_int64( stmt_loc_overlap_lookup, ":loc_id" , loc_id );
  
  while ( sql.step( stmt_loc_overlap_lookup ) )
    {

      uint64_t id1 = sql.get_int64( stmt_loc_overlap_lookup , 0 ) ;
      uint64_t id2 = sql.get_int64( stmt_loc_overlap_lookup , 1 ) ;

      if ( id1 == loc_id ) 
	id1 = id2;
      
      sql.bind_int64( stmt_loc_lookup_id, ":loc_id" , id1 );

      if ( sql.step( stmt_loc_lookup_id ) )
	regions.insert( construct_region( stmt_loc_lookup_id ) );
      
      sql.reset( stmt_loc_lookup_id );

    }
  
  sql.reset( stmt_loc_overlap_lookup );

  return regions;
  
}


bool LocDBase::get_regions_and_overlap( void (*f)( Region&,Region&, int, int, void * ) ,
					void * data )
{

  while ( 1 ) 
    {
            
      // loc_id    0   
      // name      1   
      // group_id  2   
      // chr       3   
      // bp1       4   
      // bp2       5   
      
      // loc_id1   8
      // loc_id2   9
      // val_inter 10
      // val_union 11
      
      bool next1 = sql.step( stmt_loc_lookup_group_with_overlap_p1 );
      bool next2 = sql.step( stmt_loc_lookup_group_with_overlap_p2 );
      
      if ( ! ( next1 && next2 ) ) break;
      
      Region r1 = construct_region( stmt_loc_lookup_group_with_overlap_p1 );
      Region r2 = construct_region( stmt_loc_lookup_group_with_overlap_p2 );
      
      int v_int = sql.get_int( stmt_loc_lookup_group_with_overlap_p1 , 8 );
      int v_union = sql.get_int( stmt_loc_lookup_group_with_overlap_p1 , 9 );
      
      // Call user-defined function
      
      f( r1, r2, v_int, v_union , data );
      
    }
  
  sql.reset( stmt_loc_lookup_group_with_overlap_p1 );
  sql.reset( stmt_loc_lookup_group_with_overlap_p2 );
  
  return true;
}


Region LocDBase::get_region(const std::string & g, const std::string & genename )
{
  return get_region( lookup_group_id( g ) , genename );
}

Region LocDBase::get_region(const int group_id , const std::string & genename )
{
  Region r;
  if ( group_id == 0 ) return r;
  sql.bind_int64( stmt_loc_lookup_group_and_name, ":group_id" , group_id );
  sql.bind_text( stmt_loc_lookup_group_and_name, ":name" , genename );
  
  if ( sql.step( stmt_loc_lookup_group_and_name ) )
    r = construct_region( stmt_loc_lookup_group_and_name );
  sql.reset( stmt_loc_lookup_group_and_name );  
  return r;
}
 

Region LocDBase::get_region( const uint64_t loc_id ) 
{
  // lookup from a single ID
  
  sql.bind_int64( stmt_loc_lookup_id, ":loc_id" , loc_id );
  
  Region r;

  if ( sql.step( stmt_loc_lookup_id ) )
    r = construct_region( stmt_loc_lookup_id ) ;
    
  sql.reset( stmt_loc_lookup_id );

  return r;
  
}

std::set<Region> LocDBase::get_regions( const std::string & grp )
{
  std::set<Region> s;
  if ( ! attached() ) return s;
  int gid = lookup_group_id( grp );
  if ( gid == 0 ) return s;
  return get_regions(gid);
}


std::set<Region> LocDBase::get_regions(uint64_t gid)
{

  std::set<Region> regions;
  
  sql.begin();

  sql.bind_int64( stmt_loc_lookup_group, ":group_id" , gid );
  while ( sql.step( stmt_loc_lookup_group ) )
      {
	  regions.insert( construct_region( stmt_loc_lookup_group ) );
      }

  sql.commit();

  sql.reset( stmt_loc_lookup_group );

  return regions;
}


uint64_t LocDBase::insert_indiv( const std::string & indiv_id )
{
  if ( indmap.find( indiv_id ) != indmap.end() ) return indmap[ indiv_id ];
  sql.bind_text( stmt_insert_indiv , ":name", indiv_id );
  sql.step( stmt_insert_indiv );
  sql.reset( stmt_insert_indiv );
  uint64_t i = lookup_indiv_id( indiv_id ); 
  indmap[ indiv_id ] = i;
  return i;
}

uint64_t LocDBase::lookup_indiv_id( const std::string & indiv_id ) 
{
  
  if ( indmap.find( indiv_id ) != indmap.end() ) return indmap[ indiv_id ];  
  uint64_t i = 0;
  sql.bind_text( stmt_lookup_indiv_id , ":name" , indiv_id );
  if ( sql.step( stmt_lookup_indiv_id ) ) 
    i = sql.get_int64( stmt_lookup_indiv_id, 0 );
  sql.reset( stmt_lookup_indiv_id );
  indmap[ indiv_id ] = i;
  return i;
}


void LocDBase::insert_segment( const std::string & indiv_id , Region & segment )
{
  
}


std::set<Region> LocDBase::get_indiv_regions( const std::string & group , const std::string & person )
{
  uint64_t grp_id = lookup_group_id( group );
  return get_indiv_regions( grp_id , person );
}

std::set<Region> LocDBase::get_indiv_regions( const uint64_t grp_id , const std::string & person )
{
  std::set<Region> segs;
  if ( grp_id == 0 ) return segs;

  uint64_t ind_id = lookup_indiv_id( person );
  if ( ind_id == 0 ) return segs;

  return get_indiv_regions( grp_id , ind_id );

}

std::set<Region> LocDBase::get_indiv_regions( uint64_t group_id , uint64_t indiv_id )
{
  
  std::set<Region> segments;
  
  sql.bind_int64( stmt_fetch_segment , ":group_id" , group_id );
  sql.bind_int64( stmt_fetch_segment , ":indiv_id" , indiv_id );

  while ( sql.step( stmt_fetch_segment ) )
    {
      segments.insert( construct_region( stmt_fetch_segment ) );
    }
  sql.reset( stmt_fetch_segment );
  
  return segments;
}


std::vector<uint64_t> LocDBase::get_region_ids(uint64_t gid, int chr, int bp1, int bp2)
{
  std::vector<uint64_t> regions;
  sql.begin();
  sql.bind_int64( stmt_loc_lookup_id_group_and_range, ":group_id" , gid );
  sql.bind_int( stmt_loc_lookup_id_group_and_range, ":chr" , chr );
  sql.bind_int( stmt_loc_lookup_id_group_and_range, ":start" , bp1 );
  sql.bind_int( stmt_loc_lookup_id_group_and_range, ":end" , bp2 );
  
  while ( sql.step( stmt_loc_lookup_id_group_and_range ) )
    {
      regions.push_back( sql.get_int64( stmt_loc_lookup_id_group_and_range , 0  ) );
    }
  sql.reset( stmt_loc_lookup_id_group_and_range );
  
  sql.commit();
  return regions;
}

bool LocDBase::contains( const std::string & group , const int chr , const int bp1 , const int bp2 )
{
  if ( ! attached() ) return false;
  uint64_t id = lookup_group_id( group );
  if ( id == 0 ) return false;
  return contains( id , chr , bp1 , bp2 );
  
}

bool LocDBase::contains( const int gid , const int chr , const int bp1 , const int bp2 )
{
  sql.bind_int64( stmt_loc_lookup_group_and_range, ":group_id" , gid );
  sql.bind_int( stmt_loc_lookup_group_and_range, ":chr" , chr );
  sql.bind_int( stmt_loc_lookup_group_and_range, ":start" , bp1 );
  sql.bind_int( stmt_loc_lookup_group_and_range, ":end" , bp2 );
  bool found = false;
  while ( sql.step( stmt_loc_lookup_group_and_range ) )
    {
      found = true;
      break;
    }
  sql.reset( stmt_loc_lookup_group_and_range );
  return found; 
}

std::set<Region> LocDBase::get_regions( const std::string & group , const Region & r )
{
  std::set<Region> reg;
  uint64_t id = lookup_group_id( group );
  if ( id == 0 ) return reg;
  return get_regions( id , r.chromosome() , r.start.position() , r.stop.position() ) ;
}


std::set<Region> LocDBase::get_regions( const std::string & group , const Variant & v )
{
  std::set<Region> r;
  uint64_t id = lookup_group_id( group );
  if ( id == 0 ) return r;
  return get_regions( id , v.chromosome() , v.position() , v.stop() == 0 ? v.position() : v.stop() );
}


std::set<Region> LocDBase::get_regions(uint64_t gid, int chr, int bp1, int bp2)
{
  // Get all regions that span the specified region
  
  std::set<Region> regions;
  
  sql.begin();
  
  sql.bind_int64( stmt_loc_lookup_group_and_range, ":group_id" , gid );
  sql.bind_int( stmt_loc_lookup_group_and_range, ":chr" , chr );
  sql.bind_int( stmt_loc_lookup_group_and_range, ":start" , bp1 );
  sql.bind_int( stmt_loc_lookup_group_and_range, ":end" , bp2 );
  
  while ( sql.step( stmt_loc_lookup_group_and_range ) )
    {
      regions.insert( construct_region( stmt_loc_lookup_group_and_range ) );
    }
  
  sql.commit();
  sql.reset( stmt_loc_lookup_group_and_range );
  return regions;
}


int LocDBase::count(uint64_t grp)
{
  std::vector<int> r;

  if ( grp == -1 )
    r = sql.intTable( "SELECT count(*) FROM loci;" , 1 );
  else
    r = sql.intTable( "SELECT count(*) FROM loci WHERE group_id == " + Helper::int2str(grp) +" ;" , 1 );

  if ( r.size() == 1 ) 
    return r[0];
  return -1;
}



uint64_t LocDBase::span(uint64_t grp)
{
  std::vector<uint64_t> r;
  
  if ( grp == -1 )
    r = sql.int64Table( "SELECT sum(bp2-bp1) FROM loci;" , 1 );
  else
    r = sql.int64Table( "SELECT sum(bp2-bp1) FROM loci WHERE group_id == " + Helper::int2str(grp) +" ;" , 1 );
  
  if ( r.size() == 1 ) 
    return r[0];

  // TODO: used to be -1, but was clearly incorrect. Check this function is not redundant; check that
  // 0 return code is okay as NA.
  return 0;
}



MetaInformation<LocMeta> LocDBase::meta( uint64_t loc_id )
{

  MetaInformation<LocMeta> m;
  
  sql.bind_int64(stmt_loc_get_meta, ":loc_id" , loc_id ); 
  while ( sql.step( stmt_loc_get_meta ) )
    {
      std::string value = sql.get_text( stmt_loc_get_meta , 0 );
      m.parse( value ); 
    }
  
  sql.reset( stmt_loc_get_meta );

  return m;
}


MetaInformation<LocMeta> LocDBase::submeta( uint64_t sub_id )
{

  MetaInformation<LocMeta> m;
  
  sql.bind_int64(stmt_loc_get_submeta, ":sub_id" , sub_id ); 
	
  while ( sql.step( stmt_loc_get_submeta ) )
    {      
      std::string value = sql.get_text( stmt_loc_get_submeta , 0 );     
      m.parse( value ); 
    }
  
  sql.reset( stmt_loc_get_submeta );

  return m;
}


uint64_t LocDBase::merge( const std::string & grp_name, const std::string & name, const std::string & mergeField)
{

  if ( ! attached() ) return 0;

  //
  // Does this group exist?
  //

  uint64_t grp_id = lookup_group_id( grp_name );

  if ( grp_id == 0 ) return 0;
  

  //
  // Merge regions with the same "name" into a new group, where each
  // region is now a subregion. Specific meta-information will be
  // lost.
  //
  
  // We are assuming that this function is being applied only to
  // features read from a GTF2.2 file as specific elements (CDS/exon
  // distinctions) will be stored/handled differently.

  // If mergeField != "" then use that meta-information field. If not
  // present, skip that region
  
  bool useName = mergeField == "";


  // Assign a new group ID, or append to an existing one

  uint64_t new_group_id = set_group_id( name );

  sql.bind_int64( stmt_loc_lookup_group, ":group_id" , grp_id );

  get_subregions( false );
  get_meta( true ); // as any meta-info will be copied to the subregion
  
  // Hold new regions in memory

  sql.begin();

  std::map<NameAndChr,Region> merged;
  int cnt = 0;
  while ( sql.step( stmt_loc_lookup_group ) )
    {
      
      Region r = construct_region( stmt_loc_lookup_group );
      
      std::string name;
      
      if ( useName ) 
	name = r.name;
      else
	{
	  if ( ! r.meta.has_field( mergeField ) )
	    continue;
	  name = r.meta.get1_string( mergeField );	  
	}
      
      // Have we seen a region of this name before?
      NameAndChr nameChr(name, r.chromosome());
      std::map<NameAndChr,Region>::iterator i = merged.find( nameChr );

      if ( i != merged.end() )
	{
	  Region & par = i->second;
	  
	  if ( r.chromosome() != par.chromosome() )
	    continue;
	  
	  if ( r.start.position() < par.start.position() )
	    par.start.position( r.start.position() );
	  
	  if ( r.stop.position() > par.stop.position() )
	    par.stop.position( r.stop.position() );
	  
	  par.addSubRegion( r );
	  
	}
      else
	{
	    Region par( r.chromosome(), 
			r.start.position() , 
			r.stop.position() , 
			name , 
			r.altname, 
			new_group_id );	  	    	    

	    // hmm, this is quite messy, but live with for now
	    // will end up with _STRAND tags in the Region for no
	    // good reason, but the subregion strand/frame members
	    // will be okay
	    
	    par.meta.set( PLINKSeq::TRANSCRIPT_FRAME() , r.meta.get1_int( PLINKSeq::TRANSCRIPT_FRAME() )  );
	    par.meta.set( PLINKSeq::TRANSCRIPT_STRAND() , r.meta.get1_int( PLINKSeq::TRANSCRIPT_STRAND() )  );

	    par.addSubRegion( r );
	    
	    merged.insert( std::make_pair( nameChr , par ) );
	}
      
    }
      
  sql.reset( stmt_loc_lookup_group );

  sql.commit();


  /////////////////////////////////////////
  //                                     //
  // Insert new, merged segments into DB //
  //                                     //
  /////////////////////////////////////////

  sql.begin();

  std::map<NameAndChr,Region>::iterator i = merged.begin();
  
  while ( i != merged.end() )
    {
      range_insertion( i->second );
      ++i;
    }

  sql.commit();
  
  plog << "inserted " << merged.size() << " merged regions\n";
  
  return new_group_id;
 
}




std::set<Region> LocDBase::get_set( uint64_t set_id )
{
  std::set<Region> regions;
  sql.bind_int64( stmt_loc_lookup_set, ":set" , set_id );
  while ( sql.step( stmt_loc_lookup_set ) )
    regions.insert( construct_region( stmt_loc_lookup_set ) );
  sql.reset( stmt_loc_lookup_set );
  return regions;
}


void LocDBase::load_alias( const std::string & filename )
{

  if ( ! attached() ) return;
  
  Helper::checkFileExists( filename );
  
  InFile f( filename );
   
  // We require a header row with names, tab-delimited

  std::vector<std::string> header = f.tokenizeLine("\t");
  if ( header.size() < 2 ) Helper::halt("empty/incomplete alias table");
  if ( header[0][0] != '#' ) Helper::halt("expecting #header in alias table");
  if ( header[0] == "#" ) Helper::halt("expecting #header in alias table (not # header) ");  
  header[0] = header[0].substr(1);
  
  int coln = header.size();  
  std::vector<int> ids(coln);
  
  for (int i=0; i<coln; i++)
    ids[i] = insert_alias_group( header[i] );
  
 
  //
  // Now read table, and insert all pairs all non-missing entries (with double-entry). Nulls 
  // are allowed in the second column only (i.e. to explicitly show there is no equivalent)
  //
  
  sql.begin();
  
  int inserted = 0;
  long int inserted_pairs = 0;
  
  while ( ! f.eof () )
    {

      // require a tab-delimited set
      std::vector<std::string> buffer = f.tokenizeLine("\t");
      
      if( buffer.size() != coln ) continue;
      
      // convert 'n/a' to missing

      for (int i=0; i<coln; i++)	  
	if ( buffer[i] == "." || buffer[i] == "n/a" || buffer[i] == "N/A" ) 
	  buffer[i] = "";
      
      ++inserted;
      
      for (int i=0; i < (coln-1); i++)
	for (int j=i+1; j<coln; j++)	  
	  {
	    int grpi = ids[i];
	    int grpj = ids[j];
	    
	    if ( buffer[i] != "" )
	      {		
		sql.bind_int64( stmt_loc_alias_insert , ":group_id1" , ids[i] );
		sql.bind_text( stmt_loc_alias_insert , ":name1" , buffer[i] );
		
		sql.bind_int64( stmt_loc_alias_insert , ":group_id2" , ids[j] );
		sql.bind_text( stmt_loc_alias_insert , ":name2" , buffer[j] );
		
		sql.step( stmt_loc_alias_insert );
		sql.reset( stmt_loc_alias_insert );
		++inserted_pairs;
	      }
	    
	    // double-entry

	    if ( buffer[j] != "" )
	      {		
		sql.bind_int64( stmt_loc_alias_insert , ":group_id1" , ids[j] );
		sql.bind_text( stmt_loc_alias_insert , ":name1" , buffer[j] );
		
		sql.bind_int64( stmt_loc_alias_insert , ":group_id2" , ids[i] );
		sql.bind_text( stmt_loc_alias_insert , ":name2" , buffer[i] );
		
		sql.step( stmt_loc_alias_insert );
		sql.reset( stmt_loc_alias_insert );
		++inserted_pairs;
	      }
	  }

      // next row of aliases
    }
  
  sql.commit();	
  
  f.close();
  
  plog << "inserted " 
       << inserted_pairs 
       << " alias pairs, from " 
       << inserted 
       << " table rows\n";
  
}

void LocDBase::delete_aliases()
{
  sql.query(" DELETE FROM aliases; " );
  sql.query(" DELETE FROM alias_groups; ");

}


uint64_t LocDBase::alias_id( const std::string & g ) 
{
  return alias_group_table[ g ];
}

std::string LocDBase::alias( const std::string & query , uint64_t query_grp_id , uint64_t alias_grp_id )
{
  return Helper::stringize( targetted_lookup_alias( query , query_grp_id , alias_grp_id ) );
}

std::string LocDBase::alias( const std::string & query , bool show_keys )
{
  return attached() ? Helper::stringizeKeyPairList( lookup_alias( query , 0 ) , show_keys ) : "." ;
}

std::map<std::string,std::string> LocDBase::lookup_alias( const std::string & query , const std::string & alias_group  )
{
  
  if ( alias_group_table.find( alias_group ) == alias_group_table.end() ) 
    {
      std::map<std::string,std::string> s;
      return s;
    }
  else
    return lookup_alias( query , alias_group_table[ alias_group ] );
}

std::map<std::string,std::string> LocDBase::lookup_alias( const std::string & name , const uint64_t alias_group_id )
{
  std::map<std::string,std::string> m;
  if ( ! attached() ) return m;

  sqlite3_stmt * s = alias_group_id == 0 ? stmt_loc_alias_lookup : stmt_loc_group_alias_lookup ;
  if ( alias_group_id )  sql.bind_int( s , ":group_id" , alias_group_id );
  sql.bind_text( s , ":name" , name );  

  while ( sql.step( s ) )
    {
      std::string k = alias_group_reverse_table[ sql.get_int( s , 0 ) ];
      if ( m.find(k) == m.end() ) 
	m[k] = sql.get_text( s , 1 ) ;  
      else
	m[k] += "," + sql.get_text( s , 1 ) ;
    }
  sql.reset( s );
  return m;
}

std::set<std::string> LocDBase::targetted_lookup_alias( const std::string & query , 
							const std::string & query_group , 
							const std::string & alias_group )
{
  return targetted_lookup_alias( query , alias_group_table[ query_group ] , alias_group_table[ alias_group ]);
}

std::set<std::string> LocDBase::targetted_lookup_alias( const std::string & query , 
							const uint64_t gquery , 
							const uint64_t gtarget )
{
  std::set<std::string> s;
  if ( ! attached() ) return s;
  if ( gquery == 0 || gtarget == 0 ) return s;

  sql.bind_text( stmt_loc_targetted_group_alias_lookup , ":name" , query );
  sql.bind_int(  stmt_loc_targetted_group_alias_lookup , ":group_id1" , gquery );
  sql.bind_int(  stmt_loc_targetted_group_alias_lookup , ":group_id2" , gtarget );

  while ( sql.step( stmt_loc_targetted_group_alias_lookup ) )
    {      
      s.insert( sql.get_text( stmt_loc_targetted_group_alias_lookup , 0 ) ) ;  
    }
  sql.reset( stmt_loc_targetted_group_alias_lookup );

  return s;  
}


uint64_t LocDBase::merge_overlap(uint64_t grp_id, const std::string & name, bool storeSubregions)
{
  
  // Get a new ID code, which we will return also
  
  uint64_t group_id = set_group_id( name );
  
  Region current;
  
  int cnt = 0;
  
  sql.begin();

  sql.bind_int64( stmt_loc_iterate_group , ":group_id" , group_id ); 

  // Iterate through position-ordered list of regions

  while ( sql.step( stmt_loc_iterate_group ) )
    {
      
      Region r = construct_region( stmt_loc_iterate_group );
      
      // Does this overlap with current?
      
      if ( current.overlaps( r ) )
	{
	  current.addSubRegion(r);
	}
      else
	{
	  
	  // Create a single region that represents the span of all
	  // sub-regions (optionally keeping the subregions too)
	  
	  current.collapse( storeSubregions );	  
	  range_insertion( current );
	  current = r;
	  
	}
    }
  

  // Insert final region
  
  current.collapse( storeSubregions );
  range_insertion( current );
  
  
  // Clean up & return

  sql.reset( stmt_loc_iterate_group );
  sql.commit();
  return group_id;
  
}



uint64_t LocDBase::extract(uint64_t group1_id , uint64_t group2_id, const std::string & name )
{
    
  // Take the regions from group1 that have a name-matching entry in group2 and 
  // place in a new group, 
  
  uint64_t new_group_id = set_group_id( name );
  
  // Insert all names into a temporary table
  
  sql.begin();
  
  sql.bind_int64( stmt_loc_intersect , ":group1_id" , group1_id ); 
  sql.bind_int64( stmt_loc_intersect , ":group2_id" , group2_id ); 
  
  while ( sql.step( stmt_loc_intersect ) )
    {
      
      Region r = construct_region( stmt_loc_intersect );
      
      // Update group ID of region
      
      r.group = new_group_id;
      
      // Insert back into database
      range_insertion( r );
      
    }
  
  sql.reset( stmt_loc_intersect );
  
  sql.commit();
  
  return new_group_id;
}


uint64_t LocDBase::merge(uint64_t grp1_id, uint64_t grp2_id , const std::string & name )
{
    //
}
   
uint64_t LocDBase::rename( uint64_t grp_id, uint64_t alias_id , const std::string & name )
{

//     // Take the regions from group1 that have a name-matching entry in group2 and 
//     // place in a new group, 

//     uint64_t new_group_id = set_group_id( name );
    
//     // Insert all names into a temporary table
    
//     sql.begin();
    
//     sql.bind_int64( stmt_loc_lookup_group , ":group_id" , grp_id ); 

//     while ( sql.step( stmt_loc_lookup_group ) )
//     {
// 	Region r = construct_region( stmt_loc_lookup_group );
	
// 	// Update region
// 	r.group = new_group_id;
// 	vector<string> alias = get_alias( r.name, alias_id );	
// 	for (int i=0; i<alias.size(); i++)
// 	{	    
// 	    r.name = alias[i];
// 	    range_insertion( r );
// 	}
//     }

//     sql.reset( stmt_loc_lookup_group );

//     sql.commit();

//    return new_group_id;
}

std::set<GroupInfo> LocDBase::group_information()
{

  std::set<GroupInfo> gset;
  
  while ( sql.step( stmt_loc_group_list ) ) 
    {
      GroupInfo g;
      g.idx = sql.get_int64( stmt_loc_group_list, 0 );
      g.name = sql.get_text( stmt_loc_group_list , 1 );
      g.temp = sql.get_int( stmt_loc_group_list , 2 );
      g.description = sql.get_text( stmt_loc_group_list , 3 );
      gset.insert( g ) ;
    }
  sql.reset( stmt_loc_group_list );
  return gset;    
}

std::set<GroupInfo> LocDBase::set_information()
{
  
  // Use 'idx' slot to represent the group_id
  
  std::set<GroupInfo> gset;
  
  while ( sql.step( stmt_dump_all_sets ) ) 
    {
      GroupInfo g;
      g.idx = sql.get_int64( stmt_dump_all_sets, 0 );
      g.name = sql.get_text( stmt_dump_all_sets , 1 );
      g.description = sql.get_text( stmt_dump_all_sets , 2 );
      gset.insert( g ) ;
    }
  sql.reset( stmt_dump_all_sets );
  return gset;      
}

std::string LocDBase::summary( bool ugly )
{
  
  std::stringstream ss;
  
  if ( ! ugly ) ss << "---Locus DB summary---\n\n";

  std::set<GroupInfo> g = group_information();
  std::set<GroupInfo>::iterator i = g.begin();
  

  while ( i != g.end() )
    {
      if ( ugly ) 
	ss << "LOCDB\t"
	   << "NAME=" << i->name << "\t"
	   << "N=" << count( i->idx ) << "\t"
	   << "DESC=" << i->description
	   << "\n";
      else
	ss << "Group : " << i->name << " (" << count( i->idx ) << " entries) " << i->description << "\n";
      
      ++i;
    }
  

  // Sets
  
  std::set<GroupInfo> h = set_information();

  if ( g.size() == 0 && h.size() == 0 ) ss << "(empty)\n";

  i = h.begin();
  while ( i != h.end() )
    {
      std::string grp = lookup_group_id( i->idx );
      
      std::vector<std::string> t = fetch_set_names( grp , i->name );
      
      if ( ugly ) 
	ss << "LOCDB_SET\t"
	   << "GROUP=" << grp << "\t"
	   << "NAME=" << i->name << "\t"
	   << "N=" << t.size() << "\t"
	   << "DESC=" << i->description << "\n";
      else
	ss << "Locus set : " << i->name << " (" << t.size() << " entries) " << i->description << "\n";
      
      ++i;
    }
  
  return ss.str();
}


uint64_t LocDBase::load_set(const std::string & filename, 
			    const std::string & label, 
			    const std::string & set_group , 
			    bool use_altname ) 
{
  
  if ( ! attached() ) return 0;

  if ( ! Helper::fileExists( filename ) ) return 0;
  
  InFile f( filename );

  // Specify which group of loci to look at, e.g. CCDS
  
  uint64_t sgroup_id = lookup_group_id( set_group );

  if ( sgroup_id == 0 ) return 0;
  
  // Exect file that contains two columns of information only
  //  GENE-NAME  PATHWAY-NAME

  // Assign a -set group ID

  uint64_t group_id = set_set_id( label , sgroup_id , 1 , "n/a" );
  
  // Track gene and pathway IDs below

  std::map<std::string,std::vector<int> > gmap;
  std::map<std::string,int> pmap;
  
  //
  // Begin SQL transaction               
  // 
  
  sql.begin();
  
  int inserted = 0;

  // 
  // Process each row of input
  // 

  
  while ( ! f.eof() )
    {
      
      // split of tabs
      std::vector<std::string> tok = Helper::char_split( f.readLine() , '\t' );
      
      // Should contain exactly 2 tab-delimited elements      

      if ( tok.size() == 0 ) continue;

      if ( tok.size() != 2 )
	{	 
	  plog.warn("not 2 tab-delimited columns in geneset file");
	  continue;
	}

      const std::string gene_name = tok[0];
      const std::string pathway_name = tok[1];
      
      // Look-up gene-name, pathway name
      
      std::map<std::string, std::vector<int> >::iterator i = gmap.find( gene_name );

      std::vector<int> g_id;
      int p_id = 0;

      if ( i == gmap.end() ) 
	{
	  
	  if ( use_altname ) 
	    {
		sql.bind_text( stmt_loc_lookup_real_name , ":altname" , gene_name );
		sql.bind_int( stmt_loc_lookup_real_name , ":group_id" , sgroup_id ); 
		while ( sql.step( stmt_loc_lookup_real_name ) ) 
		    g_id.push_back( sql.get_int( stmt_loc_lookup_real_name , 0 ) );	      

		sql.reset( stmt_loc_lookup_real_name );
	    }
	  else
	    {
	      sql.bind_text( stmt_loc_lookup_group_and_name , ":name" , gene_name );
	      sql.bind_int( stmt_loc_lookup_group_and_name , ":group_id" , sgroup_id );
	      while ( sql.step( stmt_loc_lookup_group_and_name ) )	    
		g_id.push_back( sql.get_int( stmt_loc_lookup_group_and_name , 0 ) );	      
	      sql.reset( stmt_loc_lookup_group_and_name );
	    }
	  gmap[ gene_name ] = g_id;
	}
      else
	g_id = i->second;
      
      // Requires that gene actually exists
      if ( g_id.size() == 0 ) 
	{
	  plog.warn("could not find gene specified in locus-set");
	  continue;
	}

      std::map<std::string,int>::iterator j = pmap.find( pathway_name );
      if ( j == pmap.end() )
	{
	  sql.bind_text( stmt_set_member_insert , ":name" , pathway_name );
	  sql.bind_int( stmt_set_member_insert , ":group_id" , group_id );
	  sql.step( stmt_set_member_insert ); 
	  p_id = sql.last_insert_rowid();
	  sql.reset( stmt_set_member_insert );
	  pmap[ pathway_name ] = p_id;
	}
      else
	p_id = j->second;
      
      
      //
      // Add to set database
      //
      

      sql.bind_int( stmt_set_data_insert , ":set_id" , p_id );
      for ( int g= 0; g < g_id.size(); g++ )
	{
	  sql.bind_int( stmt_set_data_insert , ":loc_id" , g_id[g] );
	  sql.step( stmt_set_data_insert );
	  sql.reset( stmt_set_data_insert );
	}
      
      ++inserted;

    }

  //
  // Finish transaction
  //    

  sql.commit();
  
  f.close();

  plog << "inserted " << inserted << " rows\n";
  
  return group_id;

}

std::vector<Region> LocDBase::fetch( const std::string & grp, const std::vector<std::string> & names )
{

  std::vector<Region> regions;
  
  uint64_t group_id = lookup_group_id(grp);
  if ( group_id == 0 ) return regions;
  
  sql.bind_int64( stmt_loc_lookup_group_and_name, ":group_id" , group_id );
  
  for (int i = 0 ; i < names.size(); i++)
    {
      sql.bind_text( stmt_loc_lookup_group_and_name, ":name" , names[i] );
      
      while ( sql.step( stmt_loc_lookup_group_and_name ) ) 
	{
	  regions.push_back( construct_region( stmt_loc_lookup_group_and_name ) );    
	}
      sql.reset( stmt_loc_lookup_group_and_name );
    }

  return regions;
}


void LocDBase::replace_real_names( const int grp , const std::string & name , const std::string & newname , bool search_alternate )
{
  if ( search_alternate )
    {
//	std::cout <<" search " << grp << " " << name << " " << newname << "\n";
	
	sql.bind_int64( stmt_loc_replace_real_name_alternate , ":group_id" , grp );
	sql.bind_text( stmt_loc_replace_real_name_alternate , ":oldname" , name );
	sql.bind_text( stmt_loc_replace_real_name_alternate , ":newname" , newname );

	while ( sql.step( stmt_loc_replace_real_name_alternate ) ) 
	{ 
	    // std::cout << " done\t" << sql.get_text( stmt_loc_replace_real_name_alternate , 0 )  << "\n" ;
	} 
	sql.reset( stmt_loc_replace_real_name_alternate );  
    }
  else
    {
      sql.bind_int64( stmt_loc_replace_real_name , ":group_id" , grp );
      sql.bind_text( stmt_loc_replace_real_name , ":name" , name );
      sql.bind_text( stmt_loc_replace_real_name , ":altname" , newname );
      while ( sql.step( stmt_loc_replace_real_name ) ) { } 
      sql.reset( stmt_loc_replace_real_name );  
    }

}

std::vector<Region> LocDBase::fetch_real_names( const std::string & grp, const std::string & altname )
{
  std::vector<Region> regions;
  uint64_t group_id = lookup_group_id(grp);
  if ( group_id == 0 ) return regions;
  
  sql.bind_int64( stmt_loc_lookup_real_name, ":group_id" , group_id );
  sql.bind_text( stmt_loc_lookup_real_name, ":altname" , altname );

  while ( sql.step( stmt_loc_lookup_real_name ) ) 
    {
      regions.push_back( construct_region( stmt_loc_lookup_real_name ) );    
    }
  sql.reset( stmt_loc_lookup_real_name );
  return regions;
}

std::set<std::string> LocDBase::fetch_names( const std::string & loc_group )
{
  
  std::set<std::string> results;
  if ( ! attached() ) return results;

  uint64_t id = lookup_group_id( loc_group );
  if ( id == 0 ) return results;
  
  sql.bind_int64( stmt_loc_name_dump , ":group_id" ,id );
  while ( sql.step( stmt_loc_name_dump ) ) 
    results.insert( sql.get_text( stmt_loc_name_dump , 0 ) );
  sql.reset( stmt_loc_name_dump );
  return results;
}

std::vector<std::string> LocDBase::fetch_set_names( const std::string & loc_group, 
						    const std::string & set_group )
{
  std::vector<std::string> results;
  if ( ! attached() ) return results;
  uint64_t id = lookup_set_id( loc_group , set_group );
  if ( id == 0 ) return results;
  sql.bind_int64( stmt_set_names_fetch, ":group_id" , id );
  while ( sql.step( stmt_set_names_fetch ) )
    {
      results.push_back( sql.get_text( stmt_set_names_fetch , 0 ) );
    }
  sql.reset( stmt_set_names_fetch );
  return results;
}


std::vector<std::string> LocDBase::fetch_set_members( const std::string & loc_group, 
						      const std::string & set_group,
						      const std::string & set_name )
{
  std::vector<std::string> results;
  if ( ! attached() ) return results;
  
  uint64_t id = lookup_set_id( loc_group , set_group );
  if ( id == 0 ) return results;

  sql.bind_int64( stmt_set_member_lookup, ":group_id" , id );
  sql.bind_text( stmt_set_member_lookup, ":name" , set_name );
  
  uint64_t mem_id = 0;
  
  if ( sql.step( stmt_set_member_lookup ) )
    {
      mem_id = sql.get_int64( stmt_set_member_lookup , 0 );
    }
  sql.reset( stmt_set_member_lookup );
  if ( mem_id == 0 ) return results;
  
  sql.bind_int64( stmt_set_members_fetch, ":set_id" , mem_id );
  while ( sql.step( stmt_set_members_fetch ) )
    {
      results.push_back( sql.get_text( stmt_set_members_fetch , 0 ) );
    }
  sql.reset( stmt_set_members_fetch );
  
  return results;
}


std::map<std::string,std::set<Region> > LocDBase::fetch_set_regions( const std::string & loc_group , 
								     const std::string & set_group )
{

  std::map<std::string,std::set<Region> > r;
  if ( ! attached() ) return r;

  uint64_t id = lookup_set_id( loc_group , set_group );
  if ( id == 0 ) return r;

  std::vector<std::string> sets = fetch_set_names( loc_group , set_group );

  for (int s = 0 ; s < sets.size(); s++ ) 
    {
      std::set<Region> t;

      sql.bind_int64( stmt_set_member_lookup, ":group_id" , id );
      sql.bind_text( stmt_set_member_lookup, ":name" , sets[s] );
      
      uint64_t mem_id = 0;      
      if ( sql.step( stmt_set_member_lookup ) )
	mem_id = sql.get_int64( stmt_set_member_lookup , 0 );	
      sql.reset( stmt_set_member_lookup );
      int mm = mem_id;

      if ( mem_id != 0 ) 
	{

	  sql.bind_int64( stmt_set_members_fetch_regions, ":set_id" , mem_id );
	  
	  while ( sql.step( stmt_set_members_fetch_regions ) )
	    {	      
	      t.insert( construct_region( stmt_set_members_fetch_regions ) ) ;
	    }
	  sql.reset( stmt_set_members_fetch_regions );	  
	}

      // add to the list
      r[ sets[s] ] = t;
      
    }
  
  return r;

}


void LocDBase::append_metainformation( Variant & v , const std::set<int> & grp )
{
  // For this variant, based on chromosomal position, lookup all loci
  // that overlap and append as meta-information

  sql.bind_int( stmt_loc_lookup_group_and_range , ":chr" , v.chromosome() );
  sql.bind_int( stmt_loc_lookup_group_and_range , ":start" , v.position() );
  sql.bind_int( stmt_loc_lookup_group_and_range , ":end" , v.position() );

  std::set<int>::iterator i = grp.begin();
  while ( i != grp.end() )
    {
      sql.bind_int( stmt_loc_lookup_group_and_range , ":group_id" , *i );      
      while ( sql.step( stmt_loc_lookup_group_and_range ) )
	{
	  std::string name = sql.get_text( stmt_loc_lookup_group_and_range , 1 );
	  if ( v.meta.add_if_unique( PLINKSeq::META_LSET() , name ) )
	    v.meta.add( PLINKSeq::META_LGRP() , *i );	          

	}
      sql.reset( stmt_loc_lookup_group_and_range );
      ++i;
    }
    
}


std::vector<std::string> LocDBase::fetch_name_given_altname( const std::string & loc_group , 
							     const std::string & altname )
{
  std::vector<std::string> s;
  if ( ! attached() ) return s;
  
  uint64_t id = lookup_group_id( loc_group );
  if ( id == 0 ) return s;

  sql.bind_int64( stmt_loc_lookup_real_name_only, ":group_id" , id );  
  sql.bind_text( stmt_loc_lookup_real_name_only , ":altname" , altname );
  while ( sql.step( stmt_loc_lookup_real_name_only ) )
    s.push_back( sql.get_text( stmt_loc_lookup_real_name_only , 0 ) );
  sql.reset( stmt_loc_lookup_real_name_only );
  return s;
}


void LocDBase::insert_special( const std::string & key , 
			       const std::vector<std::string> & values )
{
  if ( ! attached() ) return;
  for ( int i=0; i<values.size(); i++)
    {
      sql.bind_text( stmt_insert_special , ":name" , key );
      sql.bind_text( stmt_insert_special , ":value" , values[i] );
      sql.step( stmt_insert_special );
      sql.reset( stmt_insert_special );
    }
}
 

std::vector<std::string> LocDBase::fetch_special( const std::string & key ) 
{
  std::vector<std::string> r;
  if ( ! attached() ) return r;
  sql.bind_text( stmt_fetch_special , ":name" , key );
  while ( sql.step( stmt_fetch_special ) )
    {
      std::string s = sql.get_text( stmt_fetch_special , 0 );
      r.push_back(s);
    }
  return r;
}

void LocDBase::clear_special()
{
  if ( ! attached() ) return;
  sql.query(" DELETE FROM special; "); 
}



bool LocDBase::populate_set_structures( const std::string & loc_group , 
					const std::string & set_group , 
					std::map<int,std::string> * genes ,
					std::map<int,std::string> * sets , 
					std::map<int,std::set<int> > * s2g,
					std::map<int,std::set<int> > * g2s )
{

  if ( ! attached() ) return false;
  
  // get set ID
  uint64_t id = lookup_set_id( loc_group , set_group );
  if ( id == 0 ) return false;

  // get locus ID
  uint64_t loc_id = lookup_group_id( loc_group );
  if ( loc_id == 0 ) return false;

  // Get all members of set
  sql.bind_int64( stmt_set_names_and_id_fetch , ":group_id" , id );
  while ( sql.step( stmt_set_names_and_id_fetch ) ) 
    {
      uint64_t mem_id = sql.get_int64( stmt_set_names_and_id_fetch , 1 );
      std::string name = sql.get_text( stmt_set_names_and_id_fetch , 0 );
      (*sets)[ mem_id ] = name; 
    }
  sql.reset( stmt_set_names_and_id_fetch );


  // Get locus-IDs
  sql.bind_int64( stmt_loc_lookup_group , ":group_id" , loc_id );		 
  while ( sql.step(stmt_loc_lookup_group) ) 
    {      
      uint64_t lid = sql.get_int64( stmt_loc_lookup_group , 0 );
      std::string name = sql.get_text( stmt_loc_lookup_group , 1 );      
      (*genes)[ lid ] = name;
    }
  sql.reset( stmt_loc_lookup_group );

  // Get all locus-set pairings (dumps all, manually extract relevant groups)
  while ( sql.step( stmt_set_data_dumper ) )
    {
      uint64_t set_id = sql.get_int64( stmt_set_data_dumper , 1 ) ;
      if ( sets->find( set_id ) == sets->end() ) continue;
      uint64_t loc_id = sql.get_int64( stmt_set_data_dumper , 0 ) ;
      (*s2g)[ set_id ].insert( loc_id );
      (*g2s)[ loc_id ].insert( set_id );
    }
  sql.reset( stmt_set_data_dumper );
  
  return true;
}



void LocDBase::check_version()
{

  if ( ! sql.table_exists( "dbmeta" ) )
    Helper::halt( "old database format, expecting LOCDB v"
                  + Helper::int2str( PLINKSeq::LOCDB_VERSION_NUMBER() )
                  + " : to fix, remake your project" );
  
  // expected version # is given by  PLINKSeq::LOCDB_VERSION_NUMBER()                                                                      
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
      v = PLINKSeq::LOCDB_VERSION_NUMBER();
      sql.bind_text( si , ":x" , vn );
      sql.bind_int( si , ":y" , v );
      sql.step(si);
      sql.finalise(si);
    }

  if ( v != PLINKSeq::LOCDB_VERSION_NUMBER() )
    Helper::halt("LOCDB version "
                 + Helper::int2str( v ) + " but expected "
                 + Helper::int2str( PLINKSeq::LOCDB_VERSION_NUMBER() )
                 + " : to fix, remake your LOCDB" );

  return;

}



std::string LocDBase::get_genename( const Variant & var , uint64_t group_id , const std::string &  delim )
{

  std::string s = ".";
  
  if ( var.length() == 1 ) 
    {
      sql.bind_int64( stmt_loc_fetch_altnames , ":group_id" , group_id );
      sql.bind_int( stmt_loc_fetch_altnames , ":chr" , var.chromosome() );
      sql.bind_int( stmt_loc_fetch_altnames , ":bp" , var.position() );
      while ( sql.step( stmt_loc_fetch_altnames ) )
	{
	  if ( s == "." ) s = sql.get_text( stmt_loc_fetch_altnames , 0 );
	  else s += delim + sql.get_text( stmt_loc_fetch_altnames , 0 );
	}
      sql.reset( stmt_loc_fetch_altnames );
    }
  else
    {
      sql.bind_int64( stmt_loc_fetch_altnames_indel , ":group_id" , group_id );
      sql.bind_int( stmt_loc_fetch_altnames_indel , ":chr" , var.chromosome() );
      sql.bind_int( stmt_loc_fetch_altnames_indel , ":start" , var.position() );
      sql.bind_int( stmt_loc_fetch_altnames_indel , ":stop" , var.stop() );
      while ( sql.step( stmt_loc_fetch_altnames_indel ) )
	{
	  if ( s == "." ) s = sql.get_text( stmt_loc_fetch_altnames_indel , 0 );
	  else s += delim + sql.get_text( stmt_loc_fetch_altnames_indel , 0 );
	}
      sql.reset( stmt_loc_fetch_altnames_indel );
    }
  return s;
}

uint64_t LocDBase::get_region_id( uint64_t gid , const std::string & name )
{
  sql.bind_int( stmt_loc_fetch_id_given_name , ":group_id" , gid );
  sql.bind_text( stmt_loc_fetch_id_given_name , ":name" , name );
  uint64_t t = 0;
  if ( sql.step( stmt_loc_fetch_id_given_name ) )
    t = sql.get_int( stmt_loc_fetch_id_given_name , 0 );
  sql.reset( stmt_loc_fetch_id_given_name );
  return t;
}

bool LocDBase::check_GTF(const std::string & filename , bool use_geneid , bool add_frame )
{
  
  // as well as checking the basic syntax/format of the GTF, this function
  // will check the integrity of the implied transcript models (see below)


  // could add option to add frame, if that is not present

  //
  // Expect GTF2.2 format, as documented 
  //  http://mblab.wustl.edu/GTF22.html

  if ( ! Helper::fileExists( filename ) ) return false;
  
  InFile f( filename );
    
  // keep track of transcripts we've seen

  int inserted = 0;  
  std::map<NameAndChr,std::vector<Region> > regions;
  
  /////////////////////////////////////////
  //                                     //
  // Process each row of input           //
  //                                     //
  /////////////////////////////////////////

  int errs = 0;
  int warns = 0;

  while ( ! f.eof() )
    {

      // tab-delimited line

      std::string str = f.readLine();

      std::vector<std::string> tok = Helper::char_split( str , '\t' );

      // skip comments and empty lines; not errors
      if ( tok.size() == 0 ) continue;
      if ( tok[0][0] == '#' ) continue;

      // Should contain exactly 9 tab-delimited elements      
      if ( tok.size() != 9 ) 
	{
	  plog << "!!!\tfile_format\t"
	       << "other than 9 tab-delimited fields: " << str << "\n\n";
	  ++errs;
	  continue;
	}
      
      
      // <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [attributes] [comments]
      // AB000381 Twinscan  CDS          380   401   .   +   0  gene_id "001"; transcript_id "001.1";

            
      // Start/stop positions

      int p1,p2;
      if ( ! Helper::str2int( tok[ 3 ] , p1 ) ) 
	{
	  plog << "!!!\tfile_format\t"
	       << "non-numeric value for bp1 " << tok[3] << " for row: " << str << "\n";
	  ++errs;
	  continue;
	}

      if ( ! Helper::str2int( tok[ 4 ] , p2 ) ) 
	{
	  plog << "!!!\tfile_format\t"
	       << "non-numeric value for bp2 " << tok[4] << " for row: " << str << "\n";
	  ++errs;
	  continue;
	}      
 
     // Name (from transcript_id)

      std::vector<std::string> tok2 = Helper::char_split( tok[8] , ' ' , false );
      
      // requires atleast 4 manadatory fields: gene_id XXX transcript_id XXX
      
      if ( tok2.size() < 4 ) 
	{
	  plog << "!!!\tfile_format\t"
	       << "expecting at least [ gene_id \"ID\"; transcript_id \"ID\"; ] : " << str << "\n";
	  ++errs;
	  continue;
	}
      
      // remove quotes and last semi-colon

      if ( tok2[0] != "gene_id" ) 
	{
	  plog << "!!!\tfile_format\t"
	       << "expecting gene_id in GTF :" << str << "\n";
	  ++errs;
	  continue;
	}
 
      if ( tok2[2] != "transcript_id" ) 
	{
	  plog << "!!!\tfile_format\t"
	       << "expecting gene_id in GTF :" << str << "\n";
	  ++errs;
	  continue;
	}
      
      if ( tok2[1].substr( tok2[1].size()-1 ) == ";" ) tok2[1] = tok2[1].substr( 0, tok2[1].size()-1 );
      else 
	{	  
	  plog << "!!!\tfile_format\t"
	       << "expecting ';' character after gene_id : " << str << "\n";
	  ++errs;
	  continue;
	}

      if ( tok2[3].substr( tok2[3].size()-1 ) == ";" ) tok2[3] = tok2[3].substr( 0, tok2[3].size()-1 );
      else 
	{
	  plog << "!!!\tfile_format\t"
	       << "expecting ';' character after gene_id : " << str << "\n";
	  ++errs;
	  continue;
	}
      
      
      tok2[1] = Helper::unquote( tok2[1] );
      tok2[3] = Helper::unquote( tok2[3] );

      // tok2[0] should equal  gene_id
      // use either gene-name or transcript name as 'name' 
      
      std::string name = use_geneid ? tok2[1] : tok2[3];
      
      int chromosome = Helper::chrCode( tok[0] ) ;
      
      if ( chromosome == 0 ) 
	{
	  plog << name << "\t"
	       << tok[0] << "\t"
	       << "unknown_chrcode" << "\n";
	       
	  ++errs;
	  continue;
	}
      

      // make a region:

      Region r( chromosome ,
		p1 , p2 , 
		name , 
		0 );


      //
      // Track gene-name, if unique name is a transcript
      //
      
      if ( ! use_geneid ) r.altname = tok2[1];
      else r.altname = tok2[3]; 
      
      std::string & feature = tok[2];
      
      if ( feature != "exon" && 
	   feature != "CDS" && 
	   feature != "start_codon" && 
	   feature != "stop_codon" ) 
	{
	  plog << "??? did not recognize feature ( not 'exon', 'CDS', 'start_codon' or 'stop_codon' ) : " << str << "\n";
	  ++warns;
	  continue;
	}
      
      r.meta.set( "feature" , feature );
      

      // Track strand for CDS (-1,1,  0=missing)
      
      // Use 'strand' to encode two distinct things: a) whether
      //  transcript is on +ve or -ve strand, b) but also some
      //  meta-data about the type of feature. Not particularly
      //  elegant, but will work for now...

      //         <other>     0
      //         CDS         -1,+1
      //         exon        -2,+2
      //         start_codon -3,+3
      //         stop_codon  -4,+4

      int strand = 0;
      if      ( tok[6] == "-" ) strand = -1;
      else if ( tok[6] == "+" ) strand =  1;
      else
	{
	  plog << name << "\t"
	       << tok[0] << "\t"
	       << "invalid_strand" << "\t"
	       << tok[6] << "\n";
	  ++errs;
	  continue;
	}
      
      if      ( feature == "exon" )        strand *= 2;
      else if ( feature == "start_codon" ) strand *= 3;
      else if ( feature == "stop_codon" )  strand *= 4;
      
      r.meta.set( PLINKSeq::TRANSCRIPT_STRAND() , strand );      
      
      //
      // Frame: for CDS 0, 1 or 2 
      //
      
      int frame = 0;
      
      if ( Helper::str2int( tok[7] , frame ) ) 
	{	

	  if ( feature == "CDS" )
	    if ( frame != 0 && frame != 1 && frame != 2 ) 
	      {
		plog << "expecting frame to be 0, 1 or 2, not [" << tok[7] << "] : " << str << "\n"; 
		++errs;
		continue;
	      }
	  
	  r.meta.set( PLINKSeq::TRANSCRIPT_FRAME() , frame );
	}
      else
	{

	  if ( tok[7] != "." )
	    {
	      plog << "expecting frame to be 0, 1 or 2 or '.', not [" << tok[7] << "] : " << str << "\n"; 
	      ++errs;
	      continue;
	    }
	}

    
      ///////////////////////////
      // Add region to database
      
      //
      // store
      //

      regions[ NameAndChr(name,r.chromosome()) ].push_back( r );
      
    }

  //
  // Now validate 
  //
  
  // Specifically, will complain if any of the following (up to user
  // to then fix / remove the offending transcripts, etc)
  
  // a) that has one whole start_codon, one stop_codon (even if split across exons)
  // b) that stop_codon is after last 
  // c) that all transcripts are >= 1 base in length
  // d) that total transcript length is mod 3
  // e) that each new CDS starts at the correct frame
  // f) that start is < stop (of otherwise if -ve strand)
  // g) if has non-standard chr codes, give warning
  // h) that CDS do not overlap
  // i) that has at least one CDS (currently)
  // j) transcripts that map to different to difference chromosomes
  // k) ununusally large or small transcripts (including large introns)
  // l) that the start_cond / CDS is a same position, has frame 0, and is Met AA
  // m) that there are no nonsense codons within the transcript
  // n) that stop_codon is one base after final CDS exon

  
  bool has_seqdb = GP->seqdb.attached();
  
  if ( ! has_seqdb ) 
    plog.warn( "no SEQDB attached, unable to perform certain checks" );
  
  // Q. what about 0-based genomic encoding? 
  std::map<NameAndChr, std::vector<Region> >::iterator ii = regions.begin();
  
  // still flag when transcripts map to multiple chromosomes
  std::map<std::string,std::set<std::string> > track_chrs_per_transcript;
  
  while ( ii != regions.end() )
    {
      const std::string transcriptName = ii->first.name;
      const uint64_t transcriptChr = ii->first.chr;
      
      // track transcript->chromosome mapping
      track_chrs_per_transcript[ transcriptName ].insert( Helper::chrCode( transcriptChr ) );
      
      // get sequence
      std::string sequence = "";
      
      // implied sequence length
      int cds_length = 0;
      std::vector<Region> & ints = ii->second;
      
      
      //
      // start / stop codon
      //
      
      int start_seg = 0 , start_pos = 0;
      int stop_seg = 0 , stop_pos = 0;
      

      bool negative_strand = false;
      bool positive_strand = false;

      for (int i=0;i < ints.size(); i++)
	{
	  int s = ints[i].meta.get1_int( PLINKSeq::TRANSCRIPT_STRAND() );
	  if ( s == -1 ) negative_strand = true;
	  if ( s == +1 ) positive_strand = true;
	}

      std::set<Region> ordered;
      
      for (int i=0;i<ints.size();i++)
	{
	  
	  int s = ints[i].meta.get1_int( PLINKSeq::TRANSCRIPT_STRAND() );
	  

	  bool is_cds = false , is_start = false , is_stop = false;
	  
	  if ( s == -1 || s == +1 ) is_cds = true;
	  if ( s == -3 || s == +3 ) is_start = true;
	  if ( s == -4 || s == +4 ) is_stop = true;
	  
	  if ( is_cds )
	    {	      
	      // track ordering 
	      
	      if ( ints[i].stop.position() < ints[i].start.position() ) 
		{
		  plog << transcriptName << "\t" 
		       << Helper::chrCode(transcriptChr) << "\t"
		       << "backwards_interval" << "\t"
		       << ints[i].coordinate() << "\n";
		}
	      
	      ordered.insert( ints[i] );	     
	      cds_length += ints[i].length();

	    }
	  else if ( is_start )
	    {	      
	      if ( start_seg == 0 ) start_pos = negative_strand ? ints[i].stop.position() : ints[i].start.position();
	      start_seg += ints[i].length();
	    }
	  else if ( is_stop )
	    {
	      if ( stop_seg == 0 ) stop_pos = negative_strand ? ints[i].stop.position() : ints[i].start.position();
	      stop_seg += ints[i].length();	      
	    }
	  // otherwise, no need to check other feature types (for now only 'exon' is allowed anyway)
	}
      

      // allow for non-coding transcripts: if no CDS, start, stop, but only 'exon'
      if ( ordered.size() > 1 )
	{
	  if ( start_seg != 3 ) plog << transcriptName << "\t" << Helper::chrCode(transcriptChr) << "\t" << "bad_start_codon" << "\t" << "size " << start_seg << "\n";
	  if ( stop_seg  != 3 ) plog << transcriptName << "\t" << Helper::chrCode(transcriptChr) << "\t" << "bad_stop_codon" << "\t" << "size " << stop_seg << "\n";
	}

      //
      // is this transcript a protein-coding gene w/ CDS?
      //

      bool has_cds = ordered.size() > 0 ;
      
      
      if ( has_cds )
	{
	  if ( negative_strand == positive_strand ) 
	    {
	      plog << transcriptName << "\t" << Helper::chrCode(transcriptChr) << "\t" << "bad_strand\n";
	    }
	}
	  
      
      

      //
      // check for overlap of any CDS 
      //
      
      if ( ordered.size() > 1 ) 
	{
	  std::set<Region>::iterator jj = ordered.begin();
	  std::set<Region>::iterator kk = ordered.begin();
	  ++kk;
	  while ( kk != ordered.end() ) 
	    {
	      if ( jj->start.chromosome() == kk->start.chromosome() &&  jj->stop.position() >= kk->start.position()  ) 
		plog << transcriptName << "\t" << Helper::chrCode(transcriptChr) << "\toverlapping_CDS_exons\n";
	      else if ( jj->start.chromosome() != kk->start.chromosome() )
		plog << transcriptName << "\t" << Helper::chrCode(transcriptChr) << "\tmapped_to_multiple_chromosomes\t" << Helper::chrCode( jj->start.chromosome() ) << " and " << Helper::chrCode( kk->start.chromosome() ) << "\n";
	      else if ( kk->start.position() - jj->stop.position() + 1 > 1000000 ) 
		plog << transcriptName << "\t" << Helper::chrCode(transcriptChr) << "\tover_1mb_intronic_gap\t" << kk->start.position() - jj->stop.position() + 1 << "\n";
	      
	      ++jj;
	      ++kk;
	    }
	}
      
      //
      // Check correct frame specification of CDS
      //

      int cds_pos = 1;

      std::set<Region>::iterator jj = ordered.begin();
      while ( jj != ordered.end() ) 
	{
	  
	  int f = jj->meta.get1_int( PLINKSeq::TRANSCRIPT_FRAME() );
	  
	  // GTF has frame defined as # of bases to the next codon start 
	  // so         021 021
	  // instead of 012 012 
	  
	  if ( f == 1 ) f = 2; 
	  else if ( f == 2 ) f = 1;

	  int exp_f = ( cds_pos - 1 ) % 3 ;	      

	  if ( f != exp_f && positive_strand ) plog << transcriptName << "\t" << Helper::chrCode(transcriptChr) << "\t"
						    << "incorrect_frame" << "\t"
						    << jj->coordinate() << "\t"
						    << f << " vs " << exp_f << "\n";
	  
	  if ( has_seqdb ) sequence += GP->seqdb.lookup( *jj );
	  
	  cds_pos += jj->length();
	  
	  ++jj;
	  
	}

      // Now check for any stop codons

      if ( has_seqdb )
	{
	  
	  if ( has_cds )
	    {
	      std::vector<std::string> ref_codon;

	      if ( negative_strand ) sequence = Annotate::getrc( sequence );

	      std::string trans_ref = Annotate::translate( sequence , 0 , ref_codon );
	      

	      if ( ref_codon.size() < 20 ) plog << transcriptName << "\t" << Helper::chrCode(transcriptChr) << "\t"
						<< "small_cds" << "\t"
						<< ref_codon.size() << " amino acids\n";

	      if ( ref_codon.size() > 0 ) 
		{
		  if ( trans_ref.substr(0,1) != "M" )
		    plog << transcriptName << "\t" << Helper::chrCode(transcriptChr) << "\t"
			 << "non_Met_start" << "\t"
			 << trans_ref.substr(0,1) << "\t"
			 << Annotate::aa[ trans_ref.substr(0,1) ] << "  "<< trans_ref << "\n";
		  
		}
	      
	      for (int a=0;a<trans_ref.size(); a++)
		{
		  if ( trans_ref[a] == '*' ) plog << transcriptName << "\t" << Helper::chrCode(transcriptChr) << "\t"
						  << "stop_in_reference_CDS" << "\t"
						  << "amino acid position " << (a+1) << "\n";
		}
	    }

	}

      if ( cds_length % 3 ) plog << transcriptName << "\t" << Helper::chrCode(transcriptChr) << "\t"
				 << "CDS_not_mod3" << "\t"
				 << cds_length << " bases\n";
      
      // go to next transcript
      
      ++ii;
    }
  
  std::map<std::string,std::set<std::string> >::iterator jj = track_chrs_per_transcript.begin();
  while ( jj != track_chrs_per_transcript.end() )
    {
      if ( jj->second.size() > 1 ) 
	{
	  std::set<std::string>::iterator kk = jj->second.begin();
	  while ( kk != jj->second.end() )
	    {
	      plog << jj->first << "\t"
		   << *kk << "\t"
		   << "mapped_to_multiple_chromosomes\n";
	      ++kk;
	    }
	}
      ++jj;
    }

}

