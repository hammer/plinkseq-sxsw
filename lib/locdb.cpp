#include "locdb.h"
#include "defs.h"
#include "regions.h"

#include <iostream>
#include <set>


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
  
  if ( n == "-" ) { dettach(); return false; } 

  if ( attached() ) dettach();
  
  // 
  // If the database already exists, just open it
  //
  
  if ( Helper::fileExists(n) )
    {
      sql.open(n);
      fname = n;
      init();
      set_metatypes();
      read_alias_groups();
      return true;
    }
  
  
  //
  // Otherwise, create it
  //
  
  bool status = true;
  
  sql.open(n); 
  
  sql.synchronous(false);
  
  fname = n;
  
  // Main locus table
  
  sql.query(" CREATE TABLE IF NOT EXISTS loci("
	    "   loc_id   INTEGER PRIMARY KEY , "
	    "   name     VARCHAR(20) , "	      
	    "   group_id INTEGER NOT NULL , "
	    "   chr      INTEGER  , "
	    "   bp1      INTEGER  , "
	    "   bp2      INTEGER  , "
	    "   altname  VARCHAR(20)  ); " );
  
  
  // Sub-region table
  
  sql.query(" CREATE TABLE IF NOT EXISTS subloci("
	    "   sub_id   INTEGER PRIMARY KEY , "
	    "   loc_id   INTEGER NOT NULL , "
	    "   name     VARCHAR(20) , "
	    "   chr      INTEGER  , "
	    "   bp1      INTEGER  , "
	    "   bp2      INTEGER  ); " );    
  
  // Name alias tables
    
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
            "   indiv_id INTEGER NOT NULL , "
	    "   name     VARCHAR(20) NOT NULL , "
            " CONSTRAINT uniqID UNIQUE ( name ) ); " );


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
	     "   field_id  INTEGER , "
	     "   value     NUMERIC ); " );

  sql.query( " CREATE TABLE IF NOT EXISTS subloc_meta("
	     "   sub_id    INTEGER NOT NULL , "
	     "   field_id  INTEGER , "
	     "   value     NUMERIC ); " );


  // Group information table
  
  sql.query( " CREATE TABLE IF NOT EXISTS groups("
	     "   group_id     INTEGER PRIMARY KEY , "
	     "   name         VARCHAR(20) NOT NULL , "
             "   temp         CHAR(1) , "
	     "   description  TEXT ); " );


  // Region overlap table (i.e. all pairs of loci that overlap)
  
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
  
  if ( ! status ) 
    plog.warn( "Problem attaching LOCDB " + n );

  return status;

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
      sql.prepare("UPDATE OR IGNORE loci SET altname = :altname WHERE group_id == :group_id AND altname == :altname ; " );

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
    
    stmt_loc_name_list = 
      sql.prepare(" SELECT name FROM loci WHERE group_id == :group_id ORDER BY name;");

    stmt_loc_altname_list = 
      sql.prepare(" SELECT altname FROM loci WHERE group_id == :group_id ORDER BY altname;");

    stmt_loc_lookup_group_and_name = 
      sql.prepare(" SELECT * FROM loci WHERE group_id == :group_id AND name == :name ; " );
    
    stmt_loc_lookup_id = 
      sql.prepare(" SELECT * FROM loci WHERE loc_id == :loc_id ; " );    
        
    stmt_loc_lookup_range = 
      sql.prepare(" SELECT * FROM loci WHERE chr == :chr AND bp1 <= :end AND bp2 >= :start ; " );
    
    stmt_loc_lookup_group_and_range = 
      sql.prepare(" SELECT * FROM loci WHERE group_id == :group_id "
		  " AND chr == :chr AND bp1 <= :end AND bp2 >= :start ; " );
    
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
      sql.prepare("INSERT OR REPLACE INTO subloci ( loc_id, name, chr, bp1, bp2 ) "
		  " values ( :loc_id, :name, :chr, :bp1, :bp2 ) ; " );
    
    
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
    
    stmt_loc_meta_insert_prep = 
      sql.prepare(" SELECT field_id,type FROM metatypes WHERE name == :name ;" );
    
    stmt_loc_meta_insert_prep2 = 
      sql.prepare(" INSERT INTO metatypes (name,type,number,description) values( :name, :type, :number, :description ); " );
  


  stmt_loc_meta_insert = 
    sql.prepare(" INSERT OR REPLACE INTO loc_meta (loc_id,field_id,value)"
		"  values( :reg_id, :field_id, :value ); ");

  stmt_loc_submeta_insert = 
    sql.prepare(" INSERT OR REPLACE INTO subloc_meta (sub_id,field_id,value)"
		"  values( :reg_id, :field_id, :value ); ");


  stmt_loc_get_meta = 
    sql.prepare(" SELECT field_id, value FROM loc_meta "
		" WHERE loc_id == :loc_id " );
		
  stmt_loc_get_submeta = 
    sql.prepare(" SELECT field_id,value FROM subloc_meta "
		" WHERE sub_id == :sub_id ; " );

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
  
  stmt_set_members_fetch = 
    sql.prepare( "SELECT l.name FROM loci AS l , set_data AS sd "
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
    sql.prepare( "SELECT * FROM loci WHERE loc_id IN "
		 " ( SELECT loc_id FROM segments "
		 "    WHERE group_id == :group_id "
		 "      AND indiv_id == :indid_id ) ; " );
  
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
  sql.finalise(stmt_loc_meta_insert_prep );
  sql.finalise(stmt_loc_meta_insert_prep2 );
  sql.finalise(stmt_loc_meta_insert );
  sql.finalise(stmt_loc_get_meta );

  sql.finalise(stmt_loc_alias_insert);
  sql.finalise(stmt_loc_alias_lookup);
  sql.finalise(stmt_loc_group_alias_lookup);

  sql.finalise(stmt_set_names_fetch);
  sql.finalise(stmt_set_members_fetch);

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
  sql.finalise(stmt_loc_name_list);
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

	sql.bind_int64( stmt_loc_subregion_lookup, ":loc_id" , id );

	while ( sql.step( stmt_loc_subregion_lookup ) )
	{

	    uint64_t id = sql.get_int64( stmt_loc_subregion_lookup , 0 );
	    std::string name = sql.get_text( stmt_loc_subregion_lookup , 2 ) ;
	    int bp1 = sql.get_int( stmt_loc_subregion_lookup , 4 ) ;
	    int bp2 = sql.get_int( stmt_loc_subregion_lookup , 5 ) ;
	    
	    r.addSubRegion( id, name, chr, bp1, bp2 );	    
	    
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


void LocDBase::add_overlap_table(uint64_t group1_id, uint64_t group2_id )
{

    clear_overlaps();

//     if ( group2_id != 0 ) 
// 	clearOverlaps(group1_id,group2_id);
//     else if ( group1_id == 0 )
// 	clearOverlaps();
//     else
// 	clearOverlaps(group1_id);    
	

  // For each overlapping pair of regions, give intersection and union
  // distances From these, any subsequent definition of overlap can be
  // computed Need to handle issue of subregions properly also
  
  // Iterate over all regions, sorted by position
  

    std::set<Region> current;

  int cnt = 0;

  sql.begin();

  while ( sql.step( stmt_loc_iterate ) )
    {
      
      Region r = construct_region( stmt_loc_iterate );

      // Are we imposing constraints on which groups 
      // to consider
      
      if ( group2_id == 0 ) 
      {
	  if ( group1_id != 0 && group1_id != r.group ) 
	      continue;
      }
      else
      {
	  if ( group1_id != 0 && 
	       ! ( group1_id == r.group || group2_id == r.group ) ) 
	      continue;
      }
      
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


  sql.reset( stmt_loc_iterate );

  sql.commit();
  
}



bool LocDBase::range_insertion(const Region & region , uint64_t indiv_id )
{

  sql.bind_text( stmt_loc_insert , ":name" , region.name );
  sql.bind_int64( stmt_loc_insert , ":group_id" , (uint64_t)region.group );
  sql.bind_int( stmt_loc_insert , ":chr" , region.start.chromosome() );
  sql.bind_int( stmt_loc_insert , ":bp1" , region.start.position() );
  sql.bind_int( stmt_loc_insert , ":bp2" , region.stop.position() );
  sql.bind_text( stmt_loc_insert , ":altname" , region.altname );
  
  sql.step( stmt_loc_insert );
  sql.reset( stmt_loc_insert );
  
  // Subregions?
  uint64_t loc_id = sql.last_insert_rowid();
  int chr = region.start.chromosome();

  for (int s = 0 ; s < region.subregion.size(); s++)
    {

       sql.bind_int64( stmt_loc_subregion_insert , ":loc_id" , loc_id );
       sql.bind_text( stmt_loc_subregion_insert , ":name" , region.subregion[s].name );
       sql.bind_int( stmt_loc_subregion_insert , ":chr" , chr );
       sql.bind_int( stmt_loc_subregion_insert , ":bp1" , region.subregion[s].start.position() );
       sql.bind_int( stmt_loc_subregion_insert , ":bp2" , region.subregion[s].stop.position() );
       sql.step( stmt_loc_subregion_insert );
       sql.reset( stmt_loc_subregion_insert );
       uint64_t sub_id = sql.last_insert_rowid();
       
       // Subregion meta-information?
       insertMeta( stmt_loc_submeta_insert , region.subregion[s].meta , sub_id, true );
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


void LocDBase::insertMeta( sqlite3_stmt * s , const MetaInformation<LocMeta> & meta, int id, bool subregion )
{
    
    std::vector<std::string> keys = meta.keys();
    
    for (int j=0; j<keys.size(); j++)
    {
	
	sql.bind_text( stmt_loc_meta_insert_prep , ":name" , keys[j] );
	
	uint64_t field_id = 0;
	
	if ( sql.step( stmt_loc_meta_insert_prep ) )
	{
	    field_id = sql.get_int64( stmt_loc_meta_insert_prep , 0 );
	}
	else // insert this as a new field
	{
	    meta_index_t midx = MetaInformation<LocMeta>::field( keys[j] );
	    sql.bind_text( stmt_loc_meta_insert_prep2 , ":name" , keys[j] );
	    sql.bind_int( stmt_loc_meta_insert_prep2 , ":type" , midx.mt );
	    sql.bind_int( stmt_loc_meta_insert_prep2 , ":number" , midx.len );
	    sql.bind_text( stmt_loc_meta_insert_prep2 , ":description" , midx.description );
	    sql.step( stmt_loc_meta_insert_prep2 );		
	    //mtmap[ i->first ] = sql.last_insert_rowid();
	    field_id = sql.last_insert_rowid();
	    sql.reset(stmt_loc_meta_insert_prep2);

	}
	
	sql.reset( stmt_loc_meta_insert_prep );
	
	// Meta-value
	// NOTE: for now, only let a single value here, not a vector...
	
	mType mt = MetaInformation<LocMeta>::type( keys[j] );
	
	
	// Now we've got field_id, insert actual value

	sql.bind_int64( s , ":reg_id" , id );
	sql.bind_int64( s , ":field_id" , field_id );
      
	//
	// Note -- should be able to swap in parse_set() here
	//

	switch ( mt ) {
	    case META_INT :
	    {
		sql.bind_int( s , ":value" , meta.get1_int( keys[j] ) );  
		break;
	    }
	    case META_BOOL :
	  {
	      sql.bind_int( s , ":value" , meta.get1_int( keys[j] ) );  
	      break;
	  }
	  case META_FLOAT :
	  {
	      sql.bind_double( s , ":value" , meta.get1_double( keys[j] ) );      
	      break;
	  }
	  default :
	  {
	      sql.bind_text( s , ":value" , meta.get1_string( keys[j] ) );  
	  }
      }
      
      sql.step( s );
      sql.reset( s );	  
    }
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
  
  registerMetatype( PLINKSeq::TRANSCRIPT_FRAME() , 
		    META_INT, 1 , META_GROUP_LOC , "CDS Frame" );
  
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
      
      r.meta.set( "feature" , tok[2] );
      

      // Track strand and frame

      r.meta.set( PLINKSeq::TRANSCRIPT_STRAND() , tok[6] );
      
      int frame = 0;

      // implies "." --> 0  (i.e. for stop_codon )

      if ( Helper::str2int( tok[7] , frame ) ) 
	r.meta.set( PLINKSeq::TRANSCRIPT_FRAME() , frame );
      

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



uint64_t LocDBase::load_regions( const std::string & filename,
				 const std::string & grp, 
				 int col_pos,
				 int col_chr, 
				 int col_bp1, 
				 int col_bp2,
				 int col_name, 
				 int col_sub,
				 int col_indiv,
				 std::map<std::string,int> * meta )
{

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
	  chromosome = t.chromosome();
	  p1 = t.start.position();
	  p2 = t.stop.position();
	}
      else
	{
	  // CHR, BP1 and BP2 specified separately
	  if ( chr ) chromosome = Helper::chrCode( buffer[ col_chr ] ) ;
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
      
      if ( meta )
	{
	  std::map<std::string,int>::iterator m = meta->begin();
	  while ( m != meta->end() )
	    {
	      r.meta.set( m->first , buffer[ m->second ] );
	      ++m;
	    }

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
      
      
      if ( ! ( next1 && next2 ) )
	break;
      
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
  sql.bind_text( stmt_insert_indiv , ":name", indiv_id );
  sql.step( stmt_insert_indiv );
  sql.reset( stmt_insert_indiv );
  return 0;
}

uint64_t LocDBase::lookup_indiv_id( const std::string & name ) 
{
  uint64_t i = 0;
  sql.bind_text( stmt_lookup_indiv_id , ":name" , name );
  while ( sql.step( stmt_lookup_indiv_id ) )
    i = sql.get_int64( stmt_lookup_indiv_id, 0 );
  sql.reset( stmt_lookup_indiv_id );
  return i;
}


void LocDBase::insert_segment( const std::string & indiv_id , Region & segment )
{
  
}



std::set<Region> LocDBase::get_indiv_regions( const std::string & group , const std::string & person )
{
  std::set<Region> segs;
  
  uint64_t grp_id = lookup_group_id( group );
  if ( grp_id == 0 ) return segs;
  
  uint64_t ind_id = lookup_indiv_id( person );
  if ( ind_id == 0 ) return segs;
  
  return get_indiv_regions( grp_id , ind_id );

}

std::set<Region> LocDBase::get_indiv_regions( uint64_t group_id , uint64_t indiv_id )
{
  
  std::set<Region> segments;
  
  sql.begin();
  
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
  return get_regions( id , v.chromosome() , v.position() , v.position() );
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
      int field = sql.get_int( stmt_loc_get_meta , 0 );
      std::string value = sql.get_text( stmt_loc_get_meta , 1 );           
      m.parse_set( meta_fields[field] , value ); 
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
	int field = sql.get_int( stmt_loc_get_submeta , 0 );
	std::string value = sql.get_text( stmt_loc_get_submeta , 1 );           
	m.parse_set( meta_fields[field] , value ); 
    }
  
  sql.reset( stmt_loc_get_submeta );

  return m;
}


uint64_t LocDBase::merge( const std::string & grp_name, const std::string & name, const std::string & mergeField)
{
  
  if ( ! attached() ) return 0;
  uint64_t grp_id = lookup_group_id( grp_name );
  if ( grp_id == 0 ) return 0;
  
  
  // Merge regions with the same "name" into a new group, where each
  // region is now a subregion. Specific meta-information will be
  // lost.
  
  // e.g. make one-region == exon into a gene with exons as subregions

  // if mergeField != "" then use that meta-information field. If not present, 
  // skip that region
  
  bool useName = mergeField == "";


  // Assign a new group ID, or append to an existing one

  uint64_t new_group_id = set_group_id( name );

  sql.bind_int64( stmt_loc_lookup_group, ":group_id" , grp_id );

  get_subregions( false );
  get_meta( true ); // as any meta-info will be copied to the subregion
  
  // Hold new regions in memory

  sql.begin();

  std::map<std::string,Region> merged;
  int cnt = 0;
  while ( sql.step( stmt_loc_lookup_group ) )
    {
      
      Region r = construct_region( stmt_loc_lookup_group );
      
      std::string name;
      
      if ( useName ) 
	name = r.name;
      else
	{
	  if ( ! r.meta.hasField( mergeField ) )
	    continue;
	  name = r.meta.get1_string( mergeField );	  
	}
      
      // Have we seen a region of this name before?
      
      std::map<std::string,Region>::iterator i = merged.find( name );

      if ( i != merged.end() )
	{
	  Region & par = i->second;
	  
	  if ( r.chromosome() != par.chromosome() )
	    continue;
	  
	  if ( r.start.position() < par.start.position() )
	    par.start.position( r.start.position() );
	  
	  if ( r.stop.position() > par.start.position() )
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

	    par.addSubRegion( r );
	    
	    merged.insert( make_pair( name , par ) );
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

  std::map<std::string,Region>::iterator i = merged.begin();
  
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
	if ( buffer[i] == "n/a" || buffer[i] == "N/A" ) 
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

std::string LocDBase::alias( const std::string & query , bool show_keys ) 
{
  return attached() ? Helper::stringizeKeyPairList( lookup_alias( query , 0 ) , show_keys ) : "." ;
}

std::map<std::string,std::string> LocDBase::lookup_alias( const std::string & query , const std::string & alias_group  )
{
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

  std::set<std::string> s;
  if ( ! attached() ) return s;
  int gquery = alias_group_table[ query_group ] ;
  int gtarget = alias_group_table[ alias_group ] ;
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

std::string LocDBase::summary()
{
  std::stringstream ss;
  std::set<GroupInfo> g = group_information();
  std::set<GroupInfo>::iterator i = g.begin();
  while ( i != g.end() )
    {

      ss << "LOCDB\t"
      //	 << "GROUP_N=" << i->idx << "\t"
	 << "NAME=" << i->name << "\t"
	 << "N=" << count( i->idx ) << "\t"
	 << "DESC=" << i->description
	 << "\n";
      ++i;
    }
  
  // Sets
  
  g = set_information();
  i = g.begin();
  while ( i != g.end() )
    {
      std::string grp = lookup_group_id( i->idx );
      
      std::vector<std::string> t = fetch_set_names( grp , i->name );
      
      ss << "LOCDB_SET\t"
	 << "GROUP=" << grp << "\t"
	 << "NAME=" << i->name << "\t"
	 << "N=" << t.size() << "\t"
	 << "DESC=" << i->description << "\n";
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
      sql.bind_int64( stmt_loc_replace_real_name_alternate , ":group_id" , grp );
      sql.bind_text( stmt_loc_replace_real_name_alternate , ":altname" , name );
      sql.bind_text( stmt_loc_replace_real_name_alternate , ":altname" , newname );
      while ( sql.step( stmt_loc_replace_real_name_alternate ) ) { } 
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

std::vector<std::string> LocDBase::fetch_names( const std::string & loc_group , bool alternate )
{
  
  std::vector<std::string> results;
  if ( ! attached() ) return results;
  uint64_t id = lookup_group_id( loc_group );
  if ( id == 0 ) return results;
  
  sqlite3_stmt * s = alternate ? stmt_loc_altname_list : stmt_loc_name_list ;
  
  sql.bind_int64( s , ":group_id" ,id );
  while ( sql.step( s ) ) 
    {
      results.push_back( sql.get_text( s , 0 ) );
    }
  sql.reset( s );
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
