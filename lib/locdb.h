#ifndef __LOCDB_H__
#define __LOCDB_H__

#include "sqlwrap.h"
#include "regions.h"

#include <string>

class LocDBase {
  
 public:
  
  LocDBase()
    {
      // sql.version();
      vget_meta = true;
      vget_subregions = true;
    }
  
  ~LocDBase()
    {
      dettach();
      sql.close();
    }
  
  bool wipe( const std::string & name );
  bool attach( const std::string & name);
  bool dettach();
  
  bool init();
  bool release();
  
  bool attached() { return sql.is_open(); }
  
  std::string summary();

  // Queries

  // bool initIterateVariants( Region region = Region() );
  // uint64_t offset(uint64_t);
  
  void append_metainformation( Variant & , const std::set<int> & );
  
  void replace_real_names( const int group_id , 
			   const std::string & id , 
			   const std::string & real_name , 
			   bool alternate_search = false );

  std::vector<std::string> fetch_names( const std::string & loc_group , bool alternate = false);

  std::vector<std::string> fetch_name_given_altname( const std::string & loc_group , const std::string & altname );
  
  std::vector<Region> fetch_real_names( const std::string & g , const std::string & altname );

  std::vector<Region> fetch( const std::string & grp, const std::vector<std::string> & names );

  bool  contains( const std::string & grp , const int chr , const int bp1 , const int bp2 );
  bool  contains( const int grp , const int chr , const int bp1 , const int bp2 );
  
  //
  // Group level function
  //

  
  // Get a group ID by matching group name, creating a new entry if not already done
  
  uint64_t set_group_id(const std::string & grp, const bool temp = false , const std::string & desc = "n/a");

  // Only fetch an ID, returning 0 if does not exist
  uint64_t lookup_group_id(const std::string & grp);
  std::string lookup_group_id(const int );

  uint64_t set_set_id(const std::string & name, 
		      const int set_group_id , 
		      const bool temp , 
		      const std::string & desc );
  
  uint64_t lookup_set_id(const std::string & grp, const std::string & sname);
  
  // Group  (i.e. KEGG)
  // Set    (i.e. a specific KEGG pathway)
  // Member (i.e. a gene in a specific KEGG pathway)

  std::vector<std::string> fetch_set_names( const std::string & loc_group, 
					    const std::string & set_group );
  
  std::vector<std::string> fetch_set_members( const std::string & loc_group, 
					      const std::string & set_group,
					      const std::string & set_name );
  
  // Handle temporary groups
  
  void temporary(const uint64_t id, const bool b);
  bool temporary(const uint64_t id);
  

  // Delete all temporary groups

  void flush();
  
  
  // Delete a particular group
  
  void flush(const uint64_t id);
  
  void flush( const std::string & grp );
   
   
   //
   // Insertions of ranges
   //
  
   bool     range_insertion(const Region & , uint64_t indiv_id = 0 );

   void insertMeta( sqlite3_stmt * s , const MetaInformation<LocMeta> & , int id, bool subregion = false );


   //
   // Individual-based segments
   //

   uint64_t insert_indiv( const std::string & indiv_id );
   
   void insert_segment( const std::string & indiv_id , Region & segment );
   
   std::set<Region> get_indiv_regions( uint64_t , uint64_t );

   std::set<Region> get_indiv_regions( const std::string & group_id , const std::string & indiv_id );

   uint64_t lookup_indiv_id( const std::string & );


   //
   // Meta-information
   //

   MetaInformation<LocMeta> meta(uint64_t loc_id );
   MetaInformation<LocMeta> submeta(uint64_t sub_id );

   //
   // Wrapper functions to input a large number of ranges from a file
   //

   uint64_t load_GTF(const std::string & name, const std::string & grp, bool use_transcript_id = true);

   uint64_t load_set(const std::string &, const std::string &, const std::string &, bool use_altname = false );

   uint64_t load_regions(const std::string & name, 
			 const std::string & grp, 
			 int col_pos = -1,
			 int col_chr = 0,
			 int col_bp1 = 1,
			 int col_bp2 = 2, 
			 int col_name = 3,
			 int col_sub = -1, 
			 int col_indiv = -1,
			 std::map<std::string,int> * meta = NULL );
   
   
   //
   // Functions to handle sets of regions
   //

   std::set<Region> get_set( uint64_t set_id );  
   

   //
   // Name alias lookup functions
   //

   void load_alias( const std::string & filename );
   std::map<std::string,std::string> lookup_alias( const std::string & , const uint64_t alias_group_id );
   std::map<std::string,std::string> lookup_alias( const std::string & , const std::string & alias_group = "" );

   std::set<std::string> targetted_lookup_alias( const std::string & query , 
						 const std::string & query_group , 
						 const std::string & alias_group );
   std::string alias( const std::string & , bool show_key = false );
   

   void delete_aliases();

   uint64_t insert_alias_group( const std::string & alias );
   void clear_alias_groups();
   void read_alias_groups();

      

   //
   // Output functions for client
   //

   std::set<Region> get_regions( const std::string & );
   std::set<Region> get_regions(uint64_t grp_id);
   std::set<Region> get_regions(uint64_t grp_id,int,int,int);
   std::vector<uint64_t> get_region_ids(uint64_t grp_id,int,int,int);
   
   std::set<Region> get_regions( const std::string & group , const Region & );
   std::set<Region> get_regions( const std::string & group , const Variant & );

   Region get_region( const std::string &, const std::string & );
   Region get_region( const int , const std::string & );

   std::set<Region> get_overlaps(uint64_t loc_id);

   std::set<Region> get_overlaps(uint64_t loc_id, OverlapDefinition & );

   
   bool get_regions_and_overlap( void (*f)( Region&,Region&, int, int, void * ) ,
				 void * data );


   // Extract regions from a group based on name, and create a 
   // new group, returning that ID

   uint64_t extract(uint64_t , uint64_t , const std::string & );
   
   // Merge two groups into a new one

   uint64_t merge(uint64_t grp1_id, uint64_t grp2_id , const std::string & name );
   
   // Alter name based on an alias table
   
   uint64_t rename(uint64_t grp_id, uint64_t alias_id , const std::string & name );
   
   //
   // Operations
   // 
   
   // Create an overlap table of N groups
   
   bool clear_overlaps();
   bool clear_overlaps(uint64_t);
   bool clear_overlaps(uint64_t,uint64_t);
   
   void add_overlap_table(uint64_t group1_id = 0, uint64_t group2_id = 0);
   
   // Name-based merging of Regions -> Region/Subregions

   uint64_t merge( const std::string & grp_name, const std::string & , const std::string & mergeField = "");

   // Overlap-based merging of Regions -> Region (optionally, storing
   // original as subregions)

   uint64_t merge_overlap(uint64_t grp_id, const std::string & name, bool storeSubregions);
   

   //
   // Other
   //

   bool index();
   bool drop_index();
   
   int count(uint64_t grp = -1);
   uint64_t span(uint64_t grp = -1);


   //
   // Helper functions
   //

   void set_metatypes( bool clear = false);

   void begin() { sql.begin(); }
   void commit() { sql.commit(); }
   
   void get_meta(bool b) { vget_meta = b; }
   void get_subregions(bool b) { vget_subregions = b; }

   std::set<GroupInfo> group_information();
   std::set<GroupInfo> set_information();

   std::string filename() const { return fname; }
   
   void populate_meta_field_map()
       {
	   meta_fields.clear();
	   sqlite3_stmt * s = 
	       sql.prepare(" SELECT field_id, name FROM metatypes;" );
	   while ( sql.step(s) )
	   {
	       int k = sql.get_int( s , 0 );
	       std::string n = sql.get_text( s , 1 );
	       meta_fields[k] = n;
	   }
	   sql.finalise(s);
       }

 private:
  
  // Assumes a single database connection

  SQL sql;
  std::string fname;

  bool vget_meta;
  bool vget_subregions;
  
  std::map<int,std::string> meta_fields;

  std::map<std::string,int> alias_group_table;
  std::map<int,std::string> alias_group_reverse_table;


  // Prepared queries
  
  sqlite3_stmt * stmt_loc_insert; 
  sqlite3_stmt * stmt_loc_lookup_name;
  sqlite3_stmt * stmt_loc_lookup_group_id;
  sqlite3_stmt * stmt_loc_lookup_group;
  sqlite3_stmt * stmt_loc_lookup_group_and_name;
  sqlite3_stmt * stmt_loc_lookup_id_group_and_range;

  sqlite3_stmt * stmt_loc_replace_real_name;
  sqlite3_stmt * stmt_loc_replace_real_name_alternate;

  sqlite3_stmt * stmt_loc_name_list;
  sqlite3_stmt * stmt_loc_altname_list;
  sqlite3_stmt * stmt_loc_lookup_real_name;
  sqlite3_stmt * stmt_loc_lookup_real_name_only;
  
  // Individuals/segments

  sqlite3_stmt * stmt_insert_indiv;
  sqlite3_stmt * stmt_lookup_indiv_id;

  sqlite3_stmt * stmt_insert_segment;
  sqlite3_stmt * stmt_fetch_segment;

  // Locus intersection
    
  sqlite3_stmt * stmt_loc_intersect;

  sqlite3_stmt * stmt_loc_lookup_set;
  sqlite3_stmt * stmt_loc_set_insert1;
  sqlite3_stmt * stmt_loc_set_insert2;
  
  sqlite3_stmt * stmt_loc_lookup_group_with_overlap;

  sqlite3_stmt * stmt_loc_lookup_group_with_overlap_p1;
  sqlite3_stmt * stmt_loc_lookup_group_with_overlap_p2;

  sqlite3_stmt * stmt_loc_lookup_id;
  sqlite3_stmt * stmt_loc_lookup_range;
  sqlite3_stmt * stmt_loc_lookup_group_and_range;

  sqlite3_stmt * stmt_loc_alias_lookup;
  sqlite3_stmt * stmt_loc_group_alias_lookup;
  sqlite3_stmt * stmt_loc_alias_insert;

  sqlite3_stmt * stmt_loc_iterate;
  sqlite3_stmt * stmt_loc_iterate_group;
  sqlite3_stmt * stmt_loc_iterate_overlap;

  sqlite3_stmt * stmt_loc_subregion_insert; 
  sqlite3_stmt * stmt_loc_subregion_lookup;

  // Group information

  sqlite3_stmt * stmt_loc_group_list;
  sqlite3_stmt * stmt_loc_insert_group_name; 
  sqlite3_stmt * stmt_loc_lookup_group_name; 
  sqlite3_stmt * stmt_loc_remove_group1; 
  sqlite3_stmt * stmt_loc_remove_group2; 
  sqlite3_stmt * stmt_loc_update_temp_status;
  sqlite3_stmt * stmt_loc_lookup_temp_status;

  // Aliases

  sqlite3_stmt * stmt_loc_alias_group_insert;
  sqlite3_stmt * stmt_loc_alias_group_dump;
  sqlite3_stmt * stmt_loc_targetted_group_alias_lookup;

  // Sets

  sqlite3_stmt * stmt_set_group_insert;
  sqlite3_stmt * stmt_set_group_lookup;
  sqlite3_stmt * stmt_set_member_lookup;
  sqlite3_stmt * stmt_set_member_insert;
  sqlite3_stmt * stmt_set_members_fetch;
  sqlite3_stmt * stmt_set_names_fetch;
  sqlite3_stmt * stmt_set_data_insert;
  sqlite3_stmt * stmt_dump_all_sets;

  // Meta-information

  sqlite3_stmt * stmt_loc_meta_insert_prep;
  sqlite3_stmt * stmt_loc_meta_insert_prep2;
  sqlite3_stmt * stmt_loc_meta_insert;
  sqlite3_stmt * stmt_loc_submeta_insert;
  sqlite3_stmt * stmt_loc_get_meta;
  sqlite3_stmt * stmt_loc_get_submeta;
  sqlite3_stmt * stmt_fetch_metatypes;

  // Overlap table

  sqlite3_stmt * stmt_loc_overlap_insert;
  sqlite3_stmt * stmt_loc_overlap_lookup;

  sqlite3_stmt * stmt_loc_insert_set_group_name;

  // 
  // Helper functions
  //

  Region construct_region(sqlite3_stmt*);

};


#endif
