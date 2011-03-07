#include "gstore.h"

#include "varfunc.h"

std::vector<bool> VarFunc::missing_genotype_mask( const Variant & v )
{
  std::vector<bool> b( v.size() , false );
  for (int i=0; i<v.size(); i++) if ( v(i).null() ) b[i] = true;
  return b;
}

std::vector<bool> VarFunc::missing_genotype_mask( const VariantGroup & vars )
{
  // case-wise missingness
  const int n = vars.n_individuals();
  std::vector<bool> b( n , false );
  for (int v=0; v<vars.size(); v++)
    {
      const Variant & var = vars(v);
      for (int i=0; i<n; i++) 
	if ( var(i).null() ) b[i] = true;
    }
  return b;
}


struct D1 {
  std::set<Region> * regions;
  LocDBase *    locdb;
  uint64_t      grp_id;
};



void f_add_set( Variant & v , void * p )
{

//     D1 * d = (D1*)p;    
    
//     Region r( v.chromosome() , v.position() , v.position() );    
    
//     r.group = d->grp_id;
    
//     // Add the variant ID to this region, so we can keep track of 
//     // which variant is which 
    
//     r.meta.set( "index" , uint64_t2str( v.index() ) );
    
//     // If this variant overlaps at least 1 region, insert
//     // into locus-database
    
//     if ( r.within( *(d->regions) ) )
//       {
//       d->locdb->range_insertion(r);
//     }
}


struct D2 {
  uint64_t                 var_grp_id;
  VarDBase *               vardb;
  std::map<uint64_t,uint64_t> * region_table;
};

void f_get_overlap(Region & r1, Region & r2, 
		   int r_intersection, int r_union , 
		   void * p)
{
    
  
    // Two regions overlap -- figure out which is the variant, then add that
    // to the loc-database
  
//   D2 * d = (D2*)p;
  
//   if ( r1.group == d->var_grp_id && r2.group != d->var_grp_id ) 
//     {
//       uint64_t region_id = (*d->region_table)[r2.id];	
//       uint64_t var_id;
//       str2uint64_t( r1.meta.get1_string( "index" ) , var_id );
//       d->vardb->set_add_variant( region_id , var_id );
//     }
//   else if ( r1.group != d->var_grp_id && r2.group == d->var_grp_id ) 
//     {
//       uint64_t region_id = (*d->region_table)[r1.id];	
//       uint64_t var_id;
//       str2uint64_t( r2.meta.get1_string( "index" ) , var_id );
//       d->vardb->set_add_variant( region_id , var_id );
//     }
  
}


void GStore::vardb_make_set( const std::string & loc_grp_name, const std::string & name )
{
 
//   if ( ! locdb.attached() ) return;

//   if ( ! vardb.attached() ) return;
  
//   uint64_t loc_grp_id = locdb.lookup_group_id( loc_grp_name );

//   if ( loc_grp_id == 0 ) return;

  
//   // 1) Extract all regions from loc-db, as specified by loc_grp_id
  
//   // 2) Create a new set in var-db, and insert all regions as members
  
//   // 3) Iterate over all variants, and get a list of those that overlap
//   //    at least 1 region in the group

//   //   This is designed to work on the level of 1 SampleVariant at a time
  
//   // 4) Insert those into loc-db as a new, temporary set, then 
//   //    construct an overlap table between those two groups
  
//   // 5) Use that table to populate the set group in var-db
  
//   // 6) Remove temporary overlap list and variant set in loc-db
  
  
//   //
//   // 1) Extract all regions from loc-db, as specified by loc_grp_id
//   //
  
  
//   set<Region> reg = locdb.get_regions( loc_grp_id );
  

//   //
//   // 2) Create a new set in var-db, and insert all regions as members
//   //
  
//   uint64_t set_grp_id = vardb.set_group_id( name );
  
//   map<uint64_t,uint64_t> region_table;
  
//   vardb.begin();
//   set<Region>::iterator i = reg.begin();
//   while ( i != reg.end() )
//     {
//       uint64_t id = vardb.set_member_id( set_grp_id , i->name );
//       region_table[ i->id ] = id;
//       ++i;
//     }
//   vardb.commit();
  
    
//   //
//   // 3) Iterate over all variants, and get a list of those that overlap
//   //    at least 1 region in the group
//   //
  
//   uint64_t temp_var_id = locdb.set_group_id( "TMP" );
  
//   D1 d;
//   d.regions = &reg;
//   d.locdb   = &locdb;
//   d.grp_id  = temp_var_id;
  
//   locdb.begin();
//   vardb.iterate( f_add_set , &d );
//   locdb.commit();
  
//   //
//   // 4) Having inserted those into loc-db as a new, temporary set, then
//   //    construct an overlap table between those two groups
//   //
  
//   locdb.get_meta(false);
//   locdb.get_subregions(false);
  
//   locdb.add_overlap_table( temp_var_id , loc_grp_id );
  
  
//   //
//   // 5) Use that table to populate the set group in var-db
//   //
  
//   D2 d2;
//   d2.var_grp_id   = temp_var_id;
//   d2.vardb        = &vardb;
//   d2.region_table = &region_table;
  

//   //
//   // We need to access the meta-information (to tie the VAR-region
//   // back to a VAR)
//   //
  
//   locdb.get_meta(true);
  

//   //
//   // Perform inside a VAR transaction, as we will potentially be
//   // performing a large number of insertions
//   //

//   vardb.begin();
//   locdb.get_regions_and_overlap( f_get_overlap , &d2 );
//   vardb.commit();
  
//   locdb.get_subregions(true);
  
  
//   //
//   // 6) Remove temporary overlap list and variant set in loc-db
//   //
  
//   locdb.flush( temp_var_id );
//   locdb.clear_overlaps();

}


void g_consolidate( Variant & var , void * p )
{

//   VarDBase * db = (VarDBase*)p;
  
//   //  vars.collapse();
  
//   // Should just be a single variant now
  
//   //   plog >> vars << "\n"
//   //        << "--------------------\n\n";
  
//   plog >> var.chromosome() << "\t" << var.position() <<"\n";

//   // TODO: although, when writing a variant, we probably want to do 
//   // better than this??   Or, no perhas this is the point...

//   db->insert_consensus( 1 , var );  
  
}

void GStore::vardb_consolidate( const std::string & name , Mask & mask )
{
  
//   // Create new database
  
//   VarDBase newdb;
  
//   newdb.attach( name );
  
  
//   // Add headers and meta-information, etc

// // vardb->insert_header( 1, tok[0] , "" ); 
// // vardb->insert_header( file_id , tok[0] , tok[1] ); 
// // vardb->insert_metatype( file_id , name , mt, num, mgrp, desc );
// // vardb->insert_header( file_id , s , "" ); 
  
//   // Ensure iteration will group by genomic variant position
  
//   mask.group_file();

  
//   // Add individuals

//   IndividualAlignment alignment;
//   int n_uniq_indiv = vardb.populate_individual_alignment( alignment , mask );

//   map<int,string> imap = alignment.map_slot_to_id();
//   map<int,string>::iterator i = imap.begin();
  
//   while ( i != imap.end() )
//     {
//       Individual person( i->second );

//       newdb.insert( 1 , person );

//       ++i;
//     }


//   // Add variants
  
//   vardb.iterate( g_consolidate , &newdb , mask );


//   return;	 

//   // Update file-index

  
//  exit(0);
}
  

