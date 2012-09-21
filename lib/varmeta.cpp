#include "plinkseq/vardb.h"
#include "plinkseq/helper.h"


int VarDBase::loader_indep_meta( const std::string & filename , int f , const std::string & meta_name )
{
  
  std::map<int,std::string> files = fetch_files();
  
  if ( f != 0 && files.find( f ) == files.end() ) 
    Helper::halt("could not find file specified by --id; implied file_id " + Helper::int2str(f) );
  
  if ( ! Helper::fileExists( filename ) ) 
    Helper::halt("could not find file " + filename ) ;
  
  InFile file( filename );
  

  //
  // Add group info tag
  //

  sql.bind_text(stmt_insert_indep_meta_group , ":name" , meta_name );
  sql.bind_int(stmt_insert_indep_meta_group , ":file_id" , f );
  sql.step( stmt_insert_indep_meta_group );
  sql.reset( stmt_insert_indep_meta_group );
  
  //
  // Retrive group ID
  //

  sql.bind_text( stmt_fetch_indep_meta_group , ":name" , meta_name );
  sql.bind_int( stmt_fetch_indep_meta_group , ":file_id" , f );
  sql.step( stmt_fetch_indep_meta_group );
  int g = sql.get_int( stmt_fetch_indep_meta_group , 0 );
  sql.reset( stmt_fetch_indep_meta_group );
  

  // Read in variant information.  If we are given a clue about the
  // meta-field, insert it from the header.  Otherwise, assume it is
  // text


  sql.begin();
  
  while ( ! file.eof() )
    {

      std::string line = file.readLine();
      
      if ( line == "" ) continue;
      if ( file.eof() ) break;
      
      // Is this a header field? 
      // ##Name,Len,Type,"Description"
      // ##DB,1,Flag,"In dbSNP"
      
      if ( line.substr(0,2) == "##" ) 
	{
	  process_indep_meta_header( line , f );
	  continue;
	}
      
      // Otherwise, assume contains data, tab-delimited
      //  Variant   Key  Value
      
      std::vector<std::string> tok = Helper::parse( line , "\t" );
      
      if ( tok.size() != 3 ) 
	{
	  plog.warn("skipping row of input, not 3 tab-delimited fields",line);
	  continue;
	}
      
      
      // Resolve variant specifier -- positional or ID-based?

      bool is_region = false;
      
      Region r( tok[0] , is_region );
      
      std::set<uint64_t> var_ids;
      
      // Lookup based on region; just allow start-base lookup for now, positionally

      if ( is_region )
	{
	  sql.bind_int( stmt_fetch_var_from_position , ":chr" , r.chromosome() );
	  sql.bind_int( stmt_fetch_var_from_position , ":bp1" , r.start.position() );

	  while ( sql.step( stmt_fetch_var_from_position ) )
	    {
	      int var_id = sql.get_int64( stmt_fetch_var_from_position , 0 );
	      int file_id = sql.get_int64( stmt_fetch_var_from_position , 1 );
	      if ( f == 0 || f == file_id ) var_ids.insert( var_id );
	    }
	  
	  sql.reset( stmt_fetch_var_from_position );	  
	}
      else // or based on ID
	{

	  sql.bind_text( stmt_fetch_var_from_name , ":name" , tok[0] );

	  while ( sql.step( stmt_fetch_var_from_name ) )
	    {
	      int var_id = sql.get_int64( stmt_fetch_var_from_name , 0 );
	      int file_id = sql.get_int64( stmt_fetch_var_from_name , 1 );
	      if ( f == 0 || f == file_id ) var_ids.insert( var_id );
	    }

	  sql.reset( stmt_fetch_var_from_name );
	}      


      //
      // Insert for each specified variant/file combination
      //

      std::string skey = tok[1];
      std::string value = tok[2];
      
      //
      // Do we recognise this key?  If not, default to single textual entry
      //
      
      int nkey = 0;
      
      if ( indep_metamap.find( skey ) == indep_metamap.end() ) 
	{
	  process_indep_meta_header( "##" + skey + ",1,String," , f  );	      
	}

      nkey = indep_metamap[ skey ];
      

      //
      // A valid value?
      //
      
  //     meta_index_t midx = MetaInformation<VarMeta>::field( skey );	  
      
//       if ( midx.mt == META_INT ) 
// 	{
	  
// 	}

// 	|| midx.mt == META_BOOL )
// 	target.meta.set( key , sql.get_int( stmt_fetch_indep_meta_value , 1 ) );
// 	  else if ( midx.mt == META_FLOAT ) 
// 	    target.meta.set( key , sql.get_double( stmt_fetch_indep_meta_value , 1 ) );
// 	  else if ( midx.mt == META_FLAG && sql.get_int( stmt_fetch_indep_meta_value , 1 ) != 0 ) 
// 	    target.meta.set( key );
// 	  else // META_TEXT as default
// 	    target.meta.set( key , sql.get_text( stmt_fetch_indep_meta_value , 1 ) );
  
      
      
      //
      // Insert actual value (as text for now)
      //
      
      std::set<uint64_t>::iterator i = var_ids.begin();      
      while ( i != var_ids.end() )
	{
	  
	  sql.bind_int64( stmt_insert_indep_meta_value , ":group_id" , g );
	  sql.bind_int64( stmt_insert_indep_meta_value , ":var_id" , *i );
	  sql.bind_int64( stmt_insert_indep_meta_value , ":meta_id" , nkey );
	  sql.bind_text( stmt_insert_indep_meta_value , ":value" , value );
	  sql.step( stmt_insert_indep_meta_value );
	  sql.reset( stmt_insert_indep_meta_value );
	  
	  // next value/var-id
	  ++i;
	}
                  
    }
  
  sql.commit();

}


void VarDBase::populate_indep_metadata_map()
{

  // meta-groups
  
  indep_group_metamap.clear();
  while ( sql.step( stmt_dump_indep_meta_group ) )
    {
      int g = sql.get_int( stmt_dump_indep_meta_group , 0 );
      std::string n = sql.get_text( stmt_dump_indep_meta_group , 1 );
      int f = sql.get_int( stmt_dump_indep_meta_group , 2 );
      indep_group_metamap[n] = int2(g,f);      
    }
  sql.reset( stmt_dump_indep_meta_group ) ;
  

  // meta-types
  indep_metamap.clear();
  reverse_indep_metamap.clear();
    
  while ( sql.step( stmt_fetch_indep_meta_type ) )
    {

      int id = sql.get_int( stmt_fetch_indep_meta_type , 0 ) ;
      std::string name = sql.get_text( stmt_fetch_indep_meta_type , 1 ) ;
      int len = sql.get_int( stmt_fetch_indep_meta_type , 2 ) ;
      mType mt = (mType)sql.get_int( stmt_fetch_indep_meta_type , 3 ) ;
      std::string desc = sql.get_text( stmt_fetch_indep_meta_type , 4 ) ;      
      
      MetaInformation<VarMeta>::field( name , mt , len , desc );
      
      indep_metamap[ name ] = id;
      reverse_indep_metamap[ id ] = name;

    }

  sql.reset( stmt_fetch_indep_meta_type );  

}



bool VarDBase::attach_indep_metadata( const uint64_t & svar_id , 
				      SampleVariant & target , 
				      Variant & parent , 
				      const std::set<std::string> * grps )
{

  

  sql.bind_int64( stmt_fetch_indep_meta_value , ":var_id" , svar_id );
  
  while ( sql.step( stmt_fetch_indep_meta_value ) )
    {

      int meta_id = sql.get_int( stmt_fetch_indep_meta_value , 0 );
      
      if ( reverse_indep_metamap.find( meta_id ) != reverse_indep_metamap.end() )
	{
	  
	  std::string key = reverse_indep_metamap[ meta_id ] ;
	  
	  // if only a subset of meta-fields specified
	  if ( grps && grps->find( key ) == grps->end() ) continue;

	  meta_index_t midx = MetaInformation<VarMeta>::field( key );	  
	  
	  MetaInformation<VarMeta> & m = MetaMeta::static_variant( key ) 
	    ? parent.meta : target.meta ; 
	  
	  if ( midx.mt == META_INT || midx.mt == META_BOOL )
	    m.set( key , sql.get_int( stmt_fetch_indep_meta_value , 1 ) );
	  else if ( midx.mt == META_FLOAT ) 
	    m.set( key , sql.get_double( stmt_fetch_indep_meta_value , 1 ) );
	  else if ( midx.mt == META_FLAG && sql.get_int( stmt_fetch_indep_meta_value , 1 ) != 0 ) 
	    {
	      m.set( key );
	    }
	  else // META_TEXT as default
	    m.set( key , sql.get_text( stmt_fetch_indep_meta_value , 1 ) );
	  
	}

    }
  sql.reset( stmt_fetch_indep_meta_value );
  return true;
}


int VarDBase::process_indep_meta_header( const std::string & line , const int f )
{

  //
  // strip leading ## and tokenize
  //

  std::vector<std::string> tok = Helper::quoted_parse( line.substr(2) );
   
  std::string name = "";
  int num = -9;
  std::string type = "";
  mType mt = META_UNDEFINED;
  std::string desc;
  
  if ( tok.size() != 4 ) return 0;

  name = tok[0];
  if ( ! Helper::str2int( tok[1] , num ) ) num = -1;
  type = tok[2];
  desc = tok[3];
  
  if ( Helper::is_int( type ) ) mt = META_INT;
  else if ( Helper::is_float( type ) ) mt = META_FLOAT;
  else if ( Helper::is_text( type ) ) mt = META_TEXT;
  else if ( Helper::is_flag( type ) ) { mt = META_FLAG; }
  

  //
  // Does this contain valid information?
  //

  if ( name == "" || mt == META_UNDEFINED || num < -1 ) return 0;


  //
  // Insert in DB as meta-field and register with internal meta-field
  // store (for either variant, or genotype)
  //

  sql.bind_text( stmt_insert_indep_meta_type , ":name" , name );
  sql.bind_int( stmt_insert_indep_meta_type , ":length" , num );
  sql.bind_int( stmt_insert_indep_meta_type , ":type" , mt );
  sql.bind_text( stmt_insert_indep_meta_type , ":desc" , desc );
  sql.step( stmt_insert_indep_meta_type );
  sql.reset( stmt_insert_indep_meta_type );


  MetaInformation<VarMeta>::field( name , mt , num , desc );
  
  
  // 
  // Return meta ID 
  //
  
  populate_indep_metadata_map();

  return 0;

}  



int VarDBase::flush_indep_meta( const std::string & name )
{
  
  sqlite3_stmt * s1 = sql.prepare( "SELECT group_id FROM indep_meta_groups WHERE name == :name ;" ) ;
  sql.bind_text( s1 , ":name" , name );
  int grp_id = 0;
  if ( sql.step( s1 ) ) 
    grp_id = sql.get_int( s1 , 0 );
  sql.finalise(s1);
  
  if ( grp_id == 0 ) return 0;
  sql.query(" DELETE FROM indep_meta_data WHERE group_id == " + Helper::int2str( grp_id ) + " ; " );
  sql.query(" DELETE FROM indep_meta_groups WHERE group_id == " + Helper::int2str( grp_id ) + " ; " );
  
  return 1;
}

void VarDBase::flush_indep_meta( )
{
  sql.query( " DELETE FROM indep_meta_groups; " );
  sql.query( " DELETE FROM indep_meta_types; " );
  sql.query( " DELETE FROM indep_meta_data; " );
}
