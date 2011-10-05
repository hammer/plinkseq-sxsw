#ifndef __PSEQ_FUNC_H__
#define __PSEQ_FUNC_H__

#include "pseq.h"

#include <string>

class OptGMatrix;
class OptUniq;

namespace Pseq 
{
  
  void finished();

  namespace Util {
    class Options;
    class Commands;
    class ArgMap;
  }

  class VStat;    
  
  bool new_project( std::string project , Pseq::Util::ArgMap & args );

  bool set_project( std::string project );

  namespace VarDB
    {
      
      bool summary( Mask & );

      bool attach( std::string db );

      bool load_VCF( Mask & );
      
      bool header_VCF( const bool meta , const bool header , Mask & );
      
      bool load_PLINK( const std::vector<std::string> & , const Pseq::Util::Options & , const std::string & );
      
      bool flush( const std::vector<int> & );
      
      bool vacuum();

      bool write_VCF(Mask & m);

      bool write_BCF(Mask & m, const std::string & );

      bool write_PED(Mask & m, std::string, bool);

      bool write_lik(Mask & m);

      bool write_matrix(Mask & m);

      bool write_meta_matrix(Mask & m);

      bool write_var_meta_matrix(Mask & m, std::string & name);

      bool write_gene_matrix(Mask & m, OptGMatrix & opt);
      
      bool dump_indiv();

      bool consolidate(Mask & m, std::string label);

      bool uniq_report( std::vector<std::string> & , Mask & m , OptUniq & );
      
      bool gene_stats_header( VStat & );

      bool variant_stats(Mask & m);

      bool annotate_loc( const std::string & , Mask & m);
      
      bool write_vardb( const std::string & , const std::string & , Mask & m);
      
      bool check_concordance(Mask &m);
      
      bool cluster_scan(Mask &m);

      bool proximity_scan(Mask &m);
      
      bool make_counts_file(Mask &m , const std::string & );

      bool simple_counts(Mask &m, bool);

      bool lookup_list(const std::string & filename, Mask &m);

      bool simple_sim();

      bool vdist_summary(Mask & m, long int );

    }    
  
  namespace PPH2DB
  {
    bool load( const std::string & dbname , const std::string & filename );
    bool score( Mask & m , const std::string & dbname ); 
  }

  namespace NetDB
  {

    bool loader( const std::string & db , const std::string & file );

    bool lookup( const std::string & db , const std::string & gene , const std::string & grp );

  }


  namespace LocDB
    {

      bool summary( LocDBase * );      

      bool attach( std::string db );

      bool load_GTF( std::string gtf , std::string label , bool );

      bool load_generic_regions( std::string & filename , const std::string & label , Pseq::Util::Options & , bool );

      bool load_segments( std::string filename , std::string label , Pseq::Util::Options & opt );

      bool load_set( std::string file , std::string set_label , std::string group_label , bool , bool );

      bool merge( std::string label1 , std::string label2 , bool );

      bool index();

      bool swap_alternate_names( const std::string & group , const std::string & filename );

      bool overlap_analysis( std::string label1, std::string label2 );

      bool intersection( std::string filename , std::string group , bool );
      
      bool load_pathway( std::string file1 , std::string file2 , std::string group , bool , bool use_altname = false );

      bool loc_view( const std::string & , const std::vector<std::string> & );
    }


  namespace RefDB
    {
      
      bool summary();      
      
      bool attach( std::string );
      
      bool load_refvar( const std::string & filename , const std::string & name , Pseq::Util::Options & opt );
      
      bool load_VCF( const std::string & filename , const std::string & name );
    }
  
  namespace IndDB
    {
      
      bool attach( std::string );
      
      bool summary();
      
      bool set_phenotype( std::string );

      bool make_phenotype( std::string );
      
      bool load_ped_info( std::string );
      
      bool load_phenotypes( std::string );
      
      bool dump_table(Mask &);
      
    }
  
  

  namespace Util
    {
      
      Mask construct_mask( std::string );
      
      void set_default( VStat & );

      bool file_summary();
      
      bool meta_summary();
      
      // helper class for Options and ArgMap
      // (to let Commands know which options/arguments are appropriate for a 
      //  given command)

      struct opt_t { 	
        opt_t( const std::string & name , const std::string & type ) : name(name) , type(type) { } 
	std::string name;  // e.g. file 
	std::string type;   //      str
	bool operator<( const opt_t & rhs ) const { return name < rhs.name; } 
      };
      
      
      class Options {
	
      public:
	
	void set(  const std::vector<std::string > & s ) 
	{
	  // asusme all key=value,key,pairs
	  for (int i=0; i<s.size(); i++)
	    {
	      if ( s[i].find("=") != std::string::npos )
		{
		  const std::string key = s[i].substr(0, s[i].find("=") );
		  const std::string vallist = s[i].substr( s[i].find("=") + 1 );
		  optdata[ key ] = Helper::parseCommaList( vallist );
		  simpledata[ key ] = vallist;
		}
	      else
		{
		  // evoke default set constructor
		  optdata[ s[i] ];		  
		}
	    }
	}
	
	void dump() const
	{
	  std::map< std::string, std::set< std::string> >::const_iterator i = optdata.begin();
	  while ( i != optdata.end() )
	    {
	      plog << i->first << "\t";
	      std::set<std::string>::iterator j = i->second.begin();
	      while ( j != i->second.end() )
		{
		  plog << *j << " ";
		  ++j;
		}
	      plog << "\n";
	      ++i;
	    }
	}


	bool key( const std::string & k ) const
	{
	  std::map< std::string, std::set< std::string > >::const_iterator i = optdata.find(k);
	  return  i != optdata.end();
	}
	
	bool value( const std::string & k , const std::string & val ) const
	{
	  
	  std::map< std::string, std::set< std::string > >::const_iterator i = optdata.find(k);
	  if ( i == optdata.end() ) return false;
	  return i->second.find(val) != i->second.end();
	}
	
	std::set<std::string> get_set( const std::string & k ) const
	  {
	    std::set<std::string> s;
 	    std::map< std::string, std::set< std::string > >::const_iterator i = optdata.find(k);
	    if ( i == optdata.end() ) return s;
	    else return i->second;
	  }
	
	std::string simple_string( const std::string & k ) const
	  {
	    
 	    std::map< std::string,std::string>::const_iterator i = simpledata.find(k);
	    if ( i == simpledata.end() ) return "";
	    return i->second;
	  }
	
	template<class T> 
        T as( const std::string & k ) const 
	  {
            T t;
            if ( ! key(k) ) return t;
	    std::map< std::string, std::set< std::string > >::const_iterator i = optdata.find(k);
	    // We only know how to handle single items here
	    if ( i->second.size() != 1 ) return t;	    
	    return Helper::lexical_cast<T>( *(i->second.begin()) );	    
	  }
	  
	void attach( const std::string & command , const std::string & opt )
	{
	  // opt will be comma-delimited list, with types
	  // ref=str-list,loc=str-list
	  std::vector<std::string> opts = Helper::char_split( opt , ',' );
	  for (int i=0; i<opts.size(); i++)
	    {
	      std::vector<std::string> o = Helper::char_split( opts[i] , '=' );
	      if ( o.size() == 2 ) 
		comm2opt[ command ].insert( opt_t( o[0] , o[1] ) );
	    }
	}
	
	std::string attached( const std::string & command ) 
	  {
	    std::map<std::string,std::set<opt_t> >::iterator i = comm2opt.find( command );
	    std::string s = "";
	    if ( i != comm2opt.end() )
	      {
		std::set<opt_t>::iterator j = i->second.begin();
		while ( j != i->second.end() )
		  {
		    s += "OPT\t" + i->first + "\t" + j->name + "\t" + j->type + "\n";
		    ++j;
		  }      
	      }
	    return s;
	  }

      private:
	
	std::map< std::string, std::set< std::string> > optdata;
	std::map< std::string,std::string> simpledata;
	
	std::map< std::string,std::set<opt_t> > comm2opt;
	
      };
      

      class ArgMap { 
	
      public:	

	static enum type_t { NONE , 
			     STRING , 
			     INT , 
			     FLOAT , 
			     STRING_VECTOR , 
			     INT_VECTOR , 
			     FLOAT_VECTOR } types ;
	
	ArgMap( int , char** );	
	
	bool has( const std::string & ) const;
	
	bool known( const std::string & ) const;

	std::string desc() const;

	std::vector<std::string> as_string_vector( const std::string & a ) const;

	std::string as_string( const std::string & a ) const;
	
	int as_int( const std::string & a ) const;

	std::vector<int> as_int_vector( const std::string & a ) const;

	double as_float( const std::string & a ) const;
	
	void reg( const std::string & s , const type_t & , const std::string & desc = "" );	

	std::string command() const { return command_str; } 
	
	std::string project_file() const { return project_str; } 
	
	void shortform( const std::string & a , const std::string & b );

	void attach( const std::string & command , const std::string & arg );
	
	std::string attached( const std::string & );
	
	std::string type( const std::string & arg )
	  {
	    std::map<std::string,type_t>::iterator i = known_type.find( arg );
	    if ( i == known_type.end() ) return ".";
	    if ( i->second == NONE ) return ".";
	    if ( i->second == STRING ) return "str";
	    if ( i->second == INT ) return "int";
	    if ( i->second == FLOAT ) return "float";
	    if ( i->second == STRING_VECTOR ) return "str-list";
	    if ( i->second == INT_VECTOR ) return "int-list";
	    if ( i->second == FLOAT_VECTOR ) return "float-list";
	  }

      private:
	
	
	// --flag    val1 val2 val3
	//   KEY    =     VALUE 
	
	std::map<std::string,std::vector<std::string> > data;
	std::map<std::string,type_t> known_type;
	std::map<std::string,std::string> known_desc;
	std::map<std::string,std::string> shortcuts;
	
	// attached args to commands
	std::map<std::string,std::set<std::string> > comm2arg;

	std::string command_str;
	std::string project_str;
      };
      
      
      class Commands {
	
      public:
	
	Commands & operator<<( const std::string & s )
	  {

	    // command|group|description|GRP|VCF|NOGENO
	    // aleays 3 |-delimited fields
	    // optionally, GRP meaning group-iteration
	    //             VCF meaning applicable to individual VCFs
	    //             NOGENO means that, by default, we do not need to look at genotype data
	    //             ARG:arg1,arg2
	    //             OPT:opt1=str-list,opt2=str   (we need format for opt)

	    std::vector<std::string> d = Helper::char_split( s, '|' );
	    
	    if ( d.size() < 3 ) 
	      Helper::halt( "internal error: malformed command spec:" + s );
	    
	    // first 3 standard arguments

	    // description
	    comm_desc[ d[0] ] = d[2];

	    // 1 or more groups (comma-delimited)
	    std::vector<std::string> grps = Helper::char_split( d[1] , ',' );
	    for (int g = 0; g < grps.size(); g++)
	      {
		check_defined( grps[g] );
		add_command( grps[g] , d[0] );
	      }	    
	    
	    // process optional arguments
	    for (int i=3; i<d.size(); i++)
	      {
		if ( d[i] == "GRP" ) comm_group_iteration.insert( d[0] );
		if ( d[i] == "VCF" ) comm_single_vcf.insert( d[0] );
		if ( d[i] == "NOGENO" ) comm_no_genotypes.insert( d[0] );
		if ( d[i].substr(0,4) == "ARG:" ) pargs->attach( d[0] , d[i].substr(4) );
		if ( d[i].substr(0,4) == "OPT:" ) popt->attach( d[0] , d[i].substr(4) );
	      }
	    
	    return *this;
	  }
	
	void attach( ArgMap * p ) { pargs = p; }

	void attach( Options * p ) { popt = p; }

	std::string group_description( const std::string & group ) const
	  {
	    std::map<std::string,std::string>::const_iterator i = group_desc.find( group );
	    return i == group_desc.end() ? "." : i->second;
	  }
	
	void new_group( const std::string & group , const std::string & desc )
	{
	  group_desc[ group ] = desc;
	}

	void check_defined( const std::string & g ) const
	{
	  if ( group_desc.find(g) == group_desc.end() ) 
	    Helper::halt( "internal error, command-group not defined: " + g );
	}

	void add_to_group( const std::string & group1 , const std::string & group2 )
	{
	  check_defined( group1 );
	  check_defined( group2 );

	  // 1-to-many relationship between groups 
	  group_group[ group1 ].push_back( group2 );
	  group_group_rmap[ group2 ] = group1 ;
	}
	
	void add_command( const std::string & g , const std::string & c )
	{
	  comm_group[ g ].push_back( c );
	  comm_group_rmap[ c ].insert( g );
	}

	bool known( const std::string & c ) const 
	{
	  return comm_desc.find(c) != comm_desc.end() ;
	}
	
	std::string description( const std::string & c ) 
	  {
	    return known(c) ? comm_desc[c] : "?" ;
	  }

	bool groups( const std::string & c ) 
	{
	  return comm_group_iteration.find(c) != comm_group_iteration.end();
	}

	bool single_VCF_mode( const std::string & c )
	{
	  return comm_single_vcf.find(c) != comm_single_vcf.end();
	}
	
	bool need_genotypes( const std::string & c ) const 
	{
	  return comm_no_genotypes.find(c) == comm_no_genotypes.end();
	}

	std::set<std::string> command_belongs_to() const
	  {
	    std::set<std::string> s;
	    return s;
	  }
	
	void display( const std::string & n ) const
	{

	  if ( group_group.find(n) != group_group.end() )
	    {
	      const std::vector<std::string> & s = group_group.find(n)->second;
	      std::vector<std::string>::const_iterator i = s.begin();
	      while ( i != s.end())
		{
		  plog << "GROUP" << "\t"
		       << *i << "\t"
		       << group_desc.find( *i )->second << "\n";
		  ++i;
		}
	    }
	  
	  if ( comm_group.find(n) != comm_group.end() )
	    {
	      const std::vector<std::string> & s = comm_group.find(n)->second;
	      std::vector<std::string>::const_iterator i = s.begin();
	      while ( i != s.end())
		{
		  plog << "COMM" << "\t"
		       << *i << "\t"
		       << comm_desc.find( *i )->second << "\n";
		  ++i;
		}
	    }

	  // if a command (rather than a group) the list all args/opts
	  
	  if ( known(n) ) 
	    {
	      plog << pargs->attached( n ) 
		   << popt->attached( n );
	    }
	  
	}


      private:
	
	std::map<std::string, std::string> comm_desc;
	std::set<std::string> comm_group_iteration;
	std::set<std::string> comm_single_vcf;
	std::set<std::string> comm_no_genotypes;

	std::map<std::string,std::vector<std::string> > comm_group;
	std::map<std::string,std::set<std::string> > comm_group_rmap;
	std::map<std::string,std::string> group_desc;
	std::map<std::string,std::vector<std::string> > group_group;
	std::map<std::string,std::string> group_group_rmap;

	ArgMap * pargs;
	Options * popt;
	
      };


      template<class T>
	std::vector<T> n_arguments( Pseq::Util::ArgMap & args , const std::string & a , int n = 0 )
	{
	  if ( ! args.has( a ) )
	    Helper::halt("no --" + a + " specified");
	  
	  std::vector<std::string> v = args.as_string_vector(a);
	  
	  // expand out comma-separated values to a new list of type T
	  std::vector<std::string> v2;
	  for (int i=0; i<v.size(); i++)
	    {
	      std::vector<std::string> t = Helper::parse( v[i] , "," );
	      for (int j=0; j<t.size(); j++)
		v2.push_back(t[j]);
	    }

	  std::vector<T> fv; // final values
	  for (int i=0; i<v2.size(); i++)
	    {
	      T a;
	      try { a = Helper::lexical_cast<T>( v2[i] ); }
	      catch ( std::exception& e) { Helper::halt("trouble parsing " + v2[i] ); }
	      fv.push_back(a);	      
	    }
	  if ( n > 0 && fv.size() != n ) 
	    Helper::halt("need " + Helper::int2str(n) + " argument(s) for --" + a );
	  return fv;
	}
      
      template<class T> 
 	T single_argument( Pseq::Util::ArgMap & args , const std::string & a ) 
 	{ 
	  return n_arguments<T>(args,a,1)[0];
 	} 



    }
  
  namespace SeqDB
    {
      
      bool load_FASTA( const std::string & );
      
      bool summary();
      
      bool lookup( Region & );
      
      bool load_transcripts( std::string );

      bool loc_stats( const std::string & , const std::string & );

      bool loc_translate( const std::string & );
    }

}


#endif
