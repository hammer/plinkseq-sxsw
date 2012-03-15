#ifndef __PSEQ_UTIL_H__
#define __PSEQ_UTIL_H__

#include "pseq.h"

#include <string>


namespace Pseq 
{
  
  namespace Util
    {
	
	class Commands;
	
	void populate_commands( Pseq::Util::Commands & );	


	class Options { 
	    
	public:	
	    
	    // --arg keyword
	    
	    // keywords can have their own type too, e.g. INT
	    // keywords start with = and are ',' delimited
	    
	    //   arg  KEYWORD_LIST  is_keyword = F
	    //   val1 NONE          is_keyword = T
	    //   val2 INT_VECTOR    is_keyword = T
	    
	    //  --arg val val2=1,2,3
	    
	    // (this is similar to Mask specification, but as the Mask is
	    // more complex and not limited to PSEQ, we let the Mask class
	    // take care of that itself)
	    
	    static enum type_t { NONE , 
				 STRING , 
				 KEYWORD ,
				 INT , 
				 INT_RANGE ,
				 FLOAT , 
				 FLOAT_RANGE ,
				 STRING_VECTOR , 
				 INT_VECTOR , 
				 INT_RANGE_VECTOR ,
				 FLOAT_VECTOR , 
				 FLOAT_RANGE_VECTOR ,
				 VARIABLE_TYPE } types;
	    
	    
	    void load( int , char** );	
	    
	    bool has( const std::string & ) const;
	    
	    bool known( const std::string & ) const;
	    
	    
	    // for KEYWORD  style args
	    
	    bool has( const std::string & , const std::string & ) const;

	    bool has( const std::string & , const std::string & , const std::string & ) const;

	    
	    bool known( const std::string & , const std::string & ) const;
	    
	    std::set<std::string> keywords( const std::string & ) const;
	    
	    bool has_args( const std::string & , const std::string & ) const;
	    
	    std::string key_type( const std::string & a , const std::string & b ) const;
	    
	    std::string key_description( const std::string & a , const std::string & b ) const;


	    // misc
	    
	    bool help() const;
	    std::string desc( const std::string & ) const;
	    std::string command() const { return command_str; } 	
	    std::string project_file() const { return project_str; } 
	    void shortform( const std::string & a , const std::string & b );

	    
	    // extraction/converstion functions
	    std::vector<std::string> as_string_vector( const std::string & a ) const;
	    std::string as_string( const std::string & a ) const;	
	    int as_int( const std::string & a ) const;
	    std::vector<int> as_int_vector( const std::string & a ) const;
	    double as_float( const std::string & a ) const;	
	    std::vector<double> as_float_vector( const std::string & a ) const;

	    // as above, but for keywords ( --arg value=key value=key) -- still space delimited
	    
	    std::vector<std::string> as_string_vector( const std::string & , const std::string & a ) const;
	    std::string as_string( const std::string & , const std::string & a ) const;	
	    int as_int( const std::string & , const std::string & a ) const;
	    std::vector<int> as_int_vector( const std::string & , const std::string & a ) const;
	    double as_float( const std::string & , const std::string & a ) const;	
	    
	    
	    
	    // register known args
	    
	    void reg( const std::string & s , const type_t & , const std::string & desc = "" );
	    void keyword( const std::string & s , const std::string & k , const type_t & , const std::string & desc = "" );
	    
	    void attach( const std::string & command , const std::string & arg );	
	    std::string attached( const std::string & );
	    
	    
	    // helpers for help function
	    
	    std::string arg_description( const std::string & arg ) const;
	    
	    std::string type( const std::string & arg ) const;


	    // custom extraction functions 

	    std::set<std::string> get_set( const std::string & k ) const;

	    std::set<std::string> get_set( const std::string & a , const std::string & b ) const;
	    	    
	    std::string comma_string( const std::string & k ) const;
	
	    
	
	private:
	    
	    // primary data-stores
	    
	    std::map<std::string, std::vector<std::string> > data;
	    
	    // for each arg, store actual given keywords upto the '=';
	    std::map<std::string, std::map<std::string,std::vector<std::string> > > data_kw;
	    
	    
	    // some meta-information on each argument, from reg()
	    
	    std::map<std::string,type_t>                    arg_type;
	    std::map<std::string,std::string>               arg_desc;
	    
	    // map arg --> keywords, e.g. --arg keyword keyword=value
	    
	    std::map<std::string,std::set<std::string> >              keyword_map;
	    std::map<std::string,std::map<std::string,type_t> >       keyword_type;
	    std::map<std::string,std::map<std::string,std::string> >  keyword_desc;
	    
	    // e.g. allow --mask to be -m 
	    
	    std::map<std::string,std::string>               shortcuts;
	    
	    
	    // map Commands to their arguments, and arg/key pairs
	    
	    std::map<std::string,std::set<std::string> >    comm2arg;
	    std::map<std::string,std::map<std::string,std::set<std::string> > >   comm2key;
	    
	    
	    // the first two entries always expected to be project and command
	    // (unless 'help' mode)
	    
	    std::string command_str;
	    
	    std::string project_str;
	    
	    
	    // 'help' mode 
	    
	    bool needs_help;
	    
	    std::string help_str;
	    
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
		    //             OPT:opt1
		    
		    // formats for options are registered separately
		    
		    std::vector<std::string> d = Helper::char_split( s, '|' );
		    
		    if ( d.size() < 3 ) 
			Helper::halt( "internal error: malformed command spec:" + s );
		    
		    // first 3 standard arguments

		    // command under development? 
		    
		    if ( d[0][0] == '*' ) 
		      {
			d[0] = d[0].substr(1);
			devel.insert(d[0]);
		      }

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
		    }
		    
		    return *this;
		}
	    
	    void attach( Options * p ) { pargs = p; }
	    
	    std::string group_description( const std::string & group ) const
		{
		    std::map<std::string,std::string>::const_iterator i = group_desc.find( group );
		    return i == group_desc.end() ? "." : i->second;
		}
	    
	    bool has_group( const std::string & g ) 
		{
		    return group_desc.find( g ) != group_desc.end();
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

	    bool stable( const std::string & c ) const
	    {
	      return devel.find( c ) == devel.end();
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
			plog << pargs->attached( n );		   
		    }
		    
		}
	    
	    
	    // simple list of command groups
	    
	    std::vector<std::string> groups( const std::string & g = "root" ) const 
		{
		    std::vector<std::string> s;	    
		    std::map<std::string,std::vector<std::string> >::const_iterator i = group_group.find(g);
		    if ( i == group_group.end() ) return s;
		    return i->second;
		}
	    
	    
	    // given a group, give a list of commands
	    std::vector<std::string> commands( const std::string g ) const 
		{	    
		    if ( comm_group.find( g ) != comm_group.end() )
			return comm_group.find( g )->second;
		    std::vector<std::string> s;
		    return s;
		}
	    
	    // list all commands (only once, alphabetically)
	    std::vector<std::string> all_commands() const
		{
		    std::vector<std::string> s;	    
		    std::map<std::string,std::set<std::string> >::const_iterator i = comm_group_rmap.begin();
		    while ( i != comm_group_rmap.end() )
		    {
			s.push_back( i->first );
			++i;
		    }
		    return s;
		}
	    
	    // given a command, give the description and options       
	    
	    std::string command_description( const std::string & c , bool show_args = false ) const; 
	    
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

	    // command under development, should not be shown in help
	    // registered as "*command" when entering
	    std::set<std::string> devel;
	    
	    Options * pargs;
	    
	};
	
	
	template<class T>
	    std::vector<T> n_arguments( Pseq::Util::Options & args , const std::string & a , int n = 0 )
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
	    T single_argument( Pseq::Util::Options & args , const std::string & a ) 
	    { 
		return n_arguments<T>(args,a,1)[0];
	    } 
	
    }
  
}


#endif
