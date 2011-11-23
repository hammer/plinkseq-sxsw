
#include "helper.h"
#include "func.h"

extern Pseq::Util::Commands pcomm;

enum Pseq::Util::ArgMap::type_t types;

void Pseq::Util::ArgMap::shortform( const std::string & sht , const std::string & lng ) 
{
  shortcuts[ sht ] = lng;
}

Pseq::Util::ArgMap::ArgMap( int n , char ** argv )
{
  
  reg( "help" , NONE , "produce help message" );
  reg( "vcf" , STRING_VECTOR , "VCF file locations" );
  reg( "bcf" , STRING_VECTOR , "BCF file locations" );
  reg( "resources" , STRING , "central resource folder" );
  reg( "scratch" , STRING , "scratch folder" );
  reg( "metameta" , STRING , "meta-information meta-information" );

  reg( "history" , STRING_VECTOR , "use a .history file with GSEQ" );
  
  reg( "vardb", STRING, "variant database location" );
  reg( "inddb", STRING, "individual database location" );
  reg( "refdb", STRING, "reference database location" );
  reg( "seqdb", STRING, "sequence database location" );
  reg( "segdb", STRING, "segent database location" );
  reg( "locdb", STRING, "locus database location" );
  reg( "netdb", STRING, "network database location" );  
  reg( "ibddb", STRING, "IBD segment database location" );  
  
  reg( "file" , STRING_VECTOR , "generic input file(s)" );
  reg( "group" , STRING_VECTOR , "generic group label(s)" );
  reg( "region" , STRING_VECTOR , "region(s) ");
  reg( "alias" , STRING_VECTOR , "locus alias group(s)" );
  reg( "name" , STRING_VECTOR , "generic name(s) variable" );
  reg( "key" , STRING , "key of key-value pair" );
  reg( "value" , STRING_VECTOR , "value(s) of key-value pair" );
  reg( "type", STRING , "type of project entry");
  reg( "id" , INT_VECTOR , "generic numeric IDs" );
  reg( "options" , STRING_VECTOR, "context-specific options\n");
  reg( "output", STRING, "output folder\n" );
  reg( "whitespace", NONE , "allow whitespace delimited input" );

  reg( "new-project" , STRING , "new project specification filename" );
  reg( "new-vardb" , STRING , "new VARDB name, for write-vardb" );
  
  reg( "debug", NONE , "set debug mode");
  reg( "silent", NONE , "set silent mode");
  reg( "ignore-warnings" , NONE , "turn off warnings");
  reg( "early-warnings" , NONE , "display warning as soon as happens");

  reg( "out-file", STRING , "set main output file");
  reg( "debug-file", STRING , "debug file name");
  reg( "prolix-file", STRING, "prolix output filename");
  reg( "long" , NONE , "set long output mode");
  reg( "long-header" , NONE , "set header/long output mode");
  
  reg( "mask", STRING_VECTOR, "mask specification");
  reg( "include", STRING, "filter specification");
  reg( "exclude", STRING , "filter specification");
  reg( "gene", STRING_VECTOR, "gene-group gene1 gene2 ...");
  reg( "em", FLOAT , "EM calculation of P(geno|data) from GL or PL");

  reg( "assume-ref" , NONE , "convert null genotypes to reference homozygote");

  reg( "hide", STRING_VECTOR,"hide specific meta-fields");
  reg( "show", STRING_VECTOR,"show specific meta-fields" );
  reg( "force-consensus" , NONE , "set all tags to consensus" );

  reg( "vmeta", NONE, "show variant meta-information" );
  reg( "samples", NONE, "show each specific sample variant");
  reg( "verbose", NONE, "verbose output");
  reg( "geno", NONE, "show genotypes");
  reg( "gmeta", NONE, "show genotype meta-information" );
  reg( "transpose", NONE, "transposed g-view output");
  reg( "simple" , NONE , "simple variant format, POS RET ALT" );
  
  reg( "variant", STRING , "show specific variant (v-view)");
  reg( "indiv", STRING_VECTOR , "specify individual(s)");
    
  reg( "annot" , STRING_VECTOR , "transcript(s) group for annotation" );
	
  reg( "phenotype" , STRING_VECTOR, "phenotype specification");
  reg( "make-phenotype" , STRING, "dichotomise factor");
  reg( "strata" , STRING,"stratifier variable");
  reg( "covar" , STRING_VECTOR , "covariate(s)");

  reg( "perm" , INT, "number of permutations");
  reg( "aperm" , INT_VECTOR , "adaptive perm min, max");

  reg( "weights" , STRING , "name of variant weights tag");
  
  
  //
  // Also, add in major descriptors for options
  //

  Pseq::Util::Options::reg_option( "opt1" , "str" , "this is my description" );


  //
  // Register some short-cuts
  //
  
  shortform( "-o" , "--options" );
  shortform( "-m" , "--mask" );
  shortform( "-f" , "--file" );
  shortform( "-g" , "--group" );
  shortform( "-p" , "--phenotype" );
  shortform( "-h" , "--help" );  
  shortform( "-n" , "--name" );
  shortform( "-w" , "--weights" );

  shortform( "help" , "--help" );

  shortform( "--out" , "--output" );
  shortform( "--pheno" , "--phenotype" );
  shortform( "--phe" , "--phenotype" );
  shortform( "--weight" , "--weights" );
  
  
  needs_help = n > 1 && ( std::string( argv[1] ) == "--help" 
			  || std::string( argv[1] ) == "help" 
			  || std::string( argv[1] ) == "-h" 
			  || std::string( argv[1] ) == "masks" 
			  || std::string( argv[1] ) == "mask" ) ;


  if ( needs_help ) // will be called later
    {      
      if ( n == 2 ) help_str = std::string(argv[1]).substr(0,1) == "m" ? "mask" : "." ;
      else if ( n > 2 ) help_str = argv[2];
      return;
    }

  
  // arg 0 is filename
  // position 1 should be project
  // position 2 should be command
  // then commands and flags.
  
  if ( n < 2  ) 
    Helper::halt("no project or command specified; try 'pseq help'");
  else if ( n < 3 ) 
    Helper::halt("no command specified; try 'pseq help'");

  project_str = argv[1];
  command_str = argv[2];
  
  for (int i=3 ; i < n ; i++ )
    {
      std::string s = argv[i];
      
      // Swap in long-form from short-form? e.g. -o to --options
      if ( shortcuts.find( s ) != shortcuts.end() ) 
	s = shortcuts[s];
      
      if ( s.substr(0,2) != "--" ) Helper::halt("unknown option: " + s ); 
      s =  s.substr(2); 
      if ( ! known(s) ) Helper::halt("unknown option: " + s );

      std::vector<std::string> a;

      while ( 1 ) 
	{
	  ++i;
	  if ( i == n ) break;
	  
	  std::string b = argv[i];

	  if ( shortcuts.find( b ) != shortcuts.end() ) 
	    b = shortcuts[b];
	  
	  if ( b.substr(0,2) == "--" ) 
	    {
	      --i; 
	      break;
	    }
	  else
	    a.push_back( b );
	}

      // store, or add in as extras
      if ( data.find(s) == data.end() )
	data[s] = a;      
      else
	Helper::append( data[s] , a );      

    }

}


bool Pseq::Util::ArgMap::known( const std::string & s ) const
{
  return known_type.find(s) != known_type.end();
}

void Pseq::Util::ArgMap::reg( const std::string & s , const type_t & t , const std::string & desc )
{
  known_type[s] = t;
  known_desc[s] = desc;
}

bool Pseq::Util::ArgMap::has( const std::string & s ) const
{
  return data.find(s) != data.end();
}

std::vector<std::string> Pseq::Util::ArgMap::as_string_vector( const std::string & a )  const
{  
  if ( ! known(a) ) Helper::halt("argument " + a + " not found" );
  return data.find(a)->second;
}

std::string Pseq::Util::ArgMap::as_string( const std::string & a )  const 
{
  if ( ! known(a) ) Helper::halt("argument --" + a + " not found" );
  if ( data.find(a)->second.size() != 1 ) Helper::halt("expecting 1 value for --" + a );
  return data.find(a)->second[0];
}

int Pseq::Util::ArgMap::as_int( const std::string & a )  const
{
  if ( ! known(a) ) Helper::halt("argument " + a + " not found" );
  if ( known_type.find(a)->second != INT ) Helper::halt("incorrect return type for " + a );
  int i;
  if ( ! Helper::str2int( data.find(a)->second[0] , i ) ) 
    Helper::halt("non-integer argument for --" + a );
  return i;
}

std::vector<int> Pseq::Util::ArgMap::as_int_vector( const std::string & a ) const
{
  if ( ! known(a) ) Helper::halt("argument " + a + " not found" );
  if ( known_type.find(a)->second != INT_VECTOR ) Helper::halt("incorrect return type for " + a );
  std::vector<int> d;
  for (int i=0; i<data.find(a)->second.size(); i++)
    {
      int x;
      if ( ! Helper::str2int( data.find(a)->second[i] , x ) ) 
	Helper::halt("non-integer argument for --" + a );
      d.push_back(x);
    }
  return d;
}

double Pseq::Util::ArgMap::as_float( const std::string & a ) const
{
  if ( ! known(a) ) Helper::halt("argument " + a + " not found" );
  if ( known_type.find(a)->second != FLOAT ) Helper::halt("incorrect return type for " + a );
  double d;
  if ( ! Helper::str2dbl( data.find(a)->second[0] , d ) ) 
    Helper::halt("non-integer argument for --" + a );
  return d;
}

bool Pseq::Util::ArgMap::help() const
{
  if ( ! needs_help ) return false;
  std::cout << desc( help_str ) ;
  return true;
}


std::string Pseq::Util::ArgMap::desc( const std::string & c ) const 
{

  std::stringstream ss;
  
  // --help         { list all commands } 
  // --help v-view  { list all arguments for v-view }
  // --help masks   { list all mask options }
  

 
  
  if ( c == "." )
    {

      ss << "\nusage:\tpseq {project-file|VCF} {command} {--options}\n\n";
      
      std::vector<std::string> groups = pcomm.groups( );
      
      ss << "\tCommand groups\n";
      ss << "\t---------------------------------------------------------\n";
      for ( int g = 0 ; g < groups.size() ; g++ )
	{	  
	  std::vector<std::string> cs = pcomm.commands( groups[g] );
	  if ( cs.size() ) // do not show hidden/empty ones
	    ss << "\t" << groups[g] << "\t\t\t" << pcomm.group_description( groups[g] ) << "\n";
	}      
      
      ss << "\n\tMask groups\n";
      ss << "\t---------------------------------------------------------\n";
      ss << Mask::list_groups();
      
      ss << "\n"
	 << "pseq help all\n"
	 << "pseq help {group}\n"
	 << "pseq help {command}\n";

   }
  else if ( c == "all" )
    {
      // verbose help mode
      ss << "\nusage:\tpseq {project-file|VCF} {command} {--options}\n\n";
      
      ss << "Commands\n";
      ss << "---------------------------------------------------------\n\n";
      
      std::vector<std::string> cs = pcomm.all_commands();

      for (int c = 0 ; c < cs.size() ; c++)
	{
	  ss << "\t" << pcomm.command_description( cs[c] );
	  ss << pcomm.command_description( cs[c] , true ) << "\n";
	}
      
      
      ss << "\n\nMask groups\n";
      ss << "---------------------------------------------------------\n\n";

      ss << Mask::list_groups( true );
      
      ss << "\n"
	 << "pseq help {group}\n"
	 << "pseq help {command}\n";      

    }
  else if ( pcomm.has_group( c ) )
    {      
      ss << "\n\t" << c << " : " << pcomm.group_description( c ) << "\n";
      ss << "\t---------------------------------------------------------\n";
      std::vector<std::string> cs = pcomm.commands( c );      
      for (int c = 0 ; c < cs.size() ; c++)
	ss << "\t" << pcomm.command_description( cs[c] );
      
    }
  else if ( pcomm.known( c ) )
    {
      std::vector<std::string> str = Helper::char_split( pcomm.command_description( c ) , '\t' , true );
      ss << "\n\t" << c << " : " << str[str.size()-1] ;
      ss << "\t---------------------------------------------------------\n";
      ss << pcomm.command_description( c , true );
    }
  else 
    {
      // is this a Mask group?
      
      std::string str = Mask::list_masks( c );
      if ( str == "" ) 	
	ss << "[ " << c << " ] not recognised as a command, command-group or mask-group\n";
      else
	{
	  // v. silly way to get desc... quicker for now.
	  std::string s = Mask::list_groups();	  
	  std::vector<std::string> s1 = Helper::char_split( s , '\n' );
	  std::string desc;
	  for ( int i = 0 ; i < s1.size() ; i++)
	    {
	      std::vector<std::string> s2 = Helper::char_split( s1[i] , '\t' , false );
	      if ( s2.size() == 2 && s2[0] == c ) { desc = s2[1] ; break; }
	    }
	  ss << "\n\t" << c << " : " << desc << "\n";
	  ss << "\t---------------------------------------------------------\n";
	  ss << str;
	}
    }
  ss << "\n";
  return ss.str();
}


void Pseq::Util::ArgMap::attach( const std::string & command , const std::string & arg )
{
  // arg will be comma-delimited list, with types
  // ref=str-list,loc=str-list
  std::vector<std::string> opts = Helper::char_split( arg , ',' );
  for (int i=0; i<opts.size(); i++)
    comm2arg[ command ].insert( opts[i] );
}


std::string Pseq::Util::ArgMap::attached( const std::string & command )
{
  std::map<std::string,std::set<std::string> >::iterator i = comm2arg.find( command );  
  if ( i == comm2arg.end() ) return "";
  std::string s = "";
  std::set<std::string>::iterator j = i->second.begin();
  while ( j != i->second.end() )
    {
      s += "ARG\t" + i->first + "\t" + *j + "\t" + type( *j ) + "\t" + arg_desc( *j ) + "\n";
      ++j;
    }      
  return s;
}
