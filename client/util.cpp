
#include "helper.h"
#include "func.h"

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

  reg("vardb", STRING, "variant database location" );
  reg("inddb", STRING, "individual database location" );
  reg("refdb", STRING, "reference database location" );
  reg("seqdb", STRING, "sequence database location" );
  reg("segdb", STRING, "segent database location" );
  reg("locdb", STRING, "locus database location" );
  
  
  reg("file" , STRING_VECTOR , "generic input file(s)" );
  reg("group" , STRING_VECTOR , "generic group label(s)" );
  reg("region" , STRING_VECTOR , "region(s) ");
  reg("alias" , STRING_VECTOR , "locus alias group(s)" );
  reg("name" , STRING_VECTOR , "generic name(s) variable" );
  reg("type", STRING , "type of project entry");
  reg("id" , INT_VECTOR , "generic numeric IDs" );
  reg("options" , STRING_VECTOR, "context-specific options\n");
  reg("output", STRING, "output folder\n" );
  reg("whitespace", NONE , "allow whitespace delimited input" );

  reg("new-project" , STRING , "new project specification filename" );
  reg("new-vardb" , STRING , "new VARDB name, for write-vardb" );
  
  reg("debug", NONE , "set debug mode");
  reg("silent", NONE , "set silent mode");
  reg("ignore-warnings" , NONE , "turn off warnings");

  reg("out-file", STRING , "set main output file");
  reg("debug-file", STRING , "debug file name");
  reg("prolix-file", STRING, "prolix output filename");
  reg("long" , NONE , "set long output mode");
  reg("long-header" , NONE , "set header/long output mode");
  
  reg("mask", STRING_VECTOR, "mask specification");
  reg("include", STRING, "filter specification");
  reg("exclude", STRING , "filter specification");
  reg("gene", STRING_VECTOR, "gene-group gene1 gene2 ...");
  reg("em", FLOAT , "EM calculation of P(geno|data) from GL or PL");

  reg("assume-ref" , NONE , "convert null genotypes to reference homozygote");

  reg("hide", STRING_VECTOR,"hide specific meta-fields");
  reg("show", STRING_VECTOR,"show specific meta-fields" );

  reg("vmeta", NONE, "show variant meta-information" );
  reg("samples", NONE, "show specific sample variants in v-view");
  reg("verbose", NONE, "verbose output");
  reg("geno", NONE, "show genotypes (g-view)");
  reg("gmeta", NONE, "show genotype meta-information" );
  reg("transpose", NONE, "transposed g-view output");
      
  reg("variant", STRING , "show specific variant (v-view)");
  reg("indiv", STRING_VECTOR , "specify individual(s)");
    
  reg("annot" , STRING_VECTOR , "transcript(s) group for annotation" );
	
  reg("phenotype" , STRING_VECTOR, "phenotype specification");
  reg("make-phenotype" , STRING, "dichotomise factor");
  reg("strata" , STRING,"stratifier variable");
  reg("covar" , STRING_VECTOR , "covariate(s)");

  reg("perm" , INT, "number of permutations");
  reg("aperm" , INT_VECTOR , "adaptive perm min, max");
  
  //
  // Register some short-cuts
  //
  
  shortform( "-o" , "--options" );
  shortform( "-m" , "--mask" );
  shortform( "-f" , "--file" );
  shortform( "-g" , "--group" );
  shortform( "-p" , "--phenotype" );
  shortform( "-h" , "--help" );
  
  
  if ( n == 2 && std::string(argv[1]) == "--help"  )
    {
      std::cout << desc() << "\n";
      exit(0);
    }
  
  // arg 0 is filename
  // position 1 should be project
  // position 2 should be command
  // then commands and flags.
  
  if ( n < 2  ) 
    Helper::halt("no project or command specified");
  else if ( n < 3 ) 
    Helper::halt("no command specified");

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


std::string Pseq::Util::ArgMap::desc() const 
{
  return "usage: pseq {project-file} {command} {--options}\n";
}
