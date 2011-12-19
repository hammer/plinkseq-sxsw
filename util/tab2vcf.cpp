#include <iostream>
#include <string>
#include <map>

#include "char_tok.h"

std::string opt_quote( const std::string & s ) 
{
  // if str contains = or ; or some such, quote
  if ( s.find( "=" ) != std::string::npos || s.find(";") != std::string::npos || s.find(":") != std::string::npos ) 
    return "\"" + s + "\"";
  return s;
}

std::string replace_space( const std::string & s , char c )
{
  std::string x = s;
  const int ls = x.size();
  
  bool inq = false;      
  for (int k = 0; k < ls; k++)
    {
      if ( x[k] == '"' ) inq = ! inq;
      else if ( x[k] == ' ' && ! inq ) x[k] = c;
    }
  return x;
}


int main( int argc , char ** argv ) 
{
  
  std::map<std::string,std::string> args;
  args["ID"] = "ID";
  args["REF"] = "REF";
  args["ALT"] = "ALT";
  args["QUAL"] = "QUAL";
  args["FILTER"] = "FILTER";
  
  std::map<std::string,std::string> types;
  std::map<std::string,std::string> number;
  std::map<std::string,std::string> desc;

  for (int i=1;i<argc;i++)
    {
      int t;
      char_tok tok( argv[i] , &t , '=' );
      if ( t == 2 ) args[ tok(0) ] = tok(1);
      else if ( t == 1 ) 
	{
	  int t2;
	  char_tok tok2( argv[i] , &t2 , ':' );
	  if ( t2 == 4 ) 
	    {
	      types[ tok2(0) ] = tok2(1);
	      number[ tok2(0) ] = tok2(2);
	      desc[ tok2(0) ] = tok2(3);
	    }
	}
    }

  // take a file with HEADERs from STDIN;
  // write as a VCF to STDOUT
  
  // Get header LINE
  
  std::string hdr;
  std::getline( std::cin , hdr );
  if ( hdr == "" ) exit(1);

  // remove a leading '#'
  if ( hdr.substr(0,1) == "#" ) 
    hdr = hdr.substr(1);

  // get fields
  int n;
  char_tok tok( hdr , &n , '\t' );
  
  // find positional fields
  // either POS, or CHR and BP (or BP1 and BP2)

  int pos = -1 , chr = -1 , bp1 = -1 , bp2 = -1;

  // other special fields
  int qual = -1 , id = -1 , ref = -1 , alt = -1 , filter = -1;

  
  std::vector<std::string> field(n);
  std::vector<int> info;

  for (int i=0;i<n;i++)
    {

      bool special_field = false;

      field[i] = tok(i);

      if      ( field[i] == "POS" ) { special_field = true; pos = i; } 
      else if ( field[i] == "VAR" ) { special_field = true; pos = i; }
      
      else if ( field[i] == "CHR" ) { special_field = true; chr = i; }
      else if ( field[i] == "CHROM" ) { special_field = true; chr = i; }
      
      else if ( field[i] == "BP" ) { special_field = true; bp1 = i; }
      else if ( field[i] == "BP1" ) { special_field = true; bp1 = i; }
      
      else if ( field[i] == "BP2" ) { special_field = true; bp2 = i; }  
      
      else if ( field[i] == "POS1" ) { special_field = true; bp1 = i;}
      else if ( field[i] == "POS2" ) { special_field = true;  bp2 = i; }
      
            
      std::map<std::string,std::string>::iterator ii = args.begin();
      while ( ii != args.end() )
	{
	  if ( ii->second == field[i] ) 
	    {
	      if ( ii->first == "QUAL" || ii->first == "qual" )
		{ special_field = true; qual = i; }
	      else if ( ii->first == "FILTER" || ii->first == "filter" )
		{ special_field = true; filter = i; }
	      else if ( ii->first == "ID" || ii->first == "id" ) 
		{ special_field = true; id = i; }
	      else if ( ii->first == "REF" || ii->first == "ref" ) 
		{ special_field = true; ref = i; }
	      else if ( ii->first == "ALT" || ii->first == "alt" ) 
		{ special_field = true; alt = i; }
	      else if ( ii->first == "SKIP" || ii->first == "skip" )
		special_field = true;		
	    }
	  ++ii;
	}

      // normal INFO field?
      
      if ( ! special_field ) 
	info.push_back(i);  

    } // next header entry

  
  // does this contain valid positional information?
  
  if ( pos == -1 && ( chr == -1 || bp1 == -1 ) ) exit(1);
  
  bool use_chr_bp = chr != -1 && bp1 != -1;
  bool use_bp2 = use_chr_bp && bp2 != -1;

  const int s = info.size();

  // write VCF header
  
  std::cout << "##fileformat=VCFv4.1\n"
	    << "##source=tab2vcf\n";
  
  for (int i=0;i<s;i++)
    {
      if ( types.find( field[ info[i] ] ) != types.end() ) 
	{
	  std::cout << "##INFO=<ID=" << field[info[i]] << ",Number=" << number[field[info[i]]]
		    << ",Type=" << types[field[info[i]]] << ",Description=\"" << desc[field[info[i]]] 
		    << "\">\n";
	}
      else 
	std::cout << "##INFO=<ID=" << field[info[i]] << ",Number=1,Type=Float,Description=\"n/a\">\n";
    }

  // second header
  std::cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

  // information
  
  while ( ! std::cin.eof() )
    {
      std::string l;
      int m;
      
      std::getline( std::cin , l );
      
      // replace spaces with ';' unless they are quoted      
      
      char_tok tok( l , &m , '\t' );
      
      if ( l == "" ) continue;

      if ( m != n ) 
	{
	  std::cout << "unequal number of fields for line\n" << l 
		    << "\n(expecting " << n << " found " << m << ")\n";
	  exit(1);
	}

      
      if ( use_chr_bp ) 
	{
	  // chromosome
	  std::cout << tok(chr) << "\t";
	  
	  // position
	  std::cout << tok(bp1) ;
	  if ( bp2 != -1 ) std::cout << ".." << tok(bp2);
	  std::cout << "\t";

	}
      else
	{
	  std::string pstr = tok( pos );

	  // expect 'chr:1234'
	  // expect 'chr:1234..4596'
	  // or     'chr:1234-5674'
	  int k;
	  char_tok stok( pstr , &k , ':' );
	  if ( k == 0 ) std::cout << ".\t.\t";
	  else if ( k == 1 ) std::cout << stok(0) << "\t.\t";
	  else if ( k >= 2 ) 
	    {
	      std::cout << stok(0) << "\t"
			<< stok(1) << "\t";
	    }	  
	}
      
      // id

      if ( id == -1 ) std::cout << ".\t";
      else std::cout << tok(id) << "\t";

      // ref
      if ( ref == -1 ) std::cout << ".\t";
      else std::cout << tok(ref) << "\t";
      
      // alt
      if ( alt == -1 ) std::cout << ".\t";
      else std::cout << tok(alt) << "\t";

      // qual
      if ( qual == -1 ) std::cout << ".\t";
      else std::cout << tok(qual) << "\t";
      
      // filter
      if ( filter == -1 ) std::cout << ".\t";
      else std::cout << replace_space( tok(filter) , ';' ) << "\t";
      

      // info
      for (int i=0;i<s;i++)
	{
	  if ( i ) std::cout << ";";
	  std::cout << field[info[i]] << "=" << opt_quote( replace_space( tok(info[i]) , ',' ) );
	}
      std::cout << "\n";
      
    }

}
