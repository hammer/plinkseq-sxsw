#include "das.h"
#include "socks.h"
#include "plinkseq.h"

#include <iostream>
#include <getopt.h>

PDAS_param cmdargs( int argc , char ** argv );

bool f( const char * t , const int n , void * data , std::stringstream & s ) ;

//
// Initiate
//

int main( int argc , char ** argv )
{
  
  PDAS_param param = cmdargs( argc, argv );
  
  std::cerr << "Starting PLINK/SEQ DAS server, version 0.01\n";
 
  if ( param.project == "" ) 
    Helper::halt( "no project specified" );
  
  if ( ! Helper::fileExists( param.project ) ) 
    Helper::halt( "project file " + param.project + " does not exist" );
      
  std::cerr << " project   " << param.project << "\n"
	    << " port #    " << param.portn   << "\n";

  
  //
  // Attach PLINK/Seq project 
  //

  GStore g;  

  g.set_project( param.project );



  //
  // Open socket on port, pointing to processing function and project
  //
  
  DAS das( &g , param.project );
  
  ServerSocket sock( param.portn , f , (void*)&das ) ; 
  
  //
  // Will never reach here...
  //
  
  exit(0);

} 




PDAS_param cmdargs( int argc , char ** argv )
{
  
  PDAS_param param;

  int c = -1;

  opterr = 0;
  
  // -d plinkseq data
  // -p port 
  
  while ( (c = getopt (argc, argv, "d:p:") ) != -1) {
    
    switch (c)
      {
      
      case 'd':
	param.project = optarg;
	break;
	
      case 'p':
	param.portn = atoi( optarg );
	break;
	
      case '?':
	if ( optopt == 'd' ||
	     optopt == 'p' 
	     )
	  fprintf (stderr, "Option -%c requires an argument.\n", optopt);
	else if (isprint (optopt))
	  fprintf (stderr, "Unknown option `-%c'.\n", optopt);
	else
	  fprintf (stderr,
		   "Unknown option character `\\x%x'.\n",
		   optopt);
	return param;

      default:
	abort ();
      }
  }
  
  return param;
}



bool f( const char * t , const int n , void * data , std::stringstream & ss ) 
{
  
  DAS * pdas = (DAS*)data;

  // expecting a valid HTTP/1.1 header 

  // add query to cahce
  
  pdas->add( t );

  

  // more input left? 
  if ( n == 255 ) return false; 

  std::cerr << "received request :\n" << n << "\n" << pdas->query() << "\n";

  pdas->attach_response_stream( & ss );


  //
  // parse DAS request 
  //

  int retc = pdas->parse_request(); 


  //
  // does requested DSN match the project name? 
  //

  if ( ! pdas->check_datasource() ) 
    {
      retc = 401; // bad data source
    }
  
  
  //
  // error out?
  //

  if ( retc != 200 )
    {
      std::cerr << "retc = " << retc << "\n";
      pdas->header( retc );
      return true; 
    }
  
  
  //
  // handle specific commands
  //
  
  if ( pdas->command( "sources" ) )
    {
      retc = pdas->sources();
    }
  
  else if ( pdas->command( "features" ) )
    {
      retc = pdas->features();
    }
  
  
  //
   // write header...
   //

   pdas->header( retc );


   //
   // ...and main body of response
   //

   pdas->write();

   return true;

 }



 //
 // -------  DAS object --------
 //


DAS::DAS( GStore * p , const std::string & n ) 
{

  // keep track of this (it will have to match the datasource)  
  project_name = n;
  
  qry = "";
  g = p;
  
  // init. error codes
  errcodes[ 200 ] = "OK";
  errcodes[ 400 ] = "Bad command";
  errcodes[ 401 ] = "Bad data source";
  errcodes[ 402 ] = "Bad command arguments";
  errcodes[ 403 ] = "Bad reference object";
  errcodes[ 404 ] = "Bad stylesheet";
  errcodes[ 405 ] = "Coordinate error";
  errcodes[ 500 ] = "Internal server error";
  errcodes[ 501 ] = "Unimplemented feature";
  
  // capabilities
  capabs.insert( "features/1.1" );
  capabs.insert( "sources/1.0" );
  capabs.insert( "types/1.1" );
  
  capabs.insert( "error-segment/1.0" );
  capabs.insert( "unknown-segment/1.0" );
  capabs.insert( "unknown-feature/1.0" );
}


void DAS::add( const char * t )
{
  qry += t;
}


//
// HTTP/DAS headers
//

void DAS::header( const int retc )
{
  
  std::cerr << "writing header...\n";
  
  *oss << "HTTP/1.1 200 OK\n"
       << "Date: Sun, 12 Mar 2012 16:13:51 GMT\n"
       << "Last-Modified: Fri, 16 Feb 2009 11:17:59 GMT\n"
       << "Content-Type: text/xml\n"
       << "Access-Control-Allow-Origin: *\n"
       << "Access-Control-Expose-Headers: X-DAS-Version, X-DAS-Status, X-DAS-Capabilities, X-DAS-Server\n"
       << "X-DAS-Version: DAS/1.6\n"
       << "X-DAS-Status: " << retc << "\n"
       << "X-DAS-Capabilities:";
  
  std::set<std::string>::iterator ii = capabs.begin();
  while ( ii != capabs.end() )
    {
      *oss << *ii << ";";
      ++ii;
    }
  *oss << "\n";
  
  *oss << "X-DAS-Server: PDAS/1.1\n\n";
  
}




 // print qq(<?xml version="1.0" standalone="yes"?>\n<!DOCTYPE DASDSN SYSTEM "http://www.biodas.org/dtd/dasdsn.dtd">\n);
 // print "<DASDSN>\n";


 //
 // DAS command 'dsn' -- list data-source names
 //

 int DAS::dsn()
 {

 //   for my $dsn (sort keys %$CONFIG) {
 //       print "$j<DSN>\n";
 //       print qq($j$j<SOURCE id="$dsn">$dsn</SOURCE>\n);
 //       print qq($j$j<MAPMASTER>$CONFIG->{$dsn}{DSN}{mapmaster}</MAPMASTER>\n);
 //       print qq($j$j<DESCRIPTION>$CONFIG->{$dsn}{DSN}{description}</DESCRIPTION>\n);
 //       print "$j</DSN>\n";
 //     }
 //   print "</DASDSN>\n";

   return 200;
 }


 //
 // DAS command 'types' -- list types
 //

 int DAS::types()
 {

   ss << "<DASTYPES>"
      << "<GFF version=\"1.2\" href=\"\">"
      << "<SEGMENT version=\"1.0\">";

   ss << "<TYPE id=\"datatype1\" category=\"other\">" 
      << "<TYPE id=\"datatype2\" category=\"other\">" ;

   ss << "</SEGMENT>"
      << "</GFF>"
      << "</DASTYPES>";

   return 200;
 }



//
// DAS command 'features'
//

void f_display_variant( Variant & , void * ) ;

int DAS::features( )
{
  
  
  //  SERVER/das/DSN/features?segment=RANGE
  //                         [;segment=RANGE]
  //                         [;type=TYPE]
  //			     [;type=TYPE]
  //			     [;category=CATEGORY]
  //			     [;category=CATEGORY]
  //			     [;feature_id=ID]
  //			     [;maxbins=BINS]  
  //
  
  
  // expects at least one 'segment' OR 'feature_id'
  // all other commands are optional
  
  // Look at 'segment'
  
  std::set<Region> s = get_segments( OPT );
  
  // for now just handle a single segment ( check SPEC for how >1
  // segments are handled (if they are)
  
  if ( s.size() == 1 ) 
    {
      
      //
      // Generate region mask
      //
      
      const Region & region = *s.begin();

      std::string mstr = "reg=" + region.coordinate() ;
      std::cout << "regstr [" << mstr << "]\n";
      
      // Make the mask...

      Mask mask( mstr ) ; 

      // ...and populate indmap
      
      g->indmap.populate( g->vardb, g->phmap, mask );
      
      
      //
      // Output header
      // 
      
      ss << "<?xml version=\"1.0\" standalone=\"no\"?>"
	 << "<DASGFF>"
	 << "<GFF href=\"" << HOST << "\">";
      
      
      //
      // Fetch variants
      //
      
      g->vardb.iterate( f_display_variant, this , mask );
      
      
      //
      // Footer
      //
      
      ss << "</GFF>"
	 << "</DASGFF>";
      
      
    }
  
  
  return 200;
  
}


void f_display_variant( Variant & var , void * p )
{

  std::cerr << "writing " << var << "\n";

  DAS * pdas = (DAS*)p;  
  
  // notes: for now, hardcode 'id' to be numeric (to match) chromosome build
  //        need some way of us tracking whether to write 'chr1' or '1' as the officual co-ordinate name
  
  std::string label_part = ":" + var.name() + ":" + var.reference() + "/" + var.alternate();
  
  pdas->ss << "<SEGMENT id=\"" << var.chromosome() 
	   << "\" start=\"" << var.position() 
	   << "\" stop=\"" << var.stop() 
	   << "\" label=\"" << var << label_part
	   << "\">" << var ;
  
  pdas->ss << "<FEATURE id=\"variant\" label=\"" << var << label_part << "\">"	    
	   << "<START>" << var.position() << "</START>"
	   << "<END>" << var.stop() << "</END>";
  
  // meta-information: variant, sample-variant meta-info and FILTER fields

  std::vector<std::string> keys = var.meta.keys();

  for (int m=0; m<keys.size(); m++)
    pdas->ss << "<NOTE>" << keys[m] << "=" << var.meta.print( keys[m] ) << "</NOTE>";
  
  SampleVariant & sample = var.consensus;

  if ( sample.quality() >= 0 ) 
    pdas->ss << "<NOTE>QUAL=" << sample.quality() << "</NOTE>";
  
  pdas->ss << "<NOTE>FILTER=" << sample.filter() << "</NOTE>";
  Genotype* geno = var.genotype(-1, 0);
  pdas->ss << "<NOTE>GENO=" << geno->allele_count() << "</NOTE>";
  
  std::vector<std::string> keys1 = sample.meta.keys();
  for (int m=0; m< keys1.size(); m++)
    pdas->ss << "<NOTE>" << keys1[m] << "=" 
	     << sample.meta.print( keys1[m] ) 
	     << "</NOTE>";
  
  pdas->ss << "</FEATURE>";
  
  pdas->ss << "</SEGMENT>";
  
}



//
// DAS command 'sequence'
//

int DAS::sequences() 
{
  // SERVER/das/DSN/sequence?segment=RANGE[;segment=RANGE...]
  
  return 200;
}



 //
 // DAS 'sources' command
 //

 int DAS::sources( )
 {

   std::cerr << "proc: 'sources' command\n";
   
   ss << "<?xml version='1.0' standalone=\"no\" ?> "
      << "<SOURCES>"
      << "<SOURCE uri=\"" << DSN << "\" "
      << " title=\"" << DSN << "\" "
      << " description=\"Variation data from PLINK/Seq project\">"
      << " <MAINTAINER email=\"shaun.purcell@mssm.edu\" />"
      << " <VERSION uri=\"" << DSN << "\" created=\"2010-06-16T11:53:29+0000\">"
      << " <COORDINATES uri=\"http://www.dasregistry.org/dasregistry/coordsys/CS_DS311\""
      << " taxid=\"9606\""
      << " source=\"Chromosome\""
      << " authority=\"GRCh\" version=\"37\""
      << " test_range=\"4:32211548,32711547\">GRCh_37,Chromosome,Homo sapiens</COORDINATES>"
      << " <CAPABILITY type=\"das1:sources\" query_uri=\"http://localhost:8888/das/" << DSN << "\" />"
      << " <CAPABILITY type=\"das1:features\" query_uri=\"http://localhost:8888/das/" << DSN << "/features\" />"
      << " </VERSION>"
      << " </SOURCE>"
      << " </SOURCES>";
   
   return 200;
   
}



//
// Main command to parse DAS command
//

int DAS::parse_request()
{

  // get first line, expected form
  
  // GET http://localhost/das/dsn1/command&options

  int nline = 0;

  Helper::char_tok ltok( qry , &nline , '\n' );
  
  if ( ltok.size() < 1 ) return 400; // bad request

  //
  // Parse all headers
  //

  for (int l=1;l<ltok.size();l++)
    {
      std::string t = ltok(l);
      
      // chomp off trailing \r
      if ( t[t.size()-1] == '\r' ) t = t.substr( 0 , t.size() - 1 );

      // Host: header
      if ( t.substr(0,6) == "Host: " || 
	   t.substr(0,6) == "host: " ) 
	HOST = t.substr(6);      
    }

  // HTTP/1.1 requires a Host: header
  if ( HOST == "" ) return 400; // bad request


  //
  // Parse first specific command line
  //

  std::vector<std::string> tok = Helper::parse( ltok(0) , " \t" );

  // TODO: might whitespace ever be in the URL/request?
  // TODO: be tolerant for multiple delimiters, etc
  
  if ( tok.size() != 3 ) return 400; // bad request
  
  //  if ( tok[0] != "GET" ) Helper::halt( "pdas only implements GET currently" );
  
  if ( tok[2] == "HTTP/1.1\r" ) tok[2] = "HTTP/1.1";
  
  for (int i=0;i<tok.size();i++)
    std::cerr << "tok " << i << " [" << tok[i] << "]\n";
  
  if ( ! ( tok[2] == "HTTP/1.1" || tok[2] == "HTTP/1.1\r" ) ) 
    Helper::halt( "pdas expecting HTTP/1.1" );
  
  
  std::string q = tok[1];
  
  int p1 = q.find( "/das/" );
  
  if ( p1 == std::string::npos ) return false;

  // Host: header should always have been defined above
  //  HOST            = q.substr( 0 , p1 ) ;
  std::string cmdstr  = q.substr( p1 + 5 ) ; 
  
  std::cerr << "host [" << HOST << "]\n";
  std::cerr << "cmd  [" << cmdstr << "]\n";
  
  // command can contain a DSN and/or a command
  // options are then everything to the right of '?'
  
  int p_slash = cmdstr.find( "/" );
  int p_quest = cmdstr.find( "?" );
  
  if ( p_quest == std::string::npos ) 
    p_quest = cmdstr.size();
  
  if ( p_slash == std::string::npos ) 
    {
      // implies request in form
      //  das/DSN
      
      // implies a 'sources' command for this one DSN
      CMD = "sources";
      DSN = cmdstr;
      std::cerr << "DSN: " << DSN << "\n";
      std::cerr << "Project name: " << project_name << "\n";
    }  
  else 
    if ( p_slash > p_quest ) return 400; // bad command
  
  
  if ( CMD != "sources" ) 
    {

      // parsing das/dsn/command
      //      or das/dsn/command?options
      
      std::cerr << "parsing full command\n";
      
      // no explicit command given?
      //   das/dsn1/?segment=1:12345,67890

      // implies command = 'features' (?)
      
      DSN = cmdstr.substr( 0 , p_slash );
      
      CMD = cmdstr.substr( p_slash + 1 , p_quest - p_slash - 1 );
      if ( CMD == "" ) CMD = "features";
      
      std::string opt = cmdstr.substr( p_quest + 1 );
      
      // either '&' or ';' separate request parameters
      OPT = Helper::parse( opt , ";&" );
      
    }

  std::cerr << "act dsn [" << DSN << "]\n";
  std::cerr << "act cmd [" << CMD << "]\n";
  for (int i=0;i<OPT.size();i++)
    std::cerr << "act opt " << i << "[" << OPT[i] << "]\n";
  
  return 200;
  
}

int DAS::entry_points()
{

  //  SERVER/das/DSN/entry_points[?rows=start-end]

  ss << "<DASEP>";
  
  //     << "<ENTRY_POINTS href=\"" << HOST << "/das/" << DSN << "/entry_points\" total=\"93\" start="21" end="29">
  
  // e.g. "<SEGMENT type=\"Chromosome\" id=\"8\" start=\"1\" stop=\"146364022\" orientation=\"+\" subparts=\"yes\">8</SEGMENT>"

  ss << "</ENTRY_POINTS>"
     << "</DASEP>";

  return 200;
}


std::set<Region> DAS::get_segments( const std::vector<std::string> & p ) const
{
  std::set<Region> regions;
  for (int i=0;i<p.size();i++)
    {
      if ( p[i].substr(0,8)== "segment=" )
	{
	  std::string t = p[i].substr(8);	  
	  // swap ',' character with -
	  for (int s=0;s<t.size();s++)
	    if ( t[s] == ',' ) t[s] = '-';
	  std::cout << "swapped seg name [" << t << "]\n";
	  bool okay = true;
	  Region r( t , okay );
	  if ( okay ) {
	    std::cerr << "adding segment " << r << "\n";
	    regions.insert( r );
	  }
	}
    }
  return regions;
}



//
// DAS co-ordinate system
//

 
//Authority, (assembly) Version, Type, and Organism
// NCBI_36

// <DASCOORDINATESYSTEM>
// <COORDINATES uri="http://www.dasregistry.org/dasregistry/coordsys/CS_DS311" taxid="9606" source="Chromosome" authority="GRCh" test_range="" version="37">GRCh_37,Chromosome,Homo sapiens</COORDINATES>
// </DASCOORDINATESYSTEM>



// <DASCOORDINATESYSTEM>
// <COORDINATES uri="http://www.dasregistry.org/dasregistry/coordsys/CS_DS313" taxid="9606" source="Supercontig" authority="GRCh" test_range="" version="37">GRCh_37,Supercontig,Homo sapiens</COORDINATES>
// </DASCOORDINATESYSTEM>



void DAS::write()
{
  *oss << ss.str(); 
  
  std::cerr << "total message.....................................................\n\n";
  std::cerr << ss.str() << "\n";

}
