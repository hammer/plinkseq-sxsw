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

  DAS das( &g );
  
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

  // parse DAS request 
  
  int retc = pdas->parse_request(); 

  std::cerr << "retc = " << retc << "\n";
  // error out?
  
  if ( retc != 200 )
    {
      pdas->header( retc );
      return true; 
    }

  pdas->header( retc );

  // perform query

  if ( pdas->command( "sources" ) )
    {
      
    }
  
  else if ( pdas->command( "features" ) )
    {
      retc = pdas->features();
    }
    
  //pdas->write();

  // return as DAS to server 
  
  return true;
  
}



//
// -------  DAS object --------
//

DAS::DAS( GStore * p) 
{
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
  
  *ss << "HTTP/1.1 200 OK\n"
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
      *ss << *ii << ";";
      ++ii;
    }
  *ss << "\n";

  *ss << "X-DAS-Server: PDAS/1.1\n\n";
  
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

  *ss << "<DASTYPES>"
      << "<GFF version=\"1.2\" href=\"\">"
      << "<SEGMENT version=\"1.0\">";
  
  *ss << "<TYPE id=\"datatype1\" category=\"other\">" 
      << "<TYPE id=\"datatype2\" category=\"other\">" ;
  
  *ss << "</SEGMENT>"
      << "</GFF>"
      << "</DASTYPES>";

  return 200;
}


//
// DAS command 'features'
//

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
  
  if ( s.size() == 1 ) 
    {
      
      const Region & region = *s.begin();
      
      std::cout << "looking...\n";
      
      if ( region.chromosome() == 22 ) 
	{
	  std::cout << "in reg 22 \n";

	  // RESPONSE: DASGFF document
	 	  	  
	  *ss << "<?xml version=\"1.0\" standalone=\"no\"?>"
	      << "<DASGFF>"
	      << "<GFF href=\"" << HOST << "\">";
	  
	  *ss << "<SEGMENT id=\"22\" start=\"2000000\" stop=\"3000000\" version=\"2.22\" label=\"label22\">S1";
	  
	  *ss << "<FEATURE id=\"id1\" label=\"my label\">"	    
	      << "<START>2000100</START>"
	      << "<END>2010000</END>"
	      << "<NOTE>Test from Shaun</NOTE>"
	      << "</FEATURE>";
	  
	  *ss << "</SEGMENT>";

	  std::cout << "returning ... \n"
	    
		    << "<?xml version=\"1.0\" standalone=\"no\"?>" 
		    << "<DASGFF>"
		    << "<GFF href=\"" << HOST << "\">"	  
		    << "<SEGMENT id=\"22\" start=\"2000000\" stop=\"20010000\" version=\"2.22\" label=\"label22\">S1"
		    << "<FEATURE id=\"id1\" label=\"my label\"></FEATURE>"
		    << "</SEGMENT>"	    
		    << "</GFF>"
		    << "</DASGFF>"
	    
		    << "\n\n";
	  

// 	     << "<TYPE d
//       <FEATURE id="id" label="label">
//     <TYPE id="mytype" category="category" reference="yes|no" cvId="SO:1234">My Type</TYPE>
//     <METHOD id="mymethod" cvId="ECO:5678">My Method</METHOD>
//     <START> start </START>
//     <END> end </END>
//          <SCORE> [X.XX|-] </SCORE>
//          <ORIENTATION> [0|-|+] </ORIENTATION>
//          <PHASE> [0|1|2|-]</PHASE>
//     <NOTE> note text </NOTE>
//     <LINK href="url"> link text </LINK>
//     <TARGET id="id" start="x" stop="y"> target name </TARGET>
//          <PARENT id="parent id1" />
//          <PART id="child id1" />
//          <PART id="child id2" />
//       </FEATURE>
//       <FEATURE id="child id1" label="child label">
//          ...
//       </FEATURE>
//       <FEATURE id="child id2" label="child label">
//          ...
//       </FEATURE>
//       ...
//       <FEATURE id="parent id1" label="parent label">
//          ...
//       </FEATURE>
//       ...
//   </SEGMENT>

	  * ss << "</GFF>"
	       << "</DASGFF>";
	}
    }  

  return 200;
  
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

  *ss << 
    "<?xml version='1.0' standalone=\"no\" ?> "
    "<SOURCES>"
    "<SOURCE uri=\"URI\" title=\"title\" doc_href=\"helpURL\" description=\"description\">"
    "<MAINTAINER email=\"email address\" />"
    "<VERSION uri=\"URI\" created=\"date\">"
    "<COORDINATES uri=\"URI\""
    " source=\"data type\""
    " authority=\"authority\""
    " taxid=\"taxonomy\""
    " version=\"version\""
    " test_range=\"id:start,stop\" >coordinate string</COORDINATES>"
    "<CAPABILITY type=\"das1:command\" query_uri=\"URL\" />"
    "<PROP name=\"key\" value=\"value\" />"
    "</VERSION>"
    "</SOURCE>"
    "</SOURCES>";

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
  
  if ( ltok.size() < 1 ) return false;

  std::vector<std::string> tok = Helper::parse( ltok(0) , " \t" );

  // TODO: might whitespace ever be in the URL/request?
  // TODO: be tolerant for multiple delimiters, etc
  
  if ( tok.size() != 3 ) return false;
  
  //  if ( tok[0] != "GET" ) Helper::halt( "pdas only implements GET currently" );
  
  for (int i=0;i<tok.size();i++)
    std::cerr << "tok " << i << " [" << tok[i] << "]\n";
  
  if ( ! ( tok[2] == "HTTP/1.1" || tok[2] == "HTTP/1.1\r" ) ) 
    Helper::halt( "pdas expecting HTTP/1.1" );
  
  std::string q = tok[1];
  
  int p1 = q.find( "/das/" );
  
  if ( p1 == std::string::npos ) return false;
  
  HOST                = q.substr( 0 , p1 ) ;
  std::string cmdstr  = q.substr( p1 + 5 ) ; 
  
  std::cerr << "Host [" << HOST << "]\n";
  std::cerr << "Cmd  [" << cmdstr << "]\n";
 
  // command can contain a DSN and/or a command
  // options are then everything to the right of '?'

  int p_slash = cmdstr.find( "/" );
  int p_quest = cmdstr.find( "?" );
  
  if ( p_slash == std::string::npos ) return false;
  if ( p_quest == std::string::npos ) p_quest = cmdstr.size();
  if ( p_slash > p_quest ) return false;

  // no explicit command given?
  //   das/dsn1/?segment=1:12345,67890

  // implies command = 'features' (?)

  DSN = cmdstr.substr( 0 , p_slash );
    
  CMD = cmdstr.substr( p_slash + 1 , p_quest - p_slash - 1 );
  if ( CMD == "" ) CMD = "features";
  
  std::string opt = cmdstr.substr( p_quest + 1 );
  
  // either '&' or ';' separate request parameters
  OPT = Helper::parse( opt , ";&" );
  
  std::cerr << "act dsn [" << DSN << "]\n";
  std::cerr << "act cmd [" << CMD << "]\n";
  for (int i=0;i<OPT.size();i++)
    std::cerr << "act opt " << i << "[" << OPT[i] << "]\n";

  return 200;
  
}

int DAS::entry_points()
{

  //  SERVER/das/DSN/entry_points[?rows=start-end]

  *ss << "<DASEP>";

  //     << "<ENTRY_POINTS href=\"" << HOST << "/das/" << DSN << "/entry_points\" total=\"93\" start="21" end="29">

  // e.g. "<SEGMENT type=\"Chromosome\" id=\"8\" start=\"1\" stop=\"146364022\" orientation=\"+\" subparts=\"yes\">8</SEGMENT>"

  *ss << "</ENTRY_POINTS>"
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
	  for (int s=8;s<t.size();s++)
	    if ( t[s] == ',' ) t[s] = '-';
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
