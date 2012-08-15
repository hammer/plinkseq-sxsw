#include "util.h"
#include "plinkseq.h"
#include <ctime>

extern GStore g;
extern Pseq::Util::Options args;
extern Pseq::Util::Commands pcomm;

extern std::string PSEQ_VERSION;

// point to atgu.mgh.harvard.edu
#define  PORT_NUM   80     
#define  IP_ADDR    "155.52.206.145"
#define  GET_STRING "GET /plinkseq/version.dat HTTP/1.1\nHost: atgu.mgh.harvard.edu\nConnection: close\n\n"

void Pseq::Util::webcheck( const std::string & cmd )
{

#ifdef SKIP_WEBCHECK
  plog << "Skipping web-based version check...\n";
  return;
#else
  
  
  //
  // First see if local .pseqversion file exists in the working
  // directory, and if it is up-to-date
  //
  
  time_t curr=time(0);
  std::string tdstamp = (std::string)ctime(&curr);      
  std::string buf; 
  std::stringstream ss(tdstamp); 
  std::vector<std::string> date_tokens; 
  while (ss >> buf) date_tokens.push_back(buf);
  std::string thisDate = date_tokens[0] + date_tokens[1] + date_tokens[2];
  
  bool hasRecord = Helper::fileExists( ".pseqversion" );
  
  
  //
  // Web-based message (but may be cached in local file)
  //

  std::vector<std::string> tokens;
  
  bool connect2web = true;

  plog << "Web version check (--noweb to skip)... ";
  

  //
  // If we have a record, are we up-to-date?
  //
  
  if ( hasRecord )
    {
      std::ifstream VER;
      VER.open( ".pseqversion" , std::ios::in );
      
      std::string oldDay, oldMonth, oldDate, webVersion;
      VER >> oldDay >> oldMonth >> oldDate;
      
      if ( thisDate == oldDay+oldMonth+oldDate )
	{
	  plog << "recently cached: ";
	  connect2web = false;
	  
	  // Read rest of cached web message
	  while ( ! VER.eof() )
	    {
	      std::string t;
	      VER >> t;
	      if (t=="") break;
	      tokens.push_back(t);
	    }	  
	}
      VER.close();
      
    }
  
  
  
  if ( connect2web ) 
    {
      //      plog << "Connecting to web to get version...\n";      
      tokens = socket_connection( IP_ADDR, PORT_NUM, GET_STRING);      
    }
  
  bool print = false;
  bool print2 = false;
  bool version_okay = true;

  for (int i=0; i<tokens.size(); i++)
    {

      if (tokens[i]=="END") break;

      if (tokens[i]=="END-MESSAGE")
	{
	  print2=false;
	  continue;
	}

      if (tokens[i]=="WARN")
	{

	  if ( i < tokens.size()-1 ) 
	    {

	      i++;

	      // specifically check the actual command being used 

	      if ( cmd == tokens[i] )
		{
		  plog.warn("A warning has been issued for the '" + cmd + "' command\n"
			    "Please see http://research.mssm.edu/statgen/plinkseq/warnings.shtml" );

		}
	    }
	  continue;
	}
      
      
      if ( tokens[i]=="FATAL" )
	{

	  if ( i < tokens.size()-1 ) 
	    {
	      i++;

	      if ( cmd == tokens[i] )
		Helper::halt( "Web-check gives a serious warning for '" + cmd + "'\n"
			      "Please see http://research.mssm.edu/statgen/plinkseq/warnings.shtml" );
	      
	    }
	  continue;
	}


      if (tokens[i]=="MESSAGE-ALL")
	{
	  print2=true;
	  continue;
	}

      // Display any other messages
      // Either conditional on old version (print)
      // or a broadcast to all users (print2)

      if ( ( print && !version_okay) || print2 ) 
	{
	  if (tokens[i]=="\\n")
	    plog << "\n";
	  else
	    plog << tokens[i] << " ";
	}

      // Check version code
      if ( tokens[i] == "PSEQ-VERSION" ) 
	{
	  print=true;
	  if ( i < tokens.size() - 1) 
	    {	      
	      if ( tokens[i+1] == PSEQ_VERSION ) 
		plog << " OK, current\n";
	      else
		{
		  plog.warn( "*** PSEQ UPDATE REQUIRED : ( " + PSEQ_VERSION + " -> " + tokens[i+1] + " )" );
		  version_okay=false;
		}

	      // Skip the version number
	      i++;
	    }
	}

    }

  // did we get the information we needed?
  if ( !print ) plog << " could not connect\n";


  ////////////////////////////////////////////////////
  // Create a record that we've checked

  // First line is date-stamp; then simply copy the 
  // whole message
  
  std::ofstream VER;

  VER.open( ".pseqversion", std::ios::out );

  VER << date_tokens[0] << " "
      << date_tokens[1] << " "
      << date_tokens[2] << "\n";
  for (int i=0; i<tokens.size(); i++)
    VER << tokens[i] << "\n";
  VER.close();


#endif  
}

