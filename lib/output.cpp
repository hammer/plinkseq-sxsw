#include "plinkseq/output.h"

bool Out::stdout = false;
bool Out::tofile = true;
std::string Out::fileroot = "pseq";
std::map<std::string,Out*> Out::streams; 

// Main Out class

Out::Out( const std::string & ext , const std::string & desc )
{
  name = fileroot + "." + ext;

  if( check_stream(ext))
    return;

  // add to list of currently open streams
  streams[ ext ] = this;
  
  plog << "Writing to file [ " << name << " ] ";  

  if ( desc != "" ) plog << " : " << desc;
  plog << "\n";
  
  if ( tofile ) 
    {
      if ( Helper::ends_with( name , ".gz" ) ) 
	{
	  compressed = true;
	  zoutfile = new OutFile( name );
	}
      else
	{
	  compressed = false;
	  outfile.open( name.c_str()  , std::ios::out );
	}
    }      

  init_data_mode();
}

void Out::init_data_mode()
{
  mode_tabular = true;
  mode_header = false;
  igrp = igrp_header = "";
}


Out::~Out()
{ 
  close();
}

void Out::close()
{
    if ( tofile )
    {
      if ( compressed ) 
	{
	  zoutfile->close();
	  delete zoutfile;
	}
      else
	{
	  outfile.close();
	}
    }
  
  std::map<std::string,Out*>::iterator ii = streams.find( name );
  if ( ii != streams.end() ) streams.erase( ii );
 
  tofile = ! streams.empty();
}

Out & Out::stream( const std::string & f ) 
{
  if ( streams.find( f ) == streams.end() ) 
    Helper::halt( "could not find stream " + fileroot + "." + f ) ;
  return *streams[f];
}

bool Out::check_stream( const std::string & f )
{
  if ( streams.find( f ) == streams.end() )
    return false;
  else
    return true;
}
