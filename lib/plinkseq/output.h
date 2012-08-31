#ifndef __PSEQ_OUTPUT_H__
#define __PSEQ_OUTPUT_H__

#include "plinkseq/filemap.h"

///////////////////////////////
// Primary output streams

class OutFile;

class Out {  

  // data-dumper modes
  std::string  igrp;              // track i-group IDs
  std::string  igrp_header;       // used in tabular-mode
  int          coln;              // and col max.
  std::vector<std::string> cols;  // track header columns
  bool mode_tabular;              // mode (long or tabular)
  bool mode_header;               // header/numbers in long mode
  

 public:
  
  Out( const std::string & n , const std::string & desc = "" ); 

  void init_data_mode();

  void close();

  ~Out();



  // write to stream 
  
  template<class T>
    Out & operator<<(const T & msg) 
    {

      // Output to file
      if ( tofile ) 
	{
	  if ( compressed ) 
	    {
	      *zoutfile << msg;
	      zoutfile->flush();
	    }
	  else
	    {
	      outfile << msg;
	      outfile.flush();
	    }
	}
      
      if ( stdout ) 
	{
#ifdef R_SHLIB 
	  rstream << msg;
#else
	  std::cout << msg;
	  std::cout.flush();
#endif
	}
      
      return *this;
    }
  

  void flush()
  {    
#ifdef R_SHLIB 
    // nothing
#else
    if ( compressed && zoutfile ) zoutfile->flush();
    else outfile.flush();
#endif    
  }


  std::string R_flush()
    {
      std::string r = rstream.str();
      rstream.str( std::string() );
      return r;
    }



  // set fileroot
  static void set_fileroot( const std::string & r ) { fileroot = r; }  
  static void set_stdout( const bool b ) { stdout = b; } 
  static void set_tofile( const bool b ) { tofile = b; } 

  // get a stream (fail if not open)
  static Out & stream( const std::string & ) ;

  // check if stream exists
  static bool check_stream( const std::string & f );


  // data-group options

  //
  // Long-format/tabular helper functions for data rows
  //
  
  void tabular() { mode_tabular = true; }
  void longmode() { mode_tabular = false; }
  
  void header() { mode_header = true; } 
  void noheader() { mode_header = false; } 
  
  std::map<std::string,std::string> databuf0;
  std::map<std::string,std::map<std::string,std::string> > databufk;
  
  void data_reset() 
  { 
    databuf0.clear();
    databufk.clear();
    igrp = "";
    coln = 0; 
    cols.clear(); 
  }

  void data_group_header( const std::string & s )
  {
    igrp_header = s;    
    if ( mode_tabular ) 
      {
	// this should really be the first column, 
	// and always be present
	(*this) << s;
      }   
  }
  
  void data_header( const std::string & s ) 
  { 
    ++coln;    
    if ( mode_tabular ) 
      {
	//we assume an igrp_header has already been written
	//so always put in a tab here
	(*this) << "\t" << s;
      }    
    cols.push_back(s);    
  } 

  
  void data_header_done()
  {
    if ( mode_tabular ) (*this) << "\n"; 
    coln = 0;
  }

  template<class T> bool data( const T & v )
    {
      return data( v , "0" );
    }

  template<class T> bool data( const T & v , int k )
    {
      return data( v , Helper::int2str(k) );
    }

  template<class T> bool data( const T & v , const std::string & k )
    {
      // start, coln will be 0; use 1-based codes for headers
      if ( coln == cols.size() ) return false;
      return data( v , cols[ coln++ ] , k );
    }
  
  template<class T> bool data( const T & v , const std::string & j, int k )
    {
      return data( v , j , Helper::int2str(k) );
    }
  
  template<class T> bool data( const T & v , const std::string & j , const std::string & k )
    {
           
      if ( mode_tabular ) 
	{
	  std::stringstream ss;
	  ss << v;
	  if ( k == "0" ) databuf0[ j ] = ss.str();
	  else databufk[ k ][ j ] = ss.str();
	}
      else
	{
	  (*this) << igrp << "\t"
		  << j 
		  << ( k == "0" ? "" : "{" + k + "}" ) 
		  << "\t" 
		  << v << "\n";
	}
      
      return true;
    }
  
  
  template<class T> void data_group( const T & i )
  {
    std::stringstream ss;
    ss << i;
    igrp = ss.str();
  }

  void print_data_group()
  { 

    if ( mode_tabular )
      {
	
	if ( databufk.size() == 0 )
	  {
	    (*this) << igrp;
	    for (int i=0; i<cols.size(); i++)
	      (*this) << "\t" << databuf0[ cols[i] ];
	    (*this) << "\n";
	  }
	else
	  {
	    
	    std::map<std::string,std::map<std::string,std::string> >::iterator ik = databufk.begin();
	    while ( ik != databufk.end() )
	      {

		(*this) << igrp;
				
		for (int i=0; i<cols.size(); i++)
		  {
		    std::map<std::string,std::string>::iterator f0 = databuf0.find( cols[i] );
		    if ( f0 != databuf0.end() )
		      {
			(*this) << "\t" << f0->second;
		      }
		    else
		      {
			std::map<std::string,std::string>::iterator f1 = ik->second.find( cols[i] );
			if ( f1 != ik->second.end() )
			  {
			    (*this) << "\t" << f1->second;
			  }
			else
			  (*this) << "\t.";
		      }
		  }
		
		(*this) << "\n";
		
		++ik;
	      }
	  }
      }

    // clear stores
    databuf0.clear();
    databufk.clear();
    
    // reset row counter
    coln = 0;
  }
  


 private:

  
  // one or more of these can be set */

  static bool stdout;   // output to STDOUT?  
  static bool tofile;   // output to file?
  static std::string fileroot;  // fileroot
  
  // Keep track of all open current streams
  static std::map<std::string,Out*> streams; 
  
  // properties of this stream
  bool        compressed;
  std::string name;

  // primary file output streams
  std::ofstream outfile;
  OutFile * zoutfile;

  // in R mode, use this as buffer 
  std::stringstream rstream;
  

};


#endif
