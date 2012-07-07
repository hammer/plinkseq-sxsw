#ifndef __PLINKSEQ_DAS_H__
#define __PLINKSEQ_DAS_H__

#include <string>
#include <map>
#include <set>
#include <vector>

#include "plinkseq.h"

class GStore;

class DAS { 
  
 public:
  
  DAS( GStore * , const std::string &  );

  bool check_datasource( ) const { 
    return DSN == project_name || DSN + ".pseq" == project_name ; 
  } 

  void add( const char * );

  bool command( const std::string & c ) { return CMD == c; } 

  std::string query() const { return qry; }

  void attach_response_stream( std::stringstream * p ) { oss = p; }

  int parse_request();


  // Write HTML/DAS header, error and main body

  void header( const int );

  void error( const int );

  void write();


  // DAS commands

  int capabilities();

  int dsn();

  int entry_points() ;

  int types();

  int sources();

  int features();

  int stylesheet();

  int sequences();
  
  // internal buffer (written to by iterator functions, so make this public)
  std::stringstream ss;
  
  
 private:

  // keep track of plink/seq project name
  GStore * g;
  std::string project_name;

  // unparsed entire command 
  std::string qry;

  // HEAD, GET, POST
  std::string request_method; 

  // (Host) URL
  std::string HOST;

  // DSN
  std::string DSN;
  
  // Operation
  std::string CMD;
  
  // Parameterss
  std::vector<std::string> OPT;

  // Error codes for DAS
  std::map<int,std::string> errcodes;

  // Capabilities for DAS
  std::set<std::string> capabs;

  // extract all 'segment' parameters
  std::set<Region> get_segments( const std::vector<std::string> & ) const;

  // All response to client goes via this std::stringstream
  std::stringstream * oss;
  
};


struct PDAS_param { 

  PDAS_param() 
  {
    project = "";
    portn = 8080;
  }
  std::string project;
  int portn;
};

#endif
