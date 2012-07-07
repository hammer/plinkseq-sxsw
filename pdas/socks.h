#ifndef __PSEQ_SOCKS_H__
#define __PSEQ_SOCKS_H__

#include <cstdio>
#include <cstdlib>
#include <string>
#include <unistd.h>
#include <sys/types.h> 
#include <sys/socket.h>
#include <sys/signal.h>  
#include <netinet/in.h>
#include <netdb.h>
#include <iostream>
#include <cerrno>
#include <string.h>

class ServerSocket { 
  
 public:

  ServerSocket( const int , 
		bool (*f)(const char *, const int, void * , std::stringstream & ) ,
		void * );
  
  ~ServerSocket();

 private:  
  
  bool (*f)( const char * , const int , std::string * response );
  
  int sockfd;
  int newsockfd;
  int portno;

  socklen_t clilen;
  char buffer[256];
  struct sockaddr_in serv_addr;
  struct sockaddr_in cli_addr;
  
  void error( const std::string & msg )
  {
    std::cerr << msg << "\n";
    exit(1);
  }
  
};


class ClientSocket {

 public:

  ClientSocket( const std::string & , 
		const int );
  
  ~ClientSocket();

  std::string send( const std::string & );
  
 private:

  int sockfd;
  int newsockfd;
  int portno;

  socklen_t clilen;
  char buffer[256];

  struct sockaddr_in serv_addr;
  struct sockaddr_in cli_addr;

  void error( const std::string & msg )
  {
    std::cerr << msg << "\n";
    std::cerr << "errno = " << errno << " : " 
	      << strerror(errno) << "\n";
    exit(1);
  }

};

#endif



