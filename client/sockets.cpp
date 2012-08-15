
#include "util.h"
#include "plinkseq.h"

#ifndef SKIP

#include <cstdio>
#include <cstdlib>
#include <string> 
#include <unistd.h>

// #ifdef WIN
// #include <windows.h>
// #endif

//#ifdef UNIX
#include <sys/types.h>    // Needed for system defined identifiers.
#include <netinet/in.h>   // Needed for internet address structure.
#include <sys/socket.h>   // Needed for socket(), bind(), etc...
#include <arpa/inet.h>    // Needed for inet_ntoa()
#include <fcntl.h>
#include <netdb.h>
#include <signal.h>
//#endif

#endif      // of SKIP 

std::vector<std::string> Pseq::Util::socket_connection( const std::string & ip_addr ,
							int port , 
							const std::string & message )
{
  
  int BUF_SIZE = 4096;
  
#ifndef SKIP

  //  plog << "Connecting to web... ";
  
  std::vector<std::string> tokens(0);

// #ifdef WIN
//   WORD wVersionRequested = MAKEWORD(1,1);    // Stuff for WSA functions
//   WSADATA wsaData;                           // Stuff for WSA functions
// #endif
  
  unsigned int         server_s;             // Server socket descriptor
  struct sockaddr_in   server_addr;          // Server Internet address
  char                 out_buf[BUF_SIZE+1];  // Output buffer for GET request
  char                 in_buf[BUF_SIZE+1];   // Input buffer for response
  unsigned int         retcode;              // Return code
  unsigned int         i;                    // Loop counter
  
// #ifdef WIN
//   WSAStartup(wVersionRequested, &wsaData);
// #endif

  // Create a socket
  server_s = socket(AF_INET, SOCK_STREAM, 0);
  
  // Fill-in the Web server socket's address information
  server_addr.sin_family = AF_INET;         // Address family to use
  server_addr.sin_port = htons(port);       // Port num to use
  
  server_addr.sin_addr.s_addr = inet_addr(ip_addr.c_str()); // IP address to use
  //server_addr.sin_addr = *((struct in_addr *)he->h_addr); 

  // Do a connect (connect() blocks)
  retcode = connect(server_s, (struct sockaddr *)&server_addr,
                    sizeof(server_addr));
  if (retcode != 0)
  {
    plog << " failed connection";
    
// #ifdef WIN
//   WSACleanup();
// #endif

    return tokens;
  }

  // Send a message to the server (appending a null terminator)
  
  send(server_s, (message+'\0').c_str(), message.length()+1, 0);

  // Receive from the Web server
  
  int echoStringLen = 100;
  std::string all_string = "";

  char echoBuffer[BUF_SIZE + 1];    // Buffer for echo string + \0

  // Receive the same string back from the server
  
  while ( 1 ) 
    {
      int retcode = recv(server_s, echoBuffer, BUF_SIZE, 0);
      
      // Give up if we encounter any problems
      
      if ( retcode < 0 )
	{
	  plog << "problem reading from pseq-server\n";
	  return tokens;
	}
      
      echoBuffer[retcode] = '\0';        // Terminate the string!    
      
      all_string += echoBuffer;
      
      // Is this the end of the input?
      
      if ( echoBuffer[retcode-1] == '\0' )
	break;
      
    }
  
  std::string buf; 
  std::stringstream ss(all_string); 
  while (ss >> buf)
    tokens.push_back(buf);
  return tokens;
  
  // Close all open sockets
// #ifdef WIN
//   closesocket(server_s);
// #endif

//#ifdef UNIX
  close(server_s);
//#endif
  
// #ifdef WIN
//   WSACleanup();
// #endif

#endif

}
