#include "socks.h"
#include <sstream>

ServerSocket::ServerSocket( const int p , 
			    bool (*f)(const char *, const int, void * , std::stringstream & ) , 
			    void * data )
{

  portno = p;
  
  // avoid zombie children
  signal(SIGCHLD,SIG_IGN);
  
  sockfd = socket(AF_INET, SOCK_STREAM, 0);
  if ( sockfd < 0 ) error( "ERROR opening socket" );
  
  bzero((char *) &serv_addr, sizeof(serv_addr));  
  serv_addr.sin_family = AF_INET;
  serv_addr.sin_addr.s_addr = INADDR_ANY;
  serv_addr.sin_port = htons(portno);
  
  if ( bind( sockfd, (struct sockaddr *) &serv_addr,
	     sizeof(serv_addr) ) < 0 ) error("ERROR on binding");


  listen( sockfd ,  5 );

  clilen = sizeof( cli_addr );


  // go into loop 

  while ( 1 ) 
    {
      
      newsockfd = accept(sockfd, 
			 (struct sockaddr *) &cli_addr, 
			 &clilen);
      
      if ( newsockfd < 0 ) error("ERROR on accept");

      const int pid = fork();

      if ( pid < 0 ) error("ERROR on fork");

      if ( pid == 0 )  
	{

	  close(sockfd);
	  
	  //
	  // Read message
	  //
	  
	  std::stringstream response;
	  
	  bool done = false; 
	  
	  while ( ! done ) 
	    {
	      
	      bzero(buffer,256);
	      
	      const int n = read(newsockfd,buffer,255);
	      
	      if (n < 0) error("ERROR reading from socket");
	      
	      //
	      // Process, and get any response
	      //
	      
	      done = f( buffer , n , data , response );
	      
	    }
	  
	  //
	  // Response
	  //
  
	  std::string s = response.str();

	  const int n_reply = write( newsockfd, s.data() , s.size() );

	  if (n_reply < 0) error("ERROR writing to socket");
	  
	  exit(0);
	}
      else 
	close(newsockfd);
    }
  

  close(sockfd);

  // should never get here  
}

ServerSocket::~ServerSocket()
{
  close(newsockfd);
  close(sockfd);
}


ClientSocket::ClientSocket( const std::string & hostname , 
			    const int p ) 
{
  
  struct hostent *server;

  portno = p;

  sockfd = socket(AF_INET, SOCK_STREAM, 0);

  if (sockfd < 0) error("ERROR opening socket");
  
  server = gethostbyname( hostname.c_str() );

  if (server == NULL) 
    {
      std::cerr << "ERROR, no such host\n";
      exit(0);
    }

  bzero((char *) &serv_addr, sizeof(serv_addr));

  serv_addr.sin_family = AF_INET;

  bcopy((char *)server->h_addr, 
	(char *)&serv_addr.sin_addr.s_addr,
	server->h_length);
 
 serv_addr.sin_port = htons(portno);

  if (connect(sockfd,(struct sockaddr *) &serv_addr,sizeof(serv_addr)) < 0) 
    error("ERROR connecting");
  
}

std::string ClientSocket::send( const std::string & m )
{

  
  for (int i = 0 ; i < m.size() ; i += 256 )
    {

      bzero(buffer,256);      
      //buffer = &m[i];
      int sz = 256;
      
      if ( m.size() - i < 256 ) sz = m.size() - i;

      std::cout << "want to send " << sz << " byes\n";
      
      const int n = write( sockfd, &m[i] , sz );
      
      if (n < 0) error("ERROR writing to socket");
    }
  
  //
  // Now get a response
  //

  std::stringstream ss;

  while ( 1 ) 
    {

      bzero(buffer,256);
      
      const int n = read(sockfd,buffer,255);

      ss << buffer;

      if ( n < 0 ) error("ERROR reading from socket");

      if ( n == 0 ) break; 

    }
  
  return ss.str();
}


ClientSocket::~ClientSocket()
{
  std::cout << "closing client socket\n";
  close(sockfd);
}

