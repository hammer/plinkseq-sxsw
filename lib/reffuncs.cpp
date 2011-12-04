
#include "gstore.h"

using namespace std;
using namespace Helper;

//
// Reference database
//

void GStore::refdb_new(string s)
{
  if ( fileExists(s) ) 
    Helper::remove_file( s );
  refdb.attach( s );
} 
  
void GStore::refdb_attach(string s)
{
    refdb.attach( s );
}

void GStore::refdb_summary( )
{
  // obsolete function
  refdb.summary( true );
}

RefVariant GStore::refdb_lookup(Variant & v, int g)
{
    return refdb.lookup(v,g);
}
