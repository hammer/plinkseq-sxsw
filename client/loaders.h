#ifndef __PSEQ_LOADERS_H__
#define __PSEQ_LOADERS_H__

#include "plinkseq.h"
#include "func.h"

namespace Pseq {

  namespace Util {
    class Options;
  }
  
  namespace RefDB {
    bool load_refvar( std::string & , const std::string & , Pseq::Util::Options & );
  }

}

#endif
