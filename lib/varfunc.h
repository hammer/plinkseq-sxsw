#ifndef __PLINKSEQ_VARFUNC_H__
#define __PLINKSEQ_VARFUNC_H__

namespace VarFunc {
 
  std::vector<bool> missing_genotype_mask( const Variant & v );
  std::vector<bool> missing_genotype_mask( const VariantGroup & v );
}


#endif
