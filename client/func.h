#ifndef __PSEQ_FUNC_H__
#define __PSEQ_FUNC_H__

#include "pseq.h"

#include <string>

class OptGMatrix;
class OptGMetaMatrix;
class OptUniq;

namespace Pseq 
{
  
  void finished();

  namespace Util {
    class Options;
    class Commands;
  }

  class VStat;    
  
  bool new_project( std::string project , Pseq::Util::Options & args );

  bool set_project( std::string project );

  namespace VarDB
    {
      
      bool summary( Mask & , bool );

      bool attach( std::string db );

      bool load_VCF( Mask & );
      
      bool header_VCF( const bool meta , const bool header , Mask & );
      
      bool load_PLINK( const std::vector<std::string> & , const Pseq::Util::Options & , const std::string & );

      bool load_dosage();

      bool swap_ids( const std::string & );

      bool flush( const std::vector<int> & );
      
      bool vacuum();

      bool write_VCF(Mask & m , bool compressed );

      bool write_BCF(Mask & m, const std::string & );

      bool write_PED(Mask & m, std::string, bool);

      bool write_lik(Mask & m);

      bool write_haps(Mask & m, const std::string & );

      bool write_matrix(Mask & m);

      bool write_lookup_matrix(Mask & , const std::string & , const std::vector<std::string> & );

      bool write_meta_matrix(Mask & m);

      bool write_var_meta_matrix(Mask & m, std::string & name);

      bool write_gene_matrix(Mask & m, OptGMatrix & opt);

      bool write_gene_meta_matrix(Mask & m, OptGMetaMatrix & opt);
      
      bool dump_indiv();

      bool consolidate(Mask & m, std::string label);

      bool uniq_report( std::vector<std::string> & , Mask & m , OptUniq & );
      
      bool denovo_scan( Mask & m );

      bool gene_stats_header( VStat & );

      bool variant_stats(Mask & m);

      bool annotate_loc( const std::string & , Mask & m);
      
      bool write_vardb( const std::string & , const std::string & , Mask & m);
      
      bool check_concordance(Mask &m);
      
      bool cluster_scan(Mask &m);

      bool proximity_scan(Mask &m);
      
      bool make_counts_file(Mask &m , const std::string & );

      bool simple_counts(Mask &m, bool);

      bool lookup_list(const std::string & filename, Mask &m, const std::vector<Region> * regs = NULL );

      bool simple_sim();

      bool vdist_summary(Mask & m, long int );

      bool add_to_varset( const std::string & , Mask & , const std::string & mtag = "" , const std::string & desc = "." );

      bool add_to_varset( const std::string & );

      bool add_superset_from_file( const std::string & );

      bool add_superset( const std::string & , const std::vector<std::string> & , const std::string & );

      bool insert_meta_from_file( const std::string & ) ;

      bool insert_meta_on_fly( const std::string & );

    }    
  
  namespace PPH2DB
  {
    bool load( const std::string & dbname , const std::string & filename );
    bool score( Mask & m , const std::string & dbname ); 
  }

  namespace NetDB
  {

    bool loader( const std::string & db , const std::string & file );

    bool lookup( const std::string & db , const std::string & gene , const std::string & grp );

  }


  namespace LocDB
    {

      bool summary( LocDBase * , bool );      

      bool attach( std::string db );

      bool load_GTF( std::string gtf , std::string label , bool );

      bool load_generic_regions( std::string & filename , const std::string & label , Pseq::Util::Options & , bool );

      bool update_searchtable( const std::string & g , bool locdb = true );

      bool load_segments( std::string filename , std::string label , Pseq::Util::Options & opt );

      bool load_set( std::string file , std::string set_label , std::string group_label , bool , bool );

      bool merge( std::string label1 , std::string label2 , bool );

      bool index();

      bool swap_alternate_names( const std::string & group , const std::string & filename );

      bool overlap_analysis( );

      bool intersection( std::string filename , std::string group , LocDBase &  );
      
      bool load_pathway( std::string file1 , std::string file2 , std::string group , bool , bool use_altname = false );

      bool loc_view( const std::string & , const std::vector<std::string> & , const bool meta = true , const bool subregions = false );

      //      bool loc_overlap( );
    }


  namespace RefDB
    {
      
      bool summary( bool );      
      
      bool attach( std::string );
      
      bool load_refvar( const std::string & filename , const std::string & name , Pseq::Util::Options & opt );
      
      bool load_VCF( const std::string & filename , const std::string & name );
    }
  
  namespace IndDB
    {
      
      bool attach( std::string );
      
      bool summary( bool );
      
      bool set_phenotype( std::string );

      bool make_phenotype( std::string );
      
      bool load_ped_info( std::string );
      
      bool load_phenotypes( std::string );
      
      bool swap_ids( const std::string & );

      bool dump_table(Mask &);

      bool make_residuals( const std::vector<std::string> & );
      
      bool dump_phenotypes( const std::vector<std::string> & );
    }
  
  

  namespace Util
    {
      
      Mask construct_mask( std::string );
      
      void set_default( VStat & );

      bool file_summary( bool );
      
      bool meta_summary( bool ) ;
      
    }
  
  namespace SeqDB
    {
      
      bool load_FASTA( const std::string & );
      
      bool summary( bool );
      
      bool lookup( Region & );
      
      bool load_transcripts( std::string );

      bool loc_stats( const std::string & , const std::string & );

      bool loc_translate( const std::string & );
    }

}


#endif
