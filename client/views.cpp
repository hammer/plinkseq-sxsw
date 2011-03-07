
#include "views.h"
#include "em.h"
#include "func.h"

#include "pseq.h"

#include <iostream>
#include <cmath>

extern GStore g;
extern Pseq::Util::Options options;

void f_view( Variant & v , void * p )
{

  OptVView * opt = (OptVView*)p;

  plog << v << "\t";
  
  plog << "." << "\t" 
       << v.consensus << "\t"    
       << v.n_samples() << "\t";

  plog << v.print_meta_filter();

  if ( opt->vmeta )
    {
      if ( opt->vexpand ) plog << "\n" << v.meta.display();
      else plog << "\t" << v.meta ;      

      if ( opt->vexpand ) plog << "\n" << v.consensus.meta.display();
      else plog << "\t" << v.consensus.meta ;      
    }

  plog << "\n";
  
  //
  // Show separaetly specific sample information? 
  //
  
  if ( opt->show_samples )
    {
      
      v.set_first_sample();

      while ( 1 ) 
	{

	  const SampleVariant & sample = v.sample();
	  
	  plog << v << "\t"
	       << sample.file_name() << "\t"
	       << sample << "\t"
	       << sample.filter();

	  if ( opt->vmeta ) 
	    {
	      if ( opt->vexpand ) plog << "\n" << sample.meta.display();
	      else plog << "\t" << sample.meta ;
	    }

	  plog << "\n";

	  if ( ! v.next_sample() ) break;
	}      
    }


  //
  // Show genotypes and genotype meta-information? One row per individual
  //
  
  if ( opt->geno )
    {
      
      const int n = v.size();
      
      for ( int i=0; i<n; i++)
	{
	  
	  plog << v.ind( i )->id() << "\t";
	  
	  plog << v.sample_label( i, "," ) << "\t";
	  
	  //
	  // Optionally, phenotype
	  //

	  if ( g.phmap.type() == PHE_DICHOT )
	    {
	      if ( v.ind(i)->affected() == CASE ) plog << "CASE\t";
	      else if ( v.ind(i)->affected() == CONTROL ) plog << "CONTROL\t";
	      else plog << ".\t";
	    }
	  else if ( g.phmap.type() == PHE_QT )
	    plog << v.ind(i)->qt() << "\t";
	  else if ( g.phmap.type() == PHE_FACTOR )
	    plog << v.ind(i)->group_label() << "\t";

	  //
	  // Genotype
	  //

	  plog << v.label( i , "," );

	  //
	  // Gentype meta-information
	  //

	  if ( opt->gmeta )
	    plog << "\t[" << v.gmeta_label(i) << "]";
	  
	  plog << "\n";
	}
    }



}


void g_view( VariantGroup & vars , void * p )
{
  OptGView * opt = (OptGView*)p;
  
  plog << vars.dump( opt->vmeta , 
			  opt->vexpand , 
			  opt->geno , 
			  opt->gmeta , 
			  opt->transpose , 
			  opt->rarelist , 
			  opt->show_phenotype ) << "\n";
  
}

void f_view_tped( Variant & v , void * p )
{
  std::ofstream * TPED = (std::ofstream*)p;

  *TPED << v.chromosome() << "\t";
  
  if ( v.name() == "." )
    *TPED << "var_" 
	  << Helper::chrCode( v.chromosome() ) << "_" 
	  << v.position() << "\t";
  else
    *TPED << v.name() << "\t";

  *TPED << "0\t"
	<< v.position();
  
  int n = v.size();
  
  for (int i=0;i<n;i++)
    *TPED << "\t" << v.print_PED( v(i) );

  *TPED << "\n";

}

void f_view_lik( Variant & v , void * p )
{

  //
  // output genotype likelihoods from GENO_LIK field ("GL" or "PL" by default
  // in GATK, specified in lib/defs.*)
  //

  const double EPS = 1e-5;

  GStore * g = (GStore*)p;
  
  int n = v.size();
  
  plog << v << "\t" 
	    << v.reference() << "\t"
	    << v.alternate();

  for (int i = 0 ; i < n ; i++)
    {
      
      MetaInformation<GenMeta> & genotype = v(i).meta;
      
      std::vector<double> l;
      
      // Extract either PL (or GL, if no PL)

      bool using_gl = false;
      
      if ( genotype.hasField( PLINKSeq::META_GENO_PHRED() ) )
	{
	  std::vector<int> t = genotype.get_int( PLINKSeq::META_GENO_PHRED() );
	  if ( t.size() == 3 ) 
	    {
	      l.resize(3);
	      for (int i=0;i<2;i++) l[i] = t[i];
	    }
	}
      else if ( genotype.hasField( PLINKSeq::META_GENO_LIK() ) )
	{
	  l = genotype.get_double( PLINKSeq::META_GENO_LIK() );
	  using_gl = true;
	}



      // Currently, only handle bi-allelic cases


      if ( l.size() == 3 ) 
	{

	  // subtract out maximum log-likelihood, so that exp(max)==1
	  // for a decent scaling (we can set really small values to 0
	  // then and avoid over/under flow)
	  
	  double g0 = l[0];
	  double g1 = l[1];
	  double g2 = l[2];
	  
	  if ( ! using_gl ) 
	    {
	      g0 = -g0/10.0;
	      g1 = -g1/10.0;
	      g2 = -g1/10.0;
	    }

	  double mx = g0 > g1 ? g0 : g1 ;
	  mx = g2 > mx ? g2 : mx;
	  
	  g0 -= mx;
	  g1 -= mx;
	  g2 -= mx;
	  
	  g0 = exp(g0);
	  g1 = exp(g1);
	  g2 = exp(g2);
	  
	  if ( ! Helper::realnum( g0 ) ) g0 = 0;
	  if ( ! Helper::realnum( g1 ) ) g1 = 0;
	  if ( ! Helper::realnum( g2 ) ) g2 = 0;
	  
	  double s = 1 / ( g0 + g1 + g2 );
	  
	  g0 *= s;
	  g1 *= s;
	  g2 *= s;
	  
	  // For a slightly smaller file, set all trivial likelihoods
	  // to exactly 0 (?okay)
	  
	  if ( g0 < EPS ) g0 = 0;
	  if ( g1 < EPS ) g1 = 0;
	  if ( g2 < EPS ) g2 = 0;
	  
	  if ( ! ( Helper::realnum( g0 ) && Helper::realnum( g1 ) && Helper::realnum( g2 ) ) )
	    {
	      plog << "\t0.333\t0.333\t0.333";
	    }
	  else
	    plog << "\t" << g2 
		 << "\t" << g1 
		 << "\t" << g0;
	}    
      else
	{
	  // flat ratio ~ missing data point
	  plog << "\t0.333\t0.333\t0.333";
	}	    
    }
  
  plog << "\n"; 
   
}



void f_view_matrix( Variant & v , void * p )
{
  
  // output simple 0/1/2/NA calls 
  
  plog << v << "\t" 
       << v.reference() << "\t"
       << v.alternate();
  
  int n = v.size(); 
  for (int i = 0 ; i < n ; i++)
    {
      if ( v(i).null() ) 
	plog << "\tNA";
      else
	plog << "\t" << v(i).allele_count();
    }
  plog << "\n";
  
}


void f_view_meta_matrix( Variant & v , void * p )
{

  std::map<int,std::string> * samples = (std::map<int,std::string>*)p;
  
  //
  // Core variant information
  //
  
  plog << v.chromosome() << "\t" 
       << v.position() << "\t"
       << v.name();
  
  int nelem_static = v.consensus.meta.n_visible_keys_static();
  int nelem_nonstatic = v.consensus.meta.n_visible_keys_nonstatic();

  // Need to separate out 'static' (Variant) and 'non-static'
  // (SampleVariant/Consensus) information
  
  // Static, population/variant-level meta-fields (if any)
  if ( nelem_static ) plog << "\t" << v.meta.display_row_static("NA");
  
  // Consensus information

  if ( nelem_nonstatic ) plog << "\t" << v.consensus.meta.display_row_nonstatic("NA");  

  // Filters, if only a single sample
  if ( samples->size() == 1 ) 
    {
      std::string f = v.consensus.meta_filter.display_row("NA");
      if ( f != "" ) plog << "\t" << f;
    }

  // Per-sample information

  if ( samples->size() > 1 )
    {
        
      // Samples (even if only 1 sample, loop through once and dump
      // meta_filter info from SampleVariant
            
      std::map<int,std::string>::iterator i = samples->begin(); 
      
      while ( i != samples->end() )
	{
	  
	  // note; NULL is returned if unique_svar_slot() return -1 (not seen)
	  
	  SampleVariant * s = ( ! v.multi_sample() ) ? &(v.consensus) :  v.psample( v.unique_svar_slot( i->first ) );
	  
	  if ( ! s ) 
	    {
	      
	      // If sample for this variant not present, insert missing values
	      
	      int nelem = MetaInformation<VarMeta>::n_visible_keys_nonstatic();
	      for (int i=0; i<nelem; i++) plog << "\tNA";
	      
	      nelem = MetaInformation<VarFilterMeta>::n_visible_keys();
	      for (int i=0; i<nelem; i++) plog << "\tNA";
	      
	    }
	  else
	    {
	      
	      std::string m = s->meta.display_row_nonstatic( "NA" );
	      std::string f = s->meta_filter.display_row("NA"); 

	      // Variant INFO tags
	      if ( m != "" ) plog << "\t" << m;
	      
	      // Variant FILTER tags
	      if ( f != "" ) plog << "\t" << f;
	      
	    }
	  ++i;
	}
    }

  // next row
  plog << "\n";
}



void f_view_var_meta_matrix( Variant & v , void * p )
{
  std::string name = *(std::string*)p; 
  plog << v << "\t";
  for (int i=0; i<v.size(); i++)
    plog << ( v(i).meta.has_field(name) ? 
	      v(i).meta.as_string(name) : 
	      "NA" ) << "\t";
  plog << "\n";
}



void f_view_gene_matrix( VariantGroup & vars , void * p )
{
  
  OptGMatrix * opt = (OptGMatrix*)p;
  
  // for each gene
  
  const int nv = vars.size();


  std::vector<bool> altmin( vars.size() , false );
  for (int v = 0 ; v < nv ; v++)
    altmin[v] = vars(v).n_minor_allele();
  
  std::string s = vars.name() + "\t" + Helper::int2str( nv );
     
  bool variation = false;

  const int ni = vars(0).size();

  for (int i = 0 ; i < ni ; i++)
    {
      bool valid = false;
      bool flag = false;
      
      for (int v = 0 ; v < nv ; v++)
	{
	  // need to see at least 1 non-null variant
	  if ( vars(v,i).notnull() ) 
	    valid = true;
	  
	  // we are searching for minor alleles:
	  if ( vars.geno(v,i).minor_allele( altmin[v] ) )
	    {
	      flag = true;	    
	      variation = true;
	    }
	}
      
      if ( valid ) 
	s += "\t" + Helper::int2str( (int)flag );
      else
	s += "\tNA";      
    }
  
  if( ! ( opt->hide_zero_variance && ! variation ) )
    {
      plog << s;    
      plog << "\n";
    }

}


std::map<int2,std::set<int> >::iterator XQCstats::flush( std::map<int2,std::set<int> >::iterator n )
{
  
  int2 pos = n->first;
  std::vector<std::string> & d = data[pos];
  std::vector<int> & k = datak[pos];
  
  if ( d.size() > 0 )
    {

      plog.data_group( d[0] );
      for (int i=1; i<d.size(); i++)
	plog.data( d[i] , k[i]);      

      plog.data( neighbours[ pos ].size() ); 
      
      // end of data for this variant
      plog.print_data_group();    
    }
  
  data.erase( data.find( pos ) );
  datak.erase( datak.find( pos ) );
  
  neighbours.erase( n++ );
  return n;
} 


void XQCstats::flush() 
{

  std::map<int2, std::vector<std::string> >::iterator i = data.begin();
  while ( i != data.end() )
    {	
      std::vector<std::string> & d = i->second;
      std::vector<int> & k = datak[ i->first ];

      plog.data_group( d[0] );
      for (int j=1; j<d.size(); j++)
	plog.data( d[j] , k[j]);
      
      plog.data( neighbours[ i->first ].size() , 0 );
      
      plog.print_data_group();

      ++i;
    }
  data.clear();
  datak.clear();
  neighbours.clear();
}



void f_extra_qc_metrics( Variant & var , void * p )
{

  XQCstats * x = (XQCstats*)p;
  
  int2 pos = int2( var.chromosome() , var.position() );
  
  x->add( pos , var , 0 ); // first entry will be data_group()
  x->add( pos , var.chromosome() , 0 );
  x->add( pos , var.position() , 0 );
  x->add( pos , var.reference() , 0 );
  x->add( pos , var.alternate() , 1 );

  x->add( pos , var.print_meta_filter(",") , 0 );
  x->add( pos , var.consensus.quality() == -1 ? "." : Helper::dbl2str( var.consensus.quality() ) , 0 );

  x->add( pos , var.transition() , 0 );
  x->add( pos , var.transversion() , 0 );

  //
  // Frequency information for alternate allele(s)
  //

  int c     = 0; // minor allele
  int c_tot = 0; // total counts	  
  bool altmin = var.n_minor_allele( c , c_tot );
  
  int naa, nab, nbb;
  double hwe = Helper::hwe( var , &naa , &nab , &nbb ) ;
  int t = naa + nab + nbb;

  // number of obs. genotypes
  x->add( pos , (double)t / (double)var.size() , 0 );

  // MAF
  if ( c_tot == 0 ) 
    x->add( pos , "NA" , 1 ); 
  else
    x->add( pos , (double)c/(double)c_tot , 1 ); 

  // flag if alternate is major 
  x->add( pos , !altmin , 1 );
  

  // Hardy-Weinberg equilibrium, -log10(p), and heterozygote prop.
  x->add( pos ,  hwe <= 0 ? 0 : hwe , 1 );

  // Heterozygosity
  if ( t>0 ) x->add( pos ,  (double)nab/(double)(naa+nab+nbb) , 1 );
  else x->add( pos , "NA" , 1 );
  
  
  //
  // Note -- the caller for this function will have ensured that the
  // EM-caller was run prior to coming here; thus we can just look up the
  // statistics here
  //
  
  if ( x->em_stats )
    {
      double ent_all, ent_alt;  
      
      var.em.entropy( ent_all , ent_alt );
      
      x->add( pos , var.em.mean_max_posterior() , 1 );
      x->add( pos , ent_all , 1 );
      
      if ( Helper::realnum( ent_alt ) )
	x->add( pos , ent_alt , 1 );
      else
	x->add( pos , "NA" , 1 );
      
      x->add( pos , var.em.frequency() , 1 );
      
    }


  // End of basic QC info for this variant; but we cannot finish up until 
  // we get info on the number of nearby variants; thus just track and 
  // we will print later ( via flush() )

  
  //
  // Dump contents if new are on a new chromosome
  //
  
  
  if ( x->curr_chr != var.chromosome() )
    {
      if ( x->curr_chr != -1 ) 
	x->flush();
      x->curr_chr = var.chromosome();
    }
  else
    {
      //
      // Look at all other SNPs in the library -- if this SNP is within 100bases 
      // then add
      //
      
      std::map<int2,std::set<int> >::iterator n = x->neighbours.begin();
      while ( n != x->neighbours.end() ) 
	{
	  
	  if ( abs( var.position() - n->first.p2 ) <= 100 ) 
	    {
	      n->second.insert( var.position() );
	      ++n;
	    }
	  else
	    {
	      n = x->flush( n );
	    }
	}
    }


  // this will be populated later (N SNPs)
  x->neighbours[ pos ];
 

}


//
// Dump contents of group from LOCDB
//

bool Pseq::LocDB::loc_view( const std::string & group , const std::vector<std::string> & alias )
{
  LocDBase * db = g.resolve_locgroup( group ) ;
  if ( ! db ) return false;

  int gid = db->lookup_group_id( group );
  if ( gid == 0 ) return false;
  std::set<Region> loc = db->get_regions( gid ); 
  
  std::set<Region>::iterator i = loc.begin();
  while ( i != loc.end() ) 
    {
      plog << i->name << "\t"
	   << i->coordinate() << "\t"
	   << i->altname << "\t"
	   << db->alias( i->name , false ) << "\n";
      
      ++i;
    }
  return true;
}



bool Pseq::SeqDB::loc_translate( const std::string & group )
{
  
  bool verbose = options.key("verbose");
  
  LocDBase * db = g.resolve_locgroup( group ) ;
  
  if ( ! db ) return false;
  
  int gid = db->lookup_group_id( group );
  
  if ( gid == 0 ) return false;
  
  Annotate::setDB( db , &g.seqdb );
  
  std::set<Region> loc = db->get_regions( gid ); 
  
  //
  // Header row
  //

  plog << "ID\tPOS\tNEXON\tNAA\tAASEQ\n";
  
  if ( verbose )
    {
      plog << "_BSTATS\t"
	   << "GENE" << "\t"
	   << "POS" << "\t"
	   << "EXON" << "\t"
	   << "STRAND" << "\t";
      plog << "A\tC\tG\tT\tN";
      
      std::map<std::string,std::string>::iterator pp = Annotate::aa.begin();
      while ( pp != Annotate::aa.end() )
	{
	  plog << "\t" << pp->second ;
	  ++pp;
	}
      plog << "\n";
    }
  

  //
  // Process each locus
  //

  std::set<Region>::iterator i = loc.begin();
  while ( i != loc.end() ) 
    {
      
      std::string aa = Annotate::translate_reference( *i , verbose );
      
      plog << i->name << "\t"
	   << i->coordinate() << "\t"
	   << i->subregion.size() << "\t"
	   << aa.size() << "\t";
      
      if ( aa.size() == 0 ) plog << ".";
      
      for (int p=0; p<aa.size(); p++)
	{
	  if ( p > 0 && p % 10 == 0 ) plog << " ";
	  plog << aa[p];
	}
      
      plog << "\n";

      ++i;
    }
  return true;
}


void f_simple_counts( Variant & var , void * p )
{

  OptSimpleCounts * data = (OptSimpleCounts*)p;
  
  std::string gene = ".";
  std::string annot = ".";
  std::string protein = ".";
  
  //
  // Do we have a gene attached, or do we have to look it up? 
  //
  
  if ( data->apply_annot ) 
    {
      
      if ( var.meta.has_field( PLINKSeq::META_GENE() ) )
	{
	  gene = var.meta.get1_string( PLINKSeq::META_GENE() );  
	}
      
      //
      // Do we have a functional annotation attached?
      //
      
      if ( ! var.meta.has_field( PLINKSeq::META_ANNOT() ) )
	{
	  bool exonic = Annotate::annotate( var );
	  annot = var.meta.get1_string( PLINKSeq::ANNOT_TYPE() );
	  if ( exonic ) 
	    {	  
	      gene = var.meta.get1_string( PLINKSeq::ANNOT_GENE() );
	      protein = var.meta.get1_string( PLINKSeq::ANNOT_PROTEIN() );
	    }
	}
      else
	annot = var.meta.get1_string( PLINKSeq::META_ANNOT() );
    }

  int case_count = 0 , case_tot = 0; 
  int control_count = 0 , control_tot = 0;
  
  bool ma, control_ma;

  if ( data->dichot_pheno ) 
    {
      ma = var.n_minor_allele( case_count , case_tot , CASE );
      control_ma = var.n_minor_allele( control_count , control_tot , CONTROL );
      // flip alleles?
      if ( control_ma != ma ) control_count = control_tot - control_count; 
    }
  else    
    ma = var.n_minor_allele( case_count , case_tot );
  
  // Output

  plog << "chr" << var.chromosome() << ":" << var.position() ;
  plog << "\t" << ( ma ? var.alternate() + "/" + var.reference() : var.reference() + "/" + var.alternate() ) ;

  if ( data->dichot_pheno )
    plog << "\t" << case_count
	 << "\t" << control_count ;
  else
    plog << "\t" << case_count ;
  
  if ( data->apply_annot ) 
    plog << "\t" << annot
	 << "\t" << ( gene == "" ? "." : gene ) 
	 << "\t" << protein 
	 << "\t" << g.locdb.alias( gene ) ;
  
  if ( data->meta.size() )
    {
      std::set<std::string>::iterator i = data->meta.begin();
      while ( i != data->meta.end() )
	{
	  if ( var.meta.has_field( *i ) )
	    plog << "\t" << var.meta.as_string( *i );
	  else if ( var.consensus.meta.has_field( *i) )
	    plog << "\t" << var.consensus.meta.as_string( *i );
	  else
	    plog << "\t.";
	  ++i;
	}
    }

  if ( data->show_filter ) 
    plog << "\t" << var.print_meta_filter();
  
  plog << "\n";
  
}



bool Pseq::VarDB::header_VCF( const bool show_meta , 
			      const bool show_header , 
			      Mask & mask )
{
  
  // Assumes we are in single-VCF mode, or else do nothing
  if ( ! g.single_file_mode() ) return false;
  
  std::string filename = mask.external_vcf_filename();
  Helper::checkFileExists( filename );
  
  IndividualMap imap;
  VarDBase tmpdb( imap );
  tmpdb.attach( ":memory:" );
  
  File vcffile( filename , VCF );
  VCFReader v( &vcffile , "" , &tmpdb ,  NULL );

  tmpdb.begin();

  
  while ( 1 ) 
    { 
      
      VCFReader::line_t l = v.parseLine( );
      
      if ( l == VCFReader::VCF_EOF ) break;      
      if ( l == VCFReader::VCF_INVALID ) continue;
      
      if ( l == VCFReader::VCF_META ) {  } 
      
      if ( l == VCFReader::VCF_HEADER && show_header ) 
	{
	  int n = imap.populate( tmpdb , g.phmap , mask );
	  for (int i=0; i<n; i++)
	    plog << imap(i)->id() << "\n";
	}
      
      if ( l == VCFReader::VCF_VARIANT ) break; // assume no more meta-info

    }
  
  if ( show_meta ) 
    {
      plog << MetaInformation<VarMeta>::list_fields("META_VARIANT")
	   << MetaInformation<VarFilterMeta>::list_fields("META_FILTER")
	   << MetaInformation<GenMeta>::list_fields("META_GENOTYPE");
      
    }

  return true;

}

