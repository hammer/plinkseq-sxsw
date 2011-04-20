#ifndef __PARSE_H__
#define __PARSE_H__

#include <iostream>
#include <set>
#include <vector>

#include "helper.h"
#include "defs.h"
#include "options.h"

class meta_index_t;

class Variant;
class GStore;

class AlleleCode {

 public:

  std::string label;
  int code;
  int cn;
  int len;  // For indels, store length of allele

  AlleleCode() 
    { 
      label = "";
      code = -1;
      cn = -1;
      len = 1; // default is substition
    }
  
  AlleleCode(std::string l, int c, int n=1) : label(l), code(c), cn(n) 
    { }
    
  bool operator<(const AlleleCode & b) const
  {
    return code < b.code ;
  }
  
  friend std::ostream & operator<<( std::ostream & out, AlleleCode & v)
  { 
    out << v.label; // << " (" << v.code << "x" << v.cn << ")";
    return out;
  }

  friend std::ostream & operator<<( std::ostream & out, AlleleCode const & v)
  { 
    out << v.label; // << " (" << v.code << "x" << v.cn << ")";
    return out;
  }

};


class GenotypeSpec {  

 public:

  virtual ~GenotypeSpec() { }
  virtual std::ostream & Print( std::ostream & out ) = 0;

  virtual std::string print() = 0;
  virtual std::string print_reverse_order() = 0 ;

  virtual std::string print_unphased() = 0;

  virtual std::string num_print() = 0;
  virtual std::string num_print_reverse_order() = 0;

  virtual int alleleCount(const int) = 0;
  virtual int alleleCount(const std::string &) = 0;
  virtual int copyCount() = 0;
  virtual bool haploid() = 0;
  virtual bool phased() = 0;  
  virtual bool more() = 0;
  virtual bool pswitch() = 0;
  virtual bool pat() = 0;
  virtual bool mat() = 0;
  virtual std::map<std::string,int> allele_counts() const = 0;
};


inline std::ostream & operator<<( std::ostream & out, GenotypeSpec & gs )
{
  return gs.Print( out );
}


class HaploidSpec : public GenotypeSpec {  
 public:
  HaploidSpec(AlleleCode a1) : a1(a1) { }
  AlleleCode a1;

  std::ostream & Print( std::ostream & out ) 
    {
      out << a1;
      return out;
    }

  std::string print() { return a1.label; }

  std::string print_unphased() { return a1.label; }

  std::string num_print() { return Helper::int2str( a1.code ); } 

  std::string print_reverse_order() { return a1.label; }

  std::string num_print_reverse_order() { return Helper::int2str( a1.code ); }

  int alleleCount( const int a) { return a1.code == a ? ( a1.cn > 0 ? a1.cn : 0 ) : 0 ; }
  int alleleCount( const std::string & a) { return a1.label == a ? ( a1.cn > 0 ? a1.cn : 0 ) : 0 ; }

  bool pat() { return a1.code != 0; }
  bool mat() { return a1.code != 0; }
  int copyCount() { return a1.cn; } 
  bool haploid() { return true; }
  bool phased() { return false; }
  bool pswitch() { return false; }
  bool more() 
    { 
      return a1.code > 1;
    }

  std::map<std::string,int> allele_counts() const 
    {
      std::map<std::string,int> a;
      a[a1.label] += a1.cn;
    }

};


class DiploidSpec : public GenotypeSpec {  

 public:

  DiploidSpec(AlleleCode a1_, AlleleCode a2_, bool p, bool s )
    {
      a1 = a1_;
      a2 = a2_;
      phase = p;
      pswitched = s;
    }
  
  AlleleCode a1;
  AlleleCode a2;
  
  bool phase;
  bool pswitched;

  std::ostream & Print( std::ostream & out )
    {
      out << a1;
      if ( pswitched ) out << "\\";
      else if ( phase ) out << "|";
      else out << "/";
      out << a2; 
      return out; 
    }

  std::string print_unphased() 
    {
      return a2.label > a1.label ? a1.label + "/" + a2.label : a2.label + "/" + a1.label;
    }

  std::string print() 
    { 
      return 
	a1.label 
	+ std::string( pswitched ? "\\" : phase ? "|" : "/" )
	+ a2.label; 
    }

  std::string print_reverse_order() 
    { 
      return
	a2.label 
	+ std::string( pswitched ? "\\" : phase ? "|" : "/" )
	+ a1.label; 
    }

  std::string num_print() 
    { 
      return 
	Helper::int2str( a1.code ) 
	+ std::string( pswitched ? "\\" : phase ? "|" : "/" )
	+ Helper::int2str( a2.code );
    }

  std::string num_print_reverse_order() 
    { 
      return
	Helper::int2str( a2.code )
	+ std::string( pswitched ? "\\" : phase ? "|" : "/" )
	+ Helper::int2str( a1.code );
    }

  int alleleCount( const int a) 
    { 
      int c=0;
      if ( a1.code == a ) c += a1.cn > 0 ? a1.cn : 0;
      if ( a2.code == a ) c += a2.cn > 0 ? a2.cn : 0 ;
      return c;
    }

  int alleleCount( const std::string & a) 
    { 
      int c=0;
      if ( a1.label == a ) c += a1.cn > 0 ? a1.cn : 0;
      if ( a2.label == a ) c += a2.cn > 0 ? a2.cn : 0 ;
      return c;
    }

  bool pat() { return a1.code != 0; }
  bool mat() { return a2.code != 0; }

  int copyCount() { return a1.cn + a2.cn; }
  
  std::map<std::string,int> allele_counts() const 
    {
      std::map<std::string,int> a;
      a[a1.label] += a1.cn;
      a[a2.label] += a2.cn;
    }
  
  bool haploid() { return false; }

  void phased(bool b) { phase = b; }
  bool phased() { return phase; } 

  void pswitch(bool b) { pswitched = b; }
  bool pswitch() { return pswitched; }

  bool more() 
    { 
      return a1.code > 1 || a2.code > 1;
    }

};





class VariantSpec {

  
  // Map an encoding (i.e. string descriptor) to a variant
  // specification

  
  // Genotypes, composed of 1 or 2 alleles:
  // Ordered genotype labels: 0,1,2,3,4..

  std::vector<GenotypeSpec*> genotypes;

  // Map textual genotype to actual genotype
  // (Many-to-1 mapping)

  std::map<std::string,int> genotypeMap;    // 0/1, 0/1, ...
  std::map<std::string,int> genotypeMap_ACGT;  // A/A, A/C, ...

  
  // Genotype meta-information

  std::vector<std::string> slots;
  
  static std::vector<meta_index_t*> * formats;
  static int gt_field;

 public:
  
  VariantSpec(const std::vector<GenotypeSpec*> & g, 
	      const std::map<std::string,int> & gm, 
	      const std::map<std::string,int> & gm_acgt, 
	      const std::vector<std::string> & s)
    {
      genotypes = g;
      genotypeMap = gm;
      genotypeMap_ACGT = gm_acgt;
      slots = s;      
    }

  ~VariantSpec()
    {
      for (unsigned int g=0;  g < genotypes.size(); g++)
	delete genotypes[g];
    }  

  void display();

  static void set_format( int g , std::vector<meta_index_t*> * f ) 
  {
    formats = f;
    gt_field = g;
  }

  std::string printGenotype(int i , bool unphased = false )
    {
      if ( i < 0 ) return ".";
      return unphased ? 
	genotypes[i]->print_unphased() :
	genotypes[i]->print() ;
    }

  int allele_count( const int i, const int a )
  {
    if ( i < 0 || i >= genotypes.size() ) return 0;
    return genotypes[i]->alleleCount(a);
  }

  int allele_count( const int i , const std::string & a )
  {
    if ( i < 0 || i >= genotypes.size() ) return 0;
    return genotypes[i]->alleleCount(a);
  }

  int copy_count( const int i )
  {
    if ( i < 0 || i >= genotypes.size() ) return 0;
    return genotypes[i]->copyCount();
  }

  std::map<std::string,int> allele_counts(int i)
    {
      std::map<std::string,int> a;
      if ( i < 0 ) return a;
      Helper::halt("not implemented spec allele_counts()");
    }


  std::string num_printGenotype( int i )
    {
      if ( i < 0 ) return ".";
      return genotypes[i]->num_print();
    }

  Genotype callGenotype(const std::string & s , const Variant * , bool acgt = false  );

  GenotypeSpec * genotype(int i) 
    {
      return genotypes[i];
    }
  
};


class CachedSpec {
 public:
  VariantSpec * p; // Pointer to specification
  int n;           // Number of times requested
};



class specDecoder {

  // Keep a cache of commonly used variants
  
  std::map<std::string,CachedSpec> cache;
  
  std::string input;
  
  strList tokens;
  
  // Allele labels (order as ref,alt1,alt2,alt3...)
  std::set<AlleleCode> alleles;
  
  bool phased;

  bool haploid;

  GStore * parent;

  // Downcode to numeric allele codes for input
  
  bool numeric_codes;
  
  VariantSpec * decode_from_string(const std::string & );  
  
  VariantSpec * decode_from_string(const std::string & , const std::string & , const std::string & );  


 public:

  specDecoder() { }
  
  void reset()
    {
      alleles.clear();
      phased = haploid = false;
    }


  void setParent(GStore * p) { parent = p; }

  void clearCache()
    {
      std::map<std::string,CachedSpec>::iterator i = cache.begin();
      while ( i != cache.end() )
	{
	  if ( i->second.n == 1 ) 
	    delete i->second.p;	 
	  ++i;
	}
      cache.clear();
    }
  
  void displayCache();

  VariantSpec * decode(std::string s );


};




#endif
