#ifndef __REGIONS_H__
#define __REGIONS_H__

#include <string>
#include <iostream>
#include <vector>

#include "meta.h"
#include "helper.h"
#include "variant.h"
#include "refdb.h"

class RefVariant;

//////////////////////////////
// Positon and range helpers

class Position {

 public:

  int chr;
  int bp;
  
  int chromosome() const { return chr; }
  int position() const { return bp; }

  void chromosome(const int c) { chr=c; }
  void position(const int p) { bp = p; }

  Position(int chr, int bp) : chr(chr) , bp(bp) { }

  Position() 
    { 
      chr = bp = 0;
    }

  bool operator< ( const Position & b ) const 
    {
      if ( chr < b.chr ) return true;
      if ( chr > b.chr ) return false;
      return bp  < b.bp;
    }

  bool operator> ( const Position & b ) const 
    {
      if ( chr > b.chr ) return true;
      if ( chr < b.chr ) return false;
      return bp  > b.bp;
    }
  
  bool operator== ( const Position & b ) const 
    {
      return chr == b.chr && bp == b.bp;
    }
  
  bool operator<= ( const Position & b ) const 
    {      
      if ( chr < b.chr ) return true;
      if ( chr > b.chr ) return false;
      return bp <= b.bp;
    }
  
  bool operator>= ( const Position & b ) const 
    {
      if ( chr > b.chr ) return true;
      if ( chr < b.chr ) return false;
      return bp >= b.bp;
    }




};

class Region;

class Subregion {
  
 public:

  Subregion(int chr, int bp1, int bp2)
    {
      id = 0;
      start = Position(chr,bp1);
      stop = Position(chr,bp2);
      name = "-";
      strand = frame = 0;
    }

  Subregion(uint64_t id, int chr, int bp1, int bp2) : id(id)
    {
      start = Position(chr,bp1);
      stop = Position(chr,bp2);
      name = "-";
      strand = frame = 0;
    }
  
  Subregion(uint64_t id, std::string name, int chr, int bp1, int bp2) : id(id), name(name)
    {
      start = Position(chr,bp1);
      stop = Position(chr,bp2);
      strand = frame = 0;
    }

 Subregion(uint64_t id, std::string name, int chr, int bp1, int bp2, int strand, int frame ) 
   : id(id), name(name), strand(strand), frame(frame)
  {
    start = Position(chr,bp1);
    stop = Position(chr,bp2);      
  }
  

  uint64_t         id;
  
  std::string      name;

  Position         start;
  
  Position         stop;
  
  int              strand;
  
  int              frame;
  
  bool CDS() const         { return strand == -1 || strand == 1; } 
  bool exon() const        { return strand == -2 || strand == 2; }
  bool start_codon() const { return strand == -3 || strand == 3; }
  bool stop_codon() const  { return strand == -4 || strand == 4; } 

  MetaInformation<LocMeta>  meta;

  bool overlaps(const Region & b) const;
  
  bool overlaps(const Subregion & b) const
    { 
      return stop >= b.start && start <= b.stop;
    }

  friend std::ostream & operator<<( std::ostream & out, const Subregion & r)
    {
      out << r.name << ":"
	  << Helper::chrCode( r.start.chromosome() ) << ":"  
	  << r.start.position() << ".."  
       	  << r.stop.position();       
      return out;
    }
  
  std::string coordinate() const 
    {
      std::stringstream ss;
      ss << Helper::chrCode( start.chromosome() ) << ":" 
	 << start.position() << ".."
	 << stop.position();
      return ss.str();
    }  

};


class Region {

 public:

  uint64_t          id;

  Position          start;

  Position          stop;

  std::string       name;

  std::string       altname;

  int               group;
  
  std::vector<Subregion> subregion; 

  MetaInformation<LocMeta>   meta;

  bool              gene_model;
  
  Region() 
    { 
      construct(0,0,0,0,"","",0);
    }

  // Initialise from chr1:12345..67890
  Region( const std::string & , bool & );

  Region( const std::string & n )
    {
      construct(0,0,0,0,n,"",0);
    }

  Region(int chr, int bp1)
    {
      construct(0,chr,bp1,bp1,"","",0);
    }

  Region(int chr, int bp1, int bp2)
    {
      construct(0,chr,bp1,bp2,"","",0);
    }

  Region(int chr, int bp1, int bp2, std::string n, int grp)
    {
      construct(0,chr,bp1,bp2,n,"",grp);
    }

  Region(uint64_t i, int chr, int bp1, int bp2, std::string n, int grp)
    {
      construct(i,chr,bp1,bp2,n,"",grp);
    }
  
  
  Region(int chr, int bp1, int bp2, std::string n, std::string an, int grp)
    {
      construct(0,chr,bp1,bp2,n,an,grp);
    }

  Region(uint64_t i, int chr, int bp1, int bp2, std::string n, std::string an, int grp)
    {
      construct(i,chr,bp1,bp2,n,an,grp);
    }
  
  Region(const Subregion & s)
    {
      construct(0,s.start.chromosome(), s.start.position(), s.stop.position(),"","",0);
    }
  
  Region(const RefVariant & rv);


  Region(const Variant & v)
    {
      construct(0,v.chromosome(), v.position(), v.stop(),v.name(),"",0);
    }

  void construct(uint64_t i, int chr, int bp1, int bp2, std::string n, std::string an, int grp)
    {
      id = i;
      start = Position(chr,bp1);
      stop = Position(chr,bp2);
      name = n;
      altname = an;
      group = grp;
      subregion.clear();
      gene_model = false;
    }

  void addSubRegion(int chr, int bp1, int bp2)
    { subregion.push_back( Subregion( chr, bp1, bp2 ) ); }
  
  void addSubRegion(uint64_t id, int chr, int bp1, int bp2)
    { subregion.push_back( Subregion( id, chr, bp1, bp2 ) ); }

  void addSubRegion(uint64_t id, std::string name, int chr, int bp1, int bp2)
    { subregion.push_back( Subregion( id, name, chr, bp1, bp2 ) ); }

  void addSubRegion(uint64_t id, std::string name, int chr, int bp1, int bp2, int strand , int frame )
    { subregion.push_back( Subregion( id, name, chr, bp1, bp2, strand, frame ) ); }
  
  void addSubRegion( const Region & r)
      { 

	  
	  subregion.push_back( Subregion( r.chromosome(), r.start.position(), r.stop.position() ) ) ;

          // Also copy over meta-information from region -> sub-region

	  subregion.back().meta = r.meta;
	  
	  // And specially populate frame/strand information, if it exists, from Region meta-information

	  if ( r.meta.has_field( PLINKSeq::TRANSCRIPT_FRAME() ) )
	      subregion.back().frame = r.meta.get1_int( PLINKSeq::TRANSCRIPT_FRAME() ) ;
	  
	  if ( r.meta.has_field( PLINKSeq::TRANSCRIPT_STRAND() ) )
	      subregion.back().strand = r.meta.get1_int( PLINKSeq::TRANSCRIPT_STRAND() ) ;
	  
      }
  
  bool operator< (const Region& b) const
  { 
    if ( start < b.start ) return true;
    if ( start > b.start ) return false;
    if ( stop < b.stop ) return true;
    if ( stop > b.stop ) return false;
    return name < b.name;
  }
  
  int chromosome() const { return start.chromosome(); }   

  bool contains(const Variant & v) const
    {
      // Just use single point for now (bp1)
      return overlaps( Region( v.chromosome(), v.position(), v.stop() == 0 ? v.position() : v.stop() ) ); 
    }

  bool overlaps(const Region& b) const;

  bool overlaps(const Subregion & b) const
  { 
    return stop >= b.start && start <= b.stop;
  }
  
  bool operator==(const Region& b) const
  { 
    return stop >= b.start && start <= b.stop;
  }
  
  bool after(const Region & b) const
  { 
    return start > b.stop;
  }
  
  bool before(const Region & b) const
  { 
    return stop < b.start;
  }
  


  void collapse(bool storeSubregions = false) 
  {
      for (unsigned int s=0; s < subregion.size(); s++)
	{
	  if ( subregion[s].start.position() < start.position() )
	    start.position( subregion[s].start.position() );
	  if ( subregion[s].stop.position() > stop.position() )
	    stop.position( subregion[s].stop.position() );
	}

      if ( ! storeSubregions )
	subregion.clear();
    }
  
  
  std::string coordinate() const 
    {
      std::stringstream ss;
      ss << Helper::chrCode( start.chromosome() ) << ":" 
	 << start.position() << ".."
	 << stop.position();
      return ss.str();
    }

  friend std::ostream & operator<<( std::ostream & out, const Region & r) 
    { 
      out << r.name << "(" << r.group << "):"
	  << Helper::chrCode( r.start.chromosome() ) << ":"  
	  << r.start.position() << ".."  
       	  << r.stop.position();       

      if ( r.subregion.size() > 0 )
	for (unsigned int s=0; s<r.subregion.size(); s++)
	  out << r.subregion[s] << ";";

      out << "[" << r.meta << "]";

      return out;
    }
  
  bool within( std::set<Region> & s );
  
  int length() const 
  { return stop.position() - start.position() + 1 ; }

};

class NameAndChr {
public:
	NameAndChr(const std::string& nameP, uint64_t chrP) : name(nameP), chr(chrP){}

	bool operator<(const NameAndChr& that) const {
		if (this->name != that.name)
			return this->name < that.name;

		return this->chr < that.chr;
	}

	const std::string name;
	uint64_t chr;
};


namespace RegionHelper {
  
  std::set<Region> region_merge_overlap( const std::set<Region> & r );
  
  double region_span( const std::set<Region> & r);

  std::set<Region> region_remove( const std::set<Region> & , const std::set<Region> & );
  std::set<Region> region_require( const std::set<Region> & , const std::set<Region> & );

}


class OverlapDefinition {
 public:  

};


class RangeIntersector {

 public:

  RangeIntersector(  bool (*f1)( Region & , void * p ) , 
		     void (*f2)( const Region & a , const Region & b , void * p ) )
    {
      f_next = f1;
      f_report = f2;
      finished = false;
    }

  int intersect( const Region & , void * p = NULL );
  
 private:

  // next-getters (i.e assumes access to sorted lists)
  bool (*f_next)( Region & , void * p );
  
  // report when do find intersection
  void (*f_report)( const Region & a , const Region & b , void * p );

  // current pool
  std::set<Region> current;

  // state of input stream from B 
  bool finished;
};


#endif
