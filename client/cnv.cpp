#include "cnv.h"
#include "util.h"

#include <iostream>
#include <set>
#include <cmath>
#include <sstream>
using namespace std;

extern GStore g;
extern Pseq::Util::Options args;

Pseq::VarDB::SomeQualityEventOracle::SomeQualityEventOracle(const Variant& v, double SQthresh, double NQthresh)
: EventOracle(),
  // We assume the first is always the REF allele:
  _numAltAlleles(v.n_alleles() - 1),
  _altToSamples(_numAltAlleles) {

	const int numIndivs = v.size();
	vector<vector<EventStatus> > indivAltStatus(numIndivs);

	for (int indivInd = 0; indivInd < numIndivs; ++indivInd) {
		const Genotype& gt = v(indivInd);

		if (gt.meta.has_field("DSCVR") && gt.meta.get1_string("DSCVR") == "Y")
			_discoveredIndiv.insert(indivInd);

		indivAltStatus[indivInd] = vector<EventStatus>(_numAltAlleles, MISSING);
		vector<EventStatus>& indivStatuses = indivAltStatus[indivInd];

		if (!gt.meta.has_field("SQ") || !gt.meta.has_field("NQ"))
			continue;
		vector<double> SQs = gt.meta.get_double("SQ");
		vector<double> NQs = gt.meta.get_double("NQ");
		if (SQs.size() != _numAltAlleles || NQs.size() != _numAltAlleles)
			continue;

		for (int altInd = 0; altInd < _numAltAlleles; ++altInd) {
			EventStatus status = MISSING;

			bool highSQ = SQs[altInd] >= SQthresh;
			bool highNQ = NQs[altInd] >= NQthresh;
			if (highSQ && !highNQ)
				status = HAS_EVENT;
			else if (highNQ && !highSQ)
				status = DOES_NOT_HAVE_EVENT;

			indivStatuses[altInd] = status;
		}
	}

	// Copy over the calculations:
	for (int altInd = 0; altInd < _numAltAlleles; ++altInd) {
		_altToSamples[altInd] = SampleToStatus(numIndivs, MISSING);
		SampleToStatus& altStatuses = _altToSamples[altInd];
		for (int indivInd = 0; indivInd < numIndivs; ++indivInd) {
			altStatuses[indivInd] = indivAltStatus[indivInd][altInd];
		}
	}
}

Pseq::VarDB::SomeQualityEventOracle::~SomeQualityEventOracle() {
}

#define CLASS_ALLELE_HEADER \
		<< "#CLASS" \
		<< "\t" << "CHILD" \
		<< "\t" << "LOCUS" \
		<< "\t" << "CNV" \
		<< "\t" << "NUM_TARG" \
		<< "\t" << "PAR_WITH" \
		<< "\t" << "PAR_WITHOUT" \
		<< "\t" << "PAR_MISS" \
		<< "\t" << "CHILD_WITH" \
		<< "\t" << "CHILD_WITHOUT" \
		<< "\t" << "CHILD_MISS"

#define ALLELE_OUTPUT \
		<< "\t" << child->id() \
		<< "\t" << locStr \
		<< "\t" << v.allele(alleleInd).name() \
		<< "\t" << numTargets \
		<< "\t" << parentLookup.numIndivWithEvent \
		<< "\t" << parentLookup.numIndivWithoutEvent \
		<< "\t" << parentLookup.numIndivWithMissing \
		<< "\t" << childLookup.numIndivWithEvent \
		<< "\t" << childLookup.numIndivWithoutEvent \
		<< "\t" << childLookup.numIndivWithMissing

void Pseq::VarDB::f_cnv_denovo_scan(Variant& v, void* p) {
	Out& dnCNVout = Out::stream("denovo.cnv");

	AuxCNVdeNovoData* aux = static_cast<AuxCNVdeNovoData*>(p);
	const Inds& allParentInds = aux->allParentInds;
	const Inds& allChildrenInds = aux->allChildrenInds;
	map<string, ChildTransmissionSummary>& childrenSummary = aux->childrenSummary;

	EventOracle* eventOracle = new SomeQualityEventOracle(v, aux->minSQ, aux->minNQ);

	//string locStr = v.name();
	string locStr = v.coordinate();

	const int numIndivs = v.size();
	for (int childInd = 0; childInd < numIndivs; ++childInd) {
		Individual* child = v.ind(childInd);
		Individual* pat = child->pat();
		Individual* mat = child->mat();

		// only consider individuals where both parents are present:
		if (pat == NULL || mat == NULL)
			continue;

		// TODO : there will be a better way to get the actual individual slot, but use for now...
		const int patInd = g.indmap.ind_n(pat->id());
		const int matInd = g.indmap.ind_n(mat->id());

		int childSampInd = v.ind_sample(childInd);
		if (childSampInd == 0 || v.ind_sample(patInd) != childSampInd || v.ind_sample(matInd) != childSampInd) {
			stringstream str;
			str << "Skipping child " << child->id() << " since it appears in more than one sample or in a different sample than its parents";
			plog.warn(str.str());
			continue;
		}

		//const SampleVariant& childSampVariant = v.sample_metainformation(childSampInd);
		vector<SampleVariant*> fileSamples = v.fsample(childSampInd);
		if (fileSamples.size() != 1)
			continue;
		const SampleVariant& childSampVariant = *(fileSamples[0]);

		int numTargets = -1;
		if (childSampVariant.meta.has_field("NUMT"))
			numTargets = childSampVariant.meta.get1_int("NUMT");
		else if (v.consensus.meta.has_field("NUMT")) // when only a single sample, then info is only in consensus!
			numTargets = v.consensus.meta.get1_int("NUMT");
		else
			continue;

		Inds childsParents;
		childsParents.insert(patInd);
		childsParents.insert(matInd);

		Inds childIndSet;
		childIndSet.insert(childInd);

		ChildTransmissionSummary& childSummary = childrenSummary[child->id()];

		bool lookForDeNovos = !aux->allowDeNovoWithoutDiscoveryInChild || eventOracle->discoveredInIndiv(childInd);
		for (int altInd = 0; altInd < eventOracle->getNumAltAlleles(); ++altInd) {
			// Start from 2nd allele [since we assume the first is always the REF allele]:
			const int alleleInd = altInd + 1;

			const vector<EventOracle::EventStatus>& sampStatuses = *(eventOracle->getEventStatus(altInd));

			LookupEventCounts parentLookup = lookupEventInIndivs(allParentInds, childsParents, sampStatuses);
			LookupEventCounts childLookup = lookupEventInIndivs(allChildrenInds, childIndSet, sampStatuses);

			if (lookForDeNovos) {
				if (sampStatuses[childInd] == EventOracle::HAS_EVENT) {
					++(childSummary._childCNV);

					if (sampStatuses[patInd] == EventOracle::HAS_EVENT)
						++(childSummary._inPaternal);

					if (sampStatuses[matInd] == EventOracle::HAS_EVENT)
						++(childSummary._inMaternal);

					if (sampStatuses[patInd] == EventOracle::HAS_EVENT && sampStatuses[matInd] == EventOracle::HAS_EVENT)
						++(childSummary._inPaternalAndMaternal);

					if (sampStatuses[patInd] == EventOracle::MISSING)
						++(childSummary._missingPaternal);

					if (sampStatuses[matInd] == EventOracle::MISSING)
						++(childSummary._missingMaternal);

					if (sampStatuses[patInd] == EventOracle::DOES_NOT_HAVE_EVENT && sampStatuses[matInd] == EventOracle::DOES_NOT_HAVE_EVENT) {
						//
						// TODO -- If !REQUIRE_DE_NOVO_DISCOVERY_IN_CHILD, then instead of outputting here, cache DENOVO output so that can output ONLY the maximal super-segments of de novo CNVs:
						//
						dnCNVout
						<< "DENOVO"
						ALLELE_OUTPUT
						<< "\n";

						++(childSummary._denovo);
					}
				}
			}

			if (eventOracle->discoveredInIndiv(patInd)) {
				if (sampStatuses[patInd] == EventOracle::HAS_EVENT && sampStatuses[matInd] == EventOracle::DOES_NOT_HAVE_EVENT) {
					if (sampStatuses[childInd] == EventOracle::MISSING) {
						dnCNVout
						<< "PATERNAL_UNKNOWN";
						++(childSummary._paternal_unknown);
					}
					if (sampStatuses[childInd] == EventOracle::HAS_EVENT) {
						dnCNVout
						<< "PATERNAL_TRANSMITTED";
						++(childSummary._paternal_transmitted);
					}
					if (sampStatuses[childInd] == EventOracle::DOES_NOT_HAVE_EVENT) {
						dnCNVout
						<< "PATERNAL_NON_TRANSMITTED";
						++(childSummary._paternal_non_transmitted);
					}

					dnCNVout
					ALLELE_OUTPUT
					<< "\n";
				}
			}

			if (eventOracle->discoveredInIndiv(matInd)) {
				if (sampStatuses[matInd] == EventOracle::HAS_EVENT && sampStatuses[patInd] == EventOracle::DOES_NOT_HAVE_EVENT) {
					if (sampStatuses[childInd] == EventOracle::MISSING) {
						dnCNVout
						<< "MATERNAL_UNKNOWN";
						++(childSummary._maternal_unknown);
					}
					if (sampStatuses[childInd] == EventOracle::HAS_EVENT) {
						dnCNVout
						<< "MATERNAL_TRANSMITTED";
						++(childSummary._maternal_transmitted);
					}
					if (sampStatuses[childInd] == EventOracle::DOES_NOT_HAVE_EVENT) {
						dnCNVout
						<< "MATERNAL_NON_TRANSMITTED";
						++(childSummary._maternal_non_transmitted);
					}

					dnCNVout
					ALLELE_OUTPUT
					<< "\n";
				}
			}
		}
	}

	delete eventOracle;
}

Pseq::VarDB::LookupEventCounts Pseq::VarDB::lookupEventInIndivs(const Inds& lookupInds, const Inds& excludeInds, const vector<EventOracle::EventStatus>& sampStatuses) {
	LookupEventCounts lookup;

	for (Inds::const_iterator indivIt = lookupInds.begin(); indivIt != lookupInds.end(); ++indivIt) {
		const int indiv = *indivIt;
		if (excludeInds.find(indiv) == excludeInds.end()) {
			if (sampStatuses[indiv] == EventOracle::HAS_EVENT)
				lookup.numIndivWithEvent++;
			else if (sampStatuses[indiv] == EventOracle::DOES_NOT_HAVE_EVENT)
				lookup.numIndivWithoutEvent++;
			else if (sampStatuses[indiv] == EventOracle::MISSING)
				lookup.numIndivWithMissing++;
		}
	}

	return lookup;
}

Pseq::VarDB::AuxCNVdeNovoData::AuxCNVdeNovoData(const Pseq::Util::Options& args)
: minSQ(0), minNQ(0), allowDeNovoWithoutDiscoveryInChild(false) {

	SET_INT(minSQ)
	SET_INT(minNQ)
	SET_FLAG(allowDeNovoWithoutDiscoveryInChild)
}

bool Pseq::VarDB::cnv_denovo_scan(Mask& mask) {
	Out& dnCNVout = Out::stream("denovo.cnv");
	Out& dnIndivOut = Out::stream("denovo.cnv.indiv");

	AuxCNVdeNovoData* aux = new AuxCNVdeNovoData(args);

	cerr << "Starting CNV de novo scan..." << endl;
	const int n = g.indmap.size();

	// Attach parents to children:
	for (int i = 0; i < n; ++i) {
		Individual* person = g.indmap(i);
		g.inddb.fetch(person);

		Individual* p = g.indmap.ind(person->father());
		Individual* m = g.indmap.ind(person->mother());

		if (p)
			person->pat(p);
		if (m)
			person->mat(m);

		// Add inds for child and for parents ONLY if individual i has both parents (i.e., is a child), and also initialize empty summary row for each such child:
		if (p && m) {
			aux->allChildrenInds.insert(i);

			aux->allParentInds.insert(g.indmap.ind_n(p->id()));
			aux->allParentInds.insert(g.indmap.ind_n(m->id()));

			aux->childrenSummary[person->id()] = ChildTransmissionSummary();
		}
	}

	dnCNVout
	CLASS_ALLELE_HEADER
	<< "\n";

	g.vardb.iterate(Pseq::VarDB::f_cnv_denovo_scan, aux, mask);

	dnIndivOut
	<< "#SUMMARY"
	<< "\t" << "CHILD"

	<< "\t" << "CHILD_CNV"
	<< "\t" << "IN_PATERNAL"
	<< "\t" << "IN_MATERNAL"
	<< "\t" << "IN_PATERNAL_AND_MATERNAL"
	<< "\t" << "MISSING_PATERNAL"
	<< "\t" << "MISSING_MATERNAL"
	<< "\t" << "DE_NOVO"

	<< "\t" << "PATERNAL_TRANSMITTED"
	<< "\t" << "PATERNAL_NON_TRANSMITTED"
	<< "\t" << "PATERNAL_UNKNOWN"

	<< "\t" << "MATERNAL_TRANSMITTED"
	<< "\t" << "MATERNAL_NON_TRANSMITTED"
	<< "\t" << "MATERNAL_UNKNOWN"

	<< "\n";

	for (map<string, ChildTransmissionSummary>::const_iterator it = aux->childrenSummary.begin(); it != aux->childrenSummary.end(); ++it) {
		const ChildTransmissionSummary& summ = it->second;

		dnIndivOut
		<< "SUMMARY"
		<< "\t" << it->first

		<< "\t" << summ._childCNV
		<< "\t" << summ._inPaternal
		<< "\t" << summ._inMaternal
		<< "\t" << summ._inPaternalAndMaternal
		<< "\t" << summ._missingPaternal
		<< "\t" << summ._missingMaternal
		<< "\t" << summ._denovo

		<< "\t" << summ._paternal_transmitted
		<< "\t" << summ._paternal_non_transmitted
		<< "\t" << summ._paternal_unknown

		<< "\t" << summ._maternal_transmitted
		<< "\t" << summ._maternal_non_transmitted
		<< "\t" << summ._maternal_unknown

		<< "\n";
	}

	delete aux;

	return true;
}
