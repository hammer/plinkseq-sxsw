#ifndef __CNV_H__
#define __CNV_H__

#include "plinkseq.h"
#include "plinkseq/mask.h"

#include <map>
#include <string>
#include <vector>
using namespace std;

namespace Pseq {

	namespace VarDB {
		typedef unsigned int uint;
		typedef set<int> Inds;

		bool cnv_denovo_scan(Mask& m);
		void f_cnv_denovo_scan(Variant& v, void* p);

		class EventOracle {
		public:
			EventOracle() {}
			virtual ~EventOracle() {}

			enum EventStatus {
				MISSING = 0,
				HAS_EVENT,
				DOES_NOT_HAVE_EVENT
			};

			virtual EventStatus getEventStatus(int altInd, int indivInd) = 0;
			virtual const vector<EventStatus>* getEventStatus(int altInd) = 0;

			virtual int getNumAltAlleles() const = 0;

			virtual bool discoveredInIndiv(int indivInd) const = 0;
		};

		class LookupEventCounts {
		public:
			LookupEventCounts() : numIndivWithEvent(0), numIndivWithoutEvent(0), numIndivWithMissing(0) {}

			int numIndivWithEvent, numIndivWithoutEvent, numIndivWithMissing;
		};
		LookupEventCounts lookupEventInIndivs(const Inds& lookupInds, const Inds& excludeInds, const vector<EventOracle::EventStatus>& sampStatuses);

		class SomeQualityEventOracle : public EventOracle {
		public:
			SomeQualityEventOracle(const Variant& v, double SQthresh, double NQthresh);
			~SomeQualityEventOracle();

			virtual EventStatus getEventStatus(int altInd, int indivInd) { return (*(getEventStatus(altInd)))[indivInd]; }
			virtual const vector<EventStatus>* getEventStatus(int altInd) { return &(_altToSamples[altInd]); }

			virtual int getNumAltAlleles() const { return _numAltAlleles; }

			virtual bool discoveredInIndiv(int indivInd) const { return _discoveredIndiv.find(indivInd) != _discoveredIndiv.end(); }

		private:
			const int _numAltAlleles;

			typedef vector<EventStatus> SampleToStatus;
			vector<SampleToStatus> _altToSamples;

			Inds _discoveredIndiv;
		};

		class ChildTransmissionSummary {
		public:
			ChildTransmissionSummary()
			: _childCNV(0), _inPaternal(0), _inMaternal(0), _inPaternalAndMaternal(0), _missingPaternal(0), _missingMaternal(0), _denovo(0),
			  _paternal_transmitted(0), _paternal_non_transmitted(0), _paternal_unknown(0),
			  _maternal_transmitted(0), _maternal_non_transmitted(0), _maternal_unknown(0) {}
			~ChildTransmissionSummary() {}

			uint _childCNV;
			uint _inPaternal;
			uint _inMaternal;
			uint _inPaternalAndMaternal;
			uint _missingPaternal;
			uint _missingMaternal;
			uint _denovo;

			uint _paternal_transmitted;
			uint _paternal_non_transmitted;
			uint _paternal_unknown;

			uint _maternal_transmitted;
			uint _maternal_non_transmitted;
			uint _maternal_unknown;
		};

		class AuxCNVdeNovoData {
		public:
			AuxCNVdeNovoData(const vector<double>& p);

			const int _MIN_SQ;
			const int _MIN_NQ;
			const bool _REQUIRE_DE_NOVO_DISCOVERY_IN_CHILD;

			Inds allParentInds;
			Inds allChildrenInds;
			map<string, ChildTransmissionSummary> childrenSummary;
		};

	}

}

#endif
