#ifndef __CNV_H__
#define __CNV_H__

#include "pseq.h"
#include "mask.h"

#include <map>
#include <string>
using namespace std;

namespace Pseq {

	namespace VarDB {
		typedef unsigned int uint;

		bool cnv_denovo_scan(Mask& m);

		class SampleEventOracle {
		public:
			SampleEventOracle() {}
			virtual ~SampleEventOracle() {}

			enum EventStatus {
				MISSING = 0,
				HAS_EVENT,
				DOES_NOT_HAVE_EVENT
			};

			virtual EventStatus getEventStatus(string sample);
		};

		class SomeQualityEventOracle : public SampleEventOracle {

		};

		class ChildTransmissionSummary {
		public:
			ChildTransmissionSummary() {
			}
			~ChildTransmissionSummary() {
			}

			void addEvent(const SampleEventOracle* oracle, string child, string father, string mother);

			uint missing() const {
				return YES_MISS + NO_MISS + MISS_YES + MISS_NO + MISS_MISS;
			}

			uint nonMissing() const {
				return YES_YES + YES_NO + NO_YES + NO_NO;
			}

			uint atMostOneHasEvent() const {
				return YES_NO + NO_YES + NO_NO;
			}

			uint bothHaveEvent() const {
				return YES_YES;
			}

			uint neitherHaveEvent() const {
				return NO_NO;
			}

		private:
			// PATERNAL_MATERNAL
			uint YES_YES;
			uint YES_NO;
			uint YES_MISS;

			uint NO_YES;
			uint NO_NO;
			uint NO_MISS;

			uint MISS_YES;
			uint MISS_NO;
			uint MISS_MISS;
		};

	}

}

#endif
