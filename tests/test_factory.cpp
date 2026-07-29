#include <iostream>
#include <string>

#include "../src/statistic.h"
#include "../src/gap.h"
#include "../src/kabat.h"
#include "../src/wentropy.h"
#include "../src/trident.h"
#include "../src/jensen.h"
#include "../src/mvector.h"
#include "test_helpers.h"

namespace {

/* AddAllStatistics() populates Factory<Statistic>'s internal static map
 * (shared by every StatisticFactory::CreateByName call in this process).
 * Calling it more than once is harmless: it just re-registers the same
 * six name -> constructor entries. */
void register_statistics()
{
	AddAllStatistics();
}

/* This is the whole contract "-s <name>" relies on: does each of the six
 * names main.cpp knows about actually build an instance of the matching
 * class, and nothing else? A typo in AddAllStatistics() (e.g. registering
 * "gap" twice and forgetting "kabat") would compile fine and only show up
 * at run time as a wrong or missing statistic - exactly what this test
 * guards against.
 *
 * dynamic_cast is used rather than typeid comparison so each check reads
 * as "is this (at least) a GapStat", matching how the rest of the
 * program only ever uses the base Statistic* pointer. Statistic has
 * virtual functions, so RTTI is available.
 */
void test_factory_creates_the_expected_class_for_each_registered_name()
{
	register_statistics();

	Statistic * gap = StatisticFactory::CreateByName("gap");
	expect(gap != NULL, "\"gap\" should produce an instance");
	expect(dynamic_cast<GapStat *>(gap) != NULL, "\"gap\" should produce a GapStat");
	delete gap;

	Statistic * kabat = StatisticFactory::CreateByName("kabat");
	expect(kabat != NULL, "\"kabat\" should produce an instance");
	expect(dynamic_cast<KabatStat *>(kabat) != NULL, "\"kabat\" should produce a KabatStat");
	delete kabat;

	Statistic * wentropy = StatisticFactory::CreateByName("wentropy");
	expect(wentropy != NULL, "\"wentropy\" should produce an instance");
	expect(dynamic_cast<WEntStat *>(wentropy) != NULL, "\"wentropy\" should produce a WEntStat");
	delete wentropy;

	Statistic * trident = StatisticFactory::CreateByName("trident");
	expect(trident != NULL, "\"trident\" should produce an instance");
	expect(dynamic_cast<TridStat *>(trident) != NULL, "\"trident\" should produce a TridStat");
	delete trident;

	Statistic * jensen = StatisticFactory::CreateByName("jensen");
	expect(jensen != NULL, "\"jensen\" should produce an instance");
	expect(dynamic_cast<JensenStat *>(jensen) != NULL, "\"jensen\" should produce a JensenStat");
	delete jensen;

	Statistic * mvector = StatisticFactory::CreateByName("mvector");
	expect(mvector != NULL, "\"mvector\" should produce an instance");
	expect(dynamic_cast<MVectStat *>(mvector) != NULL, "\"mvector\" should produce a MVectStat");
	delete mvector;
}

/* Two distinct names must never resolve to the same C++ type: this is
 * the specific mistake a copy-pasted StatisticFactory::Add<...>(...)
 * line would cause (e.g. registering GapStat under both "gap" and
 * "kabat"), and the per-name dynamic_cast checks above wouldn't catch it
 * if the wrong class still happened to satisfy some other check. */
void test_factory_registers_six_distinct_classes()
{
	register_statistics();

	Statistic * gap      = StatisticFactory::CreateByName("gap");
	Statistic * kabat     = StatisticFactory::CreateByName("kabat");
	Statistic * wentropy  = StatisticFactory::CreateByName("wentropy");
	Statistic * jensen    = StatisticFactory::CreateByName("jensen");

	expect(dynamic_cast<KabatStat *>(gap) == NULL, "\"gap\" should not produce a KabatStat");
	expect(dynamic_cast<GapStat *>(kabat) == NULL, "\"kabat\" should not produce a GapStat");
	expect(dynamic_cast<JensenStat *>(wentropy) == NULL, "\"wentropy\" should not produce a JensenStat");
	expect(dynamic_cast<WEntStat *>(jensen) == NULL, "\"jensen\" should not produce a WEntStat");

	delete gap;
	delete kabat;
	delete wentropy;
	delete jensen;
}

/* CreateByName's contract for an unregistered name, per factory.h:
 * throw std::runtime_error rather than returning NULL or a default
 * instance. main.cpp relies on this to report "-s some_typo" as an
 * error instead of silently doing nothing. */
void test_factory_throws_for_an_unregistered_name()
{
	register_statistics();

	bool threw = false;
	try {
		StatisticFactory::CreateByName("not_a_real_statistic");
	} catch (const std::runtime_error &) {
		threw = true;
	}
	expect(threw, "creating an unregistered statistic name should throw std::runtime_error");
}

} // namespace

int main()
{
	test_factory_creates_the_expected_class_for_each_registered_name();
	test_factory_registers_six_distinct_classes();
	test_factory_throws_for_an_unregistered_name();
	std::cout << "All factory tests passed\n";
	return 0;
}
