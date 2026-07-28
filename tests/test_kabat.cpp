#include <iostream>
#include <string>
#include <vector>

#include "../src/msa.h"
#include "../src/kabat.h"
#include "../src/options.h"
#include "test_helpers.h"

namespace {

const std::string FIXTURE     = "tests/fixtures/kabat_tiny.fasta";
const std::string OUTPUT_FILE = "tests/fixtures/.kabat_test_output.txt";

/* Same pattern as test_jensen.cpp: KabatStat::calculate() stores its
 * result in the protected col_stat vector inherited from Stat1D, only
 * reachable through print(). Each test drives the real public pipeline:
 * parse options -> build Msa -> calculate() -> print() -> re-read the file. */
void parse_test_options(bool global)
{
	if (global) {
		char *argv[] = {
			const_cast<char*>("mstatx"),
			const_cast<char*>("-i"), const_cast<char*>(FIXTURE.c_str()),
			const_cast<char*>("-o"), const_cast<char*>(OUTPUT_FILE.c_str()),
			const_cast<char*>("-g")
		};
		Options::Parse(6, argv);
	} else {
		char *argv[] = {
			const_cast<char*>("mstatx"),
			const_cast<char*>("-i"), const_cast<char*>(FIXTURE.c_str()),
			const_cast<char*>("-o"), const_cast<char*>(OUTPUT_FILE.c_str())
		};
		Options::Parse(5, argv);
	}
}

/* Fixture (tests/fixtures/kabat_tiny.fasta), 4 sequences x 3 columns:
 *
 *   seq1  AAA
 *   seq2  AAA
 *   seq3  AAC
 *   seq4  AC-
 *
 *   col 0: AAAA           -> 1 type,  most frequent count = 4
 *   col 1: AAAC           -> 2 types (A,C), most frequent count = 3 (A)
 *   col 2: AAC-           -> 3 types (A,C,-), most frequent count = 2 (A)
 *
 * Wu & Kabat (1970): V(x) = k / n1, with k = number of distinct symbol
 * types in the column (nb_type, gap counts as its own type) and n1 the
 * occurrence count of the single most frequent symbol. No sequence
 * weighting, no special-casing of gaps: unlike jensen/wentropy/trident,
 * this formula is simple enough to verify entirely by hand:
 *
 *   col 0: k=1, n1=4 -> 1/4     = 0.25
 *   col 1: k=2, n1=3 -> 2/3     = 0.666667
 *   col 2: k=3, n1=2 -> 3/2     = 1.5
 *
 * Cross-checked against `./mstatx -s kabat` on this fixture: matches
 * to 6 decimal digits.
 */
void test_kabat_nominal_values_on_synthetic_alignment()
{
	parse_test_options(/* global = */ false);
	Msa msa(FIXTURE);

	KabatStat stat;
	stat.calculate(msa);
	stat.print(msa);

	std::vector<float> values = read_col_stat_file(OUTPUT_FILE);
	expect(values.size() == 3, "expected one score per column");

	expect(almost_equal(values[0], 0.25f,     1e-4f), "column 0 (fully conserved): k=1, n1=4 -> 0.25");
	expect(almost_equal(values[1], 0.666667f, 1e-4f), "column 1 (3 A / 1 C): k=2, n1=3 -> 0.666667");
	expect(almost_equal(values[2], 1.5f,      1e-4f), "column 2 (2 A / 1 C / 1 gap): k=3, n1=2 -> 1.5");

	/* This fixture's three columns were deliberately built in increasing
	 * order of divergence (fully conserved -> 3/1 split -> 3-way split
	 * with a gap). It's not a general law of V = k/n1 for arbitrary
	 * alignments, but for this specific progression it's a cheap,
	 * formula-independent sanity check: a change that broke the ordering
	 * would be worth noticing even before comparing exact magnitudes. */
	expect(values[0] < values[1] && values[1] < values[2],
	       "V(x) = k/n1 should increase as the column gets less conserved");
}

/* Same fixture, through the --global switch: a single score, the
 * unweighted mean of the per-column scores above. */
void test_kabat_global_mode_is_the_mean_of_column_scores()
{
	parse_test_options(/* global = */ true);
	Msa msa(FIXTURE);

	KabatStat stat;
	stat.calculate(msa);
	stat.print(msa);

	float global_value = read_global_stat_file(OUTPUT_FILE);
	float expected_mean = (0.25f + 0.666667f + 1.5f) / 3.0f;
	expect(almost_equal(global_value, expected_mean, 1e-4f),
	       "--global should output the mean of the per-column scores");
}

} // namespace

int main()
{
	test_kabat_nominal_values_on_synthetic_alignment();
	test_kabat_global_mode_is_the_mean_of_column_scores();
	std::cout << "All kabat tests passed\n";
	return 0;
}
