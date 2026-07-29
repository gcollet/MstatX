#include <iostream>
#include <string>
#include <vector>

#include "../src/msa.h"
#include "../src/gap.h"
#include "../src/options.h"
#include "test_helpers.h"

namespace {

const std::string FIXTURE     = "tests/fixtures/gap_tiny.fasta";
const std::string OUTPUT_FILE = "tests/fixtures/.gap_test_output.txt";

/* Same pattern as the other Stat1D tests: GapStat::calculate() stores
 * its result in the protected col_stat vector, only reachable through
 * print(). Each test drives the real public pipeline: parse options ->
 * build Msa -> calculate() -> print() -> re-read the file. */
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

/* Fixture (tests/fixtures/gap_tiny.fasta), 4 sequences x 4 columns,
 * built to sweep the full [0, 1] range of the formula in one shot:
 *
 *   seq1  AAA-
 *   seq2  AAA-
 *   seq3  AA--
 *   seq4  A---
 *
 *   col 0: AAAA  -> 0 gap / 4  = 0
 *   col 1: AAA-  -> 1 gap / 4  = 0.25
 *   col 2: AA--  -> 2 gaps / 4 = 0.5
 *   col 3: ----  -> 4 gaps / 4 = 1.0  (all-gap column, an edge case
 *                                      worth its own line: no residue
 *                                      type at all in this column)
 *
 * GapStat::calculate() is the simplest of the six statistics: V(x) =
 * nb_gaps(x) / N, no sequence weighting, no scoring matrix. Simple
 * enough that these values were verified by hand (just counting '-' per
 * column) rather than needing an independent re-implementation, then
 * confirmed against `./mstatx -s gap` on this fixture.
 */
void test_gap_nominal_values_on_synthetic_alignment()
{
	parse_test_options(/* global = */ false);
	Msa msa(FIXTURE);

	GapStat stat;
	stat.calculate(msa);
	stat.print(msa);

	std::vector<float> values = read_col_stat_file(OUTPUT_FILE);
	expect(values.size() == 4, "expected one score per column");

	expect(almost_equal(values[0], 0.0f,  1e-6f), "column 0 (no gap)");
	expect(almost_equal(values[1], 0.25f, 1e-6f), "column 1 (1 gap out of 4)");
	expect(almost_equal(values[2], 0.5f,  1e-6f), "column 2 (2 gaps out of 4)");
	expect(almost_equal(values[3], 1.0f,  1e-6f), "column 3 (all gaps)");
}

/* Same fixture, through the --global switch: a single score, the
 * unweighted mean of the per-column scores above. */
void test_gap_global_mode_is_the_mean_of_column_scores()
{
	parse_test_options(/* global = */ true);
	Msa msa(FIXTURE);

	GapStat stat;
	stat.calculate(msa);
	stat.print(msa);

	float global_value = read_global_stat_file(OUTPUT_FILE);
	float expected_mean = (0.0f + 0.25f + 0.5f + 1.0f) / 4.0f;
	expect(almost_equal(global_value, expected_mean, 1e-6f),
	       "--global should output the mean of the per-column scores");
}

} // namespace

int main()
{
	test_gap_nominal_values_on_synthetic_alignment();
	test_gap_global_mode_is_the_mean_of_column_scores();
	std::cout << "All gap tests passed\n";
	return 0;
}
