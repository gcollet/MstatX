#include <iostream>
#include <string>
#include <vector>

#include "../src/msa.h"
#include "../src/wentropy.h"
#include "../src/options.h"
#include "test_helpers.h"

namespace {

/* Reuses the alignment from test_jensen.cpp on purpose: getSeqWeights()
 * is shared by wentropy/trident/jensen, and its values for this fixture
 * are already cross-checked there (Python re-implementation of Henikoff &
 * Henikoff, matched against the real binary). Reusing it here means the
 * only new arithmetic to verify is wentropy's own formula. */
const std::string FIXTURE     = "tests/fixtures/jensen_tiny.fasta";
const std::string OUTPUT_FILE = "tests/fixtures/.wentropy_test_output.txt";

/* Same pattern as test_jensen.cpp / test_kabat.cpp: WEntStat::calculate()
 * stores its result in the protected col_stat vector inherited from
 * Stat1D, only reachable through print(). Each test drives the real
 * public pipeline: parse options -> build Msa -> calculate() -> print()
 * -> re-read the file. */
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

/* Fixture (tests/fixtures/jensen_tiny.fasta), 4 sequences x 3 columns:
 *
 *   seq1  AAA
 *   seq2  AAA
 *   seq3  ACA
 *   seq4  AC-
 *
 *   col 0: AAAA  -> fully conserved, no gap
 *   col 1: AACC  -> split 50/50 between A and C, no gap
 *   col 2: AAA-  -> mostly conserved, one gap out of 4
 *
 * Valdar (2002), eq. (50)-(52): t(x) = lambda * (- sum_a p_a log(p_a)),
 * with p_a the Henikoff-weighted frequency of symbol a in column x
 * (gap included as its own symbol, no special-casing), and
 * lambda = 1 / log(min(K, N)). Here K = |alphabet| = 3 ("AC-"),
 * N = 4 sequences, so lambda = 1 / log(3).
 *
 * Expected values come from an independent Python re-implementation of
 * this exact formula (reusing the already cross-checked sequence
 * weights), then verified against `./mstatx -s wentropy` on this
 * fixture: both agree to 5-6 decimal digits: 0 / 0.625299 / 0.579380.
 */
void test_wentropy_nominal_values_on_synthetic_alignment()
{
	parse_test_options(/* global = */ false);
	Msa msa(FIXTURE);

	WEntStat stat;
	stat.calculate(msa);
	stat.print(msa);

	std::vector<float> values = read_col_stat_file(OUTPUT_FILE);
	expect(values.size() == 3, "expected one score per column");

	/* A fully conserved column (single symbol, weight sums to ~1) has
	 * zero Shannon entropy up to floating-point noise from the weight
	 * sum not landing on exactly 1.0 - hence the looser tolerance here
	 * compared to the other two columns. */
	expect(almost_equal(values[0], 0.0f,      1e-3f), "column 0 (fully conserved) should have ~zero weighted entropy");
	expect(almost_equal(values[1], 0.625299f, 1e-4f), "column 1 (50/50 split, no gap)");
	expect(almost_equal(values[2], 0.579380f, 1e-4f), "column 2 (conserved but one gap out of 4)");

	/* The fully conserved column must score strictly lower than either
	 * of the two more divergent columns: entropy is non-negative and
	 * only zero at full conservation, regardless of the exact formula
	 * constants (lambda, weights). */
	expect(values[0] < values[1] && values[0] < values[2],
	       "a fully conserved column should have the lowest weighted entropy");
}

/* Same fixture, through the --global switch: a single score, the
 * unweighted mean of the per-column scores above. */
void test_wentropy_global_mode_is_the_mean_of_column_scores()
{
	parse_test_options(/* global = */ true);
	Msa msa(FIXTURE);

	WEntStat stat;
	stat.calculate(msa);
	stat.print(msa);

	float global_value = read_global_stat_file(OUTPUT_FILE);
	float expected_mean = (0.0f + 0.625299f + 0.579380f) / 3.0f;
	expect(almost_equal(global_value, expected_mean, 1e-4f),
	       "--global should output the mean of the per-column scores");
}

} // namespace

int main()
{
	test_wentropy_nominal_values_on_synthetic_alignment();
	test_wentropy_global_mode_is_the_mean_of_column_scores();
	std::cout << "All wentropy tests passed\n";
	return 0;
}
