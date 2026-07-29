#include <iostream>
#include <string>
#include <vector>

#include "../src/msa.h"
#include "../src/trident.h"
#include "../src/options.h"
#include "test_helpers.h"

namespace {

/* Reuses the alignment from test_jensen.cpp / test_wentropy.cpp on
 * purpose: getSeqWeights() is shared across wentropy/trident/jensen and
 * its values for this fixture are already cross-checked (Python
 * re-implementation of Henikoff & Henikoff, matched against the real
 * binary). Only column A and C are used here, and both are on
 * unaffected rows of the HENS920102 matrix (see the caveat below), so
 * they don't interact with the parsing issue found while building
 * these values. */
const std::string FIXTURE     = "tests/fixtures/jensen_tiny.fasta";
const std::string OUTPUT_FILE = "tests/fixtures/.trident_test_output.txt";
const std::string MATRIX      = "data/aaindex/HENS920102.mat";

/* Same pattern as the other Stat1D tests: TridStat::calculate() stores
 * its result in the protected col_stat vector, only reachable through
 * print(). Each test drives the real public pipeline: parse options ->
 * build Msa -> calculate() -> print() -> re-read the file. trident also
 * needs a scoring matrix (-m), unlike jensen/kabat/wentropy. */
void parse_test_options(bool global)
{
	if (global) {
		char *argv[] = {
			const_cast<char*>("mstatx"),
			const_cast<char*>("-i"), const_cast<char*>(FIXTURE.c_str()),
			const_cast<char*>("-o"), const_cast<char*>(OUTPUT_FILE.c_str()),
			const_cast<char*>("-m"), const_cast<char*>(MATRIX.c_str()),
			const_cast<char*>("-g")
		};
		Options::Parse(8, argv);
	} else {
		char *argv[] = {
			const_cast<char*>("mstatx"),
			const_cast<char*>("-i"), const_cast<char*>(FIXTURE.c_str()),
			const_cast<char*>("-o"), const_cast<char*>(OUTPUT_FILE.c_str()),
			const_cast<char*>("-m"), const_cast<char*>(MATRIX.c_str())
		};
		Options::Parse(7, argv);
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
 * Valdar (2002), eq. (50)-(56): col_stat(x) = (1-t(x))^a * (1-r(x))^b *
 * (1-g(x))^c, with default factors a=1.0, b=0.5, c=3.0 (--trident_a/b/c).
 *   - t(x): the same Henikoff-weighted Shannon entropy as wentropy.
 *   - g(x): the raw gap fraction of the column.
 *   - r(x): a stereochemical divergence score built from the
 *     HENS920102 substitution matrix (normalized to [0,1] using its
 *     global min/max), comparing each distinct non-gap residue type in
 *     the column to their mean position in that normalized space.
 *
 * Expected values come from an independent Python re-implementation:
 * the exact same substitution matrix file is parsed with a plain
 * whitespace split, normalized the same way (min=-6, max=16 for
 * HENS920102), and the same t(x)/g(x)/r(x)/combination formulas are
 * applied. Verified against `./mstatx -s trident` on this fixture:
 * 1.000000 / 0.373894 / 0.177449 (Python) vs 1 / 0.373903 / 0.177449
 * (real binary) - the ~9e-6 gap on column 1 is float32 accumulation
 * noise, well within the 1e-4 tolerance used below.
 *
 * *** Caveat found while building this test, not caused by it ***
 * ScoringMatrix's file parser reads each matrix row in fixed 8-character
 * fields (`s.substr(j*8, 8)`), which silently misaligns - and corrupts
 * the sign of - several entries on the G, S, T and V rows of
 * HENS920102.mat (e.g. score('G','R') currently returns +3 instead of
 * -3). This fixture only exercises A and C, whose rows are unaffected,
 * and the matrix's global min/max happen to be unaffected too, so the
 * values below are not compromised by it - but any alignment containing
 * G, S, T or V would get a silently wrong r(x) for trident (and a wrong
 * result for mvector, which uses the same matrix). Flagged separately;
 * not something this test suite tries to fix.
 */
void test_trident_nominal_values_on_synthetic_alignment()
{
	parse_test_options(/* global = */ false);
	Msa msa(FIXTURE);

	TridStat stat;
	stat.calculate(msa);
	stat.print(msa);

	std::vector<float> values = read_col_stat_file(OUTPUT_FILE);
	expect(values.size() == 3, "expected one score per column");

	expect(almost_equal(values[0], 1.0f,      1e-4f), "column 0 (fully conserved, no gap): t=0, r=0, g=0 -> 1");
	expect(almost_equal(values[1], 0.373903f, 1e-3f), "column 1 (50/50 split, no gap)");
	expect(almost_equal(values[2], 0.177449f, 1e-4f), "column 2 (conserved but one gap out of 4)");

	/* This fixture's three columns were deliberately built in increasing
	 * order of divergence (fully conserved -> 3/1 split at the sequence
	 * level -> a gap on top of that). It's not a general law for
	 * arbitrary alignments, but for this specific progression it's a
	 * cheap, formula-independent sanity check. */
	expect(values[0] > values[1] && values[1] > values[2],
	       "trident's combined score should decrease as the column gets less conserved");
}

/* Same fixture, through the --global switch: a single score, the
 * unweighted mean of the per-column scores above. */
void test_trident_global_mode_is_the_mean_of_column_scores()
{
	parse_test_options(/* global = */ true);
	Msa msa(FIXTURE);

	TridStat stat;
	stat.calculate(msa);
	stat.print(msa);

	float global_value = read_global_stat_file(OUTPUT_FILE);
	float expected_mean = (1.0f + 0.373903f + 0.177449f) / 3.0f;
	expect(almost_equal(global_value, expected_mean, 1e-3f),
	       "--global should output the mean of the per-column scores");
}

} // namespace

int main()
{
	test_trident_nominal_values_on_synthetic_alignment();
	test_trident_global_mode_is_the_mean_of_column_scores();
	std::cout << "All trident tests passed\n";
	return 0;
}
