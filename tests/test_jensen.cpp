#include <iostream>
#include <string>
#include <vector>

#include "../src/msa.h"
#include "../src/jensen.h"
#include "../src/options.h"
#include "test_helpers.h"

namespace {

const std::string FIXTURE     = "tests/fixtures/jensen_tiny.fasta";
const std::string OUTPUT_FILE = "tests/fixtures/.jensen_test_output.txt";

/* JensenStat::calculate() stores its per-column result in a protected
 * vector (col_stat, inherited from Stat1D) - the only public way to read
 * it back is Stat1D::print(), which writes it to Options::Get().output_fname.
 * So each test drives the exact same public pipeline a real run would use:
 * parse options -> build Msa -> calculate() -> print() -> re-read the file.
 * This exercises the class as a black box, through its real public
 * contract, rather than reaching into its internals. */
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
 *   col 0: fully conserved (AAAA), no gap
 *   col 1: split 50/50 between A and C, no gap
 *   col 2: mostly conserved (AAA-), one gap out of 4 sequences
 *
 * Expected values are NOT hand-derived from the formula alone. They come
 * from an independent Python re-implementation of the exact Capra & Singh
 * (2007) computation in jensen.cpp (same pseudo-count handling, same
 * hardcoded background frequency table, same lambda = 0.5), then
 * cross-checked against `./mstatx -s jensen` on this fixture - both agree
 * to 6 decimal digits: 0.752800 / 0.776281 / 0.640615.
 *
 * This test's job is to catch any accidental change in that computation
 * as the surrounding code keeps getting refactored - not to judge whether
 * the formula itself is the "right" one.
 */
void test_jensen_nominal_values_on_synthetic_alignment()
{
	parse_test_options(/* global = */ false);
	Msa msa(FIXTURE);

	JensenStat stat;
	stat.calculate(msa);
	stat.print(msa);

	std::vector<float> values = read_col_stat_file(OUTPUT_FILE);
	expect(values.size() == 3, "expected one score per column");

	expect(almost_equal(values[0], 0.752800f, 1e-4f), "column 0 (fully conserved, no gap)");
	expect(almost_equal(values[1], 0.776281f, 1e-4f), "column 1 (50/50 split, no gap)");
	expect(almost_equal(values[2], 0.640615f, 1e-4f), "column 2 (conserved but one gap out of 4)");
}

/* Same fixture, but through the --global switch (Stat1D::print's other
 * branch): a single score, the unweighted mean of the per-column scores
 * above. This path is currently exercised by no other test. */
void test_jensen_global_mode_is_the_mean_of_column_scores()
{
	parse_test_options(/* global = */ true);
	Msa msa(FIXTURE);

	JensenStat stat;
	stat.calculate(msa);
	stat.print(msa);

	float global_value = read_global_stat_file(OUTPUT_FILE);
	float expected_mean = (0.752800f + 0.776281f + 0.640615f) / 3.0f;
	expect(almost_equal(global_value, expected_mean, 1e-4f),
	       "--global should output the mean of the per-column scores");
}

/* Same fixture, with -k uniform instead of the default "legacy"
 * background: values come from an independent Python re-computation of
 * the same formula with q(a) = 1/20 for every standard amino acid,
 * cross-checked against `./mstatx -s jensen -k uniform` on this
 * fixture: 0.719283 / 0.786037 / 0.618728. Different numbers from the
 * default on purpose - the whole point of -k/--background is that the
 * choice measurably changes the score. */
void test_jensen_uniform_background_gives_different_values()
{
	char *argv[] = {
		const_cast<char*>("mstatx"),
		const_cast<char*>("-i"), const_cast<char*>(FIXTURE.c_str()),
		const_cast<char*>("-o"), const_cast<char*>(OUTPUT_FILE.c_str()),
		const_cast<char*>("-k"), const_cast<char*>("uniform")
	};
	Options::Parse(sizeof(argv) / sizeof(argv[0]), argv);
	Msa msa(FIXTURE);

	JensenStat stat;
	stat.calculate(msa);
	stat.print(msa);

	std::vector<float> values = read_col_stat_file(OUTPUT_FILE);
	expect(values.size() == 3, "expected one score per column");

	expect(almost_equal(values[0], 0.719283f, 1e-4f), "column 0, uniform background");
	expect(almost_equal(values[1], 0.786037f, 1e-4f), "column 1, uniform background");
	expect(almost_equal(values[2], 0.618728f, 1e-4f), "column 2, uniform background");
}

} // namespace

int main()
{
	test_jensen_nominal_values_on_synthetic_alignment();
	test_jensen_global_mode_is_the_mean_of_column_scores();
	test_jensen_uniform_background_gives_different_values();
	std::cout << "All jensen tests passed\n";
	return 0;
}
