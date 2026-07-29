#include <iostream>
#include <string>
#include <vector>

#include "../src/msa.h"
#include "../src/mvector.h"
#include "../src/options.h"
#include "test_helpers.h"

namespace {

/* Reuses the alignment from test_jensen.cpp / test_wentropy.cpp /
 * test_trident.cpp: an alignment with a fully conserved column, a
 * 50/50 split column and a column with one gap is exactly what's
 * needed here too, and it keeps the number of fixture files down. */
const std::string FIXTURE     = "tests/fixtures/jensen_tiny.fasta";
const std::string OUTPUT_FILE = "tests/fixtures/.mvector_test_output.txt";
const std::string MATRIX      = "data/aaindex/HENS920102.mat";

/* HENS920102.mat's alphabet, in the file's own order (see its header
 * line: "M rows = ARNDCQEGHILKMFPSTWYV, cols = ARNDCQEGHILKMFPSTWYV").
 * MVectStat::print() writes one output column per symbol in this exact
 * order, which is also the order Options::Get().matrix_fname's
 * ScoringMatrix reports through getAlphabet() - not something this test
 * needs to query, since it's a fixed property of the matrix file. */
const int A = 0, C = 4, G = 7;
const int ALPHABET_SIZE = 20;

void parse_test_options()
{
	char *argv[] = {
		const_cast<char*>("mstatx"),
		const_cast<char*>("-i"), const_cast<char*>(FIXTURE.c_str()),
		const_cast<char*>("-o"), const_cast<char*>(OUTPUT_FILE.c_str()),
		const_cast<char*>("-m"), const_cast<char*>(MATRIX.c_str())
	};
	Options::Parse(7, argv);
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
 * MVectStat::calculate() computes, for each alignment column, the mean
 * normalized substitution score against every symbol of the scoring
 * matrix's own alphabet: mean_col[a] = (1/N) * sum over non-gap
 * sequences of normScore(a, observed_symbol). Note the divisor is
 * always N (the total sequence count), not the number of non-gap
 * sequences - a gapped position dilutes the mean rather than being
 * excluded from it, which is worth pinning down explicitly (column 2
 * below).
 *
 * Expected values come from the same independent Python matrix parser
 * used for the trident tests (whitespace-tokenized, not the old
 * fixed-width one), applied to this simpler formula (no sequence
 * weighting here, unlike trident/wentropy/jensen). Verified against
 * `./mstatx -s mvector` on this fixture. Three representative symbols
 * are checked per column: A and C (the only residues actually present
 * in this alignment) and G (one of the four rows previously corrupted
 * by the scoring-matrix parsing bug - a cheap extra regression guard
 * for that fix, from mvector's side).
 */
void test_mvector_nominal_values_on_synthetic_alignment()
{
	parse_test_options();
	Msa msa(FIXTURE);

	MVectStat stat;
	stat.calculate(msa);
	stat.print(msa);

	std::vector<std::vector<float> > table = read_mvector_file(OUTPUT_FILE, ALPHABET_SIZE);
	expect(table.size() == 3, "expected one row per column");

	/* Column 0: all 4 sequences carry 'A' -> mean_col[a] = normScore(a, 'A') exactly. */
	expect(almost_equal(table[0][A], 0.545455f, 1e-3f), "col 0, vs A");
	expect(almost_equal(table[0][C], 0.227273f, 1e-3f), "col 0, vs C");
	expect(almost_equal(table[0][G], 0.272727f, 1e-3f), "col 0, vs G");

	/* Column 1: 2 sequences 'A', 2 sequences 'C', no gap -> a straight
	 * average of the two, still divided by N=4. */
	expect(almost_equal(table[1][A], 0.386364f, 1e-3f), "col 1, vs A");
	expect(almost_equal(table[1][C], 0.545455f, 1e-3f), "col 1, vs C");
	expect(almost_equal(table[1][G], 0.181818f, 1e-3f), "col 1, vs G");

	/* Column 2: 3 sequences 'A', 1 gap. The gapped sequence contributes
	 * nothing to the sum, but the divisor stays N=4 (not 3): this is why
	 * table[2][A] is 0.75 * normScore(A,A) rather than normScore(A,A). */
	expect(almost_equal(table[2][A], 0.409091f, 1e-3f), "col 2, vs A (diluted by the gap)");
	expect(almost_equal(table[2][C], 0.170455f, 1e-3f), "col 2, vs C (diluted by the gap)");
	expect(almost_equal(table[2][G], 0.204545f, 1e-3f), "col 2, vs G (diluted by the gap)");
}

} // namespace

int main()
{
	test_mvector_nominal_values_on_synthetic_alignment();
	std::cout << "All mvector tests passed\n";
	return 0;
}
