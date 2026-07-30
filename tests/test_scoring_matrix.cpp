#include <iostream>
#include <string>
#include <stdexcept>

#include "../src/scoring_matrix.h"
#include "../src/options.h"
#include "test_helpers.h"

namespace {

void parse_test_options()
{
	char *argv[] = {
		const_cast<char*>("mstatx"),
		const_cast<char*>("-i"), const_cast<char*>("tests/fixtures/jensen_tiny.fasta"),
		const_cast<char*>("-m"), const_cast<char*>("data/aaindex/HENS920102.mat")
	};
	Options::Parse(sizeof(argv) / sizeof(argv[0]), argv);
}

/* Nominal behavior, checked against a tiny, entirely hand-verifiable
 * 3-letter matrix (tests/fixtures/tiny_scoring_matrix.mat):
 *
 *       A     B     C
 *   A   5.0
 *   B   2.0   4.0
 *   C  -1.0   0.0   3.0
 *
 * min = -1.0, max = 5.0, so normScore(x,y) = (score(x,y) + 1) / 6.
 */
void test_scoring_matrix_nominal_values()
{
	parse_test_options();
	ScoringMatrix sm("tests/fixtures/tiny_scoring_matrix.mat");

	expect(sm.isSet(), "scoring matrix should be loaded");
	expect(sm.getAlphabetSize() == 3, "expected a 3-letter alphabet");
	expect(sm.getAlphabet() == "ABC", "alphabet should be ABC, in file order");

	expect(sm.index('A') == 0, "A should map to index 0");
	expect(sm.index('B') == 1, "B should map to index 1");
	expect(sm.index('C') == 2, "C should map to index 2");

	expect(almost_equal(sm.score('A', 'A'), 5.0f), "A-A score should be 5.0");
	expect(almost_equal(sm.score('B', 'A'), 2.0f), "B-A score should be 2.0");
	expect(almost_equal(sm.score('B', 'B'), 4.0f), "B-B score should be 4.0");
	expect(almost_equal(sm.score('C', 'A'), -1.0f), "C-A score should be -1.0");
	expect(almost_equal(sm.score('C', 'B'), 0.0f), "C-B score should be 0.0");
	expect(almost_equal(sm.score('C', 'C'), 3.0f), "C-C score should be 3.0");
	/* score() is symmetric regardless of argument order */
	expect(almost_equal(sm.score('A', 'B'), 2.0f), "score should be symmetric: A-B == B-A");

	expect(almost_equal(sm.getMin(), -1.0f), "matrix min should be -1.0");
	expect(almost_equal(sm.getMax(), 5.0f), "matrix max should be 5.0");
	expect(almost_equal(sm.normScore('A', 'A'), 1.0f, 1e-4f), "normScore(A,A) = (5+1)/6 = 1.0");
	expect(almost_equal(sm.normScore('C', 'A'), 0.0f, 1e-4f), "normScore(C,A) = (-1+1)/6 = 0.0");
	expect(almost_equal(sm.normScore('C', 'B'), 1.0f/6.0f, 1e-4f), "normScore(C,B) = (0+1)/6");
	expect(almost_equal(sm.normScore('B', 'B'), 5.0f/6.0f, 1e-4f), "normScore(B,B) = (4+1)/6");
}

/* Regression test for a parsing bug found while building the trident
 * unit tests: ScoringMatrix used to read each matrix row in fixed
 * 8-character fields (`s.substr(j*8, 8)`). In HENS920102.mat, a value
 * that prints one character shorter than its neighbours (e.g. "0." on
 * the diagonal of G, S, T, V) shifts every following fixed-width field
 * on that row, silently dropping the '-' sign of the next value. The
 * four assertions below are pairs that used to come back positive
 * instead of negative; verified by hand against the raw matrix file. */
void test_scoring_matrix_handles_rows_with_short_leading_values()
{
	parse_test_options();
	ScoringMatrix sm("data/aaindex/HENS920102.mat");

	expect(almost_equal(sm.score('G', 'R'), -3.0f), "G-R score should be -3 (row G starts with the short value '0.')");
	expect(almost_equal(sm.score('S', 'I'), -4.0f), "S-I score should be -4 (row S starts with '2.')");
	expect(almost_equal(sm.score('T', 'H'), -3.0f), "T-H score should be -3 (row T starts with '0.')");
	expect(almost_equal(sm.score('V', 'D'), -5.0f), "V-D score should be -5 (row V starts with '0.')");

	/* The global min/max used to normalize the whole matrix happen to be
	 * unaffected for this particular file, but pin them down too: a
	 * future matrix file might not be so lucky. */
	expect(almost_equal(sm.getMin(), -6.0f), "matrix min should be -6");
	expect(almost_equal(sm.getMax(), 16.0f), "matrix max should be 16");
}

/* An empty filename should throw immediately, not attempt to open "". */
void test_scoring_matrix_empty_filename_throws()
{
	parse_test_options();

	bool threw = false;
	try {
		ScoringMatrix sm("");
	} catch (const std::runtime_error &) {
		threw = true;
	}
	expect(threw, "an empty matrix filename should throw std::runtime_error");
}

/* A path that doesn't exist should throw, not silently produce a
 * zeroed-out or partially-built matrix. */
void test_scoring_matrix_nonexistent_file_throws()
{
	parse_test_options();

	bool threw = false;
	try {
		ScoringMatrix sm("tests/fixtures/does_not_exist.mat");
	} catch (const std::runtime_error &) {
		threw = true;
	}
	expect(threw, "a nonexistent matrix file should throw std::runtime_error");
}

/* Regression test for the throw added alongside the fixed-width parsing
 * fix: a row with fewer values than its position in the alphabet
 * requires (tests/fixtures/malformed_scoring_matrix.mat: row C only
 * gives one of its three expected values) must be reported, not
 * silently read as zeros or leave the rest of the matrix uninitialized. */
void test_scoring_matrix_malformed_row_throws()
{
	parse_test_options();

	bool threw = false;
	try {
		ScoringMatrix sm("tests/fixtures/malformed_scoring_matrix.mat");
	} catch (const std::runtime_error & e) {
		threw = true;
		expect(std::string(e.what()).find('C') != std::string::npos,
		       "error message should name the offending row (C)");
	}
	expect(threw, "a matrix row with too few values should throw std::runtime_error");
}

/* Regression test for a second npos-cast bug found while writing these
 * tests: ScoringMatrix::index() used to compute
 * `int pos = (int) alphabet.find(aa);` - when aa is absent,
 * std::string::find returns std::string::npos, and casting that to int
 * wraps around to -1, which is never >= alphabet.size(): the "symbol not
 * in alphabet" check silently never fired. score('Z', 'A') on the real
 * 20-letter matrix used to return a bogus value (reading matrix[0][-1],
 * undefined behavior) instead of throwing. Fixed by comparing against
 * std::string::npos before any cast - the same class of bug, and the
 * same fix, as Msa::fitToAlphabet's earlier one. */
void test_scoring_matrix_unknown_symbol_throws()
{
	parse_test_options();
	ScoringMatrix sm("tests/fixtures/tiny_scoring_matrix.mat"); // alphabet is just "ABC"

	bool threw = false;
	try {
		sm.index('Z');
	} catch (const std::runtime_error &) {
		threw = true;
	}
	expect(threw, "index() on a symbol outside the alphabet should throw");

	threw = false;
	try {
		sm.score('Z', 'A');
	} catch (const std::runtime_error &) {
		threw = true;
	}
	expect(threw, "score() on a symbol outside the alphabet should throw, not read out of bounds");

	threw = false;
	try {
		sm.normScore('A', 'Z');
	} catch (const std::runtime_error &) {
		threw = true;
	}
	expect(threw, "normScore() on a symbol outside the alphabet should throw, not read out of bounds");
}

} // namespace

int main()
{
	test_scoring_matrix_nominal_values();
	test_scoring_matrix_handles_rows_with_short_leading_values();
	test_scoring_matrix_empty_filename_throws();
	test_scoring_matrix_nonexistent_file_throws();
	test_scoring_matrix_malformed_row_throws();
	test_scoring_matrix_unknown_symbol_throws();
	std::cout << "All scoring_matrix tests passed\n";
	return 0;
}
