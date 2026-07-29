#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "../src/msa.h"
#include "../src/scoring_matrix.h"
#include "../src/options.h"
#include "test_helpers.h"

namespace {

void parse_test_options()
{
    char *argv[] = {
        const_cast<char*>("mstatx"),
        const_cast<char*>("-i"),
        const_cast<char*>("tests/fixtures/simple_alignment.fasta"),
        const_cast<char*>("-m"),
        const_cast<char*>("data/aaindex/HENS920102.mat")
    };
    Options::Parse(5, argv);
}

void test_msa_basic_properties()
{
    parse_test_options();
    Msa msa("tests/fixtures/simple_alignment.fasta");

    expect(msa.getNseq() == 3, "expected 3 sequences");
    expect(msa.getNcol() == 4, "expected 4 columns");
    const std::string alphabet = msa.getAlphabet();
    expect(alphabet.size() == 6, "alphabet should contain 6 symbols");
    expect(alphabet.find('A') != std::string::npos, "alphabet should include A");
    expect(alphabet.find('C') != std::string::npos, "alphabet should include C");
    expect(alphabet.find('-') != std::string::npos, "alphabet should include gaps");
    expect(msa.getAaPos('A') == 0, "A should index to 0");
    expect(msa.getAaPos('C') == 1, "C should index to 1");
    expect(msa.getAaPos('-') == 2, "gap should index to 2");
}

void test_msa_gap_and_frequency()
{
    parse_test_options();
    Msa msa("tests/fixtures/simple_alignment.fasta");

    expect(msa.getGap(2) == 1, "third column should have one gap");
    expect(msa.getGap(0) == 0, "first column should have no gap");
    expect(msa.getNtype(0) == 1, "first column should contain one unique symbol");
    expect(msa.getNtype(1) == 2, "second column should contain two unique symbols");
    expect(msa.getFreq('A') > 0.0f, "A frequency should be positive");
    expect(msa.getFreq('A') < 1.0f, "A frequency should be less than 1");
}

void test_msa_seq_weights()
{
    parse_test_options();
    Msa msa("tests/fixtures/simple_alignment.fasta");

    const std::vector<float> & weights = msa.getSeqWeights();
    expect(weights.size() == 3, "weights vector should contain 3 values");
    expect(weights[0] > 0.0f, "first weight should be positive");
    expect(weights[1] > 0.0f, "second weight should be positive");
    expect(weights[2] > 0.0f, "third weight should be positive");
}

void test_scoring_matrix()
{
    parse_test_options();
    ScoringMatrix sm("data/aaindex/HENS920102.mat");

    expect(sm.isSet(), "scoring matrix should be loaded");
    expect(sm.getAlphabetSize() == 20, "expected 20-letter alphabet");
    expect(sm.index('A') == 0, "A should map to index 0");
    expect(almost_equal(sm.score('A', 'A'), 6.0f), "self score should be 6");
    expect(almost_equal(sm.score('A', 'R'), -2.0f), "A-R score should be -2");
    expect(sm.normScore('A', 'A') >= 0.0f, "normalized self score should be non-negative");
    expect(sm.normScore('A', 'A') <= 1.0f, "normalized self score should be at most 1");
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

void test_fit_to_alphabet_converts_unknown_symbols_to_gaps()
{
    parse_test_options();
    Msa msa("tests/fixtures/simple_alignment.fasta");

    msa.fitToAlphabet("AC");

    expect(msa.getSymbol(0, 2) == '-', "symbol outside alphabet should be converted to a gap");
    expect(msa.getSymbol(1, 2) == '-', "symbol outside alphabet should be converted to a gap");
    expect(msa.getSymbol(0, 3) == '-', "symbol outside alphabet should be converted to a gap");
    expect(msa.getAlphabet().find('D') == std::string::npos, "removed symbol should no longer be in the alphabet");
}

} // namespace

int main()
{
    test_msa_basic_properties();
    test_msa_gap_and_frequency();
    test_msa_seq_weights();
    test_scoring_matrix();
    test_scoring_matrix_handles_rows_with_short_leading_values();
    test_fit_to_alphabet_converts_unknown_symbols_to_gaps();
    std::cout << "All tests passed\n";
    return 0;
}
