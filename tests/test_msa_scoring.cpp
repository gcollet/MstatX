#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "../src/msa.h"
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
    test_fit_to_alphabet_converts_unknown_symbols_to_gaps();
    std::cout << "All tests passed\n";
    return 0;
}
