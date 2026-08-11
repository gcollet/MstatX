#include <iostream>
#include <string>
#include <stdexcept>

#include "../src/background.h"
#include "test_helpers.h"

namespace {

/* "uniform": 1/20 for each of the 20 standard amino acids. Cheap to
 * verify by hand for any of them. */
void test_background_uniform()
{
	BackgroundDistribution bg("uniform");

	expect(almost_equal(bg.getFreq('A'), 1.0f/20.0f, 1e-6f), "A should be 1/20 under uniform");
	expect(almost_equal(bg.getFreq('Y'), 1.0f/20.0f, 1e-6f), "Y should be 1/20 under uniform");
	expect(almost_equal(bg.getFreq('W'), 1.0f/20.0f, 1e-6f), "W should be 1/20 under uniform");
}

/* "legacy": the historical Capra & Singh (2007) table jensen.cpp
 * hardcoded before this class existed. A handful of spot checks against
 * the original literal values (not all 20 - test_jensen.cpp's nominal
 * values test already exercises the full table indirectly through
 * jensen's actual output). */
void test_background_legacy()
{
	BackgroundDistribution bg("legacy");

	expect(almost_equal(bg.getFreq('A'), 0.073f), "A under legacy");
	expect(almost_equal(bg.getFreq('C'), 0.025f), "C under legacy");
	expect(almost_equal(bg.getFreq('W'), 0.013f), "W under legacy");
	expect(almost_equal(bg.getFreq('Y'), 0.033f), "Y under legacy");
}

/* A spec that isn't "uniform" or "legacy" is treated as a file path:
 * one "<letter> <frequency>" pair per line, whitespace-separated. */
void test_background_custom_file()
{
	BackgroundDistribution bg("tests/fixtures/tiny_background.txt");

	expect(almost_equal(bg.getFreq('A'), 0.5f), "A from the custom file");
	expect(almost_equal(bg.getFreq('C'), 0.3f), "C from the custom file");
	expect(almost_equal(bg.getFreq('G'), 0.2f), "G from the custom file");
}

/* A symbol the chosen distribution doesn't cover must throw, not
 * silently return 0 - the old map<char,float>::operator[] behavior this
 * class was built to replace. */
void test_background_unknown_symbol_throws()
{
	BackgroundDistribution bg("tests/fixtures/tiny_background.txt"); // only covers A, C, G

	bool threw = false;
	try {
		bg.getFreq('W');
	} catch (const std::runtime_error &) {
		threw = true;
	}
	expect(threw, "a symbol outside the custom file's coverage should throw");
}

void test_background_nonexistent_file_throws()
{
	bool threw = false;
	try {
		BackgroundDistribution bg("tests/fixtures/does_not_exist_background.txt");
	} catch (const std::runtime_error &) {
		threw = true;
	}
	expect(threw, "a nonexistent background file should throw std::runtime_error");
}

void test_background_malformed_file_throws()
{
	bool threw = false;
	try {
		BackgroundDistribution bg("tests/fixtures/malformed_background.txt");
	} catch (const std::runtime_error & e) {
		threw = true;
		expect(std::string(e.what()).find("line 2") != std::string::npos,
		       "error message should point at the offending line");
	}
	expect(threw, "a malformed background file should throw std::runtime_error");
}

} // namespace

int main()
{
	test_background_uniform();
	test_background_legacy();
	test_background_custom_file();
	test_background_unknown_symbol_throws();
	test_background_nonexistent_file_throws();
	test_background_malformed_file_throws();
	std::cout << "All background tests passed\n";
	return 0;
}
