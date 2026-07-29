#include <iostream>
#include <string>
#include <stdexcept>

#include "../src/options.h"
#include "test_helpers.h"

namespace {

/* Defaults, when only the required -i is given. Everything else should
 * fall back to the values documented in options.h. matrix_fname is only
 * checked by suffix: its default is smat_path + "/HENS920102.mat", and
 * smat_path comes from the SCORE_MAT_PATH environment variable when
 * set, "data/aaindex" otherwise - asserting the whole path would make
 * this test depend on the environment it happens to run in. */
void test_options_defaults()
{
	char *argv[] = {
		const_cast<char*>("mstatx"),
		const_cast<char*>("-i"), const_cast<char*>("tests/fixtures/jensen_tiny.fasta")
	};
	Options::Parse(sizeof(argv) / sizeof(argv[0]), argv);

	const Options & opt = Options::Get();
	expect(opt.input_fname == "tests/fixtures/jensen_tiny.fasta", "input_fname should be the given path");
	expect(opt.output_fname == "output.txt", "default output_fname should be output.txt");
	expect(opt.statistic == "wentropy", "default statistic should be wentropy");
	expect(opt.nb_seq == 500, "default nb_seq should be 500");
	expect(opt.verbose == false, "default verbose should be false");
	expect(opt.global == false, "default global should be false");
	expect(almost_equal(opt.threshold, 0.8f), "default threshold should be 0.8");
	expect(almost_equal(opt.factor_a, 1.0f), "default factor_a should be 1.0");
	expect(almost_equal(opt.factor_b, 0.5f), "default factor_b should be 0.5");
	expect(almost_equal(opt.factor_c, 3.0f), "default factor_c should be 3.0");
	expect(opt.window == 3, "default window should be 3");
	expect(opt.matrix_fname.find("HENS920102.mat") != std::string::npos,
	       "default matrix_fname should point at HENS920102.mat");
}

/* Every switch/value argument, given a non-default value, should come
 * back exactly as passed - the mirror image of the defaults test above. */
void test_options_values_can_be_overridden()
{
	char *argv[] = {
		const_cast<char*>("mstatx"),
		const_cast<char*>("-i"), const_cast<char*>("tests/fixtures/kabat_tiny.fasta"),
		const_cast<char*>("-o"), const_cast<char*>("custom_output.txt"),
		const_cast<char*>("-s"), const_cast<char*>("trident"),
		const_cast<char*>("-m"), const_cast<char*>("data/aaindex/HENS920102.mat"),
		const_cast<char*>("-n"), const_cast<char*>("42"),
		const_cast<char*>("-v"),
		const_cast<char*>("-g"),
		const_cast<char*>("-t"), const_cast<char*>("0.5"),
		const_cast<char*>("-a"), const_cast<char*>("2.0"),
		const_cast<char*>("-b"), const_cast<char*>("1.5"),
		const_cast<char*>("-c"), const_cast<char*>("0.25"),
		const_cast<char*>("-w"), const_cast<char*>("7")
	};
	Options::Parse(sizeof(argv) / sizeof(argv[0]), argv);

	const Options & opt = Options::Get();
	expect(opt.input_fname == "tests/fixtures/kabat_tiny.fasta", "input_fname override");
	expect(opt.output_fname == "custom_output.txt", "output_fname override");
	expect(opt.statistic == "trident", "statistic override");
	expect(opt.matrix_fname == "data/aaindex/HENS920102.mat", "matrix_fname override");
	expect(opt.nb_seq == 42, "nb_seq override");
	expect(opt.verbose == true, "verbose override");
	expect(opt.global == true, "global override");
	expect(almost_equal(opt.threshold, 0.5f), "threshold override");
	expect(almost_equal(opt.factor_a, 2.0f), "factor_a override");
	expect(almost_equal(opt.factor_b, 1.5f), "factor_b override");
	expect(almost_equal(opt.factor_c, 0.25f), "factor_c override");
	expect(opt.window == 7, "window override");
}

/* -i is the one argument declared "needed" with no default: omitting it
 * must throw std::runtime_error, which is what lets main.cpp report a
 * clean error instead of silently running with an empty input path. */
void test_options_missing_required_input_throws()
{
	char *argv[] = {
		const_cast<char*>("mstatx"),
		const_cast<char*>("-o"), const_cast<char*>("some_output.txt")
	};

	bool threw = false;
	try {
		Options::Parse(sizeof(argv) / sizeof(argv[0]), argv);
	} catch (const std::runtime_error & e) {
		threw = true;
		expect(std::string(e.what()).find("-i") != std::string::npos,
		       "error message should mention the missing -i argument");
	}
	expect(threw, "a missing required -i should throw std::runtime_error");
}

/* A flag Options doesn't know about must also be rejected, rather than
 * silently ignored - this is what would catch a typo like --treshold. */
void test_options_unknown_flag_throws()
{
	char *argv[] = {
		const_cast<char*>("mstatx"),
		const_cast<char*>("-i"), const_cast<char*>("tests/fixtures/jensen_tiny.fasta"),
		const_cast<char*>("--not-a-real-flag")
	};

	bool threw = false;
	try {
		Options::Parse(sizeof(argv) / sizeof(argv[0]), argv);
	} catch (const std::runtime_error & e) {
		threw = true;
		expect(std::string(e.what()).find("Unknown flag") != std::string::npos,
		       "error message should call out the unknown flag");
	}
	expect(threw, "an unrecognized flag should throw std::runtime_error");
}

/* Regression test for a real bug found while writing the tests above:
 * Options::Init() used to populate its internal command_line buffer
 * without clearing it first. If a Parse() call threw partway through -
 * e.g. on a missing required argument, before every other flag got a
 * chance to consume its own tokens - the leftover tokens stayed in that
 * buffer and silently leaked into the NEXT Parse() call in the same
 * process, which could then fail (or, worse, succeed with the wrong
 * values) for reasons that had nothing to do with what it was actually
 * given. Fixed with a `command_line.clear()` at the top of Init().
 * This test reproduces the exact scenario: a failed parse followed by a
 * legitimate one, in the same process, and checks that the second call
 * is completely unaffected by the first. */
void test_options_parse_recovers_after_a_previous_failed_parse()
{
	char *bad_argv[] = {
		const_cast<char*>("mstatx"),
		const_cast<char*>("-o"), const_cast<char*>("leftover_output.txt")
	};
	bool threw = false;
	try {
		Options::Parse(sizeof(bad_argv) / sizeof(bad_argv[0]), bad_argv);
	} catch (const std::runtime_error &) {
		threw = true;
	}
	expect(threw, "the deliberately invalid first call should throw");

	char *good_argv[] = {
		const_cast<char*>("mstatx"),
		const_cast<char*>("-i"), const_cast<char*>("tests/fixtures/gap_tiny.fasta"),
		const_cast<char*>("-o"), const_cast<char*>("real_output.txt")
	};
	Options::Parse(sizeof(good_argv) / sizeof(good_argv[0]), good_argv);

	const Options & opt = Options::Get();
	expect(opt.input_fname == "tests/fixtures/gap_tiny.fasta", "input_fname after recovery");
	expect(opt.output_fname == "real_output.txt",
	       "output_fname after recovery should not be polluted by the previous failed call");
}

} // namespace

int main()
{
	test_options_defaults();
	test_options_values_can_be_overridden();
	test_options_missing_required_input_throws();
	test_options_unknown_flag_throws();
	test_options_parse_recovers_after_a_previous_failed_parse();
	std::cout << "All options tests passed\n";
	return 0;
}
