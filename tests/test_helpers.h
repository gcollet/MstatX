/* Shared helpers for tests/test_*.cpp.
 *
 * Extracted from tests/test_msa_scoring.cpp so every new test file
 * (test_jensen.cpp, test_gap.cpp, ...) can reuse the same expect()/
 * almost_equal() instead of redefining them.
 */

#ifndef __TEST_HELPERS_H__
#define __TEST_HELPERS_H__

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace test_helpers {

inline bool almost_equal(float a, float b, float eps = 1e-5f)
{
	return std::fabs(a - b) <= eps;
}

inline void expect(bool condition, const std::string & message)
{
	if (!condition) {
		std::cerr << "FAIL: " << message << "\n";
		std::exit(1);
	}
}

/* Reads back a file written by Stat1D::print() in its default
 * (non-global) mode, i.e. lines of "<column>\t<value>". Returns the
 * values in column order, ignoring the column index itself (the caller
 * already knows it, since it's just 1..ncol). */
inline std::vector<float> read_col_stat_file(const std::string & path)
{
	std::vector<float> values;
	std::ifstream file(path.c_str());
	expect(file.is_open(), "could not open statistic output file: " + path);
	int col;
	float value;
	while (file >> col >> value) {
		values.push_back(value);
	}
	return values;
}

/* Reads back a file written by Stat1D::print() in --global mode, i.e. a
 * single float on its own line. */
inline float read_global_stat_file(const std::string & path)
{
	std::ifstream file(path.c_str());
	expect(file.is_open(), "could not open statistic output file: " + path);
	float value;
	file >> value;
	return value;
}

/* Reads back a file written by MVectStat::print(): a header line (the
 * scoring matrix alphabet, one letter per column, ignored here since the
 * caller already knows the order it asked the matrix for) followed by
 * one line per alignment column: the column index (1-based, also
 * ignored) then one mean score per symbol in the alphabet. Returns the
 * values as [column][symbol_index], in alignment-column order. */
inline std::vector<std::vector<float> > read_mvector_file(const std::string & path, int alphabet_size)
{
	std::vector<std::vector<float> > table;
	std::ifstream file(path.c_str());
	expect(file.is_open(), "could not open statistic output file: " + path);

	std::string header_line;
	std::getline(file, header_line);

	int col;
	while (file >> col) {
		std::vector<float> row(alphabet_size);
		for (int a = 0; a < alphabet_size; ++a) {
			expect(bool(file >> row[a]), "unexpected end of mvector output file");
		}
		table.push_back(row);
	}
	return table;
}

} // namespace test_helpers

using test_helpers::almost_equal;
using test_helpers::expect;
using test_helpers::read_col_stat_file;
using test_helpers::read_global_stat_file;
using test_helpers::read_mvector_file;

#endif
