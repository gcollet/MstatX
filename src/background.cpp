#include "background.h"

#include <fstream>
#include <sstream>
#include <stdexcept>

namespace {

/* The 20 standard amino acids, used by the "uniform" preset. Order
 * doesn't matter (freq is a map), alphabetical for readability. */
const std::string STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY";

void load_uniform(std::map<char, float> & freq)
{
	float p = 1.0f / static_cast<float>(STANDARD_AA.size());
	for (char aa : STANDARD_AA){
		freq[aa] = p;
	}
}

/* The background distribution jensen.cpp hardcoded before this class
 * existed (Capra & Singh, 2007). Kept verbatim as the default so
 * -k/--background with no explicit value reproduces the exact same
 * numbers as before this feature was added. */
void load_legacy(std::map<char, float> & freq)
{
	freq['A'] = 0.073f;
	freq['C'] = 0.025f;
	freq['D'] = 0.050f;
	freq['E'] = 0.061f;
	freq['F'] = 0.042f;
	freq['G'] = 0.072f;
	freq['H'] = 0.023f;
	freq['I'] = 0.053f;
	freq['K'] = 0.064f;
	freq['L'] = 0.089f;
	freq['M'] = 0.023f;
	freq['N'] = 0.043f;
	freq['P'] = 0.052f;
	freq['Q'] = 0.040f;
	freq['R'] = 0.052f;
	freq['S'] = 0.073f;
	freq['T'] = 0.056f;
	freq['V'] = 0.063f;
	freq['W'] = 0.013f;
	freq['Y'] = 0.033f;
}

/* One "<letter> <frequency>" pair per line, whitespace-separated (not
 * fixed-width: see the ScoringMatrix parsing bug this project already
 * hit once with fixed-width fields). Blank lines are skipped. Any
 * other malformed line throws, naming the file and line number. */
void load_from_file(std::map<char, float> & freq, const std::string & path)
{
	std::ifstream file(path.c_str());
	if (!file.is_open()){
		throw std::runtime_error("Cannot open file " + path);
	}

	std::string line;
	int line_no = 0;
	while (std::getline(file, line)){
		++line_no;
		std::istringstream iss(line);
		std::string symbol;
		if (!(iss >> symbol)){
			continue; // blank line
		}
		if (symbol.size() != 1){
			throw std::runtime_error("Malformed background file " + path +
				" at line " + std::to_string(line_no) +
				": expected a single-letter symbol, got \"" + symbol + "\"");
		}
		float value;
		if (!(iss >> value)){
			throw std::runtime_error("Malformed background file " + path +
				" at line " + std::to_string(line_no) + ": missing frequency value");
		}
		if (value < 0.0f){
			throw std::runtime_error("Malformed background file " + path +
				" at line " + std::to_string(line_no) + ": negative frequency");
		}
		freq[symbol[0]] = value;
	}

	if (freq.empty()){
		throw std::runtime_error("Background file " + path + " contains no frequency data");
	}
}

} // namespace

BackgroundDistribution :: BackgroundDistribution(const std::string & spec)
{
	if (spec == "uniform"){
		load_uniform(freq);
	} else if (spec == "legacy"){
		load_legacy(freq);
	} else {
		load_from_file(freq, spec);
	}
}

float
BackgroundDistribution :: getFreq(char aa) const
{
	auto it = freq.find(aa);
	if (it == freq.end()){
		throw std::runtime_error(std::string("symbol ") + aa + " is not covered by the background distribution");
	}
	return it->second;
}
