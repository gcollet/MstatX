#pragma once

#include <map>
#include <string>

/**
 * BackgroundDistribution holds a background (reference) frequency for
 * each of the 20 standard amino acids, used by statistics that compare
 * a column's observed composition to what would be expected "by
 * chance" (currently jensen; see jensen.cpp).
 *
 * The choice of background distribution is not a neutral implementation
 * detail: the literature shows it measurably affects conservation
 * scores (Capra & Singh, 2007, and follow-up work on background
 * sensitivity). This class exists to make that choice explicit and
 * swappable via -k/--background, instead of a single distribution
 * hardcoded into jensen.cpp.
 *
 * Phase 1 (this class) supports two ways to build one:
 *   - "uniform": 1/20 for each of the 20 standard amino acids.
 *   - a file path: one "<letter> <frequency>" pair per line.
 * The historical Capra & Singh (2007) table jensen.cpp shipped with
 * before this class existed is preserved as the "legacy" spec, so the
 * default --background behaves exactly as before.
 *
 * Named, literature-sourced presets (e.g. a BLOSUM62-derived
 * background) are intentionally left for a later phase: getting their
 * numbers right requires a citable source, not a guess.
 */
class BackgroundDistribution
{
protected:
	std::map<char, float> freq;

public:
	/**
	 * spec is one of:
	 *   "uniform" - 1/20 for each of the 20 standard amino acids
	 *   "legacy"  - the historical Capra & Singh (2007) table
	 *   otherwise - treated as a path to a background frequency file
	 *               (one "<letter> <frequency>" pair per line)
	 * Throws std::runtime_error if a file path can't be opened or is
	 * malformed.
	 */
	explicit BackgroundDistribution(const std::string & spec);

	/**
	 * Returns the background frequency of amino acid aa.
	 * Throws std::runtime_error if aa isn't covered by this
	 * distribution, rather than silently returning 0 - a symbol a
	 * statistic actually needs and can't find in the chosen background
	 * is a configuration error worth surfacing, not a value to guess.
	 */
	float getFreq(char aa) const;
};
