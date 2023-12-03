#pragma once
#include "Algorithm.h"
#include "LocalMutation.h"

enum class EquiEnergyMovesType {
	ORIGINAL = 0,
	FREQUENT_FALLBACK = 1
};

template<uint32_t TDimension>
class EquiEnergyMovesAlgorithm : public Algorithm<TDimension> {
	struct Chain {
		LocalMutation<TDimension>* mutator;
		Float invTemperature;
		Pcg random;
		int swapAttempts, swaps;
		std::vector<Float> levels;
	};
	std::vector<Float> levels;
	std::vector<Chain> chains;
	Float mLargeStepProb;
	EquiEnergyMovesType mType;
	int mRingCount;
	int movesPossible, movesAttempts;
	Pcg random;
public:
	INLINE EquiEnergyMovesAlgorithm(const Integrand<TDimension>& integrand, std::vector<Float> temperatures, const int ringCount, const EquiEnergyMovesType type) : 
		Algorithm<TDimension>(integrand), mType(type), mRingCount(ringCount) {
		mLargeStepProb = 0.3f;
		std::sort(temperatures.begin(), temperatures.end());
		assert(temperatures.size() > 0 && temperatures[0] == 1.f);
		assert(mRingCount > 1);
		chains.resize(temperatures.size());
		for (uint32_t i = 0; i < chains.size(); ++i) {
			chains[i].invTemperature = 1.f / temperatures[i];
			chains[i].random.reset(1337 * i + 235, 1579 * i + 454);
			chains[i].mutator = new LocalMutation<TDimension>(chains[i].random);
			chains[i].swapAttempts = 0;
			chains[i].swaps = 0;
			computeLevels(chains[i].levels, integrand.maxValue(chains[i].invTemperature));
		}
		random.reset(0xDEAD, 0xDEAD);
		movesPossible = 0;
		movesAttempts = 0;
		computeLevels(levels, integrand.maxValue(Float(1)));
	}

	virtual ~EquiEnergyMovesAlgorithm() override {
		for (uint32_t i = 0; i < chains.size(); ++i) {
			delete chains[i].mutator;
		}
	}

	virtual void run(std::vector<SampleAndPdf<TDimension>>& samples) override {
		for (uint32_t i = 0; i < chains.size(); ++i) {
			do {
				chains[i].mutator->setState(randomVector<TDimension>(chains[i].random));
			} while (value(chains[i].mutator->getState(), i) == Float(0));
			chains[i].mutator->startAdaptation(0.3f);
		}
		std::vector<int> ringHits(mRingCount);
		std::vector<int> chainHits(chains.size());
		for (uint32_t i = 0; i < uint32_t(samples.size()) + MCMC_BURN_PERIOD; ++i) {
			for (int c = int(chains.size()) - 1; c >= 0; --c) {
				const bool localStep = i % 3 != 0;
				Vector<TDimension> proposed;
				if (localStep) {
					proposed = chains[c].mutator->mutateState();
				}
				else {
					proposed = randomVector<TDimension>(chains[c].random);
				}
				if (acceptRatio(chains[c].mutator->getState(), proposed, c) > chains[c].random()) {
					chains[c].mutator->mutationWasAccepted(proposed, localStep);
				}
				else {
					chains[c].mutator->mutationWasRejected(localStep);
				}
				if (c == 0 && i >= MCMC_BURN_PERIOD) {
					const uint32_t index = i - MCMC_BURN_PERIOD;
					samples[index].sample = chains[c].mutator->getState();
					samples[index].pdf = this->mIntegrand.value(samples[index].sample, Float(1));
				}
				if (mType == EquiEnergyMovesType::FREQUENT_FALLBACK) {
					++movesAttempts;
					chainHits.clear();
					int cRing = ringNumber(chains[c].levels, value(chains[c].mutator->getState(), c));
					for (int c2 = int(chains.size()) - 1; c2 >= 0; --c2) {
						if (c != c2) {
							int c2Ring = ringNumber(chains[c].levels, value(chains[c2].mutator->getState(), c));
							if (c2Ring == cRing) {
								chainHits.push_back(c2);
							}
						}
					}
					if (chainHits.size() > 0) {
						++movesPossible;
						const int selectedChainRelative = clamp(int(random() * chainHits.size()), 0, int(chainHits.size() - 1));
						int c2 = chainHits[selectedChainRelative];
						++chains[c].swapAttempts;
						++chains[c2].swapAttempts;
						if (swapRatio(c, c2) > random()) {
							++chains[c].swaps;
							++chains[c2].swaps;
							auto temp = chains[c].mutator->getState();
							chains[c].mutator->setState(chains[c2].mutator->getState());
							chains[c2].mutator->setState(temp);
						}
					}
					/*else {
						// Find first friendly chain
						for (int c2 = c + 1; c2 < chains.size(); ++c2) {
							int c2Ring = ringNumber(chains[c2].levels, this->mIntegrand.value(chains[c2].mutator->getState()));
							bool anyHits = false;
							for (int c3 = int(chains.size()) - 1; c3 >= 0; --c3) {
								if (c3 != c2) {
									int c3Ring = ringNumber(chains[c2].levels, this->mIntegrand.value(chains[c3].mutator->getState()));
									if (c2Ring == c3Ring) {
										anyHits = true;
									}
								}
							}
							if (!anyHits) {
								++chains[c].swapAttempts;
								++chains[c2].swapAttempts;
								if (swapRatio(c, c2) > random()) {
									++chains[c].swaps;
									++chains[c2].swaps;
									auto temp = chains[c].mutator->getState();
									chains[c].mutator->setState(chains[c2].mutator->getState());
									chains[c2].mutator->setState(temp);
								}
							}
						}
					}*/
				}
				
			}
			if (mType == EquiEnergyMovesType::ORIGINAL) {
				++movesAttempts;
				for (int & hits : ringHits) {
					hits = 0;
				}
				for (int c = int(chains.size()) - 1; c >= 0; --c) {
					++ringHits[ringNumber(levels, this->mIntegrand.value(chains[c].mutator->getState(), Float(1)))];
				}
				int usableRingsCount = 0;
				for (int hits : ringHits) {
					usableRingsCount += hits > 1;
				}
				if (usableRingsCount > 0) {
					++movesPossible;
					int selectedRingRelative = clamp(int(random() * usableRingsCount), 0, usableRingsCount - 1);
					int selectedRing = -1;
					for (unsigned int i = 0; i < ringHits.size(); ++i) {
						if (ringHits[i] > 1) {
							if (selectedRingRelative == 0) {
								selectedRing = i;
								break;
							}
							--selectedRingRelative;
						}
					}
					assert(selectedRing != -1);
					const int chainsInRing = ringHits[selectedRing];
					const int selectedChain1Relative = clamp(int(random() * chainsInRing), 0, chainsInRing - 1);
					int selectedChain2Relative = clamp(int(random() * (chainsInRing-1)), 0, chainsInRing - 2);
					if (selectedChain1Relative == selectedChain2Relative) {
						++selectedChain2Relative;
					}
					int currentChainRelative = 0, chain1 = -1, chain2 = -1;
					for (int c = int(chains.size()) - 1; c >= 0; --c) {
						const int ringIndex = ringNumber(levels, this->mIntegrand.value(chains[c].mutator->getState(), Float(1)));
						if (ringIndex == selectedRing) {
							if (currentChainRelative == selectedChain1Relative) {
								chain1 = c;
							}
							if (currentChainRelative == selectedChain2Relative) {
								chain2 = c;
							}
							++currentChainRelative;
						}
					}
					assert(chain1 != -1);
					assert(chain2 != -1);
					assert(chain1 != chain2);
					++chains[chain1].swapAttempts;
					++chains[chain2].swapAttempts;
					if (swapRatio(chain1, chain2) > random()) {
						++chains[chain1].swaps;
						++chains[chain2].swaps;
						auto temp = chains[chain1].mutator->getState();
						chains[chain1].mutator->setState(chains[chain2].mutator->getState());
						chains[chain2].mutator->setState(temp);
					}
				}
			}

		}
	}

	static std::string sName() {
		return "EEM";
	}

	virtual std::string name() const override {
		if (mType == EquiEnergyMovesType::ORIGINAL) {
			return "EEM - Original";
		}
		if (mType == EquiEnergyMovesType::FREQUENT_FALLBACK) {
			return "EEM - Frequent";
		}
		return "Unknown";
	}

	virtual bool hasNormalizedPdf() const override { return false; }

	virtual void printStats() const override {
		for (uint32_t i = 0; i < chains.size(); ++i) {
			std::cout << "#" << (i + 1) << "[ " << (1.f / chains[i].invTemperature) << " ] Acceptance rate: " << Float(100) * chains[i].mutator->acceptanceRate()
				<< " % Swap rate: " << Float(100) * (chains[i].swaps / Float(chains[i].swapAttempts)) << " %" <<  std::endl;
		}
		std::cout << "Moves rate " << Float(100) * (movesPossible / Float(movesAttempts)) << " %" << std::endl;
	}
private:
	INLINE Float acceptRatio(const Vector<TDimension>& current, const Vector<TDimension>& proposed, const uint32_t chainNo) const {
		const Float valCurr = this->mIntegrand.value(current, chains[chainNo].invTemperature);
		const Float valProp = this->mIntegrand.value(proposed, chains[chainNo].invTemperature);
		return valCurr > Float(0) ? std::min(Float(1), valProp / valCurr) : Float(1);
	}

	INLINE Float swapRatio(const uint32_t chain1, const uint32_t chain2) const {
		const auto state1 = chains[chain1].mutator->getState();
		const auto state2 = chains[chain2].mutator->getState();
		const Float val1_t1 = this->mIntegrand.value(state1, chains[chain1].invTemperature);
		const Float val1_t2 = this->mIntegrand.value(state1, chains[chain2].invTemperature);
		const Float val2_t1 = this->mIntegrand.value(state2, chains[chain1].invTemperature);
		const Float val2_t2 = this->mIntegrand.value(state2, chains[chain2].invTemperature);
		const Float factor1 = val1_t2 / val1_t1;
		const Float factor2 = val2_t1 / val2_t2;
		return std::min(Float(1), factor1 * factor2);
	}

	INLINE Float value(const Vector<TDimension>& state, const uint32_t chainNo) const {
		return this->mIntegrand.value(state, chains[chainNo].invTemperature);
	}

	INLINE void computeLevels(std::vector<Float>& levels, const Float maxValue) {
		levels.resize(mRingCount - 1);
		const Float dist = pow(maxValue, 1 / Float(mRingCount));
		Float it(0);
		for (Float & l : levels) {
			it += dist;
			l = it;
		}
	}

	INLINE uint32_t ringNumber(const std::vector<Float>& levels, Float value) {
		for (uint32_t i = 0; i < uint32_t(levels.size()); ++i) {
			if (levels[i] > value) {
				return i;
			}
		}
		return uint32_t(levels.size());
	}
	
};
