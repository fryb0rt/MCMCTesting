#pragma once
#include "Algorithm.h"
#include "LocalMutation.h"
#include "minmaxheap.h"

enum class EESType {
	ORIGINAL = 0,
	ADAPTIVE = 1
};

template<uint32_t TDimension>
class AdaptiveEESAlgorithm : public Algorithm<TDimension> {
	struct Sample {
		INLINE Sample() = default;
		Vector<TDimension> state;
		Float value;
		int chainNo;
		bool operator<(const Sample& sample) const {
			return this->value < sample.value;
		}
	};
	struct Chain {
		LocalMutation<TDimension>* mutator;
		Float invTemperature;
		Pcg random;
		int swapAttempts, swaps;
		
		std::vector<Sample> backupRings;
		std::vector<Heap<Sample>> rings;
		std::vector<Float> levels;
	};
	std::vector<Chain> chains;
	Float mLargeStepProb;
	Float mEEJProb;
	EESType mType;
public:
	INLINE AdaptiveEESAlgorithm(const Integrand<TDimension>& integrand, std::vector<Float> temperatures, const int ringCount, const Float eEJProb, const EESType type) : 
		Algorithm<TDimension>(integrand), mEEJProb(eEJProb), mType(type) {
		mLargeStepProb = 0.3f;
		std::sort(temperatures.begin(), temperatures.end());
		assert(temperatures.size() > 0 && temperatures[0] == 1.f);
		chains.resize(temperatures.size());
		for (uint32_t i = 0; i < chains.size(); ++i) {
			chains[i].invTemperature = 1.f / temperatures[i];
			chains[i].random.reset(1337 * i + 235, 1579 * i + 454);
			chains[i].mutator = new LocalMutation<TDimension>(chains[i].random);
			chains[i].swapAttempts = 0;
			chains[i].swaps = 0;
			chains[i].rings.resize(ringCount);
			computeLevels(chains[i].levels, integrand.maxValue(chains[i].invTemperature));
		}
	}

	virtual ~AdaptiveEESAlgorithm() override {
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
			for (Heap<Sample>& ring : chains[i].rings) {
				ring.clear();
			}
			chains[i].backupRings.clear();
		}

		for (int c = int(chains.size()) - 1; c >= 0; --c) {
			for (uint32_t i = 0; i < uint32_t(samples.size()) + MCMC_BURN_PERIOD; ++i) {
				
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
					/*++chains[c].swapAttempts;
					if (c != chains.size() - 1 && swapRatio(c, c + 1) > chains[c].random()) {
						++chains[c].swaps;
						auto temp = chains[c].mutator->getState();
						chains[c].mutator->setState(chains[c + 1].mutator->getState());
						chains[c + 1].mutator->setState(temp);
					}*/
					for (int c2 = int(chains.size()) - 1; c2 >= 0; --c2) {
						if (c2 != c) {
							addToRing(chains[c].mutator->getState(), c2, c);
						}
					}
				
				
					if (c == 0 && i >= MCMC_BURN_PERIOD) {
						const uint32_t index = i - MCMC_BURN_PERIOD;
						samples[index].sample = chains[c].mutator->getState();
						samples[index].pdf = this->mIntegrand.value(samples[index].sample, Float(1));
					}

					if (ringsConstructed(c) && mEEJProb > chains[c].random()) {
						++chains[c].swapAttempts;
						const Float v = value(chains[c].mutator->getState(), c);
						const int ringIndex = mType == EESType::ORIGINAL ? ringNumberLevels(chains[c].levels, v) : findRing(v, c);
						const Heap<Sample>& ring = chains[c].rings[ringIndex];
						const int sampleIndex = clamp(int(chains[c].random() * ring.size()), 0, int(ring.size() - 1));
						const Sample& s = ring.get(sampleIndex);
						if (swapRatio(c, s) > chains[c].random()) {
							++chains[c].swaps;
							chains[c].mutator->setState(s.state);
						}
					}
			}
		}
	}

	static std::string sName() {
		return "AdaptiveEES";
	}

	virtual std::string name() const override {
		return sName();
	}

	virtual bool hasNormalizedPdf() const override { return false; }

	virtual void printStats() const override {
		for (uint32_t i = 0; i < chains.size(); ++i) {
			std::cout << "#" << (i + 1) << "[ " << (1.f / chains[i].invTemperature) << " ] Acceptance rate: " << Float(100) * chains[i].mutator->acceptanceRate()
				<< " % Swap rate: " << Float(100) * (chains[i].swaps / Float(chains[i].swapAttempts)) << " %" << std::endl;
			for (int ringIndex = 0; ringIndex < chains[i].rings.size(); ++ringIndex) {
				std::cout << chains[i].rings[ringIndex].size() << " ";
			}
			std::cout << std::endl;
		}
	}
private:
	INLINE Float acceptRatio(const Vector<TDimension>& current, const Vector<TDimension>& proposed, const uint32_t chainNo) const {
		const Float valCurr = this->mIntegrand.value(current, chains[chainNo].invTemperature);
		const Float valProp = this->mIntegrand.value(proposed, chains[chainNo].invTemperature);
		return valCurr > Float(0) ? std::min(Float(1), valProp / valCurr) : Float(1);
	}

	INLINE Float swapRatio(const uint32_t chain1, const Sample& s) const {
		const auto state1 = chains[chain1].mutator->getState();
		const auto state2 = s.state;
		const Float val1_t1 = this->mIntegrand.value(state1, chains[chain1].invTemperature);
		const Float val1_t2 = this->mIntegrand.value(state1, chains[s.chainNo].invTemperature);
		const Float val2_t1 = this->mIntegrand.value(state2, chains[chain1].invTemperature);
		const Float val2_t2 = this->mIntegrand.value(state2, chains[s.chainNo].invTemperature);
		const Float factor1 = val1_t2 / val1_t1;
		const Float factor2 = val2_t1 / val2_t2;
		return std::min(Float(1), factor1 * factor2);
	}

	INLINE Float value(const Vector<TDimension>& state, const uint32_t chainNo) const {
		return this->mIntegrand.value(state, chains[chainNo].invTemperature);
	}

	INLINE bool ringsConstructed(const uint32_t chainNo) const {
		const std::vector<Heap<Sample>>& rings = chains[chainNo].rings;
		if (mType == EESType::ORIGINAL) {
			for (int ringIndex = 0; ringIndex < rings.size(); ++ringIndex) {
				if (rings[ringIndex].size() == 0) {
					return false;
				}
			}
			return true;
		}
		
		const std::vector<Sample>& backupRings = chains[chainNo].backupRings;
		return backupRings.size() == rings.size();
	}

	INLINE int findRing(const Float value, const uint32_t chainNo) {
		assert(ringsConstructed(chainNo));
		const std::vector<Heap<Sample>>& rings = chains[chainNo].rings;
		int ringIndex = 0;
		for (; ringIndex < rings.size() - 1; ++ringIndex) {
			if (rings[ringIndex].maximum().value >= value) {
				break;
			}
		}
		return ringIndex;
	}

	INLINE void computeLevels(std::vector<Float>& levels, const Float maxValue) {
		levels.resize(chains[0].rings.size()-1);
		const Float dist = pow(maxValue, 1 / Float(chains[0].rings.size()));
		Float it(0);
		for (Float & l : levels) {
			it += dist;
			l = it;
		}
	}

	INLINE uint32_t ringNumberLevels(const std::vector<Float>& levels, Float value) {
		for (uint32_t i = 0; i < uint32_t(levels.size()); ++i) {
			if (levels[i] > value) {
				return i;
			}
		}
		return uint32_t(levels.size());
	}

	INLINE void addToRing(const Vector<TDimension>& state, const uint32_t chainNo, const uint32_t sampleChainNo) {
		std::vector<Heap<Sample>>& rings = chains[chainNo].rings;
		std::vector<Sample>& backupRings = chains[chainNo].backupRings;
		Sample s;
		s.value = value(state, chainNo);
		if (s.value == Float(0)) {
			return;
		}
		s.state = state;
		s.chainNo = sampleChainNo;
		if (mType == EESType::ORIGINAL) {
			rings[ringNumberLevels(chains[chainNo].levels, s.value)].add(s);
			return;
		}
		if (!ringsConstructed(chainNo)) {
			backupRings.push_back(s);
			if (backupRings.size() == rings.size()) {
				// Copy to rings
				std::sort(backupRings.begin(), backupRings.end());
				for (int i = 0; i < rings.size(); ++i) {
					rings[i].add(backupRings[i]);
				}
			}
			return;
		}
		// Rings are already constructed and each should containt at least one sample
		int ringIndex = findRing(s.value, chainNo);
		assert(ringIndex >= 0 && ringIndex < rings.size());
		rings[ringIndex].add(s);
		// Fix the rings if necessary
		Float averageSize(0);
		for (int i = 0; i < rings.size(); ++i) {
			averageSize += rings[i].size();
		}
		averageSize /= rings.size();
		int direction = 0;
		while (rings[ringIndex].size() > (size_t(averageSize) + 5)) {
			int newIndex = -1;
			if (direction == 0) {
				if ((ringIndex + 1) < rings.size()) {
					newIndex = ringIndex + 1;
					direction = 1;
				}
				if (ringIndex > 0 && (newIndex == -1 || rings[newIndex].size() > rings[ringIndex - 1].size())) {
					newIndex = ringIndex - 1;
					direction = -1;
				}
			}
			else {
				newIndex = ringIndex + direction;
				if (newIndex < 0 || newIndex == rings.size()) {
					newIndex = -1;
				}
			}
			if (newIndex == -1 || rings[newIndex].size() >= rings[ringIndex].size()) {
				break; // Nowhere to shift
			}
			if (newIndex < ringIndex) {
				const Sample s = rings[ringIndex].minimum();
				rings[ringIndex].remove_min();
				rings[newIndex].add(s);
			}
			else {
				const Sample s = rings[ringIndex].maximum();
				rings[ringIndex].remove_max();
				rings[newIndex].add(s);
			}
			ringIndex = newIndex;
		}
		for (int i = 0; i < rings.size() - 1; ++i) {
			assert(rings[i].maximum().value <= rings[i + 1].minimum().value);
		}
	}

	INLINE Float swapRatio(const uint32_t chain1, const uint32_t chain2) const {
		const auto state1 = chains[chain1].mutator->getState();
		const auto state2 = chains[chain2].mutator->getState();
		const Float val1_t1 = this->mIntegrand.value(state1, chains[chain1].invTemperature);
		const Float val1_t2 = this->mIntegrand.value(state1, chains[chain2].invTemperature);
		const Float val2_t1 = this->mIntegrand.value(state2, chains[chain1].invTemperature);
		const Float val2_t2 = this->mIntegrand.value(state2, chains[chain2].invTemperature);
		const Float factor1 = val1_t1 / val1_t2;
		const Float factor2 = val2_t1 / val2_t2;
		return std::min(Float(1), factor1 * factor2);
	}
};