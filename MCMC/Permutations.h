#pragma once
#include "Algorithm.h"
#include "LocalMutation.h"
#include "PermutationSampler.h"

enum class PermutationsType {
	ALL = 0,
	NON_IDENTITY = 1
};

template<uint32_t TDimension>
class PermutationsAlgorithm : public Algorithm<TDimension> {
	struct Chain {
		LocalMutation<TDimension>* mutator;
		Float invTemperature;
		Pcg random;
		int swapAttempts, swaps;
	};
	std::vector<Chain> chains;
	Float mLargeStepProb;
	PermutationsType mPermutationType;
	PermutationSampler mPermutationSampler;
public:
	INLINE PermutationsAlgorithm(const Integrand<TDimension>& integrand, std::vector<Float> temperatures, const PermutationsType permutationType) : 
		Algorithm<TDimension>(integrand), mPermutationType(permutationType), mPermutationSampler(uint32_t(temperatures.size()), permutationType == PermutationsType::NON_IDENTITY) {
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
		}
	}

	virtual ~PermutationsAlgorithm() override {
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
		std::vector<Vector<TDimension>> currentStatesBackup(chains.size());
		std::vector<Float> values(chains.size() * chains.size()), permutedValues(chains.size() * chains.size());
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
			}
			for (int c1 = int(chains.size()) - 1; c1 >= 0; --c1) {
				Float average = 0;
				for (int c2 = int(chains.size()) - 1; c2 >= 0; --c2) {
					values[c1 * chains.size() + c2] = value(chains[c2].mutator->getState(), c1);
					average += values[c1 * chains.size() + c2];
				}
				average /= chains.size();
				for (int c2 = int(chains.size()) - 1; c2 >= 0; --c2) {
					values[c1 * chains.size() + c2] /= average;
				}
			}
			// DO PERMUTATION SWAP!
			for (auto& chain : chains) {
				++chain.swapAttempts;
			}
			Float nominator;
			const std::vector<uint32_t> proposal = mPermutationSampler.sample([&values, this](uint32_t chainBefore, uint32_t chainAfter) {
				return values[chainBefore * chains.size() + chainAfter];
			}, this->mRandom(), nominator);
			bool accepted = mPermutationType == PermutationsType::ALL; // Automatically accept
			if (mPermutationType == PermutationsType::NON_IDENTITY) {
				// Permute the values
				for (int c1 = int(chains.size()) - 1; c1 >= 0; --c1) {
					for (int c2 = int(chains.size()) - 1; c2 >= 0; --c2) {
						assert(c2 != proposal[c2]);
						permutedValues[c1 * chains.size() + proposal[c2]] = values[c1 * chains.size() + c2];
					}
				}
				const Float denominator = mPermutationSampler.normalization([&permutedValues, this](uint32_t chainBefore, uint32_t chainAfter) {
					return permutedValues[chainBefore * chains.size() + chainAfter];
				});
				if (nominator / denominator >= this->mRandom()) {
					accepted = true;
				}
			}
			if (accepted) {
				// Actually swap the states
				for (int c = int(chains.size()) - 1; c >= 0; --c) {
					currentStatesBackup[c] = chains[c].mutator->getState();
					if (c != proposal[c]) {
						++chains[c].swaps;
					}
				}
				for (int c = int(chains.size()) - 1; c >= 0; --c) {
					chains[c].mutator->setState(currentStatesBackup[proposal[c]]);
				}
			}
		}
	}

	static std::string sName() {
		return "Permutations";
	}

	virtual std::string name() const override {
		if (mPermutationType == PermutationsType::ALL) {
			return "Permutations - ALL";
		}
		if (mPermutationType == PermutationsType::NON_IDENTITY) {
			return "Permutations - SHUFFLE";
		}
		assert(false);
		return "Unknown";
	}

	virtual bool hasNormalizedPdf() const override { return false; }

	virtual void printStats() const override {
		for (uint32_t i = 0; i < chains.size(); ++i) {
			std::cout << "#" << (i + 1) << "[ " << (1.f / chains[i].invTemperature) << " ] Acceptance rate: " << Float(100) * chains[i].mutator->acceptanceRate()
				<< " % Swap rate: " << Float(100) * (chains[i].swaps / Float(chains[i].swapAttempts)) << " %" << std::endl;
		}
	}
private:
	INLINE Float acceptRatio(const Vector<TDimension>& current, const Vector<TDimension>& proposed, const uint32_t chainNo) const {
		const Float valCurr = this->mIntegrand.value(current, chains[chainNo].invTemperature);
		const Float valProp = this->mIntegrand.value(proposed, chains[chainNo].invTemperature);
		return valCurr > Float(0) ? std::min(Float(1), valProp / valCurr) : Float(1);
	}

	INLINE Float value(const Vector<TDimension>& state, const uint32_t chainNo) const {
		return this->mIntegrand.value(state, chains[chainNo].invTemperature);
	}
};
