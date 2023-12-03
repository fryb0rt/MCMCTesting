#pragma once
#include "Algorithm.h"
#include "LocalMutation.h"

template<uint32_t TDimension>
class SampledSwapsAlgorithm : public Algorithm<TDimension> {
	struct Chain {
		LocalMutation<TDimension>* mutator;
		Float invTemperature;
		Pcg random;
		int swapAttempts, swaps, realSwaps;
	};
	std::vector<Chain> chains;
	Float mLargeStepProb;
public:
	INLINE SampledSwapsAlgorithm(const Integrand<TDimension>& integrand, std::vector<Float> temperatures) : Algorithm<TDimension>(integrand) {
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
			chains[i].realSwaps = 0;
		}
	}

	virtual ~SampledSwapsAlgorithm() override {
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
		std::vector<Float> probabilities(chains.size());
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
				++chains[c].swapAttempts;
				Float sum = 0;
				for (int c2 = int(chains.size()) - 1; c2 >= 0; --c2) {
					/*if (c2 == c)
						continue;*/
					probabilities[c2] = sum + value(chains[c].mutator->getState(), c2) * value(chains[c2].mutator->getState(), c);
					sum = probabilities[c2];
				}
				Float rnd = chains[c].random() * sum;
				int selectedChain = int(chains.size()) - 1;
				for (int c2 = int(chains.size()) - 1; c2 >= 0; --c2) {
					/*if (c2 == c)
						continue;*/
					if (probabilities[c2] > rnd) {
						selectedChain = c2;
						break;
					}
				}
				if (selectedChain != c) {
					Float sum2 = 0;
					for (int c2 = int(chains.size()) - 1; c2 >= 0; --c2) {
						/*if (c2 == c)
							continue;*/
						if (c2 == c) {
							sum2 += value(chains[selectedChain].mutator->getState(), c2) * value(chains[selectedChain].mutator->getState(), c);
						}
						else if (c2 == selectedChain) {
							sum2 += value(chains[selectedChain].mutator->getState(), c2) * value(chains[c].mutator->getState(), c);
						}
						else {
							sum2 += value(chains[selectedChain].mutator->getState(), c2) * value(chains[c2].mutator->getState(), c);
						}
					}
					if ((sum / sum2) > chains[c].random()) {
						++chains[c].swaps;
						auto temp = chains[c].mutator->getState();
						chains[c].mutator->setState(chains[selectedChain].mutator->getState());
						chains[selectedChain].mutator->setState(temp);
					}
				}
			}
		}
	}

	static std::string sName() {
		return "SampledSwaps";
	}

	virtual std::string name() const override {
		return sName();
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
}; 
