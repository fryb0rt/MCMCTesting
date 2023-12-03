#include "TestSuite.h"

void scenario1(const std::vector<Float>& temperatures) {
	TestSuite<2> testSuite(50, 1.f, 1.f, 0.00001f, 10.f, 13370);
	//testSuite.addAlgorithm<ReferenceAlgorithm<2>>();
	//testSuite.addAlgorithm<UniformAlgorithm<4>>();
	//testSuite.addAlgorithm<HaltonAlgorithm<4>>();
	//testSuite.addAlgorithm<MetropolisHastingsAlgorithm<2>>();
	//testSuite.addAlgorithm("GlobalAdaptation");
	//testSuite.addAlgorithm<ParallelTemperingAlgorithm<2>>(temperatures);
	//testSuite.addAlgorithm<EquiEnergyMovesAlgorithm<2>>(temperatures, 8, EquiEnergyMovesType::ORIGINAL);
	//testSuite.addAlgorithm<EquiEnergyMovesAlgorithm<2>>(temperatures, 8, EquiEnergyMovesType::FREQUENT_FALLBACK);
	//testSuite.addAlgorithm<PermutationsAlgorithm<2>>(temperatures, PermutationsType::ALL);
	//testSuite.addAlgorithm<PermutationsAlgorithm<2>>(temperatures, PermutationsType::NON_IDENTITY);
	testSuite.addAlgorithm<SampledSwapsAlgorithm<2>>(temperatures);
	//testSuite.addAlgorithm<AdaptiveEESAlgorithm<2>>(temperatures, 12, Float(0.01), EESType::ADAPTIVE);
	//testSuite.addAlgorithm<AdaptiveEESAlgorithm<2>>(temperatures, 16, Float(0.1f), EESType::ADAPTIVE);
	testSuite.runAll(10, 10000);
}

void scenario2(const std::vector<Float>& temperatures) {
	TestSuite<8> testSuite(10, 1.f, 1.f, 0.0001f, 10.f, 13370);
	//testSuite.addAlgorithm<ReferenceAlgorithm<8>>();
	//testSuite.addAlgorithm<UniformAlgorithm<4>>();
	//testSuite.addAlgorithm<HaltonAlgorithm<4>>();
	testSuite.addAlgorithm<MetropolisHastingsAlgorithm<8>>();
	//testSuite.addAlgorithm("GlobalAdaptation");
	testSuite.addAlgorithm<ParallelTemperingAlgorithm<8>>(temperatures);
	testSuite.addAlgorithm<EquiEnergyMovesAlgorithm<8>>(temperatures, 8, EquiEnergyMovesType::ORIGINAL);
	testSuite.addAlgorithm<EquiEnergyMovesAlgorithm<8>>(temperatures, 8, EquiEnergyMovesType::FREQUENT_FALLBACK);
	testSuite.addAlgorithm<PermutationsAlgorithm<8>>(temperatures, PermutationsType::ALL);
	testSuite.addAlgorithm<PermutationsAlgorithm<8>>(temperatures, PermutationsType::NON_IDENTITY);
	testSuite.addAlgorithm<AdaptiveEESAlgorithm<8>>(temperatures, 16, Float(0.1f), EESType::ADAPTIVE);
	testSuite.runAll(100, 10000);
}

template<int TDim>
void scenarioVariable(const std::vector<Float>& temperatures) {
	TestSuite<TDim> testSuite(10, 1.f, 1.f, 0.0001f, 10.f, 13370);
	//testSuite.addAlgorithm<ReferenceAlgorithm<8>>();
	//testSuite.addAlgorithm<UniformAlgorithm<4>>();
	//testSuite.addAlgorithm<HaltonAlgorithm<4>>();
	//testSuite.addAlgorithm<MetropolisHastingsAlgorithm<TDim>>();
	//testSuite.addAlgorithm("GlobalAdaptation");
	//testSuite.addAlgorithm<ParallelTemperingAlgorithm<TDim>>(temperatures);
	//testSuite.addAlgorithm<EquiEnergyMovesAlgorithm<TDim>>(temperatures, 8, EquiEnergyMovesType::ORIGINAL);
	//testSuite.addAlgorithm<EquiEnergyMovesAlgorithm<TDim>>(temperatures, 8, EquiEnergyMovesType::FREQUENT_FALLBACK);
	//testSuite.addAlgorithm<PermutationsAlgorithm<TDim>>(temperatures, PermutationsType::ALL);
	//testSuite.addAlgorithm<PermutationsAlgorithm<TDim>>(temperatures, PermutationsType::NON_IDENTITY);
	testSuite.addAlgorithm<SampledSwapsAlgorithm<TDim>>(temperatures);
	//testSuite.addAlgorithm<AdaptiveEESAlgorithm<TDim>>(temperatures, 16, Float(0.1f), EESType::ADAPTIVE);
	testSuite.runAll(1, 10000);
}

int main(int argc, char ** argv) {
	const uint32_t TMP_COUNT = 8;
	const Float MAX_TEMP = Float(2500);
	const Float DIFF_TEMP = pow(MAX_TEMP, Float(1) / (TMP_COUNT - 1));
	std::vector<Float> temperatures;
	Float t(1);
	for (int i = 0; i < TMP_COUNT; ++i) {
		temperatures.push_back(t);
		t *= DIFF_TEMP;
	}
	//scenario1(temperatures);
	//scenario2(temperatures);
	//scenarioVariable<2>(temperatures);
	
	//scenarioVariable<4>(temperatures);
	
	scenarioVariable<6>(temperatures);
	
	/*scenarioVariable<8>(temperatures);
	
	scenarioVariable<10>(temperatures);

	scenarioVariable<12>(temperatures);

	scenarioVariable<14>(temperatures);*/
	return 0;
}