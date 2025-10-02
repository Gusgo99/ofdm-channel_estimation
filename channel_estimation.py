import numpy as np
import matplotlib.pyplot as plt
import Modulation
import OFDM
import Channel

N = 1024
CPSize = 16

pilotTimeSpacing = 32
pilotFreqSpacing = N // (4 * CPSize)

rng = np.random.default_rng()

ofdm = OFDM.OFDM(N, CPSize, pilotTimeSpacing, pilotFreqSpacing)
channel = Channel.Channel(rng)
modulation = Modulation.QAM16()

def get_ber(txData, rxData):
	errors = np.reshape(txData ^ rxData, (-1,))
	return sum([e.bit_count() for e in errors]) / (len(errors) * modulation.encodingRate)

def simulate_transmission(snr: float, simulatedSymbols: int, apply_dispersion: bool, eqAlgorithm: str):
	txData = rng.integers(0, modulation.symbolCount, simulatedSymbols)
	
	txQamSymbols = modulation.encode(txData)
	
	txSignal = ofdm.encode(txQamSymbols)
	
	rxSignal = channel.transmit(txSignal, CPSize if apply_dispersion else None, snr, modulation.encodingRate)
	
	rxQamSymbols = ofdm.decode(rxSignal, eqAlgorithm, simulatedSymbols)
	
	rxData = modulation.decode(rxQamSymbols)
	
	return get_ber(txData, rxData)

showBerGraphs = True

if showBerGraphs:
	snrList = np.array([i for i in range(0, 12, 2)])
	
	runCount = 100
	symbolCount = 16384

	berIdeal = np.array([np.mean([simulate_transmission(snr, symbolCount, False, "None") for _ in range(runCount)]) for snr in snrList])
	berNoEq = np.array([np.mean([simulate_transmission(snr, symbolCount, True, "None") for _ in range(runCount)]) for snr in snrList])
	berNaiveEq = np.array([np.mean([simulate_transmission(snr, symbolCount, True, "Naive") for _ in range(runCount)]) for snr in snrList])
	berLMMSE = np.array([np.mean([simulate_transmission(snr, symbolCount, True, "LMMSE") for _ in range(runCount)]) for snr in snrList])

	print(berIdeal)
	print(berNoEq)
	print(berNaiveEq)
	print(berLMMSE)

	plt.figure()
	plt.semilogy(snrList, berIdeal, label = "Ideal")
	plt.semilogy(snrList, berNoEq, label = "No equalizer")
	plt.semilogy(snrList, berNaiveEq, label = "Naive equalizer")
	plt.semilogy(snrList, berLMMSE, label = "LMMSE equalizer")
	plt.legend()
	plt.xlabel("$E_b/N0$ (dB)")
	plt.ylabel("BER")
	plt.grid()
else:
	print("Ideal:", simulate_transmission(30, 8192, False, "None"))
	print("Dispersion:", simulate_transmission(30, 8192, True, "LMMSE"))
	
plt.ioff()
plt.show()
