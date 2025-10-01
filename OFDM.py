import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import convolution_matrix

class OFDM:
	def __init__(self, subcarrierCount: int, cpSize: int, pilotTimeSpacing: int, pilotFreqSpacing: int):
		self.subcarrierCount = subcarrierCount
		self.cpSize = cpSize
		self.pilotTimeSpacing = pilotTimeSpacing
		self.pilotFreqSpacing = pilotFreqSpacing
		
		self.pilotSymbol = (1 + 1j) / np.sqrt(2)
		
		pilotDenseSymbol = np.fft.ifft(np.array([self.pilotSymbol for _ in range(self.subcarrierCount)])) * np.sqrt(self.subcarrierCount)
		doubleDenseSymbols = np.concatenate((pilotDenseSymbol, pilotDenseSymbol))
		convDenseMatrix = convolution_matrix(doubleDenseSymbols[::-1], self.cpSize)
		self.denseSymbolMatrix = convDenseMatrix[self.subcarrierCount - 1:2 * self.subcarrierCount - 1]
		
		freqPilotSparseSymbol = np.zeros_like(pilotDenseSymbol)
		freqPilotSparseSymbol[0::self.pilotFreqSpacing] = self.pilotSymbol
		freqPilotSparseSymbol[-1] = self.pilotSymbol

		pilotSparseSymbol = np.fft.ifft(freqPilotSparseSymbol) * np.sqrt(self.subcarrierCount)
		doubleSparseSymbols = np.concatenate((pilotSparseSymbol, pilotSparseSymbol))
		convSparseMatrix = convolution_matrix(doubleSparseSymbols[::-1], self.cpSize)
		self.sparseSymbolMatrix = convSparseMatrix[self.subcarrierCount - 1:2 * self.subcarrierCount - 1]
		
	def add_cp(self, symbols: np.ndarray) -> np.ndarray:
		return np.array([np.concatenate((sym[-self.cpSize:], sym)) for sym in symbols])
	
	def block_generator(self, symbols: np.ndarray):
		currentIndex = 0
		for _ in range(self.subcarrierCount):
			yield self.pilotSymbol
		while currentIndex < len(symbols):
			for _ in range(self.pilotTimeSpacing - 1):
				for _ in range(self.subcarrierCount):
					if currentIndex < len(symbols):
						yield symbols[currentIndex]
						currentIndex += 1
					else:
						yield self.pilotSymbol
				if currentIndex >= len(symbols): return
			for sc in range(self.subcarrierCount - 1):
				if (sc % self.pilotFreqSpacing) == 0:
					yield self.pilotSymbol
				else:
					if currentIndex < len(symbols):
						yield symbols[currentIndex]
						currentIndex += 1
					else:
						yield self.pilotSymbol
			yield self.pilotSymbol
	
	def generate_block(self, symbols: np.ndarray) -> np.ndarray:
		return np.array([sym for sym in self.block_generator(symbols)])
	
	def encode(self, symbols: np.ndarray) -> np.ndarray:
		block = self.generate_block(symbols)
		separatedSymbols = np.reshape(block, (-1, self.subcarrierCount))
		
		ofdmSymbols = np.fft.ifft(separatedSymbols) * np.sqrt(self.subcarrierCount)
		
		cpSymbols = self.add_cp(ofdmSymbols)
		
		return np.reshape(cpSymbols, (-1,))
	
	def separate_symbols(self, signal):
		return np.reshape(signal, (-1, self.subcarrierCount + self.cpSize))
	
	def remove_cp(self, symbols: np.ndarray):
		return np.array([sym[self.cpSize:] for sym in symbols])
	
	def non_pilot_generator(self, block: np.ndarray):
		for i, line in enumerate(block):
			if i == 0: continue
			if (i % self.pilotTimeSpacing) != 0:
				for sym in line:
					yield sym
			else:
				for j, sym in enumerate(line[:-1]):
					if (j % self.pilotFreqSpacing) != 0:
						yield sym
	
	def remove_pilots(self, block: np.ndarray) -> np.ndarray:
		return np.array(list(self.non_pilot_generator(block)))
	
	def decode(self, signal: np.array, eqAlgorithm: str, maxSize: int = None) -> np.array:
		ofdmSymbols = self.separate_symbols(signal)
		
		nonCpSymbols = self.remove_cp(ofdmSymbols)
		
		rxGroupedSymbols = np.fft.fft(nonCpSymbols) / np.sqrt(self.subcarrierCount)
		
		if eqAlgorithm == "None":
			equalizer = np.ones_like(rxGroupedSymbols[0])
		elif eqAlgorithm == "Naive":
			equalizer = rxGroupedSymbols[0] / self.pilotSymbol
		elif eqAlgorithm == "LMMSE":
			ls, *_ = np.linalg.lstsq(self.denseSymbolMatrix, nonCpSymbols[0])
			
			equalizer = np.fft.fft(ls, n = self.subcarrierCount)
			
		equalizedSymbols = np.zeros_like(rxGroupedSymbols)

		equalizedSymbols[:self.subcarrierCount] = rxGroupedSymbols[:self.subcarrierCount] / equalizer

		for i in range(self.pilotTimeSpacing, len(rxGroupedSymbols), self.pilotTimeSpacing):
			if eqAlgorithm == "LMMSE":
				receivedPilots = np.zeros_like(rxGroupedSymbols[i])
				receivedPilots[0::self.pilotFreqSpacing] = rxGroupedSymbols[i][0::self.pilotFreqSpacing]
				receivedPilots[-1] = rxGroupedSymbols[i][-1]

				pilotConstructedSymbol = np.fft.fft(receivedPilots) / np.sqrt(self.subcarrierCount)

				pls, *_ = np.linalg.lstsq(self.sparseSymbolMatrix, pilotConstructedSymbol)

				equalizerAdjust = np.fft.ifft(pls, n = self.subcarrierCount)
				equalizer = 0.7 * equalizer + 0.3 * equalizerAdjust

			j = i + self.pilotTimeSpacing
			equalizedSymbols[i:j] = rxGroupedSymbols[i:j] / equalizer

		return self.remove_pilots(equalizedSymbols)[:maxSize]
