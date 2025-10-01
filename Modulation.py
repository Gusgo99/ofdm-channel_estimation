import numpy as np

class QAM16:
	def __init__(self):
		self.symbolCount = 16 # [symbols]
		self.encodingRate = 4 # [bits per symbol]
		symbols = np.array([
			-3 + 3j,
			-3 + 1j,
			-3 - 3j,
			-3 - 1j,
			-1 + 3j,
			-1 + 1j,
			-1 - 3j,
			-1 - 1j,
			+3 + 3j,
			+3 + 1j,
			+3 - 3j,
			+3 - 1j,
			+1 + 3j,
			+1 + 1j,
			+1 - 3j,
			+1 - 1j,
		])

		avPower = np.mean(symbols * np.conj(symbols))
		self.symbolList = symbols / np.sqrt(avPower)
		
	def encode(self, data: np.ndarray) -> np.ndarray:
		return self.symbolList[data]
	
	def decode(self, symbols: np.ndarray) -> np.ndarray:
		symVec = np.reshape(symbols, (-1, 1))
		data = np.zeros(len(symVec), dtype = int)
		for i, s in enumerate(symVec):
			distances = np.abs(self.symbolList - s)
			data[i] = int(np.argmin(distances))
		return np.reshape(data, symbols.shape)
