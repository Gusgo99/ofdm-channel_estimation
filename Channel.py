import numpy as np

class Channel:
	def __init__(self, rng):
		self.rng = rng

	def transmit(self, txSignal: np.array, dispersionSize: int, snr: float, encodingRate: int) -> np.array:
		if dispersionSize is not None:
			gain = self.rng.normal(scale = 0.25)
			pathGains = np.array([gain ** i for i in range(dispersionSize)])
			
			transformSize = len(txSignal) + dispersionSize
			frequencySignal = np.fft.fft(txSignal, n = transformSize)
			channelResponse = np.fft.fft(pathGains, n = transformSize)
			
			frequencyDispersedSignal = channelResponse * frequencySignal
			
			dispersedSignal = np.fft.ifft(frequencyDispersedSignal)[:len(txSignal)]
		else:
			dispersedSignal = txSignal
			
		if snr is not None:
			normNoise = self.rng.normal(size = len(dispersedSignal)) + 1j * self.rng.normal(size = len(dispersedSignal))
			noise = ((10 ** (-snr / 20)) / (np.sqrt(2 * encodingRate))) * normNoise
		
			rxSignal = dispersedSignal + noise
		else:
			rxSignal = dispersedSignal

		return rxSignal
