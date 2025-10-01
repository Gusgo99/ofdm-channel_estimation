OFDM is a technique used to multiplex a fast signal accross many different slower channels. It is widely used for wireless telecommunications and allows for the high data-rate found today in many mainstream technologies, such as Wi-Fi, 4G and 5G.

Among the many advantages of OFDM, one which will be focused in this project is the ability of compensating the time-spreading nature of wireless channels.

Disclaimer: This project focuses on a specific part of the implementation of a digital radio. Many of the lower levels of complexity are not represented neither in the simulations nor in the following text.

# Orthogonal Frequency Division Multiplexing
Orthogonal Frequency Division Multiplexing (OFDM) is a technique that allows for high speed transmissions on time dispersive channels. It adds a preprocessing stage before the symbols are sent over the radio. This preprocessing stage is responsible for spreading a signal accross many channels, comprised of different subcarriers.

## Multiplexing
A OFDM transmission uses an Inverse Discrete Fourier Transform (IDFT) to spread a set of symbols into different frequencies. This allows for the original signal to be recovered via a Discrete Fourier Transform at the output and also allows for the techniques described in the following sections to be applied.

Each channel created by this multiplexing technique is called a subcarrier, that is, a carrier generated before the signal is modulated at the radio. The structure of the DFT guarantees complex orthogonality for the kind of channel described in the following section.

## Time-dispersive channel
A time dispersive channel, the kind of channel OFDM is designed to perform well in, is a channel in which a receiver picks up many copies of the transmitted signal, each with different gains. This is caused by reflections on the environment the system is operating in and will cause a signal to get extended.

The effect is similar to echo that can be heard in an empty room, however, since radios operate in much higher rates than human perception and with much lower signal levels, it becomes a problem even in environments where sound echos cannot be heard.

Mathematically, a time-dispersive channel can be modeled as a set of delays and gains, where the received signal is a linear combination of the `P` previous symbols sent. In signal processing terms, the channel may be described as a finite impulse response (FIR) filter. Therefore, a channel can be described, during some coherence time, as a sequence of gains. The following snippet shows how a receiver signal `y` can be calculated from an input signal `x` and channel gains `h`:

```py
y[n] = sum([x[n - k] * h[k] for k in range(P)])
```

## Cyclic prefix
To compensate for the channel, OFDM uses the property of the Discrete Fourier Transform that applying a time-domain convolution, which is the operation applied by the channel on the signal, is the same as applying a frequency-domain sample-wise multiplication. Namely, for `X = DFT(x)`, `Y = DFT(y)` and `Z = DFT(z)`:

```
Z[i] = X[i] * Y[i]
z = convolve(x, y)
```

However, the DFT is an operation applied to periodic signals, and so this property is only valid for periodic signals. This will hinder the direct application of this property in OFDM systems.

To solve that, it is possible to make the channel "appear" periodic. Since the channel is modeled as a FIR filter with size P, it is only required that the `P` previous samples be periodic. Namely, if the channel filter has size two, the signal `[1, 2, 3, 4]`
will "appear" periodic if the signal `[3, 4, 1, 2, 3, 4]` is sent through the channel.

This technique is called a cyclic prefix and it provides a simple method to revert the effect of the channel.

## Equalisation
As previously presented, the effect of the channel is a FIR filter applied over a, for this specific model, periodic signal. So we may consider a channel response as follows:

```
[1, 0.5]
```

The equaliser to be applied to the received signal can be calculated, for a system with 4 subcarriers, by

```
equaliser = DFT([1, 0.5, 0, 0])
```

The received symbols can then be recovered by dividing each symbol by the respective value in the equaliser.

## Channel estimation
Channel estimation is made by two methods, the first for the initial estimation of the channel and the second for tracking channel variations and improving previous estimations.

### Dense pilot estimation
The first method is done by placing a known symbol in every subcarrier for a time slot. This allow for the maximum amount of data to be extracted from the channel. As an example, a channel `h` sending a symbol `x` and receiving a symbol `y` will have the following behaviour:

```
h = [1, 0.5]
x = [1, 2, 3, 4]
y[0] = x[-1] * h[1] + x[0] * h[0],
y[1] = x[0] * h[1] + x[1] * h[0],
y[2] = x[1] * h[1] + x[2] * h[0],
y[3] = x[2] * h[1] + x[3] * h[0]
```

Since the values of `y` and `x` are known, this represents a linear system where `h[0]` and `h[1]` are unknown values. Since there are 4 equations and 2 unknown variables, it is possible to find the solution providing the Minimum Mean Square Error. The existance of more equations is advantegeous for reducing the effect of noise in the estimation.

This estimation coupled with the technique presented in [Equalisation] allows for the reconstruction of the sent signal.

### Scattered pilot estimation
The second method is done by placing a known symbol at a few select subcarrier of a time slot. This is more succeptible to noise than the previous method, but is a low overhead alternative that allows for tracking changes in the channel over time.

This method will require the received signal `y` to be transformed into the frequency-domain and have it's data symbols removed. Then it may be converted back to the time domain and the same technique applied previously can be performed.

The following snippet shows the values `x_pilots` and `y_pilots`, which can be directly applied to the algorithm previously shown for dense pilots.

```
Y = fft(y)

Y_pilots[pilotPositions] = Y[pilotPositions]
X_pilots[pilotPositions] = pilotValues

y_pilots = ifft(Y_pilots)
x_pilots = ifft(X_pilots)

// y_pilot can now be used for channel estimation
```

### Combining estimations
The estimation for each scaterred pilot time slot can be combined with the previous one to allow for noise better noise compensation and channel variation tracking. A simple first order filter may be used, as follows:

```
estimation[i] = (1 - a) * estimation[i] + a * newEstimation
```

The parameter `a` will determine how succeptible to noise the system will be and how fast it will respond to channel variations.

# Simulation
The code for this project will apply the techniques presented in the previous section to simulate an OFDM system with channel estimation. The channel model is time-invariant and the noise introduced in the system is Additive White Gaussian Noise (AWGN).

To run the simulation, the libraries `numpy`, `scipy` and `matplotlib` are required. Installing them to the environment and running the file `channel_estimation.py` should provide with a graph of the resulting simulation.
