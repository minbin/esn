import pickle
import numpy as np
#import cupy as np

class ESN():
  def __init__(self, n_in, n_fb, n_units, spectral_radius=1.5, noise=None):
    self.n_in = n_in
    self.n_fb = n_fb
    self.n_units = n_units
    self.random = np.random.RandomState()
    self.spectral_radius = spectral_radius
    self.noise = noise
    self._init_res()

  def _init_res(self):
    self.w = self.random.rand(self.n_units, self.n_units) - 0.5
    self.w = self.w * (self.spectral_radius / np.max(np.abs(np.linalg.eigvals(self.w))))
    print(self.w)
    self.w_in = self.random.rand(self.n_units, self.n_in) * 2 - 1
    self.w_fb = self.random.rand(self.n_units, self.n_fb) * 2 - 1

  def _next(self, state, input, feedback):
    pre = (np.dot(self.w_in, input) + np.dot(self.w, state) + np.dot(self.w_fb, feedback))
    if self.noise:
      return (np.tanh(pre) + (self.noise * (self.random.rand(self.n_units)  - 0.5)))
    else:
      return np.tanh(pre)

  def fit(self, input, feedback):
    input = np.reshape(input, (len(input), -1))
    feedback = np.reshape(feedback, (len(input), -1))
    states = np.zeros((input.shape[0], self.n_units))
    for n in range(1, input.shape[0]):
      states[n] = self._next(states[n-1], input[n], feedback[n-1])

    xstates = np.hstack((states, input))
    self.w_out = np.dot(np.linalg.pinv(xstates), feedback).T

    self.last = [input[-1], states[-1], feedback[-1]]

  def predict(self, input, cont=True):
    input = np.reshape(input, (len(input), -1))
    n_count = input.shape[0]

    last = self.last if cont else [np.zeros(self.n_units), np.zeros(self.n_in), np.zeros(self.n_fb)]

    input = np.vstack([last[0], input])
    states = np.vstack([last[1], np.zeros((n_count, self.n_units))])
    output = np.vstack([last[2], np.zeros((n_count, self.n_fb))])

    for n in range(n_count):
      states[n+1] = self._next(states[n], input[n+1], output[n])
      output[n+1] = np.dot(self.w_out, np.concatenate([states[n+1], input[n+1]]))

    return output[1:]

  def save(self):
    with open('esn.npy', 'wb') as f:
      np.save(f, self.w_in)
      np.save(f, self.w)
      np.save(f, self.w_fb)
