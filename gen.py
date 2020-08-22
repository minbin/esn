import argparse
import numpy as np

from matplotlib import pyplot as plt
from esn import ESN

parser = argparse.ArgumentParser()
parser.add_argument("-t", "--data-type", type=str, default="ig", dest="t", help="data type (mg/dnf/ig)")
parser.add_argument("-i", "--input-file", type=str, dest="i", help="spectral radius")
parser.add_argument("-s", "--spectral-radius", type=float, default=1.5, dest="s", help="spectral radius")
parser.add_argument("-u", "--units", type=int, default=4096, dest="u", help="number of reservoir units")
parser.add_argument("-tl", "--trainlen", type=int, default=1000, dest="trainlen", help="number of reservoir units")
parser.add_argument("-p", "--predlen", type=int, default=1000, dest="predlen", help="number of reservoir units")

args = parser.parse_args()

if not args.i:
  print('[ERROR] use -i FILE_NAME to set an input file')
  exit()

data = np.load(args.i)
print('loaded %d points' % len(data))

if args.t == 'dnf':
  data = np.log(data)
elif args.t == 'ig':
  data = data[0::2]/10000 - 1.1

esn = ESN(n_in=1,
          n_fb=1,
          n_units=args.u,
          spectral_radius=args.s)

trainlen = args.trainlen
predlen = args.predlen
print('using: %d, predicting: %d' % (trainlen, predlen))

print('fitting...')
fit = esn.fit(np.ones(trainlen), data[:trainlen])
print('predicting...')
pred = esn.predict(np.ones(predlen), cont=True)

if args.t == 'dnf':
  data = np.exp(data)
  pred = np.exp(pred)
elif args.t == 'ig':
  data = data+1.1 * 10000
  pred = pred+1.1 * 10000

plt.figure(figsize=(200,30))
end = trainlen + predlen + 1
datalen = min(end, len(data))
plt.plot(range(0, datalen), data[0:datalen], '--bo', linestyle="solid", linewidth=3, label="actual", alpha=0.5)
plt.plot(range(trainlen+1, end), pred, '--ro', linestyle="solid", linewidth=3, label="prediction", alpha=0.5)

fn = './%s_%d_%s.png' % (args.t, args.u, str(args.s))
print('saving %s' % fn)
plt.savefig(fn)
