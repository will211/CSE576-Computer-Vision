import sys
sys.dont_write_bytecode = True
from uwimg import *

def softmax_model(inputs, outputs):
    l = [make_layer(inputs, outputs, SOFTMAX)]
    return make_model(l)

def neural_net(inputs, outputs):
    print(inputs)
    l = [   make_layer(inputs, 64, RELU),
            make_layer(64, 32, RELU),
            make_layer(32, outputs, SOFTMAX)]
    return make_model(l)

print("loading data...")
train = load_classification_data(b"cifar.train", b"mnist.labels", 1)
test  = load_classification_data(b"cifar.test", b"mnist.labels", 1)
print("done")
print

print("training model...")
batch = 128
iters = 30000
rate = .015
momentum = .9
decay = .1

m = neural_net(train.X.cols, train.y.cols)
train_model(m, train, batch, iters, rate, momentum, decay)
print("done")
print

print("evaluating model...")
print("training accuracy: %0.2f %%"%(100*accuracy_model(m, train)))
print("test accuracy:     %0.2f %%"%(100*accuracy_model(m, test)))