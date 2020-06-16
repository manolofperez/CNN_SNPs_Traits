import numpy as np
import pandas as pd
import os
import glob
import keras
from keras.datasets import mnist
from keras.models import Sequential
from keras.layers import Input, Activation, Dense, Dropout, Flatten, concatenate
from keras.layers import Conv1D, MaxPooling1D, AveragePooling1D
from keras.models import Model
from keras import backend as K
from random import shuffle, choice
from sklearn.preprocessing import MinMaxScaler


batch_size = 250
epochs = 20
num_classes = 3

def process_structured_data(dftrain, dftest):
	"""
	Pre-processes the given dataframe by minmaxscaling the continuous features
	(fit-transforming the training data and transforming the test data)
	"""
	cs = MinMaxScaler()
	dftrain = cs.fit_transform(dftrain)
	dftest = cs.transform(dftest)
	return (dftrain, dftest)


def create_mlp(dftrain, regularizer=None):
	"""Creates a simple two-layer MLP with inputs of the given dimension"""
	model = Sequential()
	model.add(Dense(10, input_dim=dftrain.shape[1], activation="relu", kernel_regularizer=regularizer))
	model.add(Dense(6, activation="relu", kernel_regularizer=regularizer))
	return model


def create_cnn(xtest, regularizer=None):
	inputShape = (xtest.shape[1], xtest.shape[2])
	inputs = Input(shape=inputShape)
	x = inputs
	x = Conv1D(250, kernel_size=2, activation='relu',input_shape=(xtest.shape[1], xtest.shape[2]))(x)
	x = Conv1D(125, kernel_size=2, activation='relu')(x)
	x = AveragePooling1D(pool_size=2)(x)
	x = Dropout(0.75)(x)
	x = Conv1D(125, kernel_size=2, activation='relu')(x)
	x = AveragePooling1D(pool_size=2)(x)
	x = Dropout(0.75)(x)
	x = Flatten()(x)
	x = Dense(125, activation='relu')(x)
	x = Dropout(0.5)(x)
	x = Dense(125, activation='relu')(x)
	x = Dropout(0.5)(x)
	# Apply another fully-connected layer, this one to match the number of nodes coming out of the MLP
	x = Dense(6, kernel_regularizer=regularizer)(x)
	x = Activation("relu")(x)
	# Construct the CNN
	model = Model(inputs, x)
	# Return the CNN
	return model


path = r'./traits/OU' # use your path
all_files = glob.glob(path + "/*.txt")
all_files.sort(key=os.path.getmtime)

df = []
df = np.stack([np.loadtxt(f) for f in all_files])
df = np.array(df)


u1 = np.load("trainingSims/simModel1.npz")
u2 = np.load("trainingSims/simModel2.npz")
u3 = np.load("trainingSims/simModel3.npz")
x=np.concatenate((u1['simModel1'],u2['simModel2'],u3['simModel3']),axis=0)

y=[0 for i in range(len(u1['simModel1']))]
y.extend([1 for i in range(len(u2['simModel2']))])
y.extend([2 for i in range(len(u3['simModel3']))])
y = np.array(y)

print (len(x), len(y), len(df))
shf = list(range(len(x)))
shuffle(shf)

y = y[shf]
x = x[shf]
df = df[shf]

xtrain, xtest = x[100:], x[:100]
ytrain, ytest = y[100:], y[:100]
dftrain, dftest = df[100:], df[:100]


ytest = keras.utils.to_categorical(ytest, num_classes)
ytrain = keras.utils.to_categorical(ytrain, num_classes)
dftrain=dftrain.reshape((dftrain.shape[0], (dftrain.shape[1]*dftrain.shape[2])))
dftest=dftest.reshape((dftest.shape[0], (dftest.shape[1]*dftest.shape[2])))
dftrain, dftest = process_structured_data(dftrain, dftest)
# Create the MLP and CNN models
mlp = create_mlp(dftrain.shape)
cnn = create_cnn(xtest)

# Create the input to the final set of layers as the output of both the MLP and CNN
combinedInput = concatenate([mlp.output, cnn.output])

# The final fully-connected layer head will have two dense layers (one relu and one sigmoid)
x = Dense(6, activation="relu")(combinedInput)
x = Dense(num_classes, activation="softmax")(x)


# The final model accepts numerical data on the MLP input and images on the CNN input, outputting a single value
model1 = Model(inputs=[mlp.input, cnn.input], outputs=x)

model1.compile(loss=keras.losses.categorical_crossentropy,
	              optimizer=keras.optimizers.Adam(),
	              metrics=['accuracy'])

print(model1.summary())

model1.fit([dftrain, xtrain], ytrain, batch_size=batch_size,
          epochs=epochs,
          verbose=1,
          validation_data=([dftest, xtest], ytest))

model1.save(filepath='Trained_Comb_Model.acc.mod')

xCNN = Dense(num_classes, activation="softmax")(cnn.output)

model2 = Model(inputs=cnn.input, outputs=xCNN)

model2.compile(loss=keras.losses.categorical_crossentropy,
	              optimizer=keras.optimizers.Adam(),
	              metrics=['accuracy'])

print(model2.summary())

model2.fit(xtrain, ytrain, batch_size=batch_size,
          epochs=epochs,
          verbose=1,
          validation_data=(xtest, ytest))

model2.save(filepath='Trained_CNN_Model.acc.mod')


xMLP = Dense(num_classes, activation="softmax")(mlp.output)

model3 = Model(inputs=mlp.input, outputs=xMLP)

model3.compile(loss=keras.losses.categorical_crossentropy,
	              optimizer=keras.optimizers.Adam(),
	              metrics=['accuracy'])

print(model3.summary())

model3.fit(dftrain, ytrain, batch_size=batch_size,
          epochs=epochs,
          verbose=1,
          validation_data=(dftest, ytest))
          
model3.save(filepath='Trained_Traits_Model.acc.mod')


