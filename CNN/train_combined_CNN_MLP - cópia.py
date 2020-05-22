import numpy as np
import keras
from keras.datasets import mnist
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten
from keras.layers import Conv1D, MaxPooling1D, AveragePooling1D
from keras import backend as K
from random import shuffle, choice

batch_size = 250
epochs = 20
num_classes = 3

u1 = np.load("trainingSims/simModel1.npz")
u2 = np.load("trainingSims/simModel2.npz")
u3 = np.load("trainingSims/simModel3.npz")
x=np.concatenate((u1['simModel1'],u2['simModel2'],u3['simModel3']),axis=0)

y=[0 for i in xrange(len(u1['simModel1']))]
y.extend([1 for i in xrange(len(u2['simModel2']))])
y.extend([2 for i in xrange(len(u3['simModel3']))])
y = np.array(y)

print len(x), len(y)
shf = range(len(x))
shuffle(shf)

y = y[shf]
x = x[shf]

xtrain, xtest = x[10:], x[:10]
ytrain, ytest = y[10:], y[:10]

ytest = keras.utils.to_categorical(ytest, num_classes)
ytrain = keras.utils.to_categorical(ytrain, num_classes)

model = Sequential()
model.add(Conv1D(250, kernel_size=2,
                 activation='relu',
                 input_shape=(xtest.shape[1], xtest.shape[2])))
model.add(Conv1D(125, kernel_size=2, activation='relu'))
model.add(AveragePooling1D(pool_size=2))
model.add(Dropout(0.25))
model.add(Conv1D(125, kernel_size=2, activation='relu'))
model.add(AveragePooling1D(pool_size=2))
model.add(Dropout(0.5))
model.add(Flatten())
model.add(Dense(125, activation='relu'))
model.add(Dropout(0.5))
model.add(Dense(125, activation='relu'))
model.add(Dropout(0.75))
model.add(Dense(num_classes, activation='sigmoid'))
model.compile(loss=keras.losses.categorical_crossentropy,
              optimizer=keras.optimizers.Adam(),
              metrics=['accuracy'])
print(model.summary())
model.fit(xtrain, ytrain, batch_size=batch_size,
          epochs=epochs,
          verbose=1,
          validation_data=(xtest, ytest))


model.save(filepath='Trained_Model.acc.mod')
