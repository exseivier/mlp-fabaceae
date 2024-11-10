#!/usr/bin/env python3

import tensorflow as tf
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import os
import random as rnd
import subprocess as sp
import json as js
import gc
from tqdm import tqdm

def create_model():
    ################################################
    ##
    ##    Building the Clarice ML Perceptron model
    ##
    
    model = tf.keras.models.Sequential([
        tf.keras.layers.Input(shape = (1271,), sparse = True),
        tf.keras.layers.Dense(units = 500),
        tf.keras.layers.Dropout(0.2),
        tf.keras.layers.Dense(units = 250),
        tf.keras.layers.Dropout(0.2),
        tf.keras.layers.Dense(units = 100),
        tf.keras.layers.Dense(units = 50),
        tf.keras.layers.Dense(units = 1, activation=tf.keras.activations.sigmoid)
    ])
    #######################################################
    ##
    ##    Compiling the model.
    ##

    model.compile(optimizer=tf.keras.optimizers.Adam(0.0001, weight_decay=0.0000001), \
            loss = tf.keras.losses.BinaryCrossentropy(), \
            metrics=[tf.keras.metrics.BinaryAccuracy(), tf.keras.metrics.Accuracy()])

    model.summary()
    return model

def reset_history():
    history = {"accuracy" : [], \
            "loss" : [], \
            "val_accuracy" : [], \
            "val_loss" : [], \
            "binary_accuracy" : [], \
            "val_binary_accuracy" : []}
    return history

def reset_modelChackPoint():
    chkpnt = tf.keras.callbacks.ModelCheckpoint(\
                filepath="tmp.keras", \
                save_best_only=True, \
                mode = "max", \
                monitor="val_binary_accuracy")
    return chkpnt
   


def train(model, chkpnt, model_name, path_datasets, random_num, epochs_num, path_models, history):
    datasets = os.listdir(path_datasets)
    datasets = sorted(datasets)
    
    for i in range(1, len(datasets) +1):
        ds = datasets[i -1]
        print(f'Training with {path_datasets}/{ds} dataset')
        train_ds = pd.read_csv(f"{path_datasets}/{ds}")
        for r in tqdm(range(1, random_num +1)):
            train_ds = train_ds.reindex(np.random.permutation(train_ds.index))
            mask = np.random.rand(len(train_ds)) <= 0.8
            train_x = train_ds[mask]
            train_x.reset_index(inplace = True, drop = True)
            test_x = train_ds[~mask]
            test_x.reset_index(inplace = True, drop = True)
            train_label = train_x.pop("Train")
            test_label = test_x.pop("Train")
            train_x = tf.keras.utils.normalize(train_x)
            test_x = tf.keras.utils.normalize(test_x)

            filename = f"{path_models}/model-{model_name}-{i}-{r}-" + "{epoch:003d}-{val_binary_accuracy:0.4f}.keras"
            chkpnt.filepath = filename

            fitted = model.fit(\
                    train_x, \
                    train_label, \
                    batch_size=32, \
                    epochs= epochs_num, \
                    validation_split=0.2, \
                    shuffle=True, \
                    callbacks=[chkpnt], \
                    verbose=False)

            for key, value in fitted.history.items():
                history[key].extend(value)

    return history

def plot_metrics(history, start = 0, end = -1, maxy_loss_ad = 0.1, maxy_acc_ad = 0.02):

    ##########################################################
    ##
    ##    Plotting metrics
    ##
    f, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    miny_loss = min(min(history["val_loss"]), min(history["loss"]))
    maxy_loss = max(max(history["val_loss"]), max(history["loss"]))
    miny_acc = min(min(history["val_binary_accuracy"]), min(history["binary_accuracy"]))
    maxy_acc = max(max(history["val_binary_accuracy"]), max(history["binary_accuracy"]))
    ax1.plot(history["val_loss"], '-r', label="Validation loss")
    ax1.plot(history["loss"], '-b', label = "Training loss")
    ax2.plot(history["val_binary_accuracy"], '-or', label = "Validation binary accuracy")
    ax2.plot(history["binary_accuracy"], '-ob', label = "Training binary accuracy")
    ax1.set_xlabel("Iterations")
    ax1.set_ylabel("Loss values")
    ax2.set_ylabel("Accuracy percentage")
    ax1.legend(loc = "upper right")
    ax2.legend(loc = "upper left")
    ax1.set_ylim(miny_loss - 0.05, maxy_loss + maxy_loss_ad)
    ax2.set_ylim(miny_acc - 0.01, maxy_acc +maxy_acc_ad)
    plt.show()


