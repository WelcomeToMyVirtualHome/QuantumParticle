#!/usr/bin/python3
import os
from os import path
import sys
import matplotlib.pyplot as plt
import numpy as np
import csv
import pandas as pd
import glob
from scipy import stats
from matplotlib.pyplot import cm

def plot_subplots(n):
	height = 20
	width = 20
	t_label = "t[ps]"
	fig, ax = plt.subplots(2, 2, figsize=(width,height))
	ax[0,0] = plot_on_ax(ax[0,0], t_label, "E", "e_n={:d}.dat".format(n))
	ax[0,1] = plot_on_ax(ax[0,1], t_label, "X", "x_n={:d}.dat".format(n))
	ax[1,0] = plot_on_ax(ax[1,0], t_label, "N", "n_n={:d}.dat".format(n))
	ax[1,1] = plot_rho_on_ax(ax[1,1], "rho_n={:d}.dat".format(n))
	plt.savefig("plots_n={:d}.png".format(n))

def plot_on_ax(ax, x_label, y_label, file):
	file_path = path.relpath(file)
	opacity = 1
	data = pd.read_table(file, sep=" ")
	X = data[data.columns[0]].values
	Y = data[data.columns[1]].values
	name = file[:file.find(".dat")]
	ax.plot(X,Y,"-b", alpha=opacity, label=name)
	ax.set_title(name)
	ax.set_xlabel(x_label)
	ax.set_ylabel(y_label)
	return ax

def plot_rho_on_ax(ax, file):
	file_path = path.relpath(file)
	opacity = 1
	data = pd.read_table(file, sep=" ")
	X = data[data.columns[0]].values
	Y = data[data.columns[1]].values
	length = len(X)
	chunk = 100
	name = file[:file.find(".dat")]
	ax.plot(X[:chunk],Y[:chunk],"-r", alpha=opacity, label="t=0")
	ax.plot(X[int(length*0.5):int(length*0.5) + 100],Y[int(length*0.5):int(length*0.5) + 100],"-g", alpha=opacity, label="t=0.5T")
	ax.plot(X[int(length*0.8):int(length*0.8) + 100],Y[int(length*0.8):int(length*0.8) + 100],"-b", alpha=opacity, label="t=0.8T")
	ax.set_title(name)
	ax.set_ylabel("rho")
	ax.legend()
	return ax

plot_subplots(1)
plot_subplots(2)
plot_subplots(4)