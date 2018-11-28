#!/usr/bin/python3
import os
from os import path
import sys
import matplotlib.pyplot as plt
import numpy as np
import csv
import pandas as pd
from scipy import stats

def plot(x_label, y_label, file):
	file_path = path.relpath(file)
	height = 5
	width = 20
	opacity = 1
	data = pd.read_table(file, sep=" ")
	X = data[data.columns[0]].values
	Y = data[data.columns[1]].values
	name = file[:file.find(".dat")]
	plt.figure(num=None, figsize=(width, height), dpi=80, facecolor='w', edgecolor='k')
	plt.plot(X,Y,"-b", alpha=opacity, label=name)
	plt.title(name)
	plt.ylabel(y_label)
	plt.xlabel(x_label)
	plt.legend()
	plt.tight_layout()
	plt.savefig(name + ".png")

def plot_fit_line(x_label, y_label, file):
	file_path = path.relpath(file)
	height = 5
	width = 20
	opacity = 1
	data = pd.read_table(file, sep=" ")
	X = data[data.columns[0]].values
	Y = data[data.columns[1]].values
	name = file[:file.find(".dat")]
	Y_ideal = X*(3/2*125*0.00831)/(4/3*3.14*2.3**3)
	plt.figure(num=None, figsize=(width, height), dpi=80, facecolor='w', edgecolor='k')
	plt.plot(X,Y,"ob", alpha=opacity, label=name)
	plt.plot(X,Y_ideal, "-r", alpha=opacity, label='P(T)')
	plt.title(name)
	plt.ylabel(y_label)
	plt.xlabel(x_label)
	plt.legend()
	plt.tight_layout()
	plt.savefig(name + ".png")


P_label = 'P[1u nm^-1 ps^-2]'
T_label = 'T[K]'
t_label = 't[ps]'

e_label = 'E'
x_label = 'X'
n_label = 'N'

plot(t_label, e_label, "e_n=1.dat")
plot(t_label, n_label, "n_n=1.dat")
plot(t_label, x_label, "x_n=1.dat")

plot(t_label, e_label, "e_n=4.dat")
plot(t_label, n_label, "n_n=4.dat")
plot(t_label, x_label, "x_n=4.dat")

plot(t_label, e_label, "e_n=9.dat")
plot(t_label, n_label, "n_n=9.dat")
plot(t_label, x_label, "x_n=9.dat")

plot(t_label, "E_var", "e_var.dat")