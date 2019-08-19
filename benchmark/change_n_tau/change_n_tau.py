#!/usr/bin/python3
import matplotlib.pyplot as plt
import csv

csvfile=open('change_n_tau.dat','r')
f = csv.reader(csvfile,delimiter='\t')
x100 = []
y100 = []
t100 = []
x200 = []
y200 = []
t200 = []
x400 = []
y400 = []
t400 = []
x800 = []
y800 = []
t800 = []
for row in f:
    if(row[0]=='100'):
        t100.append(float(row[1]))
        x100.append(float(row[4]))
        y100.append(float(row[3]))
    elif(row[0]=='200'):
        t200.append(float(row[1]))
        x200.append(float(row[4]))
        y200.append(float(row[3]))
    elif(row[0]=='400'):
        t400.append(float(row[1]))
        x400.append(float(row[4]))
        y400.append(float(row[3]))
    elif(row[0]=='800'):
        t800.append(float(row[1]))
        x800.append(float(row[4]))
        y800.append(float(row[3]))
csvfile.close()



plt.xscale('log')
plt.xlabel('runtime [s]', fontsize=20)
plt.ylabel('${\cal S}_E$', fontsize=20)
plt.ylim([45.96,46.02])
plt.plot(x100,y100, label='n=100', marker='o')
plt.plot(x200,y200, label='n=200', marker='o')
plt.plot(x400,y400, label='n=400', marker='o')
plt.plot(x800,y800, label='n=800', marker='o')
for x,y,t in zip(x100,y100,t100):
    plt.text(x*0.81,y+0.007,'$\\tau='+str(t)+'$', rotation=90)
for x,y,t in zip(x200,y200,t200):
    plt.text(x*0.81,y+0.007,'$\\tau='+str(t)+'$', rotation=90)
for x,y,t in zip(x400,y400,t400):
    plt.text(x*0.81,y+0.007,'$\\tau='+str(t)+'$', rotation=90)
for x,y,t in zip(x800,y800,t800):
    plt.text(x*0.81,y+0.007,'$\\tau='+str(t)+'$', rotation=90)
plt.tick_params(labelsize=14)
plt.yticks([45.96, 45.97, 45.98, 45.99, 46., 46.01, 46.02],  [45.96, 45.97, 45.98, 45.99, 46., 46.01, 46.02])
plt.legend(loc='lower right')
plt.savefig('accuracy.pdf')

