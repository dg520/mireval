#!/usr/bin/python


import random

a=raw_input("Fold File: ")
e=raw_input("The fold file is from 1. pos data; 2.neg data: ")

if e=="2":

    b=raw_input("How many training sequences to be picked out: ")
    c=raw_input("How many test sequences to be picked out: ")

    d1=[]
    d2=[]

    f=open(a,"r")
    data=f.readlines()
    f.close()

    line_list=random.sample(range(0,len(data)/3),int(b)+int(c))
    line_list_train=line_list[:int(b)]
    line_list_test=line_list[int(b):]

    

    for x in line_list_train:
        d1.append(data[3*x])
        d1.append(data[3*x+1])
        d1.append(data[3*x+2])

    f=open("train_fold_neg_rep1.txt","w")
    f.writelines(d1)
    f.close()


    for x in line_list_test:
        d2.append(data[3*x])
        d2.append(data[3*x+1])
        d2.append(data[3*x+2])

    f=open("test_fold_neg_rep1.txt","w")
    f.writelines(d2)
    f.close()

if e=="1":

    b=raw_input("How many sequences to be picked out: ")
    c=raw_input("Purpose: 1. for training; 2.for test ")

    d1=[]

    f=open(a,"r")
    data=f.readlines()
    f.close()

    line_list=random.sample(range(0,len(data)/3),int(b))
    line_list_train=line_list[:int(b)]

    for x in line_list_train:
        d1.append(data[3*x])
        d1.append(data[3*x+1])
        d1.append(data[3*x+2])

    if c=="1":
        name="train_fold_pos_rep1.txt"
    if c=="2":
        name="test_fold_pos_rep1.txt"

    f=open(name,"w")
    f.writelines(d1)
    f.close()
