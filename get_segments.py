#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  7 08:25:15 2021

@author: xuchangmao
"""
import numpy as np
import scipy
from scipy.integrate import simps  # 用于计算积分
import matplotlib.pyplot as plt  # 用于画图
import pandas as pd
import random
from random import choice

raw_test_data = pd.read_excel(r'/Users/xuchangmao/Documents/工作/平台数据/6台国六-10月15-20日的数据/LFNG4HC90LTY00157（15-20.）.xlsx')
test_data = raw_test_data[:360000]


def group(L):
    first = last = L[0]
    for n in L[1:]:
        if n - 1 == last:
            last = n
        else:
            yield first, last
            first = last = n
    yield first, last


def get_segments(vehicle_speed):
    test_cycle_array = np.array(vehicle_speed)
    nonzero_index = test_cycle_array.nonzero()
    nonzero_list = list(nonzero_index[0])
    segments_group = list(group(nonzero_list))
    return segments_group


new_segments = []


def re_get_segments(raw_segments):
    raw_segments.insert(0, (0, 0))
    for i in range(len(raw_segments)):
        if i != len(raw_segments) - 1:
            new_segments.append((raw_segments[i][-1], raw_segments[i + 1][-1]))
        continue
    return new_segments


# 接下来是对不同平均速度的区间进行筛选，分到不同的列表中，列表的元素是原始片段的迄止区间点如（12，52）
urban_list = []
sub_list = []
high_list = []


def filter_segments(segments):
    for segment in new_segments:
        average_speed = test_data['车速'][segment[0]:segment[1]].mean()
        if average_speed > 75:
            high_list.append(segment)
        if average_speed > 55 and average_speed < 75:
            sub_list.append(segment)
        if average_speed < 55:
            urban_list.append(segment)


# 下面的两个sum_work是1型试验中的whtc循环的功，而max_torque是最大的扭矩，也就是扭矩百分比的计算依据
sum_work = 66
max_torque = 600

torque_percent = test_data['发动机净输出扭矩'] - test_data['摩擦扭矩']
torque_percent[torque_percent < 0] = 0
cal_w = torque_percent * max_torque * 0.01 * test_data['发动机转速'] / 9550
# 接下里需要写一个对累积功进行计算的

index = cal_w.index.tolist()
index = np.array(index)
# index相当于是时间序列，下面要转换为小时h
cal_w_np = np.array(cal_w)
y = cal_w_np
x = index / 3600  # 构造一个等差数列数组，作为横轴取值，单位为h，小时
integrals = []
for i in range(len(y)):  # 计算梯形的面积，面积其实就是所作的功，由于是累加，所以是切片"i+1"
    integrals.append(scipy.integrate.trapz(y[:i + 1], x[:i + 1]))


# integrals里面存储的是所有的累积功，而cal_w里面存储的是单点的功值


# 要创建一个函数用来计算一个含有不同segments的累积功的多少的

def get_cum_work(segments_list):
    segments_list_cum_work = 0
    for segment in segments_list:
        segments_list_cum_work = segments_list_cum_work + integrals[segment[1]] - integrals[segment[0]]
    return segments_list_cum_work


# 接下来是对3个list中的行程片段进行选择，累积功只要没超过限制那么就一直累加，采用random的形式进行选择
# 首先先定一个值，这个值就是循环功的3分之1
urban_work = 0
sub_work = 0
high_work = 0

# 函数主体运行情况
raw_segments = get_segments(test_data['车速'])
new_segments = re_get_segments(raw_segments)
filter_segments(new_segments)

pems_urban_list = []
while get_cum_work(pems_urban_list) < sum_work / 3:
    try:
        random_urban_segment = choice(urban_list)
        pems_urban_list.append(random_urban_segment)
    except IndexError:
        print('无法找到相应片段')
        break

pems_sub_list = []
while get_cum_work(pems_sub_list) < sum_work / 3:
    try:
        random_sub_segment = choice(sub_list)
        pems_sub_list.append(random_sub_segment)
    except IndexError:
        print('无法找到相应片段')
        break

pems_high_list = []
while get_cum_work(pems_high_list) < sum_work / 3:
    try:
        random_high_segment = choice(high_list)
        pems_high_list.append(random_high_segment)
    except IndexError:
        print('无法找到相应片段')
        break

print(get_cum_work(pems_urban_list))
print(get_cum_work(pems_sub_list))
print(get_cum_work(pems_high_list))

pems_urban_work = get_cum_work(pems_urban_list)
pems_sub_work = get_cum_work(pems_sub_list)
pems_high_work = get_cum_work(pems_high_list)
total_pems_work = pems_urban_work + pems_sub_work + pems_high_work

pems_urban_ratio = pems_urban_work / total_pems_work
pems_sub_ratio = pems_sub_work / total_pems_work
pems_high_ratio = pems_high_work / total_pems_work

print(pems_urban_ratio)
print(pems_sub_ratio)
print(pems_high_ratio)















