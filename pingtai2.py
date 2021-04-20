# -*- coding: utf-8 -*-
"""
Created on Thu Feb  4 14:59:19 2021

@author: CAERI
"""
import numpy as np
import scipy
from scipy.integrate import simps  # 用于计算积分
import matplotlib.pyplot as plt  # 用于画图
import pandas as pd

raw_test_data = pd.read_excel(r'C:\Users\CAERI\Desktop\6台国六-10月15-20日的数据\LFNG4HC90LTY00157（15-20.）.xlsx')
# 导入数据
test_data = raw_test_data[:14400]
torque_percent = test_data['发动机净输出扭矩'] - test_data['摩擦扭矩']
torque_percent[torque_percent < 0] = 0
cal_w = torque_percent * 600 * 0.01 * test_data['发动机转速'] / 9550
# 计算发动机的实时功率，数据中的输出扭矩减去摩擦扭矩才是实际扭矩百分比，对不同的车要乘以其标称扭矩，得到N*m单位的实际扭矩
index = cal_w.index.tolist()
index = np.array(index)
# index相当于是时间序列，下面要转换为小时h
cal_w_np = np.array(cal_w)
y = cal_w_np
x = index / 3600  # 构造一个等差数列数组，作为横轴取值，单位为h，小时
integrals = []

for i in range(len(y)):  # 计算梯形的面积，面积其实就是所作的功，由于是累加，所以是切片"i+1"
    integrals.append(scipy.integrate.trapz(y[:i + 1], x[:i + 1]))
plt.plot(x, y)
plt.show()
'''
以上就是第一部分，关于所做功的计算，
首先计算车辆的实时功率，然后计算出所做功的累积曲线
'''
plt.plot(x, integrals)
integrals_df = pd.DataFrame(integrals)
windows = []
plt.plot(x, integrals)
plt.show()

'''integrals中存储了功的累计量，接下来进行单一窗口的选取计算
选取的方法是，固定起点start，不断的往后寻找
法规中对应：

sum=integrals[i+1]-integrals[start]
sum_minus=integrals[i]-integrals[start]

'''
max_power = 120


def get_single_window(start, integrals):
    for i in range(start, len(integrals)):
        sum = integrals[i + 1] - integrals[start]
        sum_minus = integrals[i] - integrals[start]
        if (sum >= 11.64) and (sum_minus < 11.64):
            print('第i位满足窗口条件', i)
            print('当前的最终窗口与前一个窗口的累积功大小为：', str(sum) + '，' + str(sum_minus))
            if cal_w_np[start:i].mean() > 0.2 * max_power:
                windows.append((start, i))
                # print('当前的windows是',windows)
            break
    return windows


def get_average_power(window):
    average_power = cal_w[window[0]:window[1]].mean()
    return average_power


average_w_list = []


def get_eff_w(windows):
    for win in windows:
        average_w_list.append(get_average_power(win))
    average_w_pd = pd.DataFrame(average_w_list)
    judge_w = average_w_pd > (0.1 * 11.64)  # 判断有多少个超过了功率的百分之20
    judge_w[0] = judge_w[0].astype('int')  # 将true 与 false转换为1和0
    effective_w = judge_w.sum()[0]
    efficiency_w = effective_w / len(judge_w)

    return efficiency_w






