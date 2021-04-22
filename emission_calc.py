#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 15 13:31:59 2021

@author: xuchangmao
"""

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

max_power = 120
whole_windows = 0
windows = []
windows_eff = []
average_emission = []


def get_single_window(start, integrals):
    global whole_windows
    for i in range(start, len(integrals) - 1):
        sum = integrals[i + 1] - integrals[start]
        sum_minus = integrals[i] - integrals[start]
        if (sum >= 11.64) and (sum_minus < 11.64):
            whole_windows = whole_windows + 1
            if cal_w_np[start:i].mean() > 0.2 * max_power:
                windows.append((start, i))
            break
    return windows


for j in range(1, 48):
    print('-----进行第' + str(j) + '组数据处理-----')
    raw_test_data = pd.read_csv(r'/Users/xuchangmao/Documents/工作/代码/pems_cycles/' + str(j) + 'th pems_data.csv')
    # 导入数据
    test_data = raw_test_data[:]
    print('----------数据载入完成------------')
    test_data.reindex(index=test_data.index[::-1])
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

    # plt.plot(x, y)
    #
    # fig2=plt.figure(2)
    # '''
    # 以上就是第一部分，关于所做功的计算，
    # 首先计算车辆的实时功率，然后计算出所做功的累积曲线
    # '''
    # plt.plot(x, integrals)
    # integrals_df=pd.DataFrame(integrals)
    # windows=[]
    # plt.plot(x, integrals)
    # plt.title('Cumulative work')
    # plt.show()
    print('----------累积功计算完成------------')
    for i in range(0, 14400):
        get_single_window(i, integrals)
    test_data['发动机燃料流量'] = test_data['发动机燃料流量'].replace('无效', 0)
    test_data['发动机燃料流量'] = test_data['发动机燃料流量'].astype(float)
    print('----------窗口搜寻完成------------')

    fuel_mass_flow = 0.84 * test_data['发动机燃料流量'] / 3600
    air_mass_flow = test_data['进气量'] / 3600
    total_mass_flow = air_mass_flow + fuel_mass_flow
    nox_conc = test_data['SCR下游NOx传感器输出'].replace('无效', 0)
    nox_conc = nox_conc.astype(float)
    nox_emission_list = []
    for k, m in windows:
        nox_emission = 0.001587 * (nox_conc[k:m] * total_mass_flow).sum()
        nox_emission_list.append(nox_emission)
    windows_eff.append(len(windows) / whole_windows)
    average_emission.append((np.mean(nox_emission_list) * 1000) / 11.64)
    print('----------' + '第' + str(j) + '组数据计算完成------------')
info_data = {
    '窗口有效率': windows_eff,
    '平均排放（mg）': average_emission
}
info_data_pd = pd.DataFrame(info_data)
info_data_pd.to_csv('/Users/xuchangmao/Documents/工作/代码/pems_cycles/排放计算结果.csv')
print('处理完成')