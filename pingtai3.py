# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import pandas as pd
import numpy as np

data = pd.read_excel(r'C:\Users\CAERI\Desktop\data\WLTC.xlsx')
time = data['Time [s]']
speed = data['VS [km/h]']


def cal_deriv(x, y):  # x, y的类型均为列表
    diff_x = []  # 用来存储x列表中的两数之差
    for i, j in zip(x[0::], x[1::]):
        diff_x.append(j - i)

    diff_y = []  # 用来存储y列表中的两数之差
    for i, j in zip(y[0::], y[1::]):
        diff_y.append(j - i)

    slopes = []  # 用来存储斜率
    for i in range(len(diff_y)):
        slopes.append(diff_y[i] / diff_x[i])

    deriv = []  # 用来存储一阶导数
    for i, j in zip(slopes[0::], slopes[1::]):
        deriv.append((0.5 * (1 / 3.6) * (i + j)))  # 根据离散点导数的定义，计算并存储结果
    deriv.insert(0, slopes[0])  # (左)端点的导数即为与其最近点的斜率
    deriv.append(slopes[-1])  # (右)端点的导数即为与其最近点的斜率

    for i in deriv:  # 打印结果，方便检查，调用时也可注释掉
        # print(i)

        return deriv  # 返回存储一阶导数结果的列表


x = list(speed.index)
y = list(speed)
acc = cal_deriv(x, y)
acc_pd = pd.DataFrame(acc)
pos_acc = []
rpa_pd = acc_pd[acc_pd > 0.1]
speed_pd = pd.DataFrame(speed)
di = speed_pd / 3.6
di.columns = ['di']
rpa_pd.columns = ['rpa']
speed_urban = speed_pd[speed_pd < 60]
RPA_pd = rpa_pd['rpa'] * speed_pd['VS [km/h]'] / di['di']
RPA_pd = pd.DataFrame(RPA_pd)
RPA_pd.columns = ['rpa']
data_speed_rpa = RPA_pd.join(speed_pd)
data_speed_rpa_urban = data_speed_rpa[data_speed_rpa['VS [km/h]'] < 60]
rpa_urban = data_speed_rpa_urban['rpa'].mean()
data_speed_rpa_urban = data_speed_rpa[(data_speed_rpa['VS [km/h]'] > 60) & (data_speed_rpa['VS [km/h]'] < 90)]
rpa_rural = data_speed_rpa_urban['rpa'].mean()
data_speed_rpa_urban = data_speed_rpa[data_speed_rpa['VS [km/h]'] > 90]
rpa_motor = data_speed_rpa_urban['rpa'].mean()
print(rpa_urban)
print(rpa_rural)
print(rpa_motor)
