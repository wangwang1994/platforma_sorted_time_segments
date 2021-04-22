import numpy as np
import scipy
from scipy.integrate import simps  # 用于计算积分
import matplotlib.pyplot as plt  # 用于画图
import pandas as pd
import random
from random import choice
import shutil
import os
raw_test_data = pd.read_excel(
    r'/Users/xuchangmao/Documents/工作/平台数据/6台国六-10月15-20日的数据/LFNG4HC90LTY00157（15-20.）.xlsx')  # 读取文件
test_data = raw_test_data[:360000]



index = test_data.index.tolist()
x = index / 3600  # 构造一个等差数列数组，作为横轴取值，单位为h，小时
y = test_data['发动机燃料流量']
integrals = []
for i in range(len(y)):  # 计算梯形的面积，面积其实就是所作的功，由于是累加，所以是切片"i+1"
    integrals.append(scipy.integrate.trapz(y[:i + 1], x[:i + 1]))
print(integrals)
