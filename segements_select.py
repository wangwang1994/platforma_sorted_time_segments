import numpy as np
import scipy
from scipy.integrate import simps  # 用于计算积分
import matplotlib.pyplot as plt  # 用于画图
import pandas as pd

#导入原始数据
raw_test_data = pd.read_excel(r'C:\Users\CAERI\Desktop\6台国六-10月15-20日的数据\LFNG4HC90LTY00157（15-20.）.xlsx')

#原始数据拆分成为片段，涉及到对短行程的定义
test_data=raw_test_data[:14400]
#这是一个辅助函数，对区间划分使用的
def group(L):
    first = last = L[0]
    for n in L[1:]:
        if n - 1 == last:
            last = n
        else:
            yield first, last
            first = last = n
    yield first, last
#先形成一个原始的list包含有不带怠速区间的片段
def get_segments(vehicle_speed):
    test_cycle_array = np.array(vehicle_speed)
    nonzero_index = test_cycle_array.nonzero()
    nonzero_list = list(nonzero_index[0])
    segments_group = list(group(nonzero_list))
    return segments_group
new_segments=[]
#随后形成一个带有怠速区间的片段
def re_get_segments(raw_segments):
    raw_segments.insert(0,(0,0))
    for i in range(len(raw_segments)):
        if i!=len(raw_segments)-1:
            new_segments.append((raw_segments[i][-1],raw_segments[i+1][-1]))
        continue
    return new_segments
