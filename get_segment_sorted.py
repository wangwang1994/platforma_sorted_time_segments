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
from random import choice
import shutil
import os

'''
每次运行前先删除上次的结果，然后重新创建文件夹来存储本次运行结果
'''
shutil.rmtree('/Users/xuchangmao/Documents/工作/代码/urban')
shutil.rmtree('/Users/xuchangmao/Documents/工作/代码/sub')
shutil.rmtree('/Users/xuchangmao/Documents/工作/代码/high')
shutil.rmtree('/Users/xuchangmao/Documents/工作/代码/pems_urban')
shutil.rmtree('/Users/xuchangmao/Documents/工作/代码/pems_sub')
shutil.rmtree('/Users/xuchangmao/Documents/工作/代码/pems_high')
shutil.rmtree('/Users/xuchangmao/Documents/工作/代码/pems_cycles')
os.mkdir('/Users/xuchangmao/Documents/工作/代码/urban')
os.mkdir('/Users/xuchangmao/Documents/工作/代码/sub')
os.mkdir('/Users/xuchangmao/Documents/工作/代码/high')
os.mkdir('/Users/xuchangmao/Documents/工作/代码/pems_urban')
os.mkdir('/Users/xuchangmao/Documents/工作/代码/pems_sub')
os.mkdir('/Users/xuchangmao/Documents/工作/代码/pems_high')
os.mkdir('/Users/xuchangmao/Documents/工作/代码/pems_cycles')
# 读取文件
raw_test_data = pd.read_excel(
    r'/Users/xuchangmao/Documents/工作/平台数据/6台国六-10月15-20日的数据/LFNG4HC90LTY00157（15-20.）.xlsx')
# 使用前10小时数据进行分析
test_data = raw_test_data[:360000]

# 下面的3个函数都是对片段进行切割的函数
def group(L):
    '''
    这是一个生成器函数，生成非零片段的区间
    :param L: 传入的是nonzero处理后的字段，类似【1，2，3，4，5，6】
    函数的功能就是将123，456，分开提取出来，形成（1，3），（4，6）这样的效果
    :return: 返回的就是元组的形式
    '''
    first = last = L[0]
    for n in L[1:]:
        if n - 1 == last:
            last = n
        else:
            yield first, last
            first = last = n
    yield first, last


def get_segments(vehicle_speed):
    '''
    采用nonzero函数，来得到全部非零的片段
    :param vehicle_speed:输入dataframe中的车速数据
    :return:返回的是非零的片段
    '''
    test_cycle_array = np.array(vehicle_speed)
    nonzero_index = test_cycle_array.nonzero()
    nonzero_list = list(nonzero_index[0])
    segments_group = list(group(nonzero_list))
    return segments_group


def re_get_segments(raw_segments):
    '''
    通过get_segments函数得到的是不包含怠速片段的速度片段
    接下里需要对这样的raw_segments进行处理，将怠速片段纳入
    :param raw_segments:进行过简单处理的片段列表
    :return:最终的标准的短行程片段
    '''
    new_segments = []
    raw_segments.insert(0, (0, 0))
    for i in range(len(raw_segments)):
        if i != len(raw_segments) - 1:
            new_segments.append((raw_segments[i][-1], raw_segments[i + 1][-1]))
        continue
    return new_segments

'''
接下来是对不同平均速度的区间进行筛选，分到不同的列表中
形成的是3个原始的片段集合，分别代表市区，市郊和高速区间
'''
urban_list = []
sub_list = []
high_list = []

# 水温的筛选函数
def water_temp(segment_tuple):
    for i in range(segment_tuple[0], segment_tuple[1]):
        if test_data['发动机冷却液温度'][i] < 70:
            return True

# 将速度明显较低的短行程片段筛除
def urban_list_filter(segment_tuple):
    if test_data['车速'][segment_tuple[0]:segment_tuple[1]].max() < 5:
        return True

def filter_segments(segments):
    '''
    这是对短行程片段进行预处理的函数，需要从水温，车速的正常性进行判断
    '''
    # 先进行水温的判断，如果水温没有达到70摄氏度，那么就要将该区间排除在外
    segments = [x for x in segments if not water_temp(x)]
    # 按照平均车速进行分类
    for segment in segments:
        average_speed = test_data['车速'][segment[0]:segment[1]].mean()
        if average_speed > 70:
            high_list.append(segment)
        if average_speed > 45 and average_speed < 70 and test_data['车速'][segment[0]:segment[1]].max() < 85:
            sub_list.append(segment)
        if average_speed < 45 and test_data['车速'][segment[0]:segment[1]].max() < 55:
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

def get_cum_work(segments_list):
    '''
    创建一个函数用来计算一个含有不同segments的累积功的多少的
    '''
    segments_list_cum_work = 0
    for segment in segments_list:
        segments_list_cum_work = segments_list_cum_work + integrals[segment[1]] - integrals[segment[0]]
    return segments_list_cum_work

# 要创建一个函数用来计算一个含有不同的segments的时间总和是多少。
def get_cum_time(segments_list):
    '''
    创建一个函数用来计算一个含有不同的segments的时间总和是多少。
    '''
    segments_list_cum_time = 0
    for segment in segments_list:
        segments_list_cum_time = segments_list_cum_time + segment[1] - segment[0]
    return segments_list_cum_time


# 接下来是对3个list中的行程片段进行选择，累积功只要没超过限制那么就一直累加，采用random的形式进行选择
# 首先先定一个值，这个值就是循环功的3分之1
urban_work = 0
sub_work = 0
high_work = 0

# 函数主体运行情况,完成了初识的片段划分，行车了3个list
raw_segments = get_segments(test_data['车速'])
new_segments = re_get_segments(raw_segments)
filter_segments(new_segments)
urban_list = [x for x in urban_list if not urban_list_filter(x)]

# 保存筛选出来的片段
for i in range(len(urban_list)):
    plt.plot(test_data['车速'][urban_list[i][0]:urban_list[i][1]])
    plt.xlabel('time')
    plt.ylabel('vehicle speed')
    plt.title('Urban ' + str(i) + 'th segment')
    filename1 = str(i) + 'th segment' + '.png'
    filepath1 = '/Users/xuchangmao/Documents/工作/代码/urban/'
    plt.savefig(filepath1 + filename1)
    plt.close()
for i in range(len(sub_list)):
    plt.plot(test_data['车速'][sub_list[i][0]:sub_list[i][1]])
    plt.xlabel('time')
    plt.ylabel('vehicle speed')
    plt.title('Sub ' + str(i) + 'th segment')
    filename2 = str(i) + 'th segment' + '.png'
    filepath2 = '/Users/xuchangmao/Documents/工作/代码/sub/'
    plt.savefig(filepath2 + filename2)
    plt.close()
for i in range(len(high_list)):
    plt.plot(test_data['车速'][high_list[i][0]:high_list[i][1]])
    plt.xlabel('time')
    plt.ylabel('vehicle speed')
    plt.title('High ' + str(i) + 'th segment')
    filename3 = str(i) + 'th segment' + '.png'
    filepath3 = '/Users/xuchangmao/Documents/工作/代码/high/'
    plt.savefig(filepath3 + filename3)
    plt.close()
# 接下来是对选择出来的片段进行排序，按照片段时长从小到大的的顺序进行排序，将片段作为
# key,将时长作为value，组合成一个dict,然后按照value的值进行排序处理，最后形成时长
# 从小到大的片段序列。
def get_segment_time_span(segments_list):
    '''
    确定每个list中片段的时长，方便接下来按照时长进行排序。
    :param segments_list:urban list or sub list or high list
    :return: time span
    '''
    time_span = []
    for segment in segments_list:
        time_span.append(segment[1] - segment[0])
    return time_span

urban_time_span = get_segment_time_span(urban_list)
sub_time_span = get_segment_time_span(sub_list)
high_time_span = get_segment_time_span(high_list)
#对片段按照时长进行排序
urban_time_span_dict = dict(zip(urban_list, urban_time_span))
urban_list = sorted(urban_time_span_dict, key=lambda k: urban_time_span_dict[k])
sub_time_span_dict = dict(zip(sub_list, sub_time_span))
sub_list = sorted(sub_time_span_dict, key=lambda k: sub_time_span_dict[k])
high_time_span_dict = dict(zip(high_list, high_time_span))
high_list = sorted(high_time_span_dict, key=lambda k: high_time_span_dict[k])
'''
从这里开始要进行重复多次的pems拼接工作，拼接的次数由输入来定，用cishu来定义
'''
cishu=int(input('请输入需要进行多少次拼接尝试：'))
# 不同区间的比例要求明确的写出来、
urban_criteria = float(input('请输入市区比例：'))
sub_criteria = float(input('请输入市郊比例：'))
high_criteria = float(input('请输入市郊比例：'))
# 先设定总的时间长度
criteria_time =float(input('请输入总时长：'))
# 将最后选择出来参与拼接pems工况的集合定义出来
'''
开始进行多次循环
循环开始之前，将pems_urban,pems_sub等信息用列表定义出来，每进行一次拼接，就将该次拼接后的情况添加
包括，市区累积功，市郊累积功，高速累积功，总累积功，市区时长，市郊时长，高速时长，总时长
'''
pems_urban_work_list=[]
pems_sub_work_list=[]
pems_high_work_list=[]
pems_total_work_list=[]
pems_urban_time_list=[]
pems_sub_time_list=[]
pems_high_time_list=[]
pems_total_time_list=[]
pems_urban_time_ratio_list=[]
pems_sub_time_ratio_list=[]
pems_high_time_ratio_list=[]



for serial_number in range(cishu):
    pems_urban_list = []
    pems_sub_list = []
    pems_high_list = []
    total_list = []
    #先随机选择一次
    random_urban_segment = choice(urban_list)
    random_sub_segment = choice(sub_list)
    random_high_segment = choice(high_list)
    #将选择出来的片段放入工况合集中
    pems_urban_list.append(random_urban_segment)
    pems_sub_list.append(random_sub_segment)
    pems_high_list.append(random_high_segment)

    total_work = get_cum_work(pems_urban_list) + get_cum_work(pems_sub_list) + get_cum_work(pems_high_list)
    # 这种计算方式是采用，先随机选择一次，然后计算比例，比例少的就再随机选择一次，然后保证总累积功不要超过太多

    urban_criteria_time = criteria_time * urban_criteria * 3600
    sub_criteria_time = criteria_time * sub_criteria * 3600
    high_criteria_time = criteria_time * high_criteria * 3600


    # 这是原来的筛选方式，随机选择，然后补充进去。
    urban_time = get_cum_time(pems_urban_list)
    while urban_time < urban_criteria_time:
        pems_urban_list.append(choice(urban_list))
        urban_time = get_cum_time(pems_urban_list)
    sub_time = get_cum_time(pems_sub_list)
    while sub_time < sub_criteria_time:
        pems_sub_list.append(choice(sub_list))
        sub_time = get_cum_time(pems_sub_list)
    high_time = get_cum_time(pems_high_list)
    while high_time < high_criteria_time:
        pems_high_list.append(choice(high_list))
        high_time = get_cum_time(pems_high_list)

    pems_urban_work = get_cum_work(pems_urban_list)
    pems_sub_work = get_cum_work(pems_sub_list)
    pems_high_work = get_cum_work(pems_high_list)
    total_pems_work = pems_urban_work + pems_sub_work + pems_high_work

    pems_urban_ratio = pems_urban_work / total_pems_work
    pems_sub_ratio = pems_sub_work / total_pems_work
    pems_high_ratio = pems_high_work / total_pems_work

    # print('3段行程的累计功分别为：（kw*h）')
    # print(round(pems_urban_work, 2))
    # print(round(pems_sub_work, 2))
    # print(round(pems_high_work, 2))

    # print('3段行程的功率占比分别为：')
    # print('市区功率占比：' + str(round(pems_urban_ratio * 100, 2)) + '%')
    # print('市郊功率占比：' + str(round(pems_sub_ratio * 100, 2)) + '%')
    # print('高速功率占比：' + str(round(pems_high_ratio * 100, 2)) + '%')
    pems_urban_time=get_cum_time(pems_urban_list)
    pems_sub_time=get_cum_time(pems_sub_list)
    pems_high_time=get_cum_time(pems_high_list)
    total_time = pems_urban_time + pems_sub_time + pems_high_time
    pems_urban_ratio_time = pems_urban_time/ total_time
    pems_sub_ratio_time = pems_sub_time / total_time
    pems_high_ratio_time = pems_high_time / total_time
    # print('3段行程的时间占比分别为：')
    # print('市区时间占比：' + str(round(pems_urban_ratio_time * 100, 2)) + '%')
    # print('市郊时间占比：' + str(round(pems_sub_ratio_time * 100, 2)) + '%')
    # print('高速时间占比：' + str(round(pems_high_ratio_time * 100, 2)) + '%')
    pems_urban_work_list.append(round(pems_urban_work, 2))
    pems_sub_work_list.append(round(pems_sub_work, 2))
    pems_high_work_list.append(round(pems_high_work, 2))
    pems_total_work_list.append(round(pems_urban_work, 2)+round(pems_sub_work, 2)+round(pems_high_ratio_time * 100, 2))
    pems_urban_time_list.append(pems_urban_time)
    pems_sub_time_list.append(pems_sub_time)
    pems_high_time_list.append(pems_high_time)
    pems_total_time_list.append(total_time)
    pems_urban_time_ratio_list.append(round(pems_urban_ratio_time * 100, 2))
    pems_sub_time_ratio_list.append(round(pems_sub_ratio_time * 100, 2))
    pems_high_time_ratio_list.append(round(pems_high_ratio_time * 100, 2))
# 保存选择出来的片段,仅仅保留总的pems工况的情况
    pems_cycle_list = pems_urban_list + pems_sub_list + pems_high_list
    pems_cycle_pd_list = []
    for segment in pems_cycle_list:
        pems_cycle_pd_list.append(test_data[segment[0]:segment[1]])
    pems_cycle_pd = pd.concat(pems_cycle_pd_list)
    pems_cycle_pd = pems_cycle_pd.reset_index(drop=True)
    plt.plot(pems_cycle_pd['车速'])
    plt.xlabel('time')
    plt.ylabel('vehicle speed')
    plt.title(str(serial_number)+'th '+'PEMS Cycle')
    plt.savefig('/Users/xuchangmao/Documents/工作/代码/pems_cycles/'+str(serial_number)+'th '+'pems_cycle.png')
    plt.close()
    outpath = '/Users/xuchangmao/Documents/工作/代码/pems_cycles/'+str(serial_number)+'th '+'pems_data.csv'
    pems_cycle_pd.to_csv(outpath, sep=',', index=False)

# 计算一下每个区间的片段的时间长度
def get_segment_time_span(segments_list):
    time_span = []
    for segment in segments_list:
        time_span.append(segment[1] - segment[0])
    return time_span

# 计算每个区间的片段的最高车速
def get_segments_max_speed(segments_list):
    max_speed = []
    for segment in segments_list:
        max_speed.append(test_data['车速'][segment[0]:segment[1]].max())
    return max_speed

# 下面是3个区间的信息，包括时间长度和最高速度
    urban_segments_pd = pd.DataFrame(urban_list)
    urban_time_span = get_segment_time_span(urban_list)
    urban_segments_pd['time_span'] = urban_time_span
    urban_segments_pd['top speed'] = get_segments_max_speed(urban_list)
    sub_segments_pd = pd.DataFrame(sub_list)
    sub_time_span = get_segment_time_span(sub_list)
    sub_segments_pd['time_span'] = sub_time_span
    sub_segments_pd['top speed'] = get_segments_max_speed(sub_list)
    high_segments_pd = pd.DataFrame(high_list)
    high_time_span = get_segment_time_span(high_list)
    high_segments_pd['time_span'] = high_time_span
    high_segments_pd['top speed'] = get_segments_max_speed(high_list)
info_data={
    '市区累积功':pems_urban_work_list,
    '市郊累积功':pems_sub_work_list,
    '高速累积功':pems_high_work_list,
    '总累积功':pems_total_work_list,
    '市区时长':pems_urban_time_list,
    '市郊时长': pems_sub_time_list,
    '高速时长': pems_high_time_list,
    '总时长': pems_total_time_list,
}
info_data_pd=pd.DataFrame(info_data)
info_data_pd.to_csv('/Users/xuchangmao/Documents/工作/代码/pems_cycles/拼接工况数据.csv')
