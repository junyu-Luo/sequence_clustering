from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
import numpy as np
import math as m
# import matplotlib.pyplot as plt
import Levenshtein
import pandas as pd
import re

def get_dis_matrix(data):
    """
    获得邻接矩阵
    :param data: 样本集合
    :return: 邻接矩阵
    """
    nPoint = len(data)
    dis_matrix = np.zeros((nPoint, nPoint))

    for i in range(nPoint):
        for j in range(i + 1, nPoint):
            dis_matrix[i][j] = dis_matrix[j][i] = m.sqrt(np.power(data[i] - data[j], 2).sum())
    return dis_matrix


def getW(dicts):
    """
    利用 Levenshtein 距离 获得相似矩阵
    :param data: 样本集合
    :return:
    """
    Matrix = np.zeros((len(dicts), len(dicts)))
    for i, no_i in enumerate(dicts):
        for j, no_j in enumerate(dicts):
            if i == j:
                pass
            else:
                Matrix[i][j] = Levenshtein.seqratio(dicts[no_i], dicts[no_j])

    return Matrix


def getD(W):
    """
    获得度矩阵
    :param W:  相似度矩阵
    :return:   度矩阵
    """
    D = np.diag(sum(W))
    return D


def getL(D, W):
    """
    获得拉普拉斯矩阵
    :param W: 相似度矩阵
    :param D: 度矩阵
    :return: 拉普拉斯矩阵
    """
    return D - W


def getEigen(L):
    """
    从拉普拉斯矩阵获得特征(映射)矩阵
    :param L: 拉普拉斯矩阵
    :return:
    """
    eigval, eigvec = np.linalg.eig(L)
    # 取特征向量的前 4 个，但第一个都一样，可以去掉
    ix = np.argsort(eigval)[0:4]
    return eigvec[:, ix]


def plotRes(data, clusterResult, clusterNum):
    """
    结果可视化
    :param data:  样本集
    :param clusterResult: 聚类结果
    :param clusterNum:  聚类个数
    :return:
    """
    nPoints = len(data)
    scatterColors = ['black', 'blue', 'green', 'yellow', 'red', 'purple', 'orange']
    for i in range(clusterNum):
        color = scatterColors[i % len(scatterColors)]
        x1 = [];  y1 = []
        for j in range(nPoints):
            if clusterResult[j] == i:
                x1.append(data[j, 0])
                y1.append(data[j, 1])
        plt.scatter(x1, y1, c=color, alpha=1, marker='+')
    plt.show()


def find_index(lists, find):
    """
    在lists里面，找出find在lists的索引
    :param lists: 列表
    :param find: 列表中的一个元素
    :return: 一个在lists里面，找出find在lists的索引列表
    """
    flag = 0
    list_index = []
    for n in range(lists.count(find)):
        sec = flag
        flag = lists[flag:].index(find)
        list_index.append(flag + sec)
        flag = list_index[-1:][0] + 1
    return list_index


def write_file(str):
    """
    写入文件
    :param str: 字符串
    :return: 无
    """
    writefile = open("./out/前3个特征并用keams.txt", 'a+',encoding='utf-8')
    writefile.write(str + '\n')
    writefile.close()
    # return str



def get_dicts(data):
    out_list = []

    PATIENT_ID_list = []
    INPATIENT_NO_list = []
    ITEM_CODE_list = []
    DATE_BGN_list = []
    MO_NAME_list = []

    for i in range(len(data)):
        PATIENT_ID_list.append(data['PATIENT_ID'][i])
        INPATIENT_NO_list.append(data['INPATIENT_NO'][i])
        ITEM_CODE_list.append(data['ONE_HOT'][i])
        DATE_BGN_list.append(data['TIME_START'][i])
        MO_NAME_list.append(data['MO_NAME'][i])

    PATIENT_ID_seq = list(enumerate(PATIENT_ID_list))
    # ITEM_CODE_seq = list(enumerate(ITEM_CODE_list))
    DATE_BGN_seq = list(enumerate(DATE_BGN_list))

    # INPATIENT_NO_seq = list(enumerate(INPATIENT_NO_list))
    # print(PATIENT_ID_seq[0][1])
    # print(INPATIENT_NO_seq)
    # print(INPATIENT_NO_seq[0][1])

    INPATIENT_NO_seq = {}
    INPATIENT_NO_only_list = list(set(INPATIENT_NO_list))
    for Only_NO in INPATIENT_NO_only_list:
        tem_list = []
        for i, NO in enumerate(INPATIENT_NO_list):
            if Only_NO == NO:
                tem_list.append(i)
        INPATIENT_NO_seq[Only_NO] = tem_list
    # print(INPATIENT_NO_seq)

    no_i = 0
    for Only_NO in INPATIENT_NO_only_list:
        dicts = {}
        index_no = INPATIENT_NO_seq[Only_NO]
        # print(Only_NO,index_no) # 流水号：102803 [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
        time_list = []
        for index_i in index_no:
            time_list.append(DATE_BGN_seq[index_i][1])
        zip_time_index_list = list(zip(time_list, index_no))
        seq_ziped = sorted(zip_time_index_list)
        # print(seq_ziped)

        decompression = list(zip(*seq_ziped))
        decompression_time_list = list(decompression[0])
        decompression_seq_list = list(decompression[1])
        # print(Only_NO,decompression_time_list,decompression_seq_list)
        repeat_list = sorted(list(set(decompression_time_list)))
        seq_list = []
        for repeat in repeat_list:
            repeat_index = find_index(decompression_time_list, repeat)
            index_index_list = []
            for index in repeat_index:
                index_index_list.append(decompression_seq_list[index])
            seq_list.append(index_index_list)

        # print(seq_list)
        # print(ITEM_CODE_list)

        mo_name_list = []
        item_code_seq_list = []
        for seq in seq_list:
            code_tem_list = []
            name_tem_list = []
            for index_j in seq:
                code_tem_list.append(ITEM_CODE_list[index_j])
                name_tem_list.append(MO_NAME_list[index_j])
            item_code_seq_list.append("".join(sorted(list(set(code_tem_list)))))
            mo_name_list.append(name_tem_list)
        # print(item_code_seq_list)

        dicts["no"] = no_i
        dicts["INPATIENT_NO"] = Only_NO
        dicts["PATIENT_ID"] = PATIENT_ID_seq[index_no[0]][1]
        dicts["ONE_HOT"] = item_code_seq_list
        dicts['MO_NAME'] = mo_name_list
        no_i += 1
        out_list.append(dicts)

    return out_list



if __name__ == '__main__':
    # medicine_dicts = {'A': 'ICS', 'B': 'SABA', 'C': 'LABA', 'D': 'SAMA', 'E': 'LAMA', 'F': 'ICS+LABA', 'G': 'LABA+LAMA',
    #  'H': 'SABA+SAMA', 'I': '全身用糖皮质激素', 'J': '甲基黄嘌呤类', 'K': '抗过敏药', 'L': '镇咳祛痰类', 'M': '大环内酯类药', 'N': '抗感染药（不含大环内酯类）',
    #  'O': '肾素-血管紧张素系统药物', 'P': '钙通道阻滞剂', 'Q': 'β受体阻滞剂', 'R': '降血压药', 'S': '利尿剂', 'T': '防治心绞痛药', 'U': '心律失常药',
    #  'V': '治疗慢性心功能不全药物', 'W': '抗休克的血管活性药', 'X': '降糖药', 'Y': '调脂类药', 'Z': '抗凝血药', 'a': '治疗消化性溃疡和胃食管反流病药物',
    #  'b': '解热镇痛抗炎药', 'c': '抗肿瘤药', 'd': '免疫调节剂'}

    data = pd.read_excel('./excel/dataset.xlsx')
    dicts = get_dicts(data)
    print(dicts)
    test_dict = {}
    for dict in dicts:
        test_dict[dict["no"]] = dict['ONE_HOT']
    print(len(test_dict),test_dict)
    # dicts = {'NO_0': ['0', '1','24'], 'NO_1': ['lqbz', 'gaest', 'arwhht'], 'NO_2': ['2', '1', '0'],
    #          'NO_4': ['0', '1', '2', '3', '4']}
    # cluster_num = 4
    # W = getW(test_dict)
    # D = getD(W)
    # L = getL(D, W)
    # eigvec = getEigen(L)
    #
    # print(eigvec)

    for cluster_num in range(50,200):
        print((cluster_num - 50)/150,cluster_num)
        W = getW(test_dict)
        D = getD(W)
        L = getL(D, W)
        eigvec = getEigen(L)
        # print(eigvec)
#######################################################################################################################
        # 使用kmeans 聚类 选择其中一个必须注释其他
        clf = KMeans(n_clusters=cluster_num)
        s = clf.fit(eigvec)
        C = s.labels_
#######################################################################################################################
        # 使用 GMM 聚类  选择其中一个必须注释其他
        # clf = GaussianMixture(n_components=cluster_num, covariance_type='full')
        # g = clf.fit(eigvec)
        # C = g.predict(eigvec)

#######################################################################################################################
        write_file("*********【分{}类】*********".format(cluster_num))
        # print("分{}类".format(cluster_num))

        for i in range(cluster_num):
            cluster_index = find_index(list(C),i)

            tem_list = []
            for index in cluster_index:
                tem_list.append(test_dict[index])
            # print("第{}类".format(i+1),)
            write_file("…………【第{}类】…………".format(i+1),)
            # print(tem_list)
            # show_list = []
            # for elem_list in tem_list:
            #     tem_str1_list = []
            #     for string in elem_list:
            #         str_to_list = list(string)
            #         tem_str2_list = []
            #         for change in str_to_list:
            #             tem_str2_list.append(medicine_dicts[change])
            #         tem_str1_list.append(tem_str2_list)
            #     show_list.append(tem_str1_list)
            # print(tem_list)
            for print_str in tem_list:
                write_file(str(print_str))
                write_file('---------------------------------')
            write_file('==================================================================')


        # print("---------------------------------------------------------------------------------")




