# -*- coding: utf-8 -*-
# @Time    : 2023/4.1
# @Author  : Zhaike
# @project : H9N1 AIV
"""
Spyder Editor

"""

"""
####计算2511训练集测试集特征值####
import math 
import pandas as pd
from numpy import *
import time
import heapq
from multiprocessing import Pool
import os
os.chdir('/media/dulab/file/H9N2update/final')
#提取菌株对名称1+序列1+菌株对名称2+序列2
df=pd.read_csv("8.2511HA1HI.csv")
seq_dict1 = df.set_index('Virus').to_dict()['Seq1']
seq_names1 = df.Virus.to_list()
seq_dict2 = df.set_index('Serum').to_dict()['Seq2']
seq_names2 = df.Serum.to_list()

data=pd.read_csv("6.9548HA1NO2.csv")
Name = list(data['ID'])
LocN = [eval(x) for x in data['locN1']]
LocO = [eval(y) for y in data['locO1']]

#计算六个抗原表位差异
#scannet 96 6类
H9_A = [40, 45, 46, 48, 254, 267, 268, 271, 275, 276, 295]
H9_B = [176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 204, 207]
H9_C = [120, 123, 124, 125, 126, 127, 128, 143, 145, 146, 147, 148, 149, 150, 151, 152, 153, 155, 236]
H9_D = [115, 159, 161, 162, 163, 164, 165, 197, 198, 228, 229, 230, 231, 232, 234, 251]
H9_E = [54, 65, 66, 68, 69, 70, 72, 84, 86, 87, 89, 97, 139]
H9_F = [90, 91, 92, 129, 130, 131, 132, 133, 134, 135, 136, 137, 210, 211, 212, 213, 214, 215, 216, 217]

#计算抗原表位的汉明距离
def get_epitopeA(seq, epi_list=H9_A):
    '''给出一条序列，返回该序列A抗原区域的位点'''
    
    return [seq[base_idx] for base_idx in epi_list]


def get_epitopeB(seq, epi_list=H9_B):
    '''给出一条序列，返回该序列B抗原区域的位点'''
    
    return [seq[base_idx] for base_idx in epi_list]


def get_epitopeC(seq, epi_list=H9_C):
    '''给出一条序列，返回该序列C抗原区域的位点'''
    
    return [seq[base_idx] for base_idx in epi_list]

def get_epitopeD(seq, epi_list=H9_D):
    '''给出一条序列，返回该序列D抗原区域的位点'''
    
    return [seq[base_idx] for base_idx in epi_list]

def get_epitopeE(seq, epi_list=H9_E):
    '''给出一条序列，返回该序列E抗原区域的位点'''
    
    return [seq[base_idx] for base_idx in epi_list]


def get_epitopeF(seq, epi_list=H9_F):
    '''给出一条序列，返回该序列E抗原区域的位点'''
    
    return [seq[base_idx] for base_idx in epi_list]


def cal_epitope(base_set1, base_set2):
    '''给出两个抗原区域位点的列表，计算汉明距离'''

    return sum([a != b for a, b in zip(base_set1, base_set2)])
    

#计算五个氨基酸特征值
df4=pd.read_csv('0.character_index4.csv',header = None)
df4 = df4.T#转置
df4.columns= df4.iloc[0]#设置列名
df4 = df4.reindex(df4.index.drop(0))#去掉第一行
df4=df4[['Character_name','Hydrophobicity index (Fasman, 1989)', 'Volume', 'Isoelectric point', 
'Polarizability parameter', 'Average accessible surface area']]

def cal_aa(seq1, seq2):
    '''给出两条序列，5个字典，返回5个氨基酸特征平均改变量'''
    
    n=[]
    for value in dict1.values():
        f=[]
        for base_i, base_j in zip(seq1, seq2):
            f.append(abs(float(value[base_i])-float(value[base_j]))) 
        aa=sum(heapq.nlargest(3,f))/3
        n.append(aa)
    return n


#计算糖基化位点
def get_netno(name):
    '''输入菌株名，返回其糖基化位点列表'''
    
    for i in range(len(data['ID'])):
        if name==Name[i]:
            return LocN[i],LocO[i]


def cal_netno(name1,name2):
    '''给出一对菌株名，返回NO糖基化位点差异数量'''
    
    numN = len(list(set(get_netno(name1)[0])-set(get_netno(name2)[0]))+list(set(get_netno(name2)[0])-set(get_netno(name1)[0])))
    numO = len(list(set(get_netno(name1)[1])-set(get_netno(name2)[1]))+list(set(get_netno(name2)[1])-set(get_netno(name1)[1])))
    num = numN + numO
    return [numN,numO]


#计算受体结合位点

df1 = pd.read_csv("0.1jsd_CA.csv")#LOC 1-D 序列 0-D
df1 = df1[df1['E']=='A']#只要HA1 317
df1['loc'] = df1['LOC']-1 #匹配序列中残基编号和坐标编号
df3 = df1[['loc','X','Y','Z']]
df3 = df3.set_index('loc')#loc为索引对于序列残基编号
rbs = [128,129,130,131,132,91,143,145,173,180,184,185,214,215,216,217,218,219]#受体结合位点的H9成熟编号
rbs_list = [x-1 for x in rbs]#在序列列表中编号再减一，已对比无误 GTSKA YWTNVLY NGLMGR

def hamm(s1, s2):
    '''计算汉明距离'''
    
    return sum([a != b for a, b in zip(s1, s2)])


def distance(a,b):
    '''计算欧式距离'''
    dist=math.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2)
    
    return dist


def cal_rbs(list1,list2):
    '''输入菌株对的序列列表，返回最短欧式距离'''
    
    if hamm(list1,list2)==0:
        e = 0
    else:
        d=[]
        for r in range(0,317):
            if list1[r]!=list2[r]:
                d1=[]
                for o in rbs_list:
                    d1.append(distance(df3.loc[r],df3.loc[o]))#1个改变位点到18个受体结合位点的18个距离
                d.append(min(d1))#改变位点到受体结合位点的距离最小值的集合
        if len(d)<=3:
            e = mean(d)
        else:
            e = mean(heapq.nsmallest(3,d))#前三个最短欧氏距离取平均
        
    return [e]

lst=[]  
for name1,name2 in zip(seq_names1, seq_names2):
    Seq1 = seq_dict1[name1]
    Seq2 = seq_dict2[name2]

    H9_1 = [cal_epitope(get_epitopeA(Seq1),get_epitopeA(Seq2))]
    H9_2 = [cal_epitope(get_epitopeB(Seq1),get_epitopeB(Seq2))]
    H9_3 = [cal_epitope(get_epitopeC(Seq1),get_epitopeC(Seq2))]
    H9_4 = [cal_epitope(get_epitopeD(Seq1),get_epitopeD(Seq2))]
    H9_5 = [cal_epitope(get_epitopeE(Seq1),get_epitopeE(Seq2))]
    H9_6 = [cal_epitope(get_epitopeF(Seq1),get_epitopeF(Seq2))]

    line = H9_1 + H9_2 + H9_3 + H9_4 + H9_5 + H9_6 + cal_aa(Seq1,Seq2) + cal_netno(name1,name2) + cal_rbs(Seq1,Seq2)
    lst.append(line)
    
#列表转fataframe
from pandas.core.frame import DataFrame
data=DataFrame(lst)
data.columns=['H9-A','H9-B','H9-C','H9-D','H9-E','H9-F',
                  'Hydrophobicity', 'Volume', 'Charge','Polarity', 'Accessible surface area',
                  'N-Glycosylation','O-Glycosylation','Receptor binding']

df1=pd.read_csv("8.2511HA1HI.csv")
data['Y']=df1['Y']
data.to_csv('9.2511data .csv',index=None)

"""

####训练模型####随机森林和xgboost原始数据，svm和knn标准化数据（0-1）
import pandas as pd
from sklearn.tree import DecisionTreeClassifier# 决策树分类器
from sklearn.linear_model import LogisticRegression# 逻辑回归分类器
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
import xgboost as xgb
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier as KNN
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
from sklearn.metrics import auc    #AUC
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve
from xgboost import plot_importance
import os
os.chdir('/media/dulab/file/H9N2update/final')

#训练集拟合/交叉验证
data = pd.read_csv("9.2511data.csv")
from sklearn import metrics
from sklearn.model_selection import GridSearchCV, KFold, cross_val_score
from sklearn.model_selection import train_test_split
from imblearn.over_sampling import SMOTE
X = data.iloc[:,:-1]#除了最后一列
Y = data.iloc[:,-1:]#最后一列
#过采样1685+1685
x_resampled, y_resampled = SMOTE(random_state=2022).fit_resample(X,Y)
print(y_resampled['Y'].agg(['value_counts']).T)

#分割训练集和测试集2359 1011-840
train_x, test_x, train_y, test_y = train_test_split(x_resampled, y_resampled, test_size = 0.3,random_state=123)

train_x['Y']=train_y['Y']#1757
test_x['Y']=test_y['Y']#754
lt = list(train_x.columns)
#查看测试集中数据泄漏情况155条
df = test_x.merge(train_x, how='left', indicator=True, left_on=lt, right_on=lt)#1027
test = df[df['_merge']=='left_only']#599
#数据泄漏问题解决后重建训练集和测试集
train_x = train_x.iloc[:,:-1]
test_x = test.iloc[:,:-2]
test_y = test['Y']


clf = XGBClassifier(
     learning_rate= 0.1, 
     n_estimators= 230, 
     max_depth= 8,
     min_child_weight= 4,
     seed=0,
     subsample= 0.6, 
     colsample_bytree= 0.7, 
     gamma= 0.3,
     reg_alpha= 0, 
     reg_lambda= 1)

clf.fit(train_x, train_y)

# 在训练集和测试集上分布利用训练好的模型进行预测
train_predict = clf.predict(train_x)
test_predict = clf.predict(test_x)

#y_pred = clf.predict(test_x)
y_pro = clf.predict_proba(test_x)
predictions = [round(value) for value in test_predict]

cv = KFold(n_splits=5, shuffle=True, random_state=2022).split(train_x)
print('训练集accuracy:',accuracy_score(train_y,train_predict))
print('测试集accuracy:',accuracy_score(test_y,test_predict))
print('交叉验证accuracy:',cross_val_score(clf, train_x, train_y, cv=cv,
                                            scoring='accuracy').mean())
cv = KFold(n_splits=5,shuffle=True, random_state=2022).split(train_x)
print('训练集precision:',metrics.precision_score(train_y,train_predict))
print('测试集precision:',metrics.precision_score(test_y,test_predict))
print('交叉验证precision:',cross_val_score(clf, train_x, train_y, cv=cv, n_jobs=-1,
                                            scoring='precision').mean())

cv = KFold(n_splits=5,shuffle=True, random_state=2022).split(train_x)
print('训练集recall:',metrics.recall_score(train_y,train_predict))
print('测试集recall:',metrics.recall_score(test_y,test_predict))
print('交叉验证recall:',cross_val_score(clf, train_x, train_y, cv=cv, n_jobs=-1,
                                            scoring='recall').mean())

cv = KFold(n_splits=5, shuffle=True, random_state=2022).split(train_x)
print('测试集roc_auc:',metrics.roc_auc_score(test_y,test_predict))
print('交叉验证roc_auc:',cross_val_score(clf, train_x, train_y, cv=cv,
                                                scoring='roc_auc').mean())

cv = KFold(n_splits=5,shuffle=True, random_state=2022).split(train_x)
print('训练集f1:',metrics.f1_score(train_y,train_predict))
print('测试集f1:',metrics.f1_score(test_y,test_predict))
print('交叉验证f1:',cross_val_score(clf, train_x, train_y, cv=cv, n_jobs=-1,
                                           scoring='f1').mean())


print('测试集accuracy:',accuracy_score(test_y,test_predict))
print('测试集precision:',metrics.precision_score(test_y,test_predict))
print('测试集recall:',metrics.recall_score(test_y,test_predict))
print('测试集f1:',metrics.f1_score(test_y,test_predict))
print('测试集roc_auc:',metrics.roc_auc_score(test_y,test_predict))

####网格法调参####
import xgboost as xgb  #直接引用xgboost。接下来会用到其中的“cv”函数。
from xgboost.sklearn import XGBClassifier #是xgboost的sklearn包，允许我们像GBM一样使用Grid Search和并行处理
from sklearn import model_selection, metrics   #Additional     scklearn functions
from sklearn.model_selection import GridSearchCV   #Perforing grid search
from sklearn.model_selection import train_test_split
import matplotlib.pylab as plt
%matplotlib inline

data = pd.read_csv("9.2511data.csv")
X = data.iloc[:,:-1]#除了最后一列
Y = data.iloc[:,-1:]#最后一列
#过采样1685+1685
x_resampled, y_resampled = SMOTE(random_state=2022).fit_resample(X,Y)
print(y_resampled['Y'].agg(['value_counts']).T)

#分割训练集和测试集2359 1011-840
train_x, test_x, train_y, test_y = train_test_split(x_resampled, y_resampled, test_size = 0.3,random_state=123)
#910,920,930,940,950,960,970,980,990,1000,1010,1020,1030,1040,1050,1060,1070,1080,1090
#310,320,330,340,350,360,370,380,390,400,410,420,430,440,450,460,470,480,490
#100,101,102,103,104,105,106,107,108,109,110,110,120,130,140,150,160,170,180,190,200
#10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200
#先固定其他参数取值，调整n_estimators 100,200,300,400,500,600,700,800,900,1000  100 0.7816413547237077
cv_params = {'n_estimators': [110,120,130,140,150,160,170,180,190,200,210,220,230,240,250,260,270,280,290,300]}
other_params = {'learning_rate': 0.1, 'max_depth': 5, 'min_child_weight':1, 'seed': 0,
                    'subsample': 0.8, 'colsample_bytree': 0.8, 'gamma': 0, 'reg_alpha': 0, 
                'reg_lambda': 1,'scale_pos_weight' : 1}
model = xgb.XGBClassifier(**other_params)
optimized_GBM = GridSearchCV(estimator=model, param_grid=cv_params, scoring='accuracy', cv=5, verbose=1, n_jobs=4)
optimized_GBM.fit(train_x, train_y)
evalute_result = optimized_GBM.cv_results_
#print('每轮迭代运行结果:{0}'.format(evalute_result))
print('参数的最佳取值：{0}'.format(optimized_GBM.best_params_))
print('最佳模型得分:{0}'.format(optimized_GBM.best_score_))

#参数的最佳取值：{'n_estimators': 230} 最佳模型得分:0.840612292633776
#参数的最佳取值：{'n_estimators': 340} 最佳模型得分:0.8414597502608947
#参数的最佳取值：{'n_estimators': 1000} 最佳模型得分:0.8427336356112132


#调整max_depth和min_child_weight 
cv_params = {'max_depth': [1,2,3,4,5,6,7,8,9], 'min_child_weight': [1,2,3,4,5,6,7,8,9]}
other_params = {'learning_rate': 0.1, 'n_estimators': 230, 'seed': 0,
                    'subsample': 0.8, 'colsample_bytree': 0.8, 'gamma': 0, 'reg_alpha': 0, 'reg_lambda': 1}
model = xgb.XGBClassifier(**other_params)
optimized_GBM = GridSearchCV(estimator=model, param_grid=cv_params, scoring='accuracy', cv=5, verbose=1, n_jobs=4)
optimized_GBM.fit(train_x, train_y)
evalute_result = optimized_GBM.cv_results_
#print('每轮迭代运行结果:{0}'.format(evalute_result))
print('参数的最佳取值：{0}'.format(optimized_GBM.best_params_))
print('最佳模型得分:{0}'.format(optimized_GBM.best_score_))

#参数的最佳取值：{'max_depth': 8, 'min_child_weight': 4} 最佳模型得分:0.8452706106732879

#调整gamma
cv_params = {'gamma': [0,0.1,0.2,0.3,0.4,0.5,1]}
other_params = {'learning_rate': 0.1, 'n_estimators': 230, 'max_depth':8, 'min_child_weight':4, 'seed': 0,
                    'subsample': 0.8, 'colsample_bytree': 0.8,  'reg_alpha': 0, 'reg_lambda': 1}
model = xgb.XGBClassifier(**other_params)
optimized_GBM = GridSearchCV(estimator=model, param_grid=cv_params, scoring='accuracy', cv=5, verbose=1, n_jobs=4)
optimized_GBM.fit(train_x, train_y)
evalute_result = optimized_GBM.cv_results_
#print('每轮迭代运行结果:{0}'.format(evalute_result))
print('参数的最佳取值：{0}'.format(optimized_GBM.best_params_))
print('最佳模型得分:{0}'.format(optimized_GBM.best_score_))

#参数的最佳取值：{'gamma': 0.3} 最佳模型得分:0.8461198675735003


#调整subsample colsample_bytree
cv_params = {'subsample': [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9], 'colsample_bytree': [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]}
other_params = {'learning_rate': 0.1, 'n_estimators':230, 'max_depth': 8, 'min_child_weight':4, 'seed': 0,
                    'gamma': 0.3, 'reg_alpha': 0, 'reg_lambda': 1}
model = xgb.XGBClassifier(**other_params)
optimized_GBM = GridSearchCV(estimator=model, param_grid=cv_params, scoring='accuracy', cv=5, verbose=1, n_jobs=4)
optimized_GBM.fit(train_x, train_y)
evalute_result = optimized_GBM.cv_results_
#print('每轮迭代运行结果:{0}'.format(evalute_result))
print('参数的最佳取值：{0}'.format(optimized_GBM.best_params_))
print('最佳模型得分:{0}'.format(optimized_GBM.best_score_))

#参数的最佳取值：{'colsample_bytree': 0.7, 'subsample': 0.6} 最佳模型得分:0.8486667386375905


#调整reg_alpha reg_lambda
cv_params = {'reg_alpha': [0, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 3, 4, 5], 'reg_lambda': [0,0.01,0.05, 0.1,0.5, 1, 2, 3,4,5]}
other_params = {'learning_rate': 0.1, 'n_estimators':230, 'max_depth': 8, 'min_child_weight': 4, 'seed': 0,
                    'subsample': 0.6, 'colsample_bytree': 0.7, 'gamma': 0.3}
model = xgb.XGBClassifier(**other_params)
optimized_GBM = GridSearchCV(estimator=model, param_grid=cv_params, scoring='accuracy', cv=5, verbose=1, n_jobs=4)
optimized_GBM.fit(train_x, train_y)
evalute_result = optimized_GBM.cv_results_
#print('每轮迭代运行结果:{0}'.format(evalute_result))
print('参数的最佳取值：{0}'.format(optimized_GBM.best_params_))
print('最佳模型得分:{0}'.format(optimized_GBM.best_score_))

#参数的最佳取值：{'reg_alpha': 0, 'reg_lambda': 1} 最佳模型得分:0.8486667386375905

#调整learning_rate
cv_params = {'learning_rate': [0.01,0.05,0.1,0.15,0.2,0.3,0.4,0.5]}
other_params = { 'n_estimators': 230, 'max_depth': 8, 'min_child_weight':4, 'seed': 0,
                    'subsample': 0.6, 'colsample_bytree': 0.7, 'gamma': 0.3, 'reg_alpha': 0, 'reg_lambda': 1}
model = xgb.XGBClassifier(**other_params)
optimized_GBM = GridSearchCV(estimator=model, param_grid=cv_params, scoring='accuracy', cv=5, verbose=1, n_jobs=4)
optimized_GBM.fit(train_x, train_y)
evalute_result = optimized_GBM.cv_results_
#print('每轮迭代运行结果:{0}'.format(evalute_result))
print('参数的最佳取值：{0}'.format(optimized_GBM.best_params_))
print('最佳模型得分:{0}'.format(optimized_GBM.best_score_))

#参数的最佳取值：{'learning_rate': 0.1} 最佳模型得分:0.8486667386375905

#最终模型参数
#保存模型

data = pd.read_csv("9.2511data.csv")
X = data.iloc[:,:-1]#除了最后一列
Y = data.iloc[:,-1:]#最后一列
#过采样1685+1685
x_resampled, y_resampled = SMOTE(random_state=2022).fit_resample(X,Y)
print(y_resampled['Y'].agg(['value_counts']).T)

#分割训练集和测试集2359 1011-840
train_x, test_x, train_y, test_y = train_test_split(x_resampled, y_resampled, test_size = 0.3,random_state=123)

xgb1 = XGBClassifier(
     learning_rate= 0.1, 
     n_estimators= 230, 
     max_depth= 8,
     min_child_weight= 4,
     seed=0,
     subsample= 0.6, 
     colsample_bytree= 0.7, 
     gamma= 0.3,
     reg_alpha= 0, 
     reg_lambda= 1)
xgb1.fit(train_x,train_y)
#保存模型
import pickle
pickle.dump(xgb1, open("finalmodel230.pickle.dat", "wb"))



####SHAP可解释性概要图####分开输出
import shap  
import matplotlib.pyplot as plt
from cycler import cycler
explainer = shap.TreeExplainer(xgb1)
shap_values = explainer.shap_values(train_x)#计算每个样本的每个特征的SHAP值
# 特征统计值
#汇总图
fig= shap.summary_plot(shap_values, train_x,show=False, plot_type="bar")
plt.savefig('SHAPf1.pdf',bbox_inches='tight')
#天际线图
fig=shap.summary_plot(shap_values, train_x,show=False)
plt.savefig('SHAPf2.pdf',bbox_inches='tight')


# 计算每个特征的绝对值平均 SHAP 值
mean_abs_shap = np.mean(np.abs(shap_values), axis=0)
# 对 SHAP 值进行排序
sorted_index = np.argsort(mean_abs_shap)[::-1]
sorted_shap = mean_abs_shap[sorted_index]
# 绘制条形图
fig, ax = plt.subplots()
y_pos = np.arange(train_x.shape[1])
ax.barh(y_pos, sorted_shap, align='center')
ax.set_yticks(y_pos)
#ax.set_yticklabels(sorted_names)
for i, v in enumerate(sorted_shap):
    ax.text(v + 0.01, i, "{:.2f}".format(v))
plt.savefig('SHAPf3.pdf',bbox_inches='tight')   
plt.show()


# 获取特征的SHAP值和特征名称2359
shap_df = pd.DataFrame(shap_values, columns=train_x.columns)
shap_df.to_csv('SHAPvalue2.csv',index=None)

# 计算每个特征的SHAP值占总SHAP值的比例
shap_sum = shap_df.abs().sum(axis=0)
feature_shap_ratio = shap_sum / shap_sum.sum()

# 生成饼图
# 定义经典的lancet配色方案
color_cycle = (cycler(color=['#B2DF8A', '#99CC00','#33A02C', '#FB9A99', '#E31A1C',
                            '#FDBF6F', '#FF7F00', '#CAB2D6', '#6A3D9A', '#FFFF99', '#B15928',
                            '#CCCCCC','#1F78B4','#A6CEE3']))
font_title={'family':'Arial', 'size':7,'weight':'bold'}
font_lable = {'family':'Arial', 'size':7}
fig, ax = plt.subplots(figsize=(6, 6))
ax.set_prop_cycle(color_cycle)
ax.pie(feature_shap_ratio, labels=feature_shap_ratio.index, autopct='%1.1f%%', startangle=90, counterclock=False,textprops=font_lable)
ax.set_title("Feature SHAP Ratio",fontdict=font_title)
# 保存饼图为PDF文件
plt.savefig("feature_shap_ratio.pdf")
plt.show()


#输出百分比贡献值以及气泡图
explainer = shap.TreeExplainer(xgb1)
shap_values = explainer.shap_values(train_x)#计算每个样本的每个特征的SHAP值
# 计算每个特征的百分比贡献值
shap_df = pd.DataFrame(shap_values, columns=train_x.columns)
# 计算每个特征的SHAP值占总SHAP值的比例
shap_sum = shap_df.abs().sum(axis=0)
feature_shap_ratio = shap_sum / shap_sum.sum()
#转dataframe
feadf = feature_shap_ratio.to_frame(name='percent').rename_axis('index')
#小到大排序
feadf.sort_values(by='percent', ascending=True,inplace=True)
feadf = feadf.reset_index()

#14个特征分组，区分颜色和气泡大小
groups = [0,1,0,1,2,3,3,3,3,0,3,0,0,0]#纵轴从下网上分组
colors = ['#BC3C29FF',
          '#0072B5FF',          
          '#E18727FF',
          '#20854EFF']
# 绘制气泡图
plt.figure(figsize=(3.6,5))
plt.scatter(feadf['percent'], range(len(feadf)), s=feadf['percent']*1000, c=[colors[group] for group in groups])
plt.yticks(range(len(feadf)), list(feadf['index']),fontsize=7)
plt.xticks(fontsize=7)
plt.xlabel("Contribution",fontsize=7)
plt.ylabel("Feature",fontsize=7)
#plt.title("Feature Contribution")
plt.legend(loc='upper right')
# 在气泡旁边标注特征贡献百分比
for i in range(len(feadf)):
    plt.text(feadf['percent'][i] + 0.005, i, f"{feadf['percent'][i]*100:.2f}%", va='center',fontsize=7)
plt.subplots_adjust(left=0.3,right=0.9)
plt.savefig("feature_qipao2.pdf")
plt.show()


# 将14个特征分成四组并图
groups = [0,2,0,2,3,1,1,1,1,0,1,0,0,0]#纵轴从下网上分组
colors = ['#BC3C29FF',
          '#20854EFF',          
          '#0072B5FF',
          '#E18727FF']
# 计算每个组的贡献之和
group_contributions = []
for group in range(len(colors)):
    contribution_sum = sum(feadf['percent'][i] for i in range(len(feadf)) if groups[i] == group)
    group_contributions.append(contribution_sum)
    
# 设置起始角度和偏移量
start_angle = 30
offset = 10
# 绘制旭日图
plt.figure(figsize=(4,4),dpi=100)
plt.rcParams.update({'font.size': 7})#设置全局字号
fig, ax = plt.subplots()
_, sunburst = ax.pie(group_contributions, colors=colors, startangle=start_angle, radius=1, wedgeprops=dict(width=0.4, edgecolor='w'))
# 添加图例
#ax.legend(sunburst, [f"Group {i+1}" for i in range(len(colors))], loc='center right', handlelength=1.5, handletextpad=0.8)
# 在每个组的中心位置添加文本标签，显示该组的贡献之和
for group, contribution in enumerate(group_contributions):
    angle = start_angle + group * 90 + offset
    x = 1.05 * np.cos(np.deg2rad(angle))
    y = 1.05 * np.sin(np.deg2rad(angle))
    ax.text(x, y, f"{contribution*100:.2f}%", ha='center', va='center')
    
# 设置标题及调整布局
#ax.set_title("Feature Contributions by Group")
plt.axis('equal')
plt.tight_layout()
plt.savefig("feature_xuritu.pdf")
plt.show()


# 绘制饼图
fig, ax = plt.subplots()
wedges, texts, autotexts = ax.pie(group_contributions, labels=[f"Group {i+1}" for i in range(len(colors))], colors=colors, autopct='%1.1f%%', startangle=270, wedgeprops={'linewidth': 2, 'edgecolor': 'white'})
ax.set_title("Feature Contributions by Group")

# 添加图例
ax.legend(wedges, [f"Group {i+1}" for i in range(len(colors))], loc='center right')

#ax.pie(group_contributions, labels=[f"Group {i+1}" for i in range(len(colors))], colors=colors, autopct='%1.1f%%', startangle=270)
#ax.set_title("Feature Contributions by Group")

plt.show()


                         

####多个模型的ROC曲线####
"""
import pandas as pd
import numpy as np
from sklearn.model_selection import GridSearchCV, KFold, cross_val_score
from imblearn.over_sampling import SMOTE
from sklearn.metrics import roc_auc_score,roc_curve,auc
from sklearn import metrics
from sklearn.model_selection import train_test_split

data = pd.read_csv("9.2511data.csv")
X = data.iloc[:,:-1]#除了最后一列
Y = data.iloc[:,-1:]#最后一列
#过采样1685+1685
x_resampled, y_resampled = SMOTE(random_state=2022).fit_resample(X,Y)
print(y_resampled['Y'].agg(['value_counts']).T)
#分割训练集和测试集2359 1011-840
train_x, test_x, train_y, test_y = train_test_split(x_resampled, y_resampled, test_size = 0.3,random_state=123)

train_x['Y']=train_y['Y']#1757
test_x['Y']=test_y['Y']#754
lt = list(train_x.columns)
#查看测试集中数据泄漏情况155条
df = test_x.merge(train_x, how='left', indicator=True, left_on=lt, right_on=lt)#1027
test = df[df['_merge']=='left_only']#599
#数据泄漏问题解决后重建训练集和测试集
train_x = train_x.iloc[:,:-1]
test_x = test.iloc[:,:-2]
test_y = test['Y']

#clf1 = DecisionTreeClassifier()
clf2 = RandomForestClassifier(n_estimators=300, max_depth=5, min_samples_split=2,min_samples_leaf=1,max_features=19,oob_score=True, random_state=10)  
clf3 = XGBClassifier(
     learning_rate= 0.1, 
     n_estimators= 230, 
     max_depth= 8,
     min_child_weight= 4,
     seed=0,
     subsample= 0.6, 
     colsample_bytree= 0.7, 
     gamma= 0.3,
     reg_alpha= 0, 
     reg_lambda= 1)

clf4 = LogisticRegression()  
clf5 = SVC(probability=True)  
clf6 = KNN()

clf3.fit(train_x, train_y)

# 在训练集和测试集上分布利用训练好的模型进行预测
train_predict = clf3.predict(train_x)
test_predict = clf3.predict(test_x)

y_pro = clf3.predict_proba(test_x)
predictions = [round(value) for value in test_predict]

print('测试集accuracy:',accuracy_score(test_y,test_predict))
print('测试集precision:',metrics.precision_score(test_y,test_predict))
print('测试集recall:',metrics.recall_score(test_y,test_predict))
print('训练集f1:',metrics.f1_score(train_y,train_predict))
print('测试集roc_auc:',metrics.roc_auc_score(test_y,test_predict))


#函数编写
def multi_models_roc(names, sampling_methods, colors, test_x, test_y, save=True, dpin=100):
        """
        将多个机器模型的roc图输出到一张图上       
        Args:
            names: list, 多个模型的名称
            sampling_methods: list, 多个模型的实例化对象
            save: 选择是否将结果保存（默认为png格式）
            dpin控制图片的信息量（其实可以理解为清晰度           
        Returns:
            返回图片对象plt
        """
        plt.figure(figsize=(4,4), dpi=100)
        # figsize控制图片大小
        for (name, method, colorname) in zip(names, sampling_methods, colors):
            method = method.fit(train_x,train_y)
            y_test_preds = method.predict(test_x)
            y_test_predprob = method.predict_proba(test_x)[:,1]
            fpr, tpr, thresholds = roc_curve(test_y, y_test_predprob, pos_label=1)
            plt.plot(fpr, tpr, lw=2, label='{} (AUC={:.3f})'.format(name, auc(fpr, tpr)),color = colorname)
            plt.plot([0, 1], [0, 1], '--', lw=2, color = 'grey')
            plt.axis('square')
            plt.xlim([0, 1])
            plt.tick_params(labelsize=7)
            plt.ylim([0, 1])
            plt.xlabel('False Positive Rate',fontsize=7,labelpad=7)
            plt.ylabel('True Positive Rate',fontsize=7,labelpad=7)
            plt.legend(loc='lower right',fontsize=7)
            '''
        if save:
            plt.savefig('multi_models_roc.pdf')
            '''
        return plt
# In[]
names = ['XGBoost',
         'Random Forest',
         'KNN', 
         'LogisticRegression',
         'SVM'               
        ]

sampling_methods = [clf3,
                    clf2,
                    clf6,
                    clf4,
                    clf5
                    ]

colors = ['#BC3C29FF',
          '#0072B5FF',          
          '#E18727FF',
          '#20854EFF',
          '#7876B1FF']
#新英格兰配色

#ROC curves
train_roc_graph = multi_models_roc(names, sampling_methods, colors, train_x, train_y,save=True)
train_roc_graph.savefig('ROC_train.pdf')
"""

names = ['XGBoost']
sampling_methods = [clf3]

colors = ['#BC3C29FF']
#新英格兰配色

#ROC curves
train_roc_graph = multi_models_roc(names, sampling_methods, colors, test_x, test_y,save=True)
train_roc_graph.savefig('ROC_test1.pdf')


'''
#模型测试集效果及ROC曲线
import numpy as np
import matplotlib.pyplot as plt
from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from xgboost import XGBClassifier
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, roc_curve, auc

# 加载数据集
data = pd.read_csv("9.2511data.csv")
X = data.iloc[:,:-1]#除了最后一列
Y = data.iloc[:,-1:]#最后一列
#过采样1685+1685
x_resampled, y_resampled = SMOTE(random_state=2022).fit_resample(X,Y)
print(y_resampled['Y'].agg(['value_counts']).T)
#分割训练集和测试集2359 1011-840
train_x, test_x, train_y, test_y = train_test_split(x_resampled, y_resampled, test_size = 0.3,random_state=123)

train_x['Y']=train_y['Y']#1757
test_x['Y']=test_y['Y']#754
lt = list(train_x.columns)
#查看测试集中数据泄漏情况155条
df = test_x.merge(train_x, how='left', indicator=True, left_on=lt, right_on=lt)#1027
test = df[df['_merge']=='left_only']#599
#数据泄漏问题解决后重建训练集和测试集
train_x = train_x.iloc[:,:-1]
test_x = test.iloc[:,:-2]
test_y = test['Y']

# 定义分类器
classifiers = {
    'XGBoost': XGBClassifier(
         learning_rate= 0.1, 
         n_estimators= 230, 
         max_depth= 8,
         min_child_weight= 4,
         seed=0,
         subsample= 0.6, 
         colsample_bytree= 0.7, 
         gamma= 0.3,
         reg_alpha= 0, 
         reg_lambda= 1),
    'Random Forest': RandomForestClassifier(n_estimators=300, max_depth=5, min_samples_split=2,min_samples_leaf=1,max_features=19,oob_score=True, random_state=10),
    'KNN': KNeighborsClassifier(),
    'Logistic Regression': LogisticRegression(),
    'SVM': SVC(probability=True)
}

# 初始化性能指标列表
accuracies = []
precisions = []
recalls = []
f1_scores = []
auc_scores = []

# 训练和评估每个分类器
for name, clf in classifiers.items():
    # 训练模型
    clf.fit(train_x, train_y)

    # 预测测试集
    y_pred = clf.predict(test_x)

    # 计算性能指标
    accuracy = accuracy_score(test_y, y_pred)
    precision = precision_score(test_y, y_pred, average='macro')
    recall = recall_score(test_y, y_pred, average='macro')
    f1 = f1_score(test_y, y_pred, average='macro')

    # 计算AUC值
    y_prob = clf.predict_proba(test_x)
    fpr, tpr, _ = roc_curve(test_y, y_prob[:, 1])
    roc_auc = auc(fpr, tpr)

    # 将性能指标添加到列表中
    accuracies.append(accuracy)
    precisions.append(precision)
    recalls.append(recall)
    f1_scores.append(f1)
    auc_scores.append(roc_auc)

    # 绘制ROC曲线
    plt.plot(fpr, tpr, label=name + ' (AUC = %0.2f)' % roc_auc)

# 绘制ROC曲线图
plt.plot([0, 1], [0, 1], 'k--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic')
plt.legend(loc="lower right")
plt.show()

# 打印性能指标
print('Accuracy:', accuracies)
print('Precision:', precisions)
print('Recall:', recalls)
print('F1 Score:', f1_scores)
print('AUC Score:', auc_scores)

'''


####预测5000条序列两两之间抗原关系####
#cdhit去掉序列完全相同的剩下4769条，4769条序列用于预测，其余的与序列相同的类别进行匹配
#cd-hit -i 5.9548filled_HA1.fasta -o 4769HA1.fasta -c 1 -n 5 -M 16000 -d 0
data = pd.read_csv("6.9548HA1NO2.csv")
data = data.drop_duplicates('Sequence')#4769

#103计算特征值=六个抗原表位+五个氨基酸特征改变量+一个受体结合位点+一个糖基化位点，预测菌株分类  
#nontype的问题是因为菌株名和netNO文件里菌株名没完全匹配
import os
os.chdir(r'/media/dulab/file/H9N2update/final')
import pandas as pd
import numpy as np
from numpy import *
import math 
import time
import heapq
import pickle
import xgboost
from pandas.core.frame import DataFrame
import itertools
from multiprocessing import Process
from multiprocessing import Pool
model = pickle.load(open("XGB822.pickle.dat", "rb"))

data = pd.read_csv("6.9548HA1NO2.csv")
data = data.drop_duplicates('Sequence')#4769

seq_dict = data.set_index('ID').to_dict()['Sequence']
seq_names = data.ID.to_list()

Name = list(data['ID'])
LocN = [eval(x) for x in data['locN1']]
LocO = [eval(y) for y in data['locO1']]

#计算六个抗原表位scannet 96 6类
H9_A = [40, 45, 46, 48, 254, 267, 268, 271, 275, 276, 295]
H9_B = [176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 204, 207]
H9_C = [120, 123, 124, 125, 126, 127, 128, 143, 145, 146, 147, 148, 149, 150, 151, 152, 153, 155, 236]
H9_D = [115, 159, 161, 162, 163, 164, 165, 197, 198, 228, 229, 230, 231, 232, 234, 251]
H9_E = [54, 65, 66, 68, 69, 70, 72, 84, 86, 87, 89, 97, 139]
H9_F = [90, 91, 92, 129, 130, 131, 132, 133, 134, 135, 136, 137, 210, 211, 212, 213, 214, 215, 216, 217]

def get_epitopeA(seq, epi_list=H9_A):
    '''给出一条序列，返回该序列A抗原区域的位点'''
    
    return [seq[base_idx] for base_idx in epi_list]


def get_epitopeB(seq, epi_list=H9_B):
    '''给出一条序列，返回该序列B抗原区域的位点'''
    
    return [seq[base_idx] for base_idx in epi_list]


def get_epitopeC(seq, epi_list=H9_C):
    '''给出一条序列，返回该序列C抗原区域的位点'''
    
    return [seq[base_idx] for base_idx in epi_list]

def get_epitopeD(seq, epi_list=H9_D):
    '''给出一条序列，返回该序列D抗原区域的位点'''
    
    return [seq[base_idx] for base_idx in epi_list]

def get_epitopeE(seq, epi_list=H9_E):
    '''给出一条序列，返回该序列E抗原区域的位点'''
    
    return [seq[base_idx] for base_idx in epi_list]


def get_epitopeF(seq, epi_list=H9_F):
    '''给出一条序列，返回该序列E抗原区域的位点'''
    
    return [seq[base_idx] for base_idx in epi_list]


def cal_epitope(base_set1, base_set2):
    '''给出两个抗原区域位点的列表，计算汉明距离'''

    return sum([a != b for a, b in zip(base_set1, base_set2)])


#计算五个氨基酸特征值
df4=pd.read_csv('character_index4.csv',header = None)
df4 = df4.T#转置
df4.columns= df4.iloc[0]#设置列名
df4 = df4.reindex(df4.index.drop(0))#去掉第一行
df4=df4[['Character_name','Hydrophobicity index (Fasman, 1989)', 'Volume', 'Isoelectric point', 
'Polarizability parameter', 'Average accessible surface area']]
#df4.columns=[['Character_name','Hydrophobicity', 'Volume', 'Charge',
       #'Polarity', 'Accessible surface area']]
dict1 = df4.set_index('Character_name').to_dict()#5个特征

def cal_aa(seq1, seq2):
    '''给出两条序列，5个字典，返回5个氨基酸特征平均改变量'''
    
    n=[]
    for value in dict1.values():
        f=[]
        for base_i, base_j in zip(seq1, seq2):
            f.append(abs(float(value[base_i])-float(value[base_j]))) 
        aa=sum(heapq.nlargest(3,f))/3
        n.append(aa)
    return n


#计算受体结合位点

df1 = pd.read_csv("1jsd_CA.csv")#LOC 1-D 序列 0-D
df1 = df1[df1['E']=='A']#只要HA1 317
df1['loc'] = df1['LOC']-1 #匹配序列中残基编号和坐标编号
df3 = df1[['loc','X','Y','Z']]
df3 = df3.set_index('loc')#loc为索引对于序列残基编号
rbs = [128,129,130,131,132,91,143,145,173,180,184,185,214,215,216,217,218,219]#受体结合位点的H9成熟编号
rbs_list = [x-1 for x in rbs]#在序列列表中编号再减一，已对比无误 GTSKA YWTNVLY NGLMGR

def hamm(s1, s2):
    '''计算汉明距离'''
    
    return sum([a != b for a, b in zip(s1, s2)])


def distance(a,b):
    '''计算欧式距离'''
    dist=math.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2+(a[2]-b[2])**2)
    
    return dist


def cal_rbs(list1,list2):
    '''输入菌株对的序列列表，返回最短欧式距离'''
    
    if hamm(list1,list2)==0:
        e = 0
    else:
        d=[]
        for r in range(0,317):
            if list1[r]!=list2[r]:
                d1=[]
                for o in rbs_list:
                    d1.append(distance(df3.loc[r],df3.loc[o]))#1个改变位点到18个受体结合位点的18个距离
                d.append(min(d1))#所有变化位点到受体结合位点的距离最小值的集合
        if len(d)<=3:
            e = mean(d)
        else:
            e = mean(heapq.nsmallest(3,d))#前三个最短欧氏距离取平均
        
    return [e]


#计算糖基化位点
def get_netno(name):
    '''输入菌株名，返回其糖基化位点列表'''
    
    for i in range(len(data['ID'])):
        if name==Name[i]:
            return LocN[i],LocO[i]


def cal_netno(name1,name2):
    '''给出一对菌株名，返回NO糖基化位点差异数量'''
    
    numN = len(list(set(get_netno(name1)[0])-set(get_netno(name2)[0]))+list(set(get_netno(name2)[0])-set(get_netno(name1)[0])))
    numO = len(list(set(get_netno(name1)[1])-set(get_netno(name2)[1]))+list(set(get_netno(name2)[1])-set(get_netno(name1)[1])))
    num = numN + numO
    return [numN,numO]


def circulation(pair):
    '''对每一个菌株对都进行以上函数操作，并汇总输出'''
    
    name1, name2 = pair
    Seq1 = seq_dict[name1]
    Seq2 = seq_dict[name2]
    
    H9_1 = [cal_epitope(get_epitopeA(Seq1),get_epitopeA(Seq2))]
    H9_2 = [cal_epitope(get_epitopeB(Seq1),get_epitopeB(Seq2))]
    H9_3 = [cal_epitope(get_epitopeC(Seq1),get_epitopeC(Seq2))]
    H9_4 = [cal_epitope(get_epitopeD(Seq1),get_epitopeD(Seq2))]
    H9_5 = [cal_epitope(get_epitopeE(Seq1),get_epitopeE(Seq2))]
    H9_6 = [cal_epitope(get_epitopeF(Seq1),get_epitopeF(Seq2))]

    line = H9_1 + H9_2 + H9_3 + H9_4 + H9_5 + H9_6 + cal_aa(Seq1,Seq2) + cal_netno(name1,name2) + cal_rbs(Seq1,Seq2)
    
    return line


def predict(featuredata):
    '''用模型对所有菌株特征值进行预测,输出分类结果和概率'''
    
    y_pro = model.predict(featuredata)
    y_pred= model.predict_proba(featuredata)#概率
    y_proab=[]#分类情况
    for i in y_pred[:,0]:#第一列数据预测为不相似0的概率 菌株对相似的概率>0.5相似
        if i > 0.5:
            y_proab.append(0)
        else:
            y_proab.append(1)
    
    return list(y_pred[:,0]),list(y_pred[:,1]), y_proab


iters = itertools.combinations(seq_names, 2)
tasks = list(iters)[0:100]

with Pool(10) as pool:

    start_time = time.time()
    res = pool.map(circulation, tasks)
    res2 = pd.DataFrame(res)
    res2.columns=['H9-A','H9-B','H9-C','H9-D','H9-E','H9-F',
                  'Hydrophobicity', 'Volume', 'Charge','Polarity', 'Accessible surface area',
                  'N-Glycosylation','O-Glycosylation','Receptor binding']
    res3 = predict(res2)
    res4 = pd.DataFrame(res3)
    res5 = res4.T
    res6 =pd.concat((pd.DataFrame(tasks),res2,res5),axis=1)
    res6.columns=['name1', 'name2', 'H9-A','H9-B','H9-C','H9-D','H9-E','H9-F',
                  'Hydrophobicity', 'Volume', 'Charge','Polarity', 'Accessible surface area',
                  'N-Glycosylation','O-Glycosylation','Receptor binding','y_proba0','y_proba1','y_pred']
    pool.close()
    pool.join()
    
    score = list()
    for i in range(len(res6['y_pred'])):
        score.append(np.log((res6['y_proba1'][i])/(res6['y_proba0'][i])))
    res6['antigenicityscore'] = score
    res7 = res6[['name1','name2','antigenicityscore']]
    res7.columns=['Source','Target','Weight']
    res6.to_csv('4769predict1.csv', index=None)
    res7.to_csv('4769weight1.csv', index=False)
    print(f"共耗时{int(time.time() - start_time)}s")


#提取抗原相似对 txt_mcl csv_cyto  mcl结果有9520菌株分了类，9条没有归类
import os
os.chdir("/media/dulab/file/H9N2update/final/10297predict")
import pandas as pd
import numpy as np
df=pd.read_csv("4769predict.csv")
print(df.head(3)) #看列名 #1是相似 y_proba1是相似概率 抗原相似性得分大于0 相似
df2=pd.read_csv("10297weight.csv")
#只要相似对
df3=df2[df2['Weight']>0]
df3.to_csv('10297similar.txt',header=None, sep='\t',index=None)

####MCL聚类####
#final 230 XGB  4769
#I=1.6 只有4761条分类 / I=2.2 25类 4761
"""
#统计每簇多少菌株，手动在节点属性文件里标记
#统计txt文本行数--cluster有几类m
import os
os.chdir("/media/dulab/file/H9N2update/final/10297predict/")

fpath = 'out.lit2020.I17'
with open(fpath,'r') as f:
    countlines = 0  #统计全文行数
    for contents in  f :
        contents = contents.strip()     #去掉前后空格
        contents = contents.split("\n") #用split根据换行分割 返回list
        countlines+=1
print(countlines)

#统计每类包含多少菌株c
import os, sys, re
n=0
count=[]
with open("out.lit2020.I17", 'r') as f:
    for line in f:
        line = line.strip()
        line = line.split('\t')
        c = len(line)
        #print(c)
        count.append(c)
        n+=c    
print(n)
print(count)

#MCL输出节点类别，整理节点属性文件
#txt按规则转excel，tab为分割，每一个一行
import datetime
import time
import os
import sys
import xlwt #需要的模块
import os
#os.chdir("/media/dulab/file/H9N2update/final/4769predict/cluster25")

def txt2xls(filename,xlsname):  #文本转换成xls的函数，filename 表示一个要被转换的txt文本，xlsname 表示转换后的文件名
    #print 'converting xls ... '
    f = open(filename)   #打开txt文本进行读取
    x = 0                #在excel开始写的位置（y）
    y = 0                #在excel开始写的位置（x）
    xls=xlwt.Workbook()
    sheet = xls.add_sheet('sheet1',cell_overwrite_ok=True) #生成excel的方法，声明excel
    while True:  #循环，读取文本里面的所有内容
        line = f.readline() #一行一行读取
        if not line:  #如果没有内容，则退出循环
            break
        for i in line.split('\t'):#读取出相应的内容写到x
            item=i.strip()
            sheet.write(x,y,item)
            x += 1 #另起一行
        y = 0  #初始成第一列
    f.close()
    xls.save(xlsname+'.xls') #保存

if __name__ == "__main__":
    filename = r'out.lit2020.I17'
    xlsname  = r'10289节点_24'
    txt2xls(filename,xlsname)
#xls手动转csv

#自动生成节点属性列
import pandas as pd
data= pd.read_csv('10289节点_24.csv')

a = [1 for _ in range(count[0])]
b = [2 for _ in range(count[1])]
c = [3 for _ in range(count[2])]
d = [4 for _ in range(count[3])]
e = [5 for _ in range(count[4])]
f = [6 for _ in range(count[5])]
g = [7 for _ in range(count[6])]
h = [8 for _ in range(count[7])]
i = [9 for _ in range(count[8])]
j = [10 for _ in range(count[9])]
k = [11 for _ in range(count[10])]
l = [12 for _ in range(count[11])]
m = [13 for _ in range(count[12])]
n = [14 for _ in range(count[13])]
o = [15 for _ in range(count[14])]
p = [16 for _ in range(count[15])]
q = [17 for _ in range(count[16])]
r = [18 for _ in range(count[17])]
s = [19 for _ in range(count[18])]
t = [20 for _ in range(count[19])]
u = [21 for _ in range(count[20])]
v = [22 for _ in range(count[21])]
w = [23 for _ in range(count[22])]
x = [24 for _ in range(count[23])]

all=a+b+c+d+e+f+g+h+i+j+k+l+m+n+o+p+q+r+s+t+u+v+w+x
data['属性']=all
data.to_csv('10289_24cnode.csv',index=None)

"""


#####抗原模式图
#通过进化树看出BJ94：1-2，1-3,2-5,2-6;G1：8-4； Y439：9-13；
####类间转换关键位点识别（计算熵和信息增益）####
"""
#10289的注释信息
import os
os.chdir("/media/dulab/file/H9N2update/final/9529HAalltree")
df1=pd.read_csv('9529HA.csv')
#修改列名
lst = list(range(1,318))
lst2 = ['ID','Sequence']+lst
df1.columns= lst2
#df1.to_csv('9529HAaanote.txt',sep='\t',index=None)

df2=pd.read_csv('/media/dulab/file/H9N2update/final/10297predict/10297_24loctime.csv')
df2= df2[['ID','Sequence','属性']]
df2.columns = ['ID','Sequence','cluster']
df=pd.merge(df1,df2,how='left',#仅使用左frame中的键，类似于SQL左外部联接；保留关键顺序
              left_on=['ID','Sequence'],right_on=['ID','Sequence'])
#未分类的八条序列用unknown填充
#df['cluster'] = df['cluster'].fillna('unknown')
#去掉未分类的序列
df=df.dropna(axis=0,how='any')
df=df.drop('Sequence',axis=1)
#df.dtypes
#列格式转换
df['cluster'] = df['cluster'].astype(int)
df['cluster'] = df['cluster'].astype(object)
#确定每类类名之后重新注释
df['Cluster'] = df['cluster'].map({1:'BJ94-1',2:'BJ94-2',3:'BJ94-3',4:'G1-1',5:'BJ94-4',
                                6:'BJ94-5',7:'WC66',8:'G1-2',9:'Y439-1',13:'Y439-2',10:'unknown',
                                11:'unknown',12:'unknown',14:'unknown',19:'unknown',
                                15:'unknown',16:'unknown',17:'unknown',18:'unknown',
                                20:'unknown',21:'unknown',22:'unknown',23:'unknown',
                                24:'unknown',
                               },na_action='ignore')
df.to_csv('10289_24aanote.txt',sep='\t',index=None)
df.to_csv('10289_24aanote.csv',sep=',',index=None)

#计算每个位点的熵
#计算每个位点的信息增益 熵 条件信息熵
#定义计算熵的函数
def ent(data):
    prob1 = pd.value_counts(data) / len(data)
    return sum(np.log2(prob1) * prob1 * (-1))
 
    
#定义计算信息增益的函数
def gain(data,str1,str2):
    e1 = data.groupby(str1).apply(lambda x:ent(x[str2]))
    p1 = pd.value_counts(data[str1]) / len(data[str1])
    e2 = sum(e1 * p1)
    return ent(data[str2]) - e2

#return ent(data[str2]), e2, (ent(data[str2]) - e2)

import os
os.chdir("/media/dulab/file/H9N2update/final/10297predict/keysite5-11")
import numpy as np
df=pd.read_csv('10289_24aanote.csv')
df1=df
df2=df1[df1['cluster']==1]
df3=df1[df1['cluster']==3]
data=pd.concat([df2,df3])    
#换数据检查这一行a=list(data.columns[2:319])[i] 1-317才对
#a=list(data.columns[1:318])[0]

#计算每个位点的熵
Entropy=[]
for i in range(317):
    a=list(data.columns[1:318])[i]      #1-317
    Entropy.append(ent(data[a]))
#计算信息增益
IG=[]
for i in range(317):
    a=list(data.columns[1:318])[i]
    IG.append(gain(data,a,'cluster'))
    
#print(gain(data,'150','cluster'))

list1 = list(range(1,318))
from pandas.core.frame import DataFrame
c={"site":list1,
  "Entropy":Entropy,
  "IG":IG}
data2=DataFrame(c)#将字典转换成为数据框
#01Max-Min标准化
#建立MinMaxScaler对象
from sklearn import preprocessing
minmax = preprocessing.MinMaxScaler()
# 'Entropy','IG'标准化处理
data_minmax = minmax.fit_transform(data2[['Entropy','IG']])
data3=pd.DataFrame(data_minmax)

c={"site":data2['site'],
  "Entropy":data3[0],
  "IG":data3[1]
  }
data4=pd.DataFrame(c)#将字典转换成为数据框
data4.eval('sum = Entropy + IG' , inplace=True)
data4 = data4.sort_values(by= "IG" , ascending=False)
data4.to_csv('3-2.csv',index=None)

"""

####中国菌株每年疫苗株推荐####

"""
####批量计算不同年份推荐疫苗株及覆盖率
import os
os.chdir("/media/dulab/file/H9N2update/final/10297predict")
import networkx as nx
import pandas as pd

# 定义函数：获取邻居节点
def get_neighbors(graph, node):
    """
    获取 node 节点的邻居节点
    """
    return set(graph.neighbors(node))

# 定义函数：获取两个节点的邻居节点的并集
def get_union_of_neighbors(graph, node1, node2):
    """
    获取两个节点的邻居节点的并集
    """
    neighbors1 = get_neighbors(graph, node1)
    neighbors2 = get_neighbors(graph, node2)
    return neighbors1.union(neighbors2)

# 定义函数：获取三个节点的邻居节点的并集
def get_union_of_neighbors_three(graph, node1, node2, node3):
    """
    获取三个节点的邻居节点的并集
    """
    neighbors1 = get_neighbors(graph, node1)
    neighbors2 = get_neighbors(graph, node2)
    neighbors3 = get_neighbors(graph, node3)
    return neighbors1.union(neighbors2, neighbors3)

# 加载数据
data = pd.read_csv('10297weight.csv')

# 创建带权重的无向图
G = nx.Graph()
for index, row in data.iterrows():
    if row['Weight'] > 0:
        G.add_edge(row['Source'], row['Target'], weight=row['Weight'])

# 选择中国菌株
df_all = pd.read_csv('10289_24loctime.csv')
df_China = df_all[df_all['Region3'] == 'China']

# 计算每年推荐的两个疫苗株以及其抗原覆盖率
results = []
#for year in range(2010, 2023):
for year in range(2019, 2022):
    # 根据年份划分组别构建多个子网络
    df_China_year = df_China[df_China['Time'] == year]
    China_nodes_year = df_China_year['ID']

    # 构建子网络
    subgraph_year = G.subgraph(China_nodes_year)

    # 计算中国该年份子网络中所有节点的度值，并匹配上节点类别
    degree_sequence_year = [subgraph_year.degree(n) for n in subgraph_year.nodes()]
    c = {"node": subgraph_year.nodes(), "degree": degree_sequence_year}
    df_degree_year = pd.DataFrame(c).sort_values('degree', axis=0, ascending=False, inplace=False)
    df_cluster_year = pd.merge(df_degree_year, df_all[['ID', '属性']], how='left', left_on=['node'], right_on=['ID']).drop('ID', axis=1)

    # 推荐疫苗株1
    strain1 = df_cluster_year['node'][0]
    strain1_cluster = df_cluster_year['属性'][0]
    
    # 计算该年份推荐的1个疫苗株在当年的抗原覆盖率
    # 针对当年毒株的抗原覆盖率
    union_of_neighbors1 = get_neighbors(subgraph_year, strain1)
    antigen_coverage_year1 = len(union_of_neighbors1) / len(df_China_year)

    # 去掉该疫苗株及其邻居节点，剩下的节点子网络再次选出度值最高的节点
    # 疫苗株及其邻居节点列表
    remove_nodes = [strain1] + list(subgraph_year.neighbors(strain1))

    # 获取剩余节点的补图
    remaining_nodes = set(China_nodes_year) - set(remove_nodes)
    complement_graph = subgraph_year.subgraph(remaining_nodes)

    # 计算子网络中所有节点的度值，并匹配上节点类别
    degree_sequence_complement = [complement_graph.degree(n) for n in complement_graph.nodes()]
    c = {"node": complement_graph.nodes(), "degree": degree_sequence_complement}
    df_degree_complement = pd.DataFrame(c).sort_values('degree', axis=0, ascending=False, inplace=False)
    df_cluster_complement = pd.merge(df_degree_complement, df_all[['ID', '属性']], how='left', left_on=['node'], right_on=['ID']).drop('ID', axis=1)

    # 推荐疫苗株2
    strain2 = df_cluster_complement['node'][0]
    strain2_cluster = df_cluster_complement['属性'][0]
    
    # 计算该年份推荐的两个疫苗株在当年和次年的抗原覆盖率
    # 针对当年毒株的抗原覆盖率
    union_of_neighbors2 = get_union_of_neighbors(subgraph_year, strain1, strain2)
    antigen_coverage_year2 = len(union_of_neighbors2) / len(df_China_year)
    
    # 针对下一年毒株的覆盖率
    #df_China_next_year = df_China[df_China['Time'] == year + 1]
    #China_nodes_next_year = list(df_China_next_year['ID']) + [strain1, strain2]
    #subgraph_next_year = G.subgraph(China_nodes_next_year)
    #union_of_neighbors = get_union_of_neighbors(subgraph_next_year, strain1, strain2)
    #antigen_coverage_next_year = len(union_of_neighbors) / len(df_China_next_year)
    
    # 去掉疫苗株12及其邻居节点，剩下的节点子网络再次选出度值最高的节点
    # 疫苗株及其邻居节点列表
    remove_nodes2 = [strain2] + list(complement_graph.neighbors(strain2))

    # 获取剩余节点的补图
    remaining_nodes2 = set(remaining_nodes) - set(remove_nodes2)
    complement_graph2 = complement_graph.subgraph(remaining_nodes2)

    # 计算子网络中所有节点的度值，并匹配上节点类别
    degree_sequence_complement2 = [complement_graph2.degree(n) for n in complement_graph2.nodes()]
    c = {"node": complement_graph2.nodes(), "degree": degree_sequence_complement2}
    df_degree_complement2 = pd.DataFrame(c).sort_values('degree', axis=0, ascending=False, inplace=False)
    df_cluster_complement2 = pd.merge(df_degree_complement2, df_all[['ID', '属性']], how='left', left_on=['node'], right_on=['ID']).drop('ID', axis=1)

    # 推荐疫苗株3
    strain3 = df_cluster_complement2['node'][0]
    strain3_cluster = df_cluster_complement2['属性'][0]
    
    # 计算该年份推荐的三个疫苗株在当年的抗原覆盖率
    # 针对当年毒株的抗原覆盖率
    union_of_neighbors3 = get_union_of_neighbors_three(subgraph_year, strain1, strain2,strain3)
    antigen_coverage_year3 = len(union_of_neighbors3) / len(df_China_year)
    
    # 将结果存储到列表中
    results.append([year, strain1, strain1_cluster, antigen_coverage_year1, strain2, strain2_cluster, antigen_coverage_year2, strain3, strain3_cluster, antigen_coverage_year3])

# 将结果转化为 DataFrame 并保存到 csv 文件中
df_results = pd.DataFrame(results, columns=['Year', 'Strain1', 'Strain1_cluster','antigen_coverage_year1', 'Strain2', 'Strain2_cluster', 'antigen_coverage_year2','Strain3', 'Strain3_cluster','antigen_coverage_year3']) #'Antigen Coverage (Next Year)'
df_results.to_csv('2019_3vaccine_recommendation_results.csv', index=False)

"""

'''
#批量计算多个疫苗株多个年份的抗原覆盖率
import os
os.chdir("/media/dulab/file/H9N2update/final/10297predict")
import networkx as nx
import pandas as pd
# 读取数据
data = pd.read_csv('10297weight.csv')

# 创建带权重的无向图
G = nx.Graph()
for index, row in data.iterrows():
    if row['Weight'] > 0:
        G.add_edge(row['Source'], row['Target'], weight=row['Weight'])

# 选择中国菌株
df_all = pd.read_csv('10289_24loctime.csv')
df_China = df_all[df_all['Region3'] == 'China']

def calculate_antigen_coverage(vaccine_strain, year):
    
    # 提取指定年份毒株
    df_China_year = df_China[df_China['Time'] == year]
    China_nodes_year = list(df_China_year['ID'])
    China_nodes_year.append(vaccine_strain)

    # 构建子网络
    subgraph = G.subgraph(China_nodes_year)

    # 计算中国该年份疫苗株在子网络中的度值
    degree1 = subgraph.degree(vaccine_strain)

    # 计算该年份疫苗株在当年的抗原覆盖率
    antigen_coverage = degree1 / len(China_nodes_year)
    
    return antigen_coverage

# 创建用于存储结果的 DataFrame
results = pd.DataFrame(columns=['vaccine_strain', 'year', 'antigen_coverage'])
# 计算多个疫苗株在2010-2022年每年的抗原覆盖率
vaccine_strains = ['a/chicken/shandong/6/1996','a/chicken/shanghai/f/1998','a/chicken/nanjing/02/2001', 'a/chicken/guangdong/ss/1994', 'a/chicken/shandong/s2/2005','a/chicken/jiangsu/wj57/2012']
years = range(2000, 2023)

for vaccine_strain in vaccine_strains:
    for year in years:
        antigen_coverage_year = calculate_antigen_coverage(vaccine_strain, year)
        results = results.append({'vaccine_strain': vaccine_strain, 'year': year, 'antigen_coverage': antigen_coverage_year}, ignore_index=True)

# 将结果保存到 CSV 文件
results.to_csv('vaccine_strain_antigen_coverage2.csv', index=False)
'''

#验证聚类效果，计算有HI试验数据的毒株分组准确率
import os
os.chdir(r'/media/dulab/file/H9N2update/final/10297predict')
import pandas as pd
df2 = pd.read_csv('7.2511HI.csv')
df1 = pd.read_csv('10289_24cnode.csv')
# 计算Y=1（826）时毒株名Virus和毒株名Serum在df1中属性相同的概率0.7832929782082324
df3 = pd.merge(df1, df2[df2['Y'] == 1], left_on='ID', right_on='Virus')
df4 = pd.merge(df3, df1, left_on='Serum', right_on='ID')
p_same_attr = sum(df4['属性_x'] == df4['属性_y']) / len(df4)
print("Y=1时毒株名Virus和毒株名Serum在df1中属性相同的概率：", p_same_attr)

# 计算Y=0时毒株名Virus和毒株名Serum在df1中属性不同的概率0.712166172106825
df5 = pd.merge(df1, df2[df2['Y'] == 0], left_on='ID', right_on='Virus')
df6 = pd.merge(df5, df1, left_on='Serum', right_on='ID')
p_diff_attr = sum(df6['属性_x'] != df6['属性_y']) / len(df6)
print("Y=0时毒株名Virus和毒株名Serum在df1中属性不同的概率：", p_diff_attr)

df2 = pd.read_csv('0.4161HIY.csv')
df1 = pd.read_csv('10289_24cnode.csv')
# 计算Y=1（1125）时毒株名Virus和毒株名Serum在df1中属性相同的概率0.7723840345199569
df3 = pd.merge(df1, df2[df2['Y'] == 1], left_on='ID', right_on='Virus')
df4 = pd.merge(df3, df1, left_on='Serum', right_on='ID')
p_same_attr = sum(df4['属性_x'] == df4['属性_y']) / len(df4)
print("Y=1时毒株名Virus和毒株名Serum在df1中属性相同的概率：", p_same_attr)

# 计算Y=0时毒株名Virus和毒株名Serum在df1中属性不同的概率0.7090142329994729
df5 = pd.merge(df1, df2[df2['Y'] == 0], left_on='ID', right_on='Virus')
df6 = pd.merge(df5, df1, left_on='Serum', right_on='ID')
p_diff_attr = sum(df6['属性_x'] != df6['属性_y']) / len(df6)
print("Y=0时毒株名Virus和毒株名Serum在df1中属性不同的概率：", p_diff_attr)






