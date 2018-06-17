import numpy as np
from sklearn.feature_extraction import DictVectorizer

def read_data(file_name):
    with open(file_name) as f:
        content = f.readlines()
    content = [x.strip() for x in content]
    data = []
    y = []
    for i in range(len(content)):
        line = content[i].split(' ')
        line_dic = {}
        for j in range(len(line)-1):
            word, freq = line[j].split(':')
            line_dic[word] = int(freq)
        if line[-1][8:] == 'negative':
            y.append(-1)
        else:
            y.append(1)
        data.append(line_dic)
    return data, y


def organise_data(src, tgt):
    XB_neg, YB_neg = read_data('data/{}/negative.review'.format(src))
    XB_pos, YB_pos = read_data('data/{}/positive.review'.format(src))
    XB = XB_pos + XB_neg
    YB = YB_pos + YB_neg
    print(src, np.sum(YB_pos), np.sum(YB_neg))
    XD_neg, YD_neg = read_data('data/{}/negative.review'.format(tgt))
    XD_pos, YD_pos = read_data('data/{}/positive.review'.format(tgt))
    XD = XD_pos + XD_neg
    YD = YD_pos + YD_neg
    print(tgt, np.sum(YD_pos), np.sum(YD_neg))
    return XB, YB, XD, YD


def vectorise_data(XB, XD):
    # Need to jointly transform becaue of the unique fields
    vectorizer = DictVectorizer(sparse=False)
    X_BD = XB + XD
    X_tot = vectorizer.fit_transform(X_BD)
    X_src = X_tot[:len(XB),:]
    X_tgt = X_tot[len(XB):,:]
    features = vectorizer.get_feature_names()
    return X_src, X_tgt, features

def select_high_freq_data(X_src, X_tgt, features, N):
    X_tot = X_src + X_tgt
    feature_count = np.sum(X_tot,0)
    sort_idx = np.argsort(feature_count)
    X_src = X_src[:, sort_idx[-N:]]
    X_tgt = X_tgt[:, sort_idx[-N:]]
    features = [features[i] for i in sort_idx[-N:]]
    return X_src, X_tgt, features




