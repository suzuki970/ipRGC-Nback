#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  5 06:53:17 2024

@author: yutasuzuki
"""



class predByRidge:
    def __init__(self, alpha=1):
        self.w_ = None
        self.alpha = alpha

    # def fit(self, X, y):
    #     X = np.insert(X, 0, 1, axis=1)
    #     I = np.eye(X.shape[1])
    #     self.w_ = np.linalg.inv(X.T @ X + self.alpha*I) @ X.T @ y

    # def predict(self, X):
    #     X = np.insert(X, 0, 1, axis=1)
    #     return X @ self.w_
    
class predByLR:
    def __init__(self, alpha=1):
        self.w_ = None
        self.alpha = alpha


    
class predBy_xgboost:
    def __init__(self, alpha=1):
        self.w_ = None
        self.alpha = alpha

    def fit(self, X, y):
        
        Xtr = np.r_[train_feat_norm,np.ones((1,train_feat_norm.shape[1]))];
        Xte = np.r_[test_feat_norm,np.ones((1,test_feat_norm.shape[1]))];
        
        Ytr = train_pupil_norm.copy()
        
        X1X1 = np.dot(Xtr.T, Xtr)
        pitr=1
        
        # tmp_res=[]
        # for pitr in np.linspace(0.01,0.1,10):
        weight,pred_score = predScore(X1X1,Xtr,Xte,Ytr,pitr)
    
    def predict(self, X):
        
    return pred_score