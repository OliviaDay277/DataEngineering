import pandas as pd
import numpy as np
from LogisticRegression import LogisticRegression
from GaussianGenerativeModel import GaussianGenerativeModel

fruits=pd.read_csv("fruit.csv")
X = fruits[['width', 'height']].values
X_lr=np.hstack((np.ones(shape=(X.shape[0],1)),X))
Y = (fruits['fruit'] - 1).values

# Logistic Regression parameters
eta = [0.0001,0.001,0.01,0.1]
lambda_parameter = [0.01,0.05,0.1,0.15]

plt.figure(figsize=(20,20))
i=1
for e in eta:
    for l in lambda_parameter:
        lr = LogisticRegression(eta=e, lambda_parameter=l)
        lr.fit(X_lr,Y)
        
        plt.subplot(4,4,i)
        plt.plot(np.arange(len(lr.neglogL_trace)),lr.neglogL_trace)
        plt.title("eta="+str(e)+" "+"Lambda="+str(l))
        i=i+1
        
lr = LogisticRegression(eta=0.0001, lambda_parameter=0.1)
lr.fit(X_lr,Y)
lr.visualize('logistic_regression_result_eta0001.png')

lr = LogisticRegression(eta=0.1, lambda_parameter=0.1)
lr.fit(X_lr,Y)
lr.visualize('0101.png')

nb1 = GaussianGenerativeModel(isSharedCovariance=False)
nb1.fit(X,Y)
nb1.visualize("generative_result_separate_covariances.png")

nb2 = GaussianGenerativeModel(isSharedCovariance=True)
nb2.fit(X,Y)
nb2.visualize("generative_result_shared_covariances.png")

X_test = np.array([[4,11],[8.5,7]])
Y_nb1 = nb1.predict(X_test)
Y_nb2 = nb2.predict(X_test)
Y_lr = lr.predict(np.hstack(np.ones(shape=(2,1)),X_test))

print("Test fruit predictions for Gaussian Model:")
print("width 4 cm and height 11 cm: " + str(Y_nb1[0]))
print("width 8.5 cm and height 7 cm: " + str(Y_nb1[1]))

print("Test fruit predictions for Shared Covariance Gaussian Model:")
print("width 4 cm and height 11 cm: " + str(Y_nb2[0]))
print("width 8.5 cm and height 7 cm: " + str(Y_nb2[1]))

print("Test fruit predictions for Linear Regression:")
print("width 4 cm and height 11 cm: " + str(Y_lr[0]))
print("width 8.5 cm and height 7 cm: " + str(Y_lr[1]))

print("The negative log-likelihood for separate covariance GGM:",nb1.neglogL)
print("The negative log-likelihood for shared covariance GGM:",nb2.neglogL)
