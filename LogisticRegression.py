import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as c
from scipy.misc import logsumexp

# Please implement the fit and predict methods of this class. You can add additional private methods
# by beginning them with two underscores. It may look like the __dummyPrivateMethod below.
# You can feel free to change any of the class attributes, as long as you do not change any of 
# the given function headers (they must take and return the same arguments), and as long as you
# don't change anything in the .visualize() method. 
class LogisticRegression:
    def __init__(self, eta, lambda_parameter):
        self.eta = eta
        self.lambda_parameter = lambda_parameter
        self.W=np.zeros((3,3),dtype=float)
        self.neglogL_trace=[]
    def softmax(self):
        z_k_mx=np.exp(np.dot(self.X,self.W.T))
        denominator_mx=np.sum(z_k_mx,axis=1)
        sm_mx=np.zeros(z_k_mx.shape)
        for i in range(z_k_mx.shape[1]):
            sm_mx[:,i]=z_k_mx[:,i]/denominator_mx
        self.sm_mx=sm_mx
        return sm_mx
    def onehot_class(self):
        C_mx=np.zeros((len(self.C),3))
        for i in range(len(self.C)):
                C_mx[i,self.C[i]]=1
        self.C_mx=C_mx
        return 
    
    def grad_desc(self):
       
        for i in range(500000): 
            
            nlogl=0
            softmax_mx=self.softmax()
            for j in range(len(self.C)):

                p_j=np.dot(softmax_mx[j,:],self.C_mx[j,:])
                nlogl=nlogl-np.log(p_j)
            self.neglogL_trace.append(nlogl)
            gradient = np.dot((self.softmax() - self.C_mx).T, self.X) + self.lambda_parameter * self.W
            #derivative of the regularized loss function
            self.W = self.W - self.eta * gradient
            if np.absolute(self.eta * gradient).sum() < 0.00001: 
                break
        return None

    # TODO: Implement this method!
    def fit(self, X, C):
        self.X = X
        self.C = C
        self.onehot_class()
        self.softmax()
        self.grad_desc()
        return

    # TODO: Implement this method!
    def predict(self, X_to_predict):
        self.X=X_to_predict
        y_pred=self.softmax()
        c_hat=y_pred.argmax(axis=1)
        return c_hat

    def visualize(self, output_file, width=2, show_charts=False):
        X = self.X

        # Create a grid of points
        x_min, x_max = min(X[:, 1] - width), max(X[:, 1] + width)
        y_min, y_max = min(X[:, 2] - width), max(X[:, 2] + width)
        xx,yy = np.meshgrid(np.arange(x_min, x_max, .05), np.arange(y_min,
            y_max, .05))

        # Flatten the grid so the values match spec for self.predict
        xx_flat = xx.flatten()
        yy_flat = yy.flatten()
        X_topredict = np.hstack((np.ones(shape=(len(xx_flat),1)),np.vstack((xx_flat,yy_flat)).T))

        # Get the class predictions
        Y_hat = self.predict(X_topredict)
        Y_hat = Y_hat.reshape((xx.shape[0], xx.shape[1]))
        
        cMap = c.ListedColormap(['r','b','g'])

        # Visualize them.
        plt.figure()
        plt.pcolormesh(xx,yy,Y_hat, cmap=cMap)
        plt.scatter(X[:, 1], X[:, 2], c=self.C, cmap=cMap,edgecolors='black')
        plt.title("eta="+str(self.eta)+", "+"lambda="+str(self.lambda_parameter))

        plt.savefig(output_file)
        if show_charts:
            plt.show()




