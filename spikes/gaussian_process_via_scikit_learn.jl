# Based on http://scikit-learn.org/stable/modules/gaussian_process.html
using PyCall
@pyimport numpy as np
@pyimport sklearn.gaussian_process as gausspr

f(x) = x .* sin(x)

xtrain = np.atleast_2d([0.9, 1.0, 3.0, 5.0, 7.0, 9.0, 9.1])'
ytrain = f(xtrain)

gp = gausspr.GaussianProcess(theta0=1e-2, thetaL=1e-4, thetaU=1e-1)
gp["fit"](xtrain, ytrain)

xtest = np.atleast_2d(np.linspace(0, 10, 1000))'
ytest = f(xtest)
y_pred, sigma2_pred = gp["predict"](xtest, eval_MSE=true)

mse(yhat, y) = mean(sum((yhat .- y).^2))
rmse(yhat, y) = sqrt(mse(yhat, y))
mape(yhat, y) = 100.0 * mean(abs( (y .- yhat)./y ))

rmse(y_pred[:], ytest[:])

# Around 1.0 it should be pretty close:
y_pred[95:105]
ytest[95:105]
mape(y_pred[95:105], ytest[95:105])
