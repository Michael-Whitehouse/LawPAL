import numpy as np
import scipy.stats as stats

def K_eta_SEIR( beta, rho, gamma):
    
    def K_eta_matrix(x_t):

        matrix      = np.zeros((4, 4))

        matrix[0,0] = np.exp(-beta*x_t[2]/(np.sum(x_t)))
        matrix[0,1] = 1 - np.exp(-beta*x_t[2]/(np.sum(x_t)))

        matrix[1,1+0] = np.exp(-rho)
        matrix[1,1+1] = 1 - np.exp(-rho)

        matrix[2,2+0] = np.exp(-gamma)
        matrix[2,2+1] = 1 - np.exp(-gamma)

        matrix[3,3]   = 1

        return matrix

    return K_eta_matrix

def K_eta_SEEIR( beta, rho_1, rho_2, gamma):
    
    def K_eta_matrix(x_t):

        matrix      = np.zeros((5, 5))

        matrix[0,0] = np.exp(-beta*x_t[3]/(np.sum(x_t)))
        matrix[0,1] = 1 - np.exp(-beta*x_t[3]/(np.sum(x_t)))

        matrix[1,1+0] = np.exp(-rho_1)
        matrix[1,1+1] = 1 - np.exp(-rho_1)

        matrix[2,2+0] = np.exp(-rho_2)
        matrix[2,2+1] = 1 - np.exp(-rho_2)

        matrix[3,3+0] = np.exp(-gamma)
        matrix[3,3+1] = 1 - np.exp(-gamma)

        matrix[4,4]   = 1

        return matrix

    return K_eta_matrix


class Compartmental_model():

    def __init__(self, pi_0, delta, K_eta, alpha, q_pars, G, kappa, n):

        self.n     = n
        self.pi_0  = pi_0
        self.delta = delta
        self.K_eta = K_eta
        self.alpha = alpha
        self.q_pars     = q_pars    
        self.G     = G    
        self.kappa = kappa

    def step_0(self):

        return np.transpose(np.random.multinomial(self.n, self.pi_0.squeeze(), size = 1))
        # return (np.random.poisson(self.n*self.pi_0))

    def step_t(self, x_tm1):

        # Latent process
        # deaths
        # deaths    = np.random.binomial(x_tm1, (1-self.delta))
        barx_tm1  = np.random.binomial(x_tm1, self.delta)

        # transitions
        K_eta_tm1     = self.K_eta(barx_tm1)
        transitions_x = np.array([np.random.multinomial(barx_tm1[i], K_eta_tm1[i,:], 1)[0] for i in range(0, len(barx_tm1))])
        tildex_t      = np.transpose(np.sum(transitions_x, axis = 0, keepdims=True))

        # births
        births = np.random.poisson(self.alpha)
        x_t    = tildex_t + births

        # Observed process
        lower = 0
        upper = 1
        mu = self.q_pars[0]
        sigma = self.q_pars[1]
        q = np.zeros(4)
        q[2] = stats.truncnorm.rvs((lower-mu)/sigma,(upper-mu)/sigma,loc=mu,scale=sigma,size=1)
        ## Next line does random under reporting for each compartment not just infectives...
        bary_t        = np.random.binomial(x_t, q)[:,2]
        #transitions_y = np.array([np.random.multinomial(bary_t, self.G[i,:], 1)[0] for i in range(0, len(bary_t))])
        #tildey_t      = np.transpose(np.sum(transitions_y, axis = 0, keepdims=True))

        # clutter
        haty_t = np.random.poisson(self.kappa)
        y_t    = bary_t

        return x_t, y_t, q[2]

    def run(self, T):

        X = self.step_0()
        Y = np.zeros(X.shape)
        q_sim = np.zeros(1)

        for t in range(0, T):

            x_t, y_t, q = self.step_t(X[:,t:t+1])

            X = np.concatenate((X, x_t), axis = 1)
            Y = np.concatenate((Y, np.expand_dims(y_t,axis =1)), axis = 1)
            q_sim = np.concatenate((q_sim,  np.array([q])), axis = 0)

        return X, Y, q_sim

    def prediction(self, barlambda_tm1):
    
        K_eta_tm1 = self.K_eta(barlambda_tm1*self.delta)

        return np.transpose(np.dot(np.transpose(barlambda_tm1*self.delta), K_eta_tm1)) + self.alpha/self.n

    def update(self, nu_t):

        return np.transpose(np.dot(np.transpose(nu_t*self.q), self.G)) + self.kappa/self.n

    def run_nu(self, T):
        
        Nu = self.pi_0
        Nu_pred = -np.ones(self.pi_0.shape)

        for t in range(0, T):

            nu_t = self.prediction(Nu[:,t:t+1])
            nu_t_pred = self.update(nu_t)

            Nu = np.concatenate((Nu, nu_t), axis = 1)
            Nu_pred = np.concatenate((Nu_pred, nu_t_pred), axis = 1)

        return Nu, Nu_pred
    
class Compartmental_model_q():

    def __init__(self, pi_0, delta, K_eta, alpha, q, G, kappa, n):

        self.n     = n
        self.pi_0  = pi_0
        self.delta = delta
        self.K_eta = K_eta
        self.alpha = alpha
        self.q     = q    
        self.G     = G    
        self.kappa = kappa
        self.t     = 0

    def step_0(self):

        return np.transpose(np.random.multinomial(self.n, self.pi_0.squeeze(), size = 1))
        # return (np.random.poisson(self.n*self.pi_0))

    def step_t(self, x_tm1):
        
        # Latent process
        # deaths
        # deaths    = np.random.binomial(x_tm1, (1-self.delta))
        barx_tm1  = np.random.binomial(x_tm1, self.delta)

        # transitions
        K_eta_tm1     = self.K_eta(barx_tm1)
        transitions_x = np.array([np.random.multinomial(barx_tm1[i], K_eta_tm1[i,:], 1)[0] for i in range(0, len(barx_tm1))])
        tildex_t      = np.transpose(np.sum(transitions_x, axis = 0, keepdims=True))

        # births
        births = np.random.poisson(self.alpha)
        x_t    = tildex_t + births

        # Observed process

        q = np.zeros(4)
        q[2] = self.q[self.t]
        ## Next line does random under reporting for each compartment not just infectives...
        bary_t        = np.random.binomial(x_t, q)[:,2]
        #transitions_y = np.array([np.random.multinomial(bary_t, self.G[i,:], 1)[0] for i in range(0, len(bary_t))])
        #tildey_t      = np.transpose(np.sum(transitions_y, axis = 0, keepdims=True))
        self.t += 1
        # clutter
        haty_t = np.random.poisson(self.kappa)
        y_t    = bary_t

        return x_t, y_t, q[2]

    def run(self, T):

        X = self.step_0()
        Y = np.zeros(X.shape)
        q_sim = np.zeros(1)

        for t in range(0, T):

            x_t, y_t, q = self.step_t(X[:,t:t+1])

            X = np.concatenate((X, x_t), axis = 1)
            Y = np.concatenate((Y, np.expand_dims(y_t,axis =1)), axis = 1)
            q_sim = np.concatenate((q_sim,  np.array([q])), axis = 0)

        return X, Y, q_sim

    def prediction(self, barlambda_tm1):
    
        K_eta_tm1 = self.K_eta(barlambda_tm1*self.delta)

        return np.transpose(np.dot(np.transpose(barlambda_tm1*self.delta), K_eta_tm1)) + self.alpha/self.n

    def update(self, nu_t):

        return np.transpose(np.dot(np.transpose(nu_t*self.q), self.G)) + self.kappa/self.n

    def run_nu(self, T):
        
        Nu = self.pi_0
        Nu_pred = -np.ones(self.pi_0.shape)

        for t in range(0, T):

            nu_t = self.prediction(Nu[:,t:t+1])
            nu_t_pred = self.update(nu_t)

            Nu = np.concatenate((Nu, nu_t), axis = 1)
            Nu_pred = np.concatenate((Nu_pred, nu_t_pred), axis = 1)

        return Nu, Nu_pred



class LaPAL_approx():

    def __init__(self,y, pi_0, delta, K_eta, alpha, q_pars, G, kappa, n):
        self.y = y
        self.lambda_0 = pi_0*n
        self.delta = delta
        self.K_eta = K_eta
        self.alpha = alpha
        self.q_pars = q_pars
    
    def prediction_step(self, bar_lambda_tm1):
        
        K = self.K_eta(bar_lambda_tm1*self.delta)
        lambda_t = np.transpose(np.dot(np.transpose(bar_lambda_tm1*self.delta), K)) + np.transpose(self.alpha)

        return lambda_t

    def update_step(self,lambda_t,y_t):
        
        q_star = 0.5*(self.q_pars[0] - lambda_t[:,2]*(self.q_pars[1]*self.q_pars[1]) + np.sqrt((lambda_t[:,2]*self.q_pars[1]**2- self.q_pars[0])*(lambda_t[:,2]*self.q_pars[1]**2- self.q_pars[0]) + 4*y_t*self.q_pars[1]*self.q_pars[1]))

        q_var = 1/(y_t/(q_star**2) + 1/(self.q_pars[1]**2))

        q = np.array((0,0,float(q_star),0))
        y = np.array((0,0,float(y_t),0))

        bar_lambda_t = (1 - q)*lambda_t + y
       # bar_lambda_t[:,2] = update
        lower = 0
        upper = 1
        lik = stats.poisson.logpmf(y_t, lambda_t[:,2]*q_star) + stats.truncnorm.logpdf(q_star,(lower-self.q_pars[0])/self.q_pars[1],(upper-self.q_pars[0])/self.q_pars[1],loc=self.q_pars[0],scale=self.q_pars[1]) + np.log(np.sqrt(q_var*2*3.14159))

        return lik, bar_lambda_t, q_star, q_var
    
    def step_t(self, barlambda_tm1, y_t):

        lambda_t = self.prediction_step(barlambda_tm1)
        lik, barlambda_t, q_star, q_var = self.update_step(lambda_t,y_t)

        return lambda_t, barlambda_t, lik, q_star, q_var

    def run(self):
        
        lambda_ = self.lambda_0
        barlambda = self.lambda_0
        qmean = np.zeros(1)
        qvar = np.zeros(1)
        logw = np.zeros(1)
        
        for t in range(0, len(self.y)):

            lambda_t, barlambda_t, lik, q_star, q_var = self.step_t(barlambda[:,t], self.y[t]) 

            lambda_ = np.concatenate((lambda_,np.transpose(lambda_t)), axis = 1)
            barlambda = np.concatenate((barlambda,np.transpose(barlambda_t)), axis = 1)
            qmean = np.concatenate((qmean, q_star), axis = 0)
            qvar = np.concatenate((qvar, q_var), axis = 0)
            logw = np.concatenate((logw,lik), axis = 0)

        return lambda_, barlambda, qmean, qvar, logw


# pi_0_true = np.transpose(np.array([[0.99, 0.0, 0.01, 0.0]]))
# 
# # transition kernel 
# beta_true  = 0.5  # transmiss rate 
# rho_true   = 0.05 # latent period rate
# gamma_true = 0.1  # recovery rate
# 
# K_eta_true = K_eta_SEIR( beta_true, rho_true, gamma_true)
# 
# # no death
# delta_true = np.ones((4, 1))
# 
# # emission distribution
# q_true = np.transpose(np.array([[0.5, 0.1]]))
# 
# alpha_true    = np.zeros((4, 1))
# kappa_true    =  np.zeros((4, 1))
# 
# # no misreporting
# G_true = np.eye(4) 
# 
# T = 200
# 
# MODEL = Compartmental_model(pi_0_true, delta_true, K_eta_true, alpha_true, q_true, G_true, kappa_true, 1000)
# X, Y, q_sim = MODEL.run(100)
# y = Y[2,1:201]
# approx = LaPAL_approx(y, pi_0_true, delta_true[0], K_eta_true, alpha_true, q_true, G_true, kappa_true, 1000)
# lambda_, barlambda, qmean, qvar, logw =  approx.run()
