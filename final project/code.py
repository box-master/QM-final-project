import numpy as np
import matplotlib.pyplot as plt

np.set_printoptions(precision=3)

step = 10000
E =  np.linspace(-0.15,-0.01,step)

def hydrogen(n,rn_2,rn_3,l):
    return 1/(8*E*n)*((1-n)*(n*(n-2)-4*l*(l+1))*rn_3-4*(2*n-1)*rn_2) # return rn_1

def hydrogen_theory(n):
    return -1/2/n**2


# recursion, n_init = 2
def spectrum(N,l,level,test):
    rn_2 = 1 
    rn_3 = -2*E
    temp = np.ones_like(E)
    temp2 =  np.array([])
    for n in range(2,N+2):
        rn_1 = hydrogen(n,rn_2,rn_3,l)
        rn_2,rn_3 = rn_1,rn_2
        temp = np.vstack([temp,rn_1])
        if n==2:
            temp2 = np.append(temp2,rn_1)
            continue
        temp2 = np.vstack([temp2,rn_1])

    temp = temp.T
    temp2 = temp2.T
    hankel_matrix = np.zeros((N//2+1,N//2+1,step))
    shifted_hankel_matrix = np.zeros((N//2,N//2,step))
    for j in range(len(hankel_matrix)):
        for k in range(len(hankel_matrix)):
            hankel_matrix[j,k,:] = temp[:,j+k]/(temp[:,j]*temp[:,k])
            if j>0 and k>0:
                shifted_hankel_matrix[j-1,k-1,:] = temp2[:,j+k-2]/(temp2[:,j-1]*temp2[:,k-1])




    a = np.linalg.eigvals(hankel_matrix.transpose(2,0,1))
    b = np.linalg.eigvals(shifted_hankel_matrix.transpose(2,0,1))
    if test==1:
        constrain = np.intersect1d(np.where((a>0).all(axis=1)),np.where((b>0).all(axis=1)))
    elif test==2:
        constrain = np.where((a>0).all(axis=1))
    else:                                    
        constrain = np.where((b>0).all(axis=1))
    ret = np.zeros_like(E)


    for i in constrain:
        ret[i] = level
        
    return ret



fig = plt.figure()
axe = plt.gca()

# for i in range(1,11):
#     axe.plot(E,spectrum(i*5,1,i,2),'.',label="N={0}".format(i*5))
for i in range(1,5):
    axe.plot(E,spectrum(30,i,i+1,2),'.',label="l={0}".format(i))
for i in range(2,9):
    plt.axvline(x=hydrogen_theory(i),color="black")
plt.xlabel("E")
plt.title("N=30")
axe.get_yaxis().set_visible(False)
plt.ylim(0.1,15)
# plt.legend()
plt.show()


