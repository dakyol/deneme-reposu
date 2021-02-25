import matplotlib.pyplot as plt
import numpy as np
import math

def unitTangent(M,x):
    if x == 0:
        t = M.shape[1]-1
    else:
        t = x-1
    if x == M.shape[1]-1:
        s = 0
    else:
        s = x+1
    unit_tangent = ((M[:,x]-M[:,t])/np.linalg.norm((M[:,x]-M[:,t])) + (M[:,s]-M[:,x])/np.linalg.norm((M[:,s]-M[:,x])))/np.linalg.norm(((M[:,x]-M[:,t])/np.linalg.norm((M[:,x]-M[:,t])) + (M[:,s]-M[:,x])/np.linalg.norm((M[:,s]-M[:,x]))))
    return(unit_tangent)

def pseudoArc(M,x,y):
    if x<y:
        perra = 0
        time = x
        while time < y:
            perra += np.linalg.norm(M[:,time]-M[:,time+1])
            time += 1
        per = perra
    elif x>y:
        perra = 0
        time = 0
        while time < y:
            perra += np.linalg.norm(M[:,time]-M[:,time+1])
            time += 1
        time = x
        sınır = M.shape[1]-1
        while time < sınır:
            perra += np.linalg.norm(M[:,time]-M[:,time+1])
            time += 1
        perra += np.linalg.norm(M[:,0]-M[:,sınır-1])
        per = perra
    else:
        per = 0
    return(per)

def intersectionSolver(M,x,y):
    a = np.array([[unitTangent(M,x)[0], -unitTangent(M,y)[0]], [unitTangent(M,x)[1], -unitTangent(M,y)[1]]])
    b = np.array([M[0,y] - M[0,x], M[1,y] - M[1,x]])
    try:
        solv = np.linalg.solve(a, b)
    except np.linalg.LinAlgError as err:
        if 'Singular matrix' in str(err):
            solv = [-1, 0]
        # error handling block
    return(solv)

def jointPerimeter(M,x,y):
    solv = intersectionSolver(M,x,y)
    if solv[0] > 0:
        per = pseudoArc(M,x,y)
        perim = per + abs(solv[0]) + abs(solv[1])
    else:
        perim = np.NaN
    return(perim)

def perList(M):
    row = M.shape[1]
    per_list = np.array([[0 for x in range(row)] for y in range(row)]).astype(float)
    for i in range(row):
        for j in range(row):
            per_list[i][j] = jointPerimeter(M,i,j)
    return(per_list)

def intersectionPoint(M,x,y):
#    a = np.array([[unitTangent(M,x)[0], -unitTangent(M,y)[0]], [unitTangent(M,x)[1], -unitTangent(M,y)[1]]])
#    b = np.array([M[0,y] - M[0,x], M[1,y] - M[1,x]])
    solv = intersectionSolver(M,x,y)
#    solv = np.linalg.solve(a, b)
    int_point = M[:,x] + solv[0]*unitTangent(M,x)

    return(int_point)

def finalPlotting(M,R):
    row = M.shape[1]

    per_list = perList(M)
    per_list = abs(per_list - R)

    final_array = np.array([[0 for x in range(2)] for y in range(row)]).astype(float)
    print(final_array.shape[0])
    print(final_array.shape[1])

    for i in range(row):
        final_array[i]=intersectionPoint(M,i,np.nanargmin(per_list[i]))
    
    return(final_array)


N = np.array([2*np.sin(np.linspace(0,2*math.pi,1000)),np.cos(np.linspace(0,2*math.pi,1000))])

#print(intersectionPoint(N,0,99))
#print(finalPlotting(N,50))

Map = finalPlotting(N,20)
plt.scatter(Map[:,0],Map[:,1])
plt.show()
