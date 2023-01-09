import matplotlib
import numpy as np
import matplotlib.pyplot as plt

#fisic dimensions
L = [8, 4] # [Lx] or [Lx, Ly]
dim = len(L)

#discretization variables
n = [18, 9] # [nx] or [nx, ny]
dist = ['trig', 'trig'] # [distx] or [distx, disty]
geo = 'tri' # quad or tri

# if discretization is polinomial
npol = [0.7, 0.7] # [polyXDistribution, polyYDistribution]

# if discretization is senoidal
k = [2, 3] # [kx] or [kx, ky]
phi = [np.pi, 0] # [phix] or [phix, phiy]

# 1D Generator
if dim == 1:
    compX  = L[0] # bar length

    # discretization
    numberNodes = n[0]
    numberElements = numberNodes-1
    grid = np.array(range(numberNodes))

    cc = [0, numberNodes]
    notcc = list(range(1,numberNodes))
    print('\n', 'cc: ' , cc, '\n', 'notcc: ', notcc)
    
    Xv = np.zeros(numberNodes, dtype='float')

    distx = dist[0] # poly or trig

    if distx == 'poly':
        polyXDistribution = npol[0]
        for i in range(numberNodes):
            Xv[i] = (compX/(numberNodes-1)**polyXDistribution)*i**polyXDistribution
    elif distx == 'trig':
        kx = k[0]
        phix = phi[0]
        for i in range(numberNodes):
            Xv[i] = (-compX/(np.sin(phix - 2*np.pi*kx) - 2*np.pi*kx))*(np.sin(2*np.pi*kx*i/(numberNodes - 1) - phix) + 
                2*np.pi*kx*i/(numberNodes - 1))
    print('\n', 'Coordenada x dos pontos: ', Xv)
    print('\n', 'number of elements ', numberElements)

    dx = np.zeros(numberElements, dtype='float')
    for i in range(numberElements):
        dx[i] = Xv[i+1]-Xv[i]
    assert(np.isclose(sum(dx), compX, rtol=1e-12))
    print('\n', 'dx: ', dx)

    # IEN Matrix
    IEN = np.zeros( (numberElements,2),dtype='int' )
    for i in range(0,numberElements):
        IEN[i] = [i,i+1]
    print('\n', 'IEN: \n', IEN)

    plt.plot(Xv,Xv*0,'ko-')
    for i in range(numberNodes):
        plt.text(Xv[i] + 0.01, 0.005, str(i), color='b')
    for e in range(numberElements):
        [v1, v2] = IEN[e]
        xm = (Xv[v1]+Xv[v2])/2.0
        plt.text(xm, 0.005, str(e))
    ax = plt.gca()
    ax.set_xlabel('x')

    plt.show()

# 2d Generator
elif dim == 2:
    # Mesh params
    Lx = L[0]
    Ly = L[1]

    # discretization
    nx = n[0]
    ny = n[1]
    npoints = nx*ny
    if geo == 'quad':
        ne = (nx-1)*(ny-1)
    elif geo == 'tri':
        ne = 2*(nx-1)*(ny-1)

    grid = np.zeros((ny, nx), dtype='int')
    for i in range(ny):
        for j in range(nx):
            pto = i*nx+j
            grid[i,j] = pto

    cc1 = grid[0,:]
    cc2 = grid[1:-1,-1]
    cc3 = grid[-1,:]
    cc4 = grid[1:-1,0]
    cc = np.sort(np.concatenate((cc1, cc2, cc3, cc4)))
    notcc = []
    for i in range(1, ny-1):
        notcc.append(grid[i, 1:-1])
    notcc = np.concatenate(notcc, axis=0)
    print('\n','cc1: ' , cc1, '\n',
    'cc2: ', cc2, '\n',
    'cc3: ', cc3, '\n',
    'cc4: ', cc4, '\n',
    'cc: ', cc, '\n',
    'notcc: ', notcc)

    Xv = np.zeros(npoints, dtype='float')
    Yv = np.zeros(npoints, dtype='float')

    distx = dist[0]
    disty = dist[1]

    if distx == 'poly':
        npolx = npol[0]
        for i in range(ny):
            for j in range(nx):
                pto = i*nx+j
                Xv[pto] = (Lx/(nx-1)**npolx)*j**npolx
    elif distx == 'trig':
        kx = k[0] # periodos
        phix = phi[0] # fase
        for i in range(ny):
            for j in range(nx):
                pto = i*nx+j
                Xv[pto] = (-Lx/(np.sin(phix - 2*np.pi*kx) - 2*np.pi*kx))*(np.sin(2*np.pi*kx*j/(nx - 1) - phix) + 
                    2*np.pi*kx*j/(nx - 1))

    if disty == 'poly':
        npoly = npol[1]
        for i in range(ny):
            for j in range(nx):
                pto = i*nx+j
                Yv[pto] = (Ly/(ny-1)**npoly)*i**npoly
    elif disty == 'trig':
        ky = k[1] # periodos
        phiy = phi[1] # fase
        for i in range(ny):
            for j in range(nx):
                pto = i*nx+j
                Yv[pto] = (-Ly/(np.sin(phiy - 2*np.pi*ky) - 2*np.pi*ky))*(np.sin(2*np.pi*ky*i/(ny - 1) - phiy) + 
                    2*np.pi*ky*i/(ny - 1))
    XY = np.stack((Xv,Yv),axis=-1)
    print('\n', 'Coordenada (x, y) dos pontos: ', XY)
    print('\n', 'number of elements ', ne)

    dx = np.zeros(nx-1, dtype='float')
    dy = np.zeros(ny-1, dtype='float')

    for i in range(nx-1):
        dx[i] = Xv[i+1]-Xv[i]
    assert(np.isclose(sum(dx), Lx, rtol=1e-12))

    for i in range(ny-1):
        dy[i] = Yv[(i+1)*nx]-Yv[i*nx]
    assert(np.isclose(sum(dy), Ly, rtol=1e-12))
    
    print('\n', 'dx: ', dx, '\n', 'dy: ', dy)

    # gerar a IEN e calculo da area
    if geo == 'quad':
        IEN = np.zeros((ne,4), dtype='int')

        for i in range(ny-1):
            for j in range(nx-1):
                v1 = nx*i + j
                v2 = nx*i + j + 1
                v3 = nx*(i + 1) + j + 1
                v4 = nx*(i + 1 ) + j
                IEN[(nx-1)*i + j] = [v1, v2, v3, v4]

        ae = np.zeros(ne)

        for e in range(0,ne):
            [v1, v2, v3, v4] = IEN[e] 
            ae[e] = (Xv[v2] - Xv[v1])*(Yv[v4] - Yv[v1]) - (Xv[v4] - Xv[v1])*(Yv[v2] - Yv[v1])

    elif geo == 'tri':
        IEN = np.zeros((ne,3), dtype='int')

        for i in range(ny-1):
            for j in range(nx-1):
                    v1p = nx*i + j
                    v2p = nx*(i + 1) + j + 1
                    v3p = nx*(i + 1) + j
                    v1i = nx*i + j
                    v2i = nx*i + j + 1
                    v3i = nx*(i + 1) + j + 1
                    ep = 2*((nx-1)*i + j)
                    ei = 2*((nx-1)*i + j) + 1
                    IEN[ep] = [v1p, v2p, v3p]
                    IEN[ei] = [v1i, v2i, v3i]

        ae = np.zeros(ne)

        for e in range(0,ne):
            [v1, v2, v3] = IEN[e] 
            ae[e] = (Xv[v1]*Yv[v2] + Xv[v2]*Yv[v3] + Xv[v3]*Yv[v1] - Xv[v1]*Yv[v3] - Xv[v2]*Yv[v1] - Xv[v3]*Yv[v2])/2.0

    print('\n', 'IEN: \n', IEN)
    assert(np.isclose(np.sum(ae), Lx*Ly, rtol=1e-12))

    # plot
    verts = XY[IEN]

    plt.plot(Xv, Yv, 'ko')
    ax = plt.gca()
    pc = matplotlib.collections.PolyCollection(verts, edgecolors='black', facecolors='pink', linewidth=0.7)
    ax.add_collection(pc)
    for i in range(npoints):
        plt.text(Xv[i]+0.01, Yv[i]+0.01, str(i), color='b')
    if geo == 'quad':
        for e in range(ne):
            [v1, v2, v3, v4] = IEN[e]
            xm = (Xv[v1] + Xv[v2] + Xv[v3] + Xv[v4])/4.0
            ym = (Yv[v1] + Yv[v2] + Yv[v3] + Yv[v4])/4.0
            plt.text(xm, ym, str(e))
    
    if geo == 'tri':
        for e in range(ne):
            [v1, v2, v3] = IEN[e]
            xm = (Xv[v1] + Xv[v2] + Xv[v3])/3.0
            ym = (Yv[v1] + Yv[v2] + Yv[v3])/3.0
            plt.text(xm, ym, str(e))

    ax.set_aspect('equal')
    ax.set_xlabel('x')
    ax.set_ylabel('y')

    plt.show()
