import meshio
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

# propriedade do fluido
#nu = 1.0 # Re=1
nu = 0.008 # Re=10 Re=rho * vx * H/mu = vx*H/nu = Re=275
dt = 0.003
nIter = 500

# leitura de malha e classificacao de contorno por nome (ccName)
mshname = 'cilindro.msh'
#mshname = 'cilindro.msh'
msh = meshio.read('./trabalhos/trabalho-2/cilindro/' + mshname)
X = np.array(msh.points[:,0])
Y = np.array(msh.points[:,1])
npoints = len(X)
IEN = msh.cells[1].data # triangles
ne = len(IEN)
IENbound = msh.cells[0].data # lines
IENboundTypeElem = list(msh.cell_data['gmsh:physical'][0] - 1)
boundNames = list(msh.field_data.keys())
IENboundElem = [boundNames[elem] for elem in IENboundTypeElem]

# cria lista de nos do contorno
cc = np.unique(IENbound.reshape(IENbound.size))
ccName = [[] for i in range( len(X) )]
# prioridade 4
for elem in range(0,len(IENbound)):
 if IENboundElem[elem] == 'inlet':
  ccName[ IENbound[elem][0] ] = IENboundElem[elem]
  ccName[ IENbound[elem][1] ] = IENboundElem[elem]
# prioridade 3
for elem in range(0,len(IENbound)):
 if IENboundElem[elem] == 'outlet':
  ccName[ IENbound[elem][0] ] = IENboundElem[elem]
  ccName[ IENbound[elem][1] ] = IENboundElem[elem]
for elem in range(0,len(IENbound)):
 if IENboundElem[elem] == 'cilindro':
  ccName[ IENbound[elem][0] ] = IENboundElem[elem]
  ccName[ IENbound[elem][1] ] = IENboundElem[elem]
# prioridade 2
for elem in range(0,len(IENbound)):
 if IENboundElem[elem] == 'paredeInf':
  ccName[ IENbound[elem][0] ] = IENboundElem[elem]
  ccName[ IENbound[elem][1] ] = IENboundElem[elem]
# prioridade 1
for elem in range(0,len(IENbound)):
 if IENboundElem[elem] == 'paredeSup':
  ccName[ IENbound[elem][0] ] = IENboundElem[elem]
  ccName[ IENbound[elem][1] ] = IENboundElem[elem]

#--------------------------------------------------
# # plot da malha de triangulos
# fig, ax = plt.subplots()
# ax.set_aspect('equal')
# ax.set_title("Malha de Triângulos")
# 
# ax = plt.triplot(X,Y,IEN[:,0:3],color='k',linewidth=0.5)
# ax = plt.plot(X,Y,'ko')
# for i in range(0,npoints):
#  plt.text(X[i]+0.02,Y[i]+0.03,str(i),color='b')
# for e in range(0,ne):
#  v = IEN[e]
#  xm = ( X[v[0]]+X[v[1]]+X[v[2]] )/3.0
#  ym = ( Y[v[0]]+Y[v[1]]+Y[v[2]] )/3.0
#  draw_circle = plt.Circle((xm+0.015, ym+0.015), 0.05, fill=False, color='r')
#  plt.text(xm,ym,str(e),color='r')
#  plt.gcf().gca().add_artist(draw_circle)
# 
# plt.gca().set_aspect('equal')
# plt.show()
#-------------------------------------------------- 

# cria vetor de listas de elementos vizinhos a um no i
vizElem = [[] for i in range(npoints)]
for e in range(0,ne):
 v1,v2,v3 = IEN[e]
 vizElem[v1].append(e) 
 vizElem[v2].append(e) 
 vizElem[v3].append(e) 

# definicao dos vetores de condicoes de contorno para vx,vy e psi
vx_cc = np.zeros( (npoints),dtype='float' )
vy_cc = np.zeros( (npoints),dtype='float' )
psi_cc = np.zeros( (npoints),dtype='float' )
for i in cc: 
 if ccName[i] == 'paredeInf':
  vx_cc[i] = 0.0
  vy_cc[i] = 0.0
  psi_cc[i] = 0.0
 if ccName[i] == 'paredeSup':
  vx_cc[i] = 0.0
  vy_cc[i] = 0.0
  psi_cc[i] = 1.0
 if ccName[i] == 'inlet':
  vx_cc[i] = 1.0
  vy_cc[i] = 0.0
  psi_cc[i] = Y[i]
 if ccName[i] == 'cilindro':
  vx_cc[i] = 0.0
  vy_cc[i] = 0.0
  psi_cc[i] = 0.5

# inicializando com zeros a matriz A (densa) e o vetor do lado direito b
K = np.zeros( (npoints,npoints), dtype='float')
M = np.zeros( (npoints,npoints), dtype='float')
Gx = np.zeros( (npoints,npoints), dtype='float')
Gy = np.zeros( (npoints,npoints), dtype='float')

for e in range(0,ne):
 v1,v2,v3 = IEN[e]

 # calcula a area do triangulo
 area = 0.5*np.linalg.det([[1.0,X[v1],Y[v1]],
                           [1.0,X[v2],Y[v2]],
                           [1.0,X[v3],Y[v3]]])

 # definir a matriz kelem
 b1 = Y[v2]-Y[v3]
 b2 = Y[v3]-Y[v1]
 b3 = Y[v1]-Y[v2]
 
 c1 = X[v3]-X[v2]
 c2 = X[v1]-X[v3]
 c3 = X[v2]-X[v1]

 kxelem = (1.0/(4*area))*np.array([[b1*b1,b1*b2,b1*b3],
                                   [b2*b1,b2*b2,b2*b3],
                                   [b3*b1,b3*b2,b3*b3]])
 kyelem = (1.0/(4*area))*np.array([[c1*c1,c1*c2,c1*c3],
                                   [c2*c1,c2*c2,c2*c3],
                                   [c3*c1,c3*c2,c3*c3]])
 kelem = kxelem+kyelem
 melem = (area/12)*np.array([[2.0,1.0,1.0],
                             [1.0,2.0,1.0],
                             [1.0,1.0,2.0]])
 
 gxelem = (1.0/6.0)*np.array([[b1,b2,b3],
                              [b1,b2,b3],
                              [b1,b2,b3]])

 gyelem = (1.0/6.0)*np.array([[c1,c2,c3],
                              [c1,c2,c3],
                              [c1,c2,c3]])

 for ilocal in range(0,3):
  iglobal = IEN[e,ilocal]
  for jlocal in range(0,3):
   jglobal = IEN[e,jlocal]

   K[iglobal,jglobal]  += kelem[ilocal,jlocal]
   M[iglobal,jglobal]  += melem[ilocal,jlocal]
   Gx[iglobal,jglobal] += gxelem[ilocal,jlocal]
   Gy[iglobal,jglobal] += gyelem[ilocal,jlocal]

# inversao da matriz de massa
Minv = np.linalg.inv(M)
# calculo da condicao de contorno para omega_z
omega_z_cc = Minv@(Gx@vy_cc - Gy@vx_cc)
# outra forma de resolver o sistema linear
#omega_z_cc = np.linalg.solve(M,(Gx@vy - Gy@vx))

#print ("... gravando em VTK passo de tempo: " + str(n))
point_data = {'psi_cc' : psi_cc}
data_vx_cc = {'vx_cc' : vx_cc}
data_vy_cc = {'vy_cc' : vy_cc}
data_omega_z_cc = {'omega_z_cc' : omega_z_cc}
point_data.update(data_vx_cc)
point_data.update(data_vy_cc)
point_data.update(data_omega_z_cc)
meshio.write_points_cells('condicaoDeContorno.vtk',
                           msh.points, 
                           msh.cells,
                           point_data=point_data,
                           )

# condicao inicial de vx,vy
vx = np.zeros( (npoints),dtype='float' )
vy = np.zeros( (npoints),dtype='float' )
psi = np.zeros( (npoints),dtype='float' )

for i in cc:
 vx[i] = vx_cc[i]
 vy[i] = vy_cc[i]

# calculo da condicao inicial de omega_z
omega_z = omega_z_cc.copy()

point_data = {'psi' : psi}
data_vx = {'vx' : vx}
data_vy = {'vy' : vy}
data_omega_z = {'omega_z' : omega_z}
point_data.update(data_vx)
point_data.update(data_vy)
point_data.update(data_omega_z)
meshio.write_points_cells('condicaoInicial.vtk',
                           msh.points, 
                           msh.cells,
                           point_data=point_data,
                           )

#--------------------------------------------------
# xp = X[86]
# yp = Y[86]
# 
# # plot da malha de triangulos
# fig, ax = plt.subplots()
# ax.set_aspect('equal')
# ax.set_title("Malha de Triângulos")
# ax = plt.triplot(X,Y,IEN[:,0:3],color='k',linewidth=0.5)
# ax = plt.plot(xp,yp,'go',markersize='2.0')
# plt.gca().set_aspect('equal')
# plt.show()
#-------------------------------------------------- 

# LOOP no TEMPO
for n in range(0,nIter):
 # montagem (assembling)
 # inicializando com zeros a matriz A (densa) e o vetor do lado direito b
 K = np.zeros( (npoints,npoints), dtype='float')
 M = np.zeros( (npoints,npoints), dtype='float')
 Gx = np.zeros( (npoints,npoints), dtype='float')
 Gy = np.zeros( (npoints,npoints), dtype='float')

 for e in range(0,ne):
  v1,v2,v3 = IEN[e]
 
  # calcula a area do triangulo
  area = 0.5*np.linalg.det([[1.0,X[v1],Y[v1]],
                            [1.0,X[v2],Y[v2]],
                            [1.0,X[v3],Y[v3]]])
 
  # definir a matriz kelem
  b1 = Y[v2]-Y[v3]
  b2 = Y[v3]-Y[v1]
  b3 = Y[v1]-Y[v2]
  
  c1 = X[v3]-X[v2]
  c2 = X[v1]-X[v3]
  c3 = X[v2]-X[v1]
 
  kxelem = (1.0/(4*area))*np.array([[b1*b1,b1*b2,b1*b3],
                                    [b2*b1,b2*b2,b2*b3],
                                    [b3*b1,b3*b2,b3*b3]])
  kyelem = (1.0/(4*area))*np.array([[c1*c1,c1*c2,c1*c3],
                                    [c2*c1,c2*c2,c2*c3],
                                    [c3*c1,c3*c2,c3*c3]])
  kelem = kxelem+kyelem
  melem = (area/12)*np.array([[2.0,1.0,1.0],
                              [1.0,2.0,1.0],
                              [1.0,1.0,2.0]])
  
  gxelem = (1.0/6.0)*np.array([[b1,b2,b3],
                               [b1,b2,b3],
                               [b1,b2,b3]])
 
  gyelem = (1.0/6.0)*np.array([[c1,c2,c3],
                               [c1,c2,c3],
                               [c1,c2,c3]])
  
  for ilocal in range(0,3):
   iglobal = IEN[e,ilocal]
   for jlocal in range(0,3):
    jglobal = IEN[e,jlocal]
 
    K[iglobal,jglobal]  += kelem[ilocal,jlocal]
    M[iglobal,jglobal]  += melem[ilocal,jlocal]
    Gx[iglobal,jglobal] += gxelem[ilocal,jlocal]
    Gy[iglobal,jglobal] += gyelem[ilocal,jlocal]
 
 # inversao da matriz de massa
 Minv = np.linalg.inv(M)
 # calculo da condicao de contorno para omega_z


 # calculo da condicao de contorno de omega_z
 omega_z_cc = Minv@(Gx@vy - Gy@vx)
 # montagem da matriz A
 vx_diag = np.diag(vx)
 vy_diag = np.diag(vy)
 # montagem da matriz A e do vetor b de transporte de vorticidade
 #A = M/dt + nu*K + vx_diag@Gx + vy_diag@Gy # implicito para conv e difusao
 A = M/dt + nu*K # implicito para conv e difusao
 b = (M/dt)@omega_z

 X = X + vx*dt
 Y = Y + vy*dt
 
 # condicao de contorno para o sistema linear Ax=b
 for i in cc:
  if ccName[i] == 'paredeSup' or \
     ccName[i] == 'paredeInf' or \
     ccName[i] == 'cilindro' or \
     ccName[i] == 'inlet':
   A[i,:] = 0.0 # zerando a linha
   A[i,i] = 1.0 # colocando 1 na diagonal
   b[i]   = omega_z_cc[i]
 
 # solucao do sistema linear para omega_z
 omega_z = np.linalg.solve(A,b)
 
 # solucao da Eq. de funcao-corrente
 # montagem da matriz A e do vetor b de funcao-corrente
 Apsi = K.copy()
 bpsi = M@omega_z
 
 # condicao de contorno para o sistema linear Ax=b
 for i in cc:
  if ccName[i] == 'paredeSup' or \
     ccName[i] == 'paredeInf' or \
     ccName[i] == 'cilindro' or \
     ccName[i] == 'inlet':
   Apsi[i,:] = 0.0 # zerando a linha
   Apsi[i,i] = 1.0 # colocando 1 na diagonal
   bpsi[i]   = psi_cc[i]
 
 
 # solucao do sistema linear para PSI
 psi = np.linalg.solve(Apsi,bpsi)
 
 # extracao dos campos de velocidade vx,vy a partir de PSI
 vx =      Minv@(Gy@psi)
 vy = -1.0*Minv@(Gx@psi)
 # outra forma de resolver o sistema linear
 #vx =  np.linalg.solve(M,(Gy@psi))
 #vy = -np.linalg.solve(M,(Gx@psi))
 
 for i in cc:
  if ccName[i] == 'paredeSup' or \
     ccName[i] == 'paredeInf' or \
     ccName[i] == 'cilindro' or \
     ccName[i] == 'inlet':
   vx[i]   = vx_cc[i]
   vy[i]   = vy_cc[i]

#--------------------------------------------------
#  # calculando a posicao da particula em n+1 usando v[677]
#  vxp = vx[86] 
#  vyp = vy[86]
#  xp = xp + vxp*10*dt
#  yp = yp + vyp*10*dt
#  print (xp,yp)
#-------------------------------------------------- 
 
#--------------------------------------------------
#  # plot da malha de triangulos
#  fig, ax = plt.subplots()
#  ax.set_aspect('equal')
#  ax.set_title("Malha de Triângulos")
#  ax = plt.triplot(X,Y,IEN[:,0:3],color='k',linewidth=0.5)
#  ax = plt.plot(xp,yp,'go',markersize='2.0')
#  plt.gca().set_aspect('equal')
#  plt.show()
#-------------------------------------------------- 

 
 print ("... gravando em VTK passo de tempo: " + str(n))
 imsh_points = np.column_stack((X,Y,0*X))
 imsh_cells = {'triangle':IEN}
 point_data = {'psi' : psi}
 data_vx = {'vx' : vx}
 data_vy = {'vy' : vy}
 data_omega_z = {'omega_z' : omega_z}
 point_data.update(data_vx)
 point_data.update(data_vy)
 point_data.update(data_omega_z)
 meshio.write_points_cells('solucao-'+str(n)+'.vtk',
                            imsh_points, 
                            imsh_cells,
                            point_data=point_data,
                            )
 
 
