import numpy as np
import xlrd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.pyplot as plt

authors = ['Giancarlo','Nívea','Matheus']

class Solver:
    def __init__(self):
        self.Ls = []
        self.deformations = []
        self.tensions = []
        self.forces = []

    ##---------------------PUBLIC METHODS-----------------##
    def import_file(self, filename):
        """Import the file to read its values and properties"""
        self.filename = filename
        arquivo = xlrd.open_workbook(filename)
        ################################################## Ler os nos
        nos = arquivo.sheet_by_name('Nos')
        # Numero de nos
        self.nn = int(nos.cell(1,3).value)
        # Matriz dos nós
        self.N_matrix = np.zeros((2,self.nn))
        for c in range(self.nn):
            self.N_matrix[0,c] = nos.cell(c+1,0).value
            self.N_matrix[1,c] = nos.cell(c+1,1).value
        ################################################## Ler a incidencia
        incid = arquivo.sheet_by_name('Incidencia')
        # Numero de membros
        self.nm = int(incid.cell(1,5).value)
        # Matriz de incidencia
        self.Inc = np.zeros((self.nm,4))
        for c in range(self.nm):
            self.Inc[c,0] = int(incid.cell(c+1,0).value)
            self.Inc[c,1] = int(incid.cell(c+1,1).value)
            self.Inc[c,2] = incid.cell(c+1,2).value
            self.Inc[c,3] = incid.cell(c+1,3).value
        ################################################## Ler as cargas
        carg = arquivo.sheet_by_name('Carregamento')
        # Numero de cargas
        self.nc = int(carg.cell(1,4).value)
        # Vetor carregamento
        self.F = np.zeros((self.nn*2,1))
        for c in range(self.nc):
            no = carg.cell(c+1,0).value
            xouy = carg.cell(c+1,1).value
            GDL = int(no*2-(2-xouy)) 
            self.F[GDL-1,0] = carg.cell(c+1,2).value
        ################################################## Ler restricoes
        restr = arquivo.sheet_by_name('Restricao')
        # Numero de restricoes
        self.nr = int(restr.cell(1,3).value)
        # Vetor com os graus de liberdade restritos
        self.R = np.zeros((self.nr,1))
        for c in range(self.nr):
            no = restr.cell(c+1,0).value
            xouy = restr.cell(c+1,1).value
            GDL = no*2-(2-xouy) 
            self.R[c,0] = GDL-1

    def solve_problem(self):
        """Main function to solve the problem due the entry"""
        try:
            self.N = np.transpose(self.N_matrix)
            self.__get_lengths_archs()
            self.__get_global_matrix()
            self.__get_displacements_reactions()
            self.__get_tension_deformation_force()
            print("--> Problem Solved")
            self.generate_outfile(self.filename.split('.')[0])
            print(f"--> A file named {self.filename.split('.')[0]}.txt has been generated with the results")
            print("--> Use 'solve.plot(e)', where 'e' is an factor of conversion, to see graphical results")
        except Exception as error:
         print(f"--> Error: {error}")


    def plot(self,e):
        """Plot the graph of the all elements of the structure"""
        # plt.show()
        fig = plt.figure()
        # Passa por todos os membros
        for i in range(self.nm):
            # encontra no inicial [n1] e final [n2] 
            n1 = int(self.Inc[i,0])
            n2 = int(self.Inc[i,1])   
            N1 = int(self.Inc[i,0])-1
            N2 = int(self.Inc[i,1])-1  
            N1x,N1y,N2x,N2y = 2*N1,2*N1+1,2*N2,2*N2+1           
            plt.plot([self.N_matrix[0,n1-1],self.N_matrix[0,n2-1]],[self.N_matrix[1,n1-1],self.N_matrix[1,n2-1]],color='r',linewidth=3)
            plt.plot([self.N_matrix[0,n1-1]+self.displacements[N1x]*e,self.N_matrix[0,n2-1]+self.displacements[N2x]*e],[self.N_matrix[1,n1-1]+self.displacements[N1y]*e,self.N_matrix[1,n2-1]+self.displacements[N2y]*e],'y:',linewidth=3)

        plt.xlabel('x [m]')
        plt.ylabel('y [m]')
        plt.grid(True)
        plt.axis('equal')
        plt.show()

    def generate_outfile(self,name):
        """Create an outfile with the solutions values"""
        name = name + '.txt'
        f = open(f"{name}","w+")
        f.write('Reacoes de apoio [N]\n')
        f.write(str(self.reactions))
        f.write('\n\nDeslocamentos [m]\n')
        f.write(str(self.displacements))
        f.write('\n\nDeformacoes []\n')
        f.write(str(self.deformations))
        f.write('\n\nForcas internas [N]\n')
        f.write(str(self.forces))
        f.write('\n\nTensoes internas [Pa]\n')
        f.write(str(self.tensions))
        f.close()

    ##---------------------PRIVATE METHODS-----------------##
    def __get_lengths_archs(self):
        print("--> Calculating lenght and archs for every element")
        self.angs = np.zeros((self.nm,self.nm))
        for i in range(self.nm):
            nos = int(self.Inc[i,0])-1,int(self.Inc[i,1])-1
            L = np.sqrt((self.N[nos[1],0]-self.N[nos[0],0])**2 +(self.N[nos[1],1]-self.N[nos[0],1])**2)
            self.Ls.append(L)
            self.angs[i,0] = (self.N[nos[1],0]-self.N[nos[0],0])/L
            self.angs[i,1] = (self.N[nos[1],1]-self.N[nos[0],1])/L
    
    def __fill_matrix(self,E,A,l,c,s):
        return E*A/l * np.array([[c**2,c*s,-c**2,-c*s],
                        [c*s,s**2,-c*s,-s**2],
                        [-c**2,-c*s,c**2,c*s],
                        [-c*s,-s**2,c*s,s**2]])
    
    def __get_element_matrix(self):
        print("--> Creating every element matrix")
        matrix_elements = []
        for i in range(self.nm):
            matrix_elements.append(self.__fill_matrix(self.Inc[i,2],self.Inc[i,3],self.Ls[i],self.angs[i,0], self.angs[i,1]))
        return matrix_elements

    def __get_global_matrix(self):
        print("--> Obtaining global matrix")
        matrix_elements = self.__get_element_matrix()
        self.global_matrix = np.zeros((2*len(self.N[:,0]),2*len(self.N[:,0])))
        for i in range(self.nm):
            N1,N2 = int(self.Inc[i,0]),int(self.Inc[i,1])

            grau1 = 2*N1 - 1
            grau2 = 2*N1 
            grau3 = 2*N2 - 1
            grau4 = 2*N2

            self.global_matrix[grau1-1:grau2, grau1-1:grau2] += matrix_elements[i][0:2, 0:2]
            self.global_matrix[grau3-1:grau4, grau1-1:grau2] += matrix_elements[i][2:4, 0:2]
            self.global_matrix[grau1-1:grau2, grau3-1:grau4] += matrix_elements[i][0:2, 2:4]
            self.global_matrix[grau3-1:grau4, grau3-1:grau4] += matrix_elements[i][2:4, 2:4]
    
    def __jacobi(self,A,b,erro):
        n = len(b[:,0])
        x = np.zeros((n,1))
        x0 = np.zeros((n,1))
        k = 0
        erroAtual = 1
        while erroAtual > erro:
            for i in range(n):
                soma = 0
                for j in range(n):
                    if j != i:
                        soma += A[i,j]*x0[j,0]
                x[i,0] = (b[i,0] - soma)/A[i,i]
            erroAtual = np.max(np.abs(x-x0))
            x0 = np.copy(x) 
            k +=1
        return x,k

    def __gauss_seidel(self,A, b, tol):
        n = len(b[:,0])
        x = np.zeros((n,1))
        x0 = np.zeros((n,1))
        count = 0
        while True:
            x_new = np.zeros_like(x)
            for i in range(A.shape[0]):
                s1 = np.dot(A[i, :i], x_new[:i])
                s2 = np.dot(A[i, i + 1:], x[i + 1:])
                x_new[i] = (b[i] - s1 - s2) / A[i, i]
            if np.allclose(x, x_new, rtol=tol):
                break
            x = x_new
            count += 1
        return x_new,count

    def __get_displacements_reactions(self):
        #Graus de liberdade
        print("--> Calculating displacements and reactions")
        restricted_degree = list(map(lambda x:int(x[0]),self.R))
        unrestricted_degree = [x for x in range(2*len(self.N[:,0])) if x not in restricted_degree]

        #Matrizes de resolução
        global_matrix_restricted = np.take(np.take(self.global_matrix,unrestricted_degree,axis=0),unrestricted_degree,axis=1) 
        force_vector_restricted = np.take(np.array(self.F),unrestricted_degree,axis=0)
        
        #solução para deslocamentos
        U1,count1 = self.__jacobi(global_matrix_restricted, force_vector_restricted,1e-9)
        U2,count2 = self.__gauss_seidel(global_matrix_restricted, force_vector_restricted,1e-9)
        if count2 < count1:
            U = U2.copy()
            print("--> Using Gauss-seidel numerical Solution")
        else:
            U = U1.copy()
            print("--> Using Jacobi  numerical Solution")

        self.displacements = np.zeros([len(self.F),1])

        j=0
        for i in range(len(self.F)):
            if i in unrestricted_degree:
                self.displacements[i,0] = U[j,0]
                j+=1

        L = np.dot(self.global_matrix,self.displacements)
        self.reactions = np.take(L,restricted_degree,axis=0)

    def __get_tension_deformation_force(self):
        print("--> Calculating Tension, deformation and internal forces")
        for i in range(self.nm):
            N1,N2 = int(self.Inc[i,0])-1,int(self.Inc[i,1])-1
            N1x,N1y,N2x,N2y = 2*N1,2*N1+1,2*N2,2*N2+1
            m = np.array([-self.angs[i,0],-self.angs[i,1],self.angs[i,0],self.angs[i,1]]).reshape(4,1)
            g = np.array([self.displacements[N1x],self.displacements[N1y],
                self.displacements[N2x],self.displacements[N2y]]).reshape(1,4)
            p = np.dot(g,m)[0][0]
            self.tensions.append(self.Inc[i,2]*p/self.Ls[i])
            self.forces.append(self.Inc[i,2]*p*self.Inc[i,3]/(self.Ls[i]))
            self.deformations.append(p/self.Ls[i])
