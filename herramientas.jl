
__precompile__() # Este comando es para que julia precompile el paquete

module herramientas

using SymPy

export metodo_newton,integracion_rectangulo,integracion_trapecio,integracion_simpson,metodo_RK4,metodo_RK4h,metodo_eulerimp,metodo_eulerimph,metodo_euler,metodo_eulerh

"""Para metodo_newton raphson se ingresan los valores: 
(función, punto inicial, *repeticiónes [50])
*Valor opcional, [#] valores predeterminados"""
function metodo_newton(g,x0,n=50)
    x=symbols("x")
    #Volvemos a g, como un simbolo:
    f=g(x)
    #Calculamos su derivada
    df=simplify(diff(f))
    #Queremos operar con f y df como función de julia:
    dg=lambdify(df,[x])
    #Metodo de Newton-Raphson para el punto x0
    for i in 1:n
        x0=x0-g(x0)/dg(x0)
    end
return x0
end

"""Para las funciones de integracion_ se pide 
la función (f)
el intervalo de integración de (a) a (b) 
el número de intervalos de dicha partición (n)

ingresandose de la forma: (f,a,b,*n[100000])
*Valor opcional, [#] valores predeterminados"""

function integracion_rectangulo(f::Function,a,b,intervalos=100000)
    #Graba el número de intervalos como n para menor escritura
    n=intervalos
    #Aquí se genera las particiones con el número brindado, el valor predeterminado será 100,000
    v=linspace(a,b,n)
    #Variable donde guardaremos nuestra aproximación
    Itot=0
    #Suma de todos los intervalos
    for i in 1:(n-1)
    #Condicion de este metodo
        y=(v[i]+v[i+1])/2
        I=(v[i+1]-v[i])*(f(y))
    #Recopilación de los datos calculados
        Itot+=I
    end    
    return Itot
end

function integracion_trapecio(f::Function,a,b,intervalos=100000)
        #Graba el número de intervalos como n para menor escritura
    n=intervalos
        #Aquí se genera las particiones con el número brindado, el valor predeterminado será 100,000
    v=linspace(a,b,n)
        #Variable donde guardaremos nuestra aproximación
    Itot=0
        #Suma de todos los intervalos
    for i in 1:(n-1)
                #Condicion de este metodo
        I=(v[i+1]-v[i])*(f(v[i])+f(v[i+1]))/2 
                #Recopilación de los datos calculados
        Itot+=I
    end     
    return Itot
end

function integracion_simpson(f::Function,a,b,intervalos=100000)
        #Graba el número de intervalos como n para menor escritura
    n=intervalos
        #Aquí se genera las particiones con el número brindado, el valor predeterminado será 100,000
    v=linspace(a,b,n)
        #Variable donde guardaremos nuestra aproximación
    Itot=0
        #Suma de todos los intervalos
    for i in 1:(n-1)
                #Condicion de este metodo
        m=(v[i+1]+v[i])/2
        I=(v[i+1]-v[i])*(f(v[i])+4*f(m)+f(v[i+1]))/6
                #Recopilación de los datos calculados
        Itot+=I
    end     
    return Itot
end

"""Resolución de ecuaciones diferenciales ordinarias.
Para el método de euler, euler implicito y Runge-Kutta se encuentran dos versiones.
Ambas funciones tienen paramentros de entrada:
f función definida como una función con parametro de evaluación x (vector o escalar) y un tiempo. i.e. f(x,t)
x0 ó xi condición inicial al tiempo inicial
ti timpo inicial
tf tiempo final

La primera versión se da designando (n) como el número de intervalos equipartidos de tiempo que tendrá la evaluación
llamandose a la función metodo_algo(f,x0,ti,tf,n)

La segunda versión se da designando (h) como el tamaño del intervalo equipartido de tiempo que tendrá la evaluación
llamandose a la función metodo_algoh(f,x0,ti,tf,h)

Para el metodo_eulerimp ó metodo_eulerimph se emplea un factor adiciónal: (f,x0,ti,tf,n ó h,*m[50])
Donde como es implicito, se manda a llamar al metodo de Newton-Raphson para resolver una ecuación
es así como m es el número de iteraciones de este método.
*Valor opcional, [#] valores predeterminados

"""

function metodo_RK4(f,x0,ti,tf,n)
    #Buscamos h
    h=(tf-ti)/n
    #Hacemos el primer intervalo
    listt=linspace(ti,tf,n+1)
    listx=[]
    #Tomamos el primer valor de x0 en la lista solución
    push!(listx,x0)
    x=x0
    #Método de Runge-kutta orden 4
    for i in 1:length(listt)-1
        k1=f(x,listt[i])
        k2=f(x+(h/2)*k1,listt[i]+(h/2))
        k3=f(x+(h/2)*k2,listt[i]+(h/2))
        k4=f(x+h*k3,listt[i+1])
        y=x+(h/6)*(k1+2*k2+2*k3+k4)
        #Se incerta en la lista
        push!(listx,y)
        #Se puede omitir
        x=y
    end
    return listt,listx
end

function metodo_RK4h(f,xi,ti,tf,h)
    #Calculo de n dado h
    n=Int(round((tf-ti)/h))
    return metodo_RK4(f,xi,ti,tf,n)
end

#Método de Euler USUAL (COPY PASTE)
#n es el número de intervalos que tendrá la partición
function metodo_euler(f,x0,ti,tf,n)
    #El comando linespace genera subintervalos con n+1 elEentos
    listt=linspace(ti,tf,n+1)
    #Con una distancia h, y para cada tk se cumple tk=t0+mh con m de 0 a n
    h=(tf-ti)/n
    listx=[]
    x=x0
    #Se pone la condición inicial para ti
    push!(listx,x)
    #La serie de recurrencia para nuestro intervalo
    for i in 2:(n+1)
        #x es el elEento inmediato inferior y y el superior
        y=x+h*f(x,listt[i])
        push!(listx,y)
        #Aquí se pasa al siguiente dato
        x=y
    end
    return listt, listx
end


#Se crea la misma función pero evaluando en h
function metodo_eulerh(f,x0,ti,tf,h)
    n=Int(round((tf-ti)/h))
    return metodo_euler(f,x0,ti,tf,n)
end

function metodo_eulerimp(f,x0,ti,tf,n,m=50)
    #Buscamos h
    h=(tf-ti)/n
    #Hacemos el primer intervalo
    listt=linspace(ti,tf,n+1)
    listx=[]
    #TOmamos el primer valor de x0 en la lista solución
    push!(listx,x0)
    #Método de Euler implicito
    for i in 2:length(listt)
        g(x)=x-f(x,listt[i])*h-listx[i-1]
        X=listx[i-1]
        #Aquí se encuentra "x" por el método de NR
        y=metodo_newton(g,X,m)
        #Se incerta en la lista
        push!(listx,y)
        #Se puede omitir
        x=y
    end
    return listt,listx
end

function metodo_eulerimph(f,xi,ti,tf,h,m=50)
    #Calculo de n dado h
    n=Int(round((tf-ti)/h))
    return metodo_eulerimp(f,xi,ti,tf,n,m)
end

end

