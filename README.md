# Actividad-6

Dana Marian Rivera Oropeza - A00830027 
Daniela Berenice Hernández de Vicente - A01735346 
Alejandro Armenta Arellano - A01734879

Antes de iniciar con el código realizamos las tranformaciones de las matrices utilizando las matrices de movimiento según sea el caso, para poder representar el movimiento que se genera y se ajusta el código según sea el tipo de articulación y movimiento que esta genera.
Después de eso iniciamos con el reset de las variables y el comando tic que nos regresa al final del código con su contrapaarte toc el tiempo de ejecución del mismo.
Creamos las variables simbólicas que varían según el código que utilizamos dependiendo de la cantidad de articulaciones y su tipo de movimiento.

``` matlab
clear all
close all
clc

tic
%Declaración de variables simbólicas
syms th1(t) th2(t) t %Angulos de cada articulación
syms th1p(t) th2p(t) %Velocidades de cada articulación
syms th1pp(t) th2pp(t) %Aceleraciones de cada articulación
syms m1 m2 m3 Ixx1 Iyy1 Izz1 Ixx2 Iyy2 Izz2  %Masas y matrices de Inercia
syms l1 l2 lc1 lc2 %l=longitud de eslabones y lc=distancia al centro de masa de cada eslabón
syms pi g a cero
```

Generamos los vectores con las velocidades y aceleraciones deribando el vector de las posiciones

```matlab
  Q= [th1; th2]; 
  Qp= [th1p; th2p];
  Qpp= [th1pp; th2pp];
 ```
Creamos las matrices de rotación y los vectores de posición para cada articulación, según el caso del robot y generamos las matrices de transformación

``` matlab
P(:,:,1)= [l1*cos(th1);l1*sin(th1); 0];
R(:,:,1)= [cos(th1) -sin(th1)  0;
           sin(th1)  cos(th1)  0;
           0         0         1];
           
A(:,:,GDL)=simplify([R(:,:,GDL) P(:,:,GDL); Vector_Zeros 1]);
T(:,:,GDL)=simplify([R(:,:,GDL) P(:,:,GDL); Vector_Zeros 1]);
PO(:,:,GDL)= P(:,:,GDL); 
RO(:,:,GDL)= R(:,:,GDL); 
```

Realizamos las derivadas para la matriz de jacobianos y de igual manera lo obtenemos de manera analítica

```matlab
Jv11= functionalDerivative(PO(1,1,GDL), th1);
Jv12= functionalDerivative(PO(1,1,GDL), th2);
%Derivadas parciales de y respecto a th1 y th2
Jv21= functionalDerivative(PO(2,1,GDL), th1);
Jv22= functionalDerivative(PO(2,1,GDL), th2);
%Derivadas parciales de z respecto a th1 y th2
Jv31= functionalDerivative(PO(3,1,GDL), th1);
Jv32= functionalDerivative(PO(3,1,GDL), th2);

%Creamos la matríz del Jacobiano lineal
jv_d=simplify([Jv11 Jv12;
              Jv21 Jv22;
              Jv31 Jv32]);
%pretty(jv_d);

%Calculamos el jacobiano lineal de forma analítica
Jv_a(:,GDL)=PO(:,:,GDL);
Jw_a(:,GDL)=PO(:,:,GDL);
```
Creamos los vectores de posicón con respecto al centro de masa y las amtrices de inercia para cada articulación y después extraemos las velocidades lineales y angulares en cada eje y en cada ángulo de Euler respectivamente

``` matlab
P01=subs(P(:,:,1), l1, lc1);
 P12=subs(P(:,:,2), l2, lc2);

I1=[Ixx1 0 0; 
    0 Iyy1 0; 
    0 0 Izz1];

I2=[Ixx2 0 0; 
    0 Iyy2 0; 
    0 0 Izz2];

V=V(t);
Vx= V(1,1);
Vy= V(2,1);
Vz= V(3,1);

W=W(t);
W_pitch= W(1,1);
W_roll= W(2,1);
W_yaw= W(3,1);
```

Se calcula la energía cinética para cada eslabón, la energía potencial
``` matlab
V1_Total= V1+cross(W1,P01);
K1= (1/2*m1*(V1_Total))'*(1/2*m1*(V1_Total)) + (1/2*W1)'*(I1*W1);
K1= simplify (K1);

h1= P01(2);
U1=m1*g*h1;
```

Finalmente para las ecuaciones que modelan el movimiento se crea un vector con las velocidades y las aceleraciones y se deriban parcialmente para definir el torque de cada articulación

``` matlab
 Qd=[th1p(t); th2p(t); th1pp(t); th2pp(t)];

 dQ1=[diff(diff(Lagrangiano,th1p), th1),diff(diff(Lagrangiano,th1p), th2),
     diff(diff(Lagrangiano,th1p), th1p),diff(diff(Lagrangiano,th1p), th2p)];

 t1= dQ1*Qd- diff(Lagrangiano, th1);
```

Extraemos los coeficientes de las aceleraciones, y realizamos derivadas parciales con respecto a cada una de las variables, realizamos los cálculos de la energía cinética en su forma matricial para obtener la fuerza centrípeta y de coriolis y finalmente realizamos los pares gravitacionales, donde se sustituyen las velocidades y aceleraciones y a la última aceleración es usada como el torque en el motor correspondiente y es este torque en cada motor correspondiente a una articulación que obtenemos los pares gravitacionales y finalmente se encuentra la línea toc que nos regresa el tiempo de ejecución del programa.

```matlab
M=[diff(t1, th1pp), diff(t1, th2pp);
   diff(t2, th1pp), diff(t2, th2pp)];
rank (M);

M=M(t);

 M11=[diff(M(1,1),th1), diff(M(1,1),th2)]*Qp; M12=[diff(M(1,2),th1), diff(M(1,2),th2)]*Qp;

 M21=[diff(M(2,1),th1), diff(M(2,1),th2)]*Qp;
 M22=[diff(M(2,2),th1), diff(M(2,2),th2)]*Qp;

  Mp=[M11, M12;...
     M21, M22];

k=1/2*transpose(Qp)*M*Qp;

dk=[diff(k, th1); diff(k, th2)];

 C= Mp*Qp-dk;
 pretty(C)

 r=cero;
 a1=subs(t1, th1p, r);
 a2=subs(a1, th2p, r);
 a3=subs(a2, th1pp,r);
 a4=subs(a3, th2pp,r);
 
 G1=a4;
 
 b1=subs(t2, th1p, r);
 b2=subs(b1, th2p, r);
 b3=subs(b2, th1pp,r);
 b4=subs(b3, th2pp,r);

 G2=b4;

G=[G1;G2];
pretty(G)

 toc
 ```


## Resultados
### Cilíndrico:
![image](https://user-images.githubusercontent.com/100874942/224981349-53cf6b6e-f24a-49ca-b458-ff5dc31fcf8a.png)
### Péndulo:
![image](https://user-images.githubusercontent.com/100874942/224981653-7daa2218-5776-4591-b6cc-fc833e9a1a0e.png)
### Angular
![image](https://user-images.githubusercontent.com/100874942/224982900-40370d12-539b-463e-8bd0-994e080e6cbc.png)

Los resultados completos eran absurdamente extensos :(
