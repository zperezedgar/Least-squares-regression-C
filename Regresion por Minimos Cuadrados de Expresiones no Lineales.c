#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int n=1, m, c, i, j, k, p, interpolar, polinomio, h, trans;
double factor, valor, elemento, suma, prom;
double dummy, big, Sm, pol, Saj, r;

int main(int argc, char *argv[]) 
{
	printf("\nREGRESION POR MINIMOS CUADRADOS DE EXPRESIONES NO LINEALES\n");
	printf("\nIntroduzca el Numero de Puntos n\n");
	scanf("%i", &m);

	double x[m];
	double y[m];

	c=1;
	for(i=0; i<m; i++) //Peticion de puntos
	{
		printf("\nIngrese X%i\n", c);
		scanf("%lf", &elemento);
		x[i]=elemento;
		
		printf("\nIngrese Y%i\n", c);
		scanf("%lf", &elemento);
		y[i]=elemento;
		c++;
	}

	system("cls");
	
	double x2[m];	//Arreglos para las transformaciones
	double y2[m];
	
	//Pedimos que transformacion se desea utilizar
	do
	{
		printf("\nElija la Transformacion que desea utilizar.\n");
		printf("\n1: Modelo Potencial y=a(x^b)");
		printf("\n2: Modelo Exponencial y=a(e^(b*x))");
		printf("\n3: Modelo Exponencial Generalizado y=a(b^x)");
		printf("\n4: Modelo Pseudo sigmoideo y=a(e^(-b/x))");
		printf("\n5: Modelo Pseudo asintotico y=a-(b/x)");
		printf("\n6: Inversa y=a(x/(b+x))\n");
		scanf("%i", &trans);
		system("cls");
		
		
		switch(trans)
		{
			case 1:
				{
					for(i=0; i<m; i++)
					{
						x2[i]=log(x[i]);
						y2[i]=log(y[i]);
					}
				}break;
			case 2:
				{
					for(i=0; i<m; i++)
					{
						y2[i]=log(y[i]);
						x2[i]=x[i];
					}
				}break;
			case 3:
				{
					for(i=0; i<m; i++)
					{
						y2[i]=log(y[i]);
						x2[i]=x[i];
					}
				}break;
			case 4:
				{
					for(i=0; i<m; i++)
					{
						x2[i]=1/x[i];
						y2[i]=log(y[i]);
					}
				}break;
			case 5:
				{
					for(i=0; i<m; i++)
					{
						x2[i]=1/x[i];
						y2[i]=y[i];
					}
				}break;
			case 6:
				{
					for(i=0; i<m; i++)
					{
						x2[i]=1/x[i];
						y2[i]=1/y[i];
					}
				}break;
			default:
				printf("\nOpcion no Valida!!.");
		}
	
		
		//Construimos la matriz de coeficientes
		double M[n+1][n+2];
	
		for(i=0; i<(n+1); i++)
		{
			for(j=0; j<(n+1); j++)
			{
				M[i][j]=0;	//Inicializamos con cero
			
				for(k=0; k<m; k++)
				{
					if((i==0)&&(j==0))
					{
						M[i][j]=m;
					}
					else
					{
						M[i][j]+=pow(x2[k],(i+j));
					}
				}
			}
		}
	
		//Creacion de la Matriz de Terminos independientes
		for(i=0; i<(n+1); i++)
		{
			M[i][n+1]=0;
		
			for(j=0; j<m; j++)
			{
				if(i==0)
				{
					M[i][n+1]+=y2[j];
				}
				else
				{
					M[i][n+1]+=y2[j]*pow(x2[j],i);
				}
			}
		}
	
	
		//Aplicar Algoritmo de Gauss-Jordan(o cualquier otro) para resolver el sistema
		//Eliminacion hacia adelante
		h=n+1;
		for(k=0; k<=(h-1); k++)
		{
			if(M[k][k]==0)
			{
		
			//****************** PIVOTEO PARCIAL ********************************//
			p=k;
			big=abs(M[k][k]);
		
			for(i=k+1; i<h; i++)	//Rutina para identificar el coeficiente mayor
			{
				dummy=abs(M[i][k]);
				if(dummy>big) 
				{
					big=dummy;
					p=i;
				}
			}
			
			if(p!=k)	//Burbuja
			{
				for(j=k; j<=(h-1); j++)
				{
					dummy=M[p][j];
					M[p][j]=M[k][j];
					M[k][j]=dummy;
				}
				dummy=M[p][h];
				M[p][h]=M[k][h];
				M[k][h]=dummy;
				//cont++;	Si necesitamos el determinante
			}
			//******************************************************************//
			}
		
		
			for(i=0; i<=(h-1); i++)
			{
				if(i!=k)
				{
					factor=M[i][k]/M[k][k];
			
					for(j=0; j<=h; j++)
					{
						M[i][j]=M[i][j]-((M[k][j])*factor);
					}	
				}
			}
		}
	
		//Normalizacion
		for(i=0; i<=(h-1); i++)
		{
			M[i][h]=M[i][h]/M[i][i];
			M[i][i]=M[i][i]/M[i][i];
		}
	
		
		printf("\nLa Linealizacion es:\n");  //Solucion del sistema
		printf("\nY=%g+%gX\n", M[0][n+1], M[1][n+1]);
		
		//Regresamos al espacio original
		printf("\nLa Ecuacion original es:\n");  
		switch(trans)
		{
			case 1:
				{
						printf("\nY=%g*(X^%g)", exp(M[0][n+1]), M[1][n+1]);

				}break;
			case 2:
				{
						printf("\nY=%g*(e^(%gX))", exp(M[0][n+1]), M[1][n+1]);
				
				}break;
			case 3:
				{
						printf("\nY=%g*(%g^X)", exp(M[0][n+1]), exp(M[1][n+1]));
				}break;
			case 4:
				{
					printf("\nY=%g*(e^(-%g/X))", exp(M[0][n+1]), -M[1][n+1]);
				}break;
			case 5:
				{
					printf("\nY=%g-(%g/X)", M[0][n+1], -M[1][n+1]);
				}break;
			case 6:
				{
					printf("\nY=%g*(X/(%g+X))", 1/(M[0][n+1]), (M[1][n+1]))/(M[0][n+1]);
				}break;
			default:
				printf("\nOpcion no Valida!!.");
		}
		
		
	
		//calculo del coeficiente de correlacion
		
		//calculo del promedio
		suma=0.0;
		
		for(i=0; i<m; i++)
		{
			suma+=y2[i];
		}
	
		prom=suma/m;
		
		//calculo de Sm
		suma=0.0;
		
		for(i=0; i<m; i++)
		{
			suma+=pow((prom-y2[i]),2.0);	
		}	
		
		Sm=suma;
		
		//Calculo de Saj
		suma=0.0;
		
		for(i=0; i<m; i++)
		{
			//Necesitamos evaluar el polinomio en la i en turno P(Xi)
			pol=M[0][n+1];
			for(j=1; j<(n+1); j++)
			{
				pol+=M[j][n+1]*pow(x2[i],j);	 
			//	printf("\nPol%i=%g", j, pol);	//DEPURANDO!!!
			}
			
			suma+=pow((y2[i]-pol),2.0);
		}
		
		Saj=suma;
		
		r=sqrt(fabs((Sm-Saj)/Sm));
	
		printf("\n\nCoeficiente de Correlacion r=%g", r);
	
	
		do
		{
			printf("\n\nIngrese valor a Interpolar o Extrapolar\n");
			scanf("%lf", &valor);
	
			
			switch(trans)
			{
				case 1:
					{
						suma=exp(M[0][n+1])*pow(valor,M[1][n+1]);
					}break;
				case 2:
					{
						suma=exp(M[0][n+1])*pow(exp(1),(M[1][n+1]*valor));
				
					}break;
				case 3:
					{
						suma=exp(M[0][n+1])*pow(exp(M[1][n+1]),valor);
					}break;
				case 4:
					{
						suma=exp(M[0][n+1])*pow(exp(1),-(-M[1][n+1]/valor));
					}break;
				case 5:
					{
						suma=M[0][n+1]-(-M[1][n+1]/valor);
					}break;
				case 6:
					{
						suma=(1/(M[0][n+1]))*((valor)/(((M[1][n+1])*(M[0][n+1]))+valor));
					}break;
				default:
					printf("\nOpcion no Valida!!.");
			}
			
	
			system("cls");	
			printf("\nPol(%g)= %g",valor, suma);
		
			printf("\n\nIngrese 1 para Interpolar o Extrapolar otro Dato.\n");
			printf("Ingrese cualquier otro numero para salir.\n");
			scanf("%i", &interpolar);
			system("cls"); 
	
		} while(interpolar==1);
		
		printf("\nIngrese 1 para cambiar la Transformacion.\n");
		printf("Ingrese cualquier otro numero para salir.\n");
		scanf("%i", &polinomio);	
		system("cls");
		
	}while(polinomio==1);
	
	
	system("PAUSE");
	return 0;
}
