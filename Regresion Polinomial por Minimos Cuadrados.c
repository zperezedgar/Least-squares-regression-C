#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int n, m, c, i, j, k, p, interpolar, polinomio, h;
double factor, valor, elemento, suma, prom;
double dummy, big, Sm, pol, Saj, r;

int main(int argc, char *argv[]) 
{
	printf("\nBIENVENIDO AL ALGORITMO PARA LA REGRESION POLINOMIAL POR MINIMOS CUADRADOS\n");
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
	
	
	do
	{
		printf("\nIntroduzca el Grado n del Polinomio Deseado\n");
		scanf("%i", &n);
		system("cls");
		
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
						M[i][j]+=pow(x[k],(i+j));
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
					M[i][n+1]+=y[j];
				}
				else
				{
					M[i][n+1]+=y[j]*pow(x[j],i);
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
	
		
		printf("\nPolinomio de Grado %i:\n", n);  //Solucion del sistema
		printf("\nP(X)= %g + %g*X", M[0][n+1], M[1][n+1]);
		for(i=2; i<=n; i++)
		{
			printf(" + %g*(X^%i)",M[i][n+1], i);
		}
		
	
		//calculo del coeficiente de correlacion
		
		//calculo del promedio
		suma=0.0;
		
		for(i=0; i<m; i++)
		{
			suma+=y[i];
		}
	
		prom=suma/m;
		
		//calculo de Sm
		suma=0.0;
		
		for(i=0; i<m; i++)
		{
			suma+=pow((prom-y[i]),2.0);	
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
				pol+=M[j][n+1]*pow(x[i],j);	 
			//	printf("\nPol%i=%g", j, pol);	//DEPURANDO!!!
			}
			
			suma+=pow((y[i]-pol),2.0);
		}
		
		Saj=suma;
		
		r=sqrt(fabs((Sm-Saj)/Sm));
	
		printf("\n\nCoeficiente de Correlacion r=%g", r);
	
	
		do
		{
			printf("\n\nIngrese valor a Interpolar o Extrapolar\n");
			scanf("%lf", &valor);
	
			suma=M[0][n+1];  //Evitamos el error de elevar algo a la cero
		
			for(i=1; i<(n+1); i++)
			{
				suma+=M[i][n+1]*pow(valor, i);
			}
	
			system("cls");	
			printf("\nPol(%g)= %g",valor, suma);
		
			printf("\n\nIngrese 1 para Interpolar o Extrapolar otro Dato.\n");
			printf("Ingrese cualquier otro numero para salir.\n");
			scanf("%i", &interpolar);
			system("cls"); 
	
		} while(interpolar==1);
		
		
		printf("\nIngrese 1 para cambiar el Grado del Polinomio.\n");
		printf("Ingrese cualquier otro numero para salir.\n");
		scanf("%i", &polinomio);	
		system("cls");

	} while(polinomio==1);
	
	
	system("PAUSE");
	return 0;
}	
