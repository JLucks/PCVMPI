#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include "mpi.h"
#include <windows.h>

#define LIM 100
#define false 0
#define true 1

typedef struct _city_coord_{
    int nome;
    int x;
    int y;
}Cidade;

void lerCidades(char *url, Cidade S[], int tam);
void imprimirCidades(Cidade S[], int tam);
float distanciaCidades(Cidade S[], int a, int b);
void matrizDistanciaCidades(float **D, Cidade S[], int tam);
void imprimirMatrizDistanciaCidades(float **D, int tam);
void salvarCaminho(int *caminhoAtual, int *caminhoFinal, int tam);
void branchAndBound(float **adj, float curr_bound, float curr_weight,int level, int *curr_path, int tam, float *final_res, int *visited, int *path);
void PCV(float **adj, int tam, int ind, float bound, float *boundTemp, int *caminho);
void imprimirResultPCV(Cidade S[], int *caminhoFinal, int tam, float bound);
float newBound(float **M, int num);
void atualizaBound(float *atual, float novo);
void atualizaCaminho(float *atual, float novo, int *caminho, int *caminhofinal, int tam);

int main(int argc, char *argv[]){
    int tam = atoi(argv[1]);
    char* url = argv[2];
	if(tam > LIM || tam <= 2 || argc != 3){
        printf("Erro, parametros invalidos\n");
        exit(1);
	}
	Cidade S[tam];
	int i,j, *caminho;
	float **D, bound, tempoT;
	D = (float **) malloc(sizeof(float*)* tam);
	caminho = (int *)malloc(sizeof(int)*(tam+1));
	if(D == NULL || caminho == NULL){
        printf("Erro, nao foi possivel criar matriz");
        exit(1);
	}
	for(i = 0; i < tam; i++){
        D[i] = (float *) malloc(sizeof(float)* tam);
        if(D[i] == NULL){
            for(j = i; j >= 0; j--)
                free(D[j]);
            free(D);
            printf("Erro, nao foi possivel criar matriz\n");
            exit(1);
        }
	}
	lerCidades(url,S,tam);
	imprimirCidades(S,tam);
	matrizDistanciaCidades(D,S,tam);
	imprimirMatrizDistanciaCidades(D,tam);
    int numtasks, rank, len, termino = false;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    MPI_Init(&argc, &argv);
    MPI_Get_processor_name(hostname, &len);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    float tempbound, my_tempo;
    int *tempcaminho;
    int my_termino = false, k = 0;
    tempcaminho = (int *)malloc(sizeof(int)*(tam+1));
    if(tempcaminho == NULL){
        printf("Erro, nao foi possivel criar matriz");
        exit(1);
    }
    MPI_Reduce(&my_termino,&termino,1,MPI_INT,MPI_LAND,0,MPI_COMM_WORLD);
    clock_t tempo;
	tempo = clock();
    while(!termino){
        if(rank == 0){
            tempbound = newBound(D,tam);
            if(tempbound > -1)
                atualizaBound(&bound,tempbound);
            MPI_Reduce(&bound,&bound, 1, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
        }
        else{
            int localn = tam / (numtasks - 1);
            int ini = (rank - 1) * localn + 1;
            int fim = ini + localn - 1;
            for(k = ini; k<= fim; k++){
                PCV(D,tam,k,bound,&tempbound,tempcaminho);
                atualizaCaminho(&bound,tempbound,tempcaminho,caminho,tam);
                MPI_Reduce(&bound,&bound, 1, MPI_FLOAT, MPI_MIN, 0, MPI_COMM_WORLD);
            }
        }
        if(tempbound <= -1 || k > 0)
            my_termino = true;
        MPI_Reduce(&my_termino,&termino,1,MPI_INT,MPI_LAND,0,MPI_COMM_WORLD);
    }
    my_tempo = (clock() - tempo) / (double)CLOCKS_PER_SEC;
    MPI_Reduce(&my_tempo,&tempoT,1,MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD);
    MPI_Finalize();
    printf("Tempo:%f\n",tempoT);
    imprimirResultPCV(S,caminho,tam,bound);
    free(tempcaminho);
    free(caminho);
    for(i = 0; i < tam; i++)
        free(D[i]);
    free(D);
	exit(1);
}

void atualizaBound(float *atual, float novo){
     if(novo < *atual)
        *atual = novo;
}

void atualizaCaminho(float *atual, float novo, int* caminho, int *caminhofinal, int tam){
    if(novo <= *atual){
        *atual = novo;
        salvarCaminho(caminho,caminhofinal,tam);
    }
}

void lerCidades(char *url, Cidade S[], int tam){
    FILE *arq;
    int i = 0;
	arq = fopen(url, "r");
	if(arq == NULL)
			printf("Erro, nao foi possivel abrir o arquivo\n");
	else{
		while( (fscanf(arq,"%d %d\n", &S[i].x, &S[i].y))!=EOF ){
            S[i].nome = i;
            i++;
            if(i == tam)
                break;
		}
	}
	fclose(arq);
}

void imprimirCidades(Cidade S[], int tam){
    int i;
    printf("\n-------------------\n");
    printf("Cidades:\n");
    printf("-------------------\n");
    for(i=0;i < tam;i++)
        printf("%d:(%d,%d)\n",S[i].nome,S[i].x,S[i].y);
    printf("-------------------\n\n");
}

float distanciaCidades(Cidade S[], int a, int b){
    return sqrt((pow(S[b].x - S[a].x, 2))+(pow(S[b].y - S[a].y, 2)));
}

void matrizDistanciaCidades(float **D, Cidade S[], int tam){
    int i, j;
    for(i = 0; i < tam; i++)
        for(j = 0; j < tam; j++)
            D[i][j] = distanciaCidades(S,i,j);
}

void imprimirMatrizDistanciaCidades(float **D, int tam){
    int i, j;
    printf("\n------------------------------------------\n");
    printf("Matriz:\n");
    printf("--------------------------------------------\n");
    for(i = 0; i < tam; i++){
        for(j = 0; j < tam; j++){
            printf("%f ",D[i][j]);
        }
        printf("\n");
    }
    printf("--------------------------------------------\n\n");
}

void salvarCaminho(int *caminhoAtual, int *caminhoFinal, int tam){
    int i;
    for (i=0; i<=tam; i++)
        caminhoFinal[i] = caminhoAtual[i];
}

float firstMin(float **adj, int i, int tam){
    float min = FLT_MAX;
    int k;
    for (k=0; k<tam; k++)
        if (adj[i][k]<min && i != k)
            min = adj[i][k];
    return min;
}

float secondMin(float **adj, int i, int tam){
    float first = FLT_MAX, second = FLT_MAX;
    int j;
    for (j=0; j<tam; j++){
        if (i == j)
            continue;
        if (adj[i][j] <= first){
            second = first;
            first = adj[i][j];
        }
        else if (adj[i][j] <= second && adj[i][j] != first)
            second = adj[i][j];
    }
    return second;
}

void branchAndBound(float **adj, float curr_bound, float curr_weight,int level, int *curr_path, int tam, float *final_res, int *visited, int *path){
    if (level==tam){
        if (adj[curr_path[level-1]][curr_path[0]] != 0){
            float curr_res = curr_weight + adj[curr_path[level-1]][curr_path[0]];
            if (curr_res < *final_res){
                curr_path[level] = curr_path[0];
                salvarCaminho(curr_path,path,tam);
                *final_res = curr_res;
            }
        }
        return;
    }
    int i;
    for (i=0; i<tam; i++){
        if (adj[curr_path[level-1]][i] != 0 && visited[i] == false){
            int temp = curr_bound;
            curr_weight += adj[curr_path[level-1]][i];
            if (level==1)
              curr_bound -= ((firstMin(adj, curr_path[level-1],tam) + firstMin(adj, i,tam))/2);
            else
              curr_bound -= ((secondMin(adj, curr_path[level-1],tam) + firstMin(adj, i, tam))/2);
            if (curr_bound + curr_weight < *final_res){
                curr_path[level] = i;
                visited[i] = true;
                branchAndBound(adj, curr_bound, curr_weight, level+1,curr_path,tam,final_res,visited,path);
            }
            curr_weight -= adj[curr_path[level-1]][i];
            curr_bound = temp;
            int j;
            for(j = 0; j < tam; j++)
                visited[j] = false;
            for (j = 0; j<=level-1; j++)
                visited[curr_path[j]] = true;
        }
    }
}

void PCV(float **adj, int tam, int ind, float bound, float *boundTemp, int *caminho){
    float curr_bound = 0;
    int *curr_path;
    int *visited;
    int i;
    curr_path = (int *)malloc(sizeof(int)*(tam+1));
    visited = (int *)malloc(sizeof(int)*tam);
    if(curr_path == NULL || visited == NULL){
        printf("Erro, nao foi possivel criar matriz");
        exit(1);
    }
    for (i=0; i<tam+1; i++)
        curr_path[i] = -1;
    for (i=0; i<tam; i++)
        visited[i] = false;
    for (i=0; i<tam; i++)
        curr_bound += (firstMin(adj, i,tam) + secondMin(adj, i,tam));
    curr_bound = curr_bound/2;
    visited[ind] = true;
    curr_path[0] = ind;
    *boundTemp = bound;
    branchAndBound(adj,curr_bound,0,1,curr_path,tam,boundTemp,visited,caminho);
    free(curr_path);
    free(visited);
}

void imprimirResultPCV(Cidade S[], int *caminhoFinal, int tam, float bound){
    int i;
    printf("\n-------------------\n");
    printf("Custo Minimo : %f\n", bound);
    printf("Caminho: \n");
    printf("Cidades:\n");
    printf("-------------------\n");
    for (i=0; i<=tam; i++)
        printf("%d:(%d,%d)\n", S[caminhoFinal[i]].nome, S[caminhoFinal[i]].x, S[caminhoFinal[i]].y);
    printf("-------------------\n\n");
}

int checkFim(int V[],int num){
    int i;
    for(i = 0; i < num; i++)
        if(!V[i])
            return false;
    return true;
}

void boundVMP(float **M, int num, int i, int ini, int visited[], float *bound){
    if(checkFim(visited,num)){
        *bound += M[i][ini];
        return;
    }
    int j, l;
    float min = FLT_MAX;
    for(j = 0; j < num; j++){
        if(j == i)
            continue;
        if(M[i][j] < min && !visited[j]){
            min = M[i][j];
            l = j;
        }
    }
    visited[l] = true;
    *bound += min;
    boundVMP(M,num,l,ini,visited,bound);
}

float newBound(float **M, int num){
    static int i = 0;
    if(i >= num)
        return -1;
    int visited[num];
    int j;
    float bound = 0;
    for(j = 0; j < num; j++)
        visited[j] = false;
    boundVMP(M,num,i,i,visited,&bound);
    i++;
    return bound;
}
