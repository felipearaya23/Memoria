#include <bits/stdc++.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <queue>
#include <ctime>

using namespace std;

//  Struct para guardar la información de un sufijo
struct suffix
{
	int index; // Guardar el índice original
	int rank[2]; // Guarda  los ranks y el rank pair siguiente
};

// Una función comparativa usada por sort () para comparar 2 sufijos
// Compara 2 pares, y retorna 1 si el primer par es más pequeño
int cmp(struct suffix a, struct suffix b)
{
	return (a.rank[0] == b.rank[0])? (a.rank[1] < b.rank[1] ?1: 0):
		(a.rank[0] < b.rank[0] ?1: 0);
}

struct subcadena
{
  string nombre;
  int veces;
};

class comparador
{
 public:
   bool operator()(const subcadena& a, const subcadena& b)
   {
        return a.veces<b.veces;
   }
};

// Esta es la principal función que toma un string 'txt' de tamaño n como entrada,
// construye y retorna el arreglo de sufijos para el string dado
vector<int> buildSuffixArray(string txt, int n)
{
	// Estructura que almacena sufijos y sus índices
	struct suffix *suffixes = new struct suffix [n]; //cambiar a memoria dinamica

    cout << n << endl;
	// Almacena los sufijos y sus índices en un arreglo de estructuras.
	// Esta estructura es necesaria para ordenar los sufijos de manera alfabética
	// y mantener sus antiguos índices mientras se ordena
	for (int i = 0; i < n; i++)
	{
		suffixes[i].index = i;
		suffixes[i].rank[0] = txt[i] - 'a';
		suffixes[i].rank[1] = ((i+1) < n)? (txt[i + 1] - 'a'): -1;
	}

	// Ordena los sufijos usando la función comparativa
	// creada arriba.
	sort(suffixes, suffixes+n, cmp);

	// A este paso, todos los sufijos son ordenados de acuerdo a sus 2 primeros
	// caracteres. Ahora se ordenarán los sufijos de acuerdo a sus primeros 4 caracteres,
	// luego con sus primeros 8 caracteres y así sucesivamente

	int *ind = new int [n]; // Este arreglo es necesario para obtener el índice original en la estructura suffixes[] (pasado a memoria dinamica)
	// Este mapeo es necesario para la creación del próximo sufijo
	for (int k = 4; k < 2*n; k = k*2)
	{
		// Se le asigna los valores de rank e índice al primer sufijo
		int rank = 0;
		int prev_rank = suffixes[0].rank[0];
		suffixes[0].rank[0] = rank;
		ind[suffixes[0].index] = 0;

		// Se le asigna los rank a los sufijos
		for (int i = 1; i < n; i++)
		{
			// Si el primer rank y los siguientes ranks son iguales a aquellos de los sufijos
			// anteriores en el arreglo, asignar el mismo nuevo rank a ese sufijo
			if (suffixes[i].rank[0] == prev_rank &&
					suffixes[i].rank[1] == suffixes[i-1].rank[1])
			{
				prev_rank = suffixes[i].rank[0];
				suffixes[i].rank[0] = rank;
			}
			else // O si no aumentar rank y asignar
			{
				prev_rank = suffixes[i].rank[0];
				suffixes[i].rank[0] = ++rank;
			}
			ind[suffixes[i].index] = i;


		}

		// Asignar próximo rank a cada sufijo
		for (int i = 0; i < n; i++)
		{
			int nextindex = suffixes[i].index + k/2;
			suffixes[i].rank[1] = (nextindex < n)?
								suffixes[ind[nextindex]].rank[0]: -1;
		}

		// Se ordenan los sufijos según los primeros k caracteres
		sort(suffixes, suffixes+n, cmp);

		cout << k << endl;
	}

    cout << "etapa final" << endl;
	// Almacenar índices de todos los sufijos ordenados en el arreglo de sufijos
	vector<int>suffixArr;
	for (int i = 0; i < n; i++){
        suffixArr.push_back(suffixes[i].index);
        //cout << i << endl;
    }

    //suffixArr.push_back(2);
	delete [] ind;
	delete [] suffixes; // se borra la memoria dinamica utilizada
	// Retornar el arreglo de sufijos
	return suffixArr;
}

// Construcción del LCP
vector<int> lcp_str(string txt, vector<int> suffixArr)
{
	int n = suffixArr.size();

	// Para almacenar el arreglo LCP
	vector<int> lcp(n, 0);

	// Un arreglo auxiliar para almacenar el inverso del arreglo de sufijos.
	// Por ejemplo si suffixArr [0] es 5, invSuff[5] debiese ser 0.
	// Estó será usado para obtener el siguiente string del arreglo de sufijos.
	vector<int> invSuff(n, 0);

	// LLenar valores en invSuff[]
	for (int i=0; i < n; i++)
		invSuff[suffixArr[i]] = i;

	// Inicializar longitud del LCP previo
	int k = 0;

	// Procesar todos los sufijos uno por uno comenzando por
	// el primer sufijo en txt[]
	for (int i=0; i<n; i++)
	{
		/* Si el sufijo actual está ubicado en la posición n-1, entonces ya no hay
		un siguiente substring a considerar, por lo tanto el LCP no es definido
		para este substring y se hace cero. */
		if (invSuff[i] == n-1)
		{
			k = 0;
			continue;
		}

		/* j contiene el índice del siguiente substring a ser considerado
		para compararlo con el actual substring, es decir,
		el siguiente string en el arreglo de sufijos */
		int j = suffixArr[invSuff[i]+1];

		// Comienza a revisar desde el k-simo índice donde
		// al menos k-q caracteres serán similares
		while (i+k<n && j+k<n && txt[i+k]==txt[j+k])
			k++;

		lcp[invSuff[i]] = k; // LCP correspondiente para el actual sufijo.

		// Borrando el cararter de partida del string.
		if (k>0)
			k--;
	}

	// Retornar el arreglo LCP construido
	return lcp;
}

// Función para imprimir un arreglo al ejecutar el programa
void printArr(vector<int>arr, int n)
{
	for (int i = 0; i < n; i++)
		cout << arr[i] << " ";
	cout << endl;
}

// Programa principal
int main()
{

	/* Método 1, leer un archivo con una determinada cadena*/

    ifstream fin("uniprot_sprot.fasta"); //30 minutos en leer archivo principal
    if(!fin)
    {
        cerr << "Couldn't open the input file!";
        return(1);
    }

    ofstream outputfile;
    outputfile.open("substrings.txt");
    string line;
    int cantidad_ss = 1;

    getline(fin, line);

    getline(fin, line);

    while(fin){
        if(line[0] == '>'){
            outputfile << "$";
            cantidad_ss = cantidad_ss + 1;
        }
        else{
            outputfile << line;
        }
    getline(fin, line);
    }

    outputfile.close();


    int inicio = cantidad_ss - 1;
    //cantidad_ss = cantidad_ss - 1;

    cout << 0 << endl;
	ifstream file("substrings.txt");
	//ifstream file("pruebas.txt");
	//cout << 0 << endl;
	fflush(stdout); //fflush(stdout) fuerza a que se imprima el buffer de salida estándar (pantalla o monitor)
	stringstream buffer;
    //cout << 1 << endl;
    //fflush(stdout);
	buffer << file.rdbuf();
	string str = buffer.str();
	//str.erase(remove_if(str.begin(), str.end(), ::isspace), str.end());
    //cout << 1 << endl;
    //fflush(stdout);

    unsigned t0, t1, t2, t3, t4, t5;

    //tinicial = clock();
	// Creación del arreglo de sufijos
	t0 = clock();
	//fflush(stdout);
	cout << 1 << endl;
	fflush(stdout);
	vector<int>suffixArr = buildSuffixArray(str, str.length()); fflush(stdout);
	cout << 2 << endl;
	t1 = clock();
    double t12 = (double(t1-t0)/CLOCKS_PER_SEC);
	//cout << suffixArr.size() << endl;
	//cout << suffixArr.max_size() << endl;
	//cout << " \n";
	//cout << suffixArr.capacity() << endl;

	//fflush(stdout);
	int n = suffixArr.size();
	//cout << n << endl;
	//int suma = 0;
    cout << 3 << endl;
	// Se imprime el arreglo de sufijos
	//cout << "Suffix Array para "<< str << ": \n";
	//printArr(suffixArr, n);

	// Creación del LCP del string ingresado
	fflush(stdout);
	t2 = clock();
	vector<int>lcp = lcp_str(str, suffixArr);
    t3 = clock();
    double t23 = (double(t3-t2)/CLOCKS_PER_SEC);
	/*for (int ejem = 0; ejem < n; ejem++){
        cout << str.substr(suffixArr[ejem], 2) << endl;
	}*/
	cout << 4 << endl;
	fflush(stdout);

	// Se imprime el LCP
	//cout << "\nLCP Array para "<< str << ": \n";
	//printArr(lcp, n);

	// Se construye la suma de diferentes substrings que posee una palabra
	/*for(int sumatemp = 0; sumatemp < n; sumatemp++){
		if (sumatemp == 0){
			suma = n - suffixArr[0];
		} else{
			suma = suma + n - suffixArr[sumatemp] - lcp[sumatemp-1];
			//suma = guarda;
		}
		//cout << suma << " \n";
	}*/
	//cout << 5 << endl;
	//fflush(stdout);


    ofstream resultados, k_utilizados;
    resultados.open("resultados_sa_lcp.txt");
    resultados << "Construcción del SA: " << t12 << " segundos" << endl;
    resultados << " \n";
    resultados << "Construcción del LCP: " << t23 << " segundos" << endl;
    resultados << " \n";
    resultados << "Cantidad total analizada: " << cantidad_ss << " proteínas" << endl;
    resultados.close();
    //resultados << "--------------------------" << endl;
    //resultados << " \n";
    cout << 5 << endl;
    for (int k1 = 1; k1 < 4; k1++){
        if (k1 == 1){
            k_utilizados.open("resultados_k1.txt");
        }else if (k1 == 2){
            k_utilizados.open("resultados_k2.txt");
        }else if (k1 == 3){
            k_utilizados.open("resultados_k3.txt");
        }
        t4 = clock();
        subcadena arr[1];
        priority_queue<subcadena, vector<subcadena>, comparador> mypq;
        int activador = 0;
        int contador;
        int ds = 0;

        for (int temp = inicio; temp < n; temp++){
            //string sc = str.substr(suffixArr[temp], k1);
            //int scnumber = sc.size();
            //if (scnumber == k1 && sc.find_first_of("$BOUXZ")==std::string::npos){
                //cout << sc << " no hay pesos" << endl;
                if (lcp[temp] >= k1 && activador == 0){
                    string sc = str.substr(suffixArr[temp], k1);
                    if (sc.find_first_of("$BOUXZ")==std::string::npos){
                        arr[0].nombre = sc;
                        contador = 2;
                        activador = 1;
                    }
                }else if (lcp[temp] >= k1 && activador == 1){
                    string sc = str.substr(suffixArr[temp], k1);
                    if (sc.find_first_of("$BOUXZ")==std::string::npos){
                        contador = contador + 1;
                    }else{
                        arr[0].veces = contador;
                        if (contador > 0){
                            mypq.push(arr[0]);
                        }
                        ds = ds + 1;
                        activador = 0;
                    }
                }else if (lcp[temp] < k1 && activador == 0){
                    string sc = str.substr(suffixArr[temp], k1);
                    int scnumber = sc.size();
                    if (scnumber == k1){
                        if (sc.find_first_of("$BOUXZ")==std::string::npos){
                            ds = ds + 1;
                            arr[0].nombre = sc;
                            arr[0].veces = 1;
                            mypq.push(arr[0]);
                        }
                    //arr[0].nombre = sc;
                    //arr[0].veces = 1;
                    //mypq.push(arr[0]);
                    }
                }else if (lcp[temp] < k1 && activador == 1){
                    arr[0].veces = contador;
                    if (contador > 0){
                        mypq.push(arr[0]);
                    }
                    ds = ds + 1;
                    activador = 0;
                }
            //}
        }

        k_utilizados << "Diferentes substrings para " << k1 << " es " << ds << endl;
        k_utilizados << " \n";

        //for (int posicion = 0; posicion < 19; posicion++){
        while (!mypq.empty()){
            k_utilizados << mypq.top().nombre << " " << mypq.top().veces;
            mypq.pop();
            k_utilizados << endl;
        }
        t5 = clock();
        double t45 = (double(t5-t4)/CLOCKS_PER_SEC);
        k_utilizados << " \n";
        k_utilizados << "Tiempo utilizado: " << t45 << " segundos" << endl;
        k_utilizados << "--------------------------" << endl;
        k_utilizados << " \n";

        cout << "6- " << k1 << endl;
        fflush(stdout);

        k_utilizados.close();


    }

    for (int k1 = 4; k1 < 6; k1++){
        vector <string> combinaciones;
        set <string> auxiliar;
        if (k1 == 4){
            k_utilizados.open("resultados_k4.txt");
            //vector <string> combinaciones;
            ifstream bases("substringsbase4.txt");
            string residuos;

            getline(bases,residuos);

            while (bases){

                combinaciones.push_back(residuos);
                getline(bases,residuos);
            }

        }else if (k1 == 5){
            k_utilizados.open("resultados_k5.txt");
            //vector <string> combinaciones;
            ifstream bases("substringsbase5.txt");
            string residuos;

            getline(bases,residuos);

            while (bases){

                combinaciones.push_back(residuos);
                getline(bases,residuos);
            }
        }
        t4 = clock();
        subcadena arr[1];
        priority_queue<subcadena, vector<subcadena>, comparador> mypq;
        int activador = 0;
        int contador;
        int ds = 0;

        for (int temp = inicio; temp < n; temp++){
            //string sc = str.substr(suffixArr[temp], k1);
            //int scnumber = sc.size();
            //if (scnumber == k1 && sc.find_first_of("$BOUXZ")==std::string::npos){
                //cout << sc << " no hay pesos" << endl;
                if (lcp[temp] >= k1 && activador == 0){
                    string sc = str.substr(suffixArr[temp], k1);
                    if (sc.find_first_of("$BOUXZ")==std::string::npos){
                        arr[0].nombre = sc;
                        auxiliar.insert(sc);
                        //combinaciones.erase(remove(combinaciones.begin(), combinaciones.end(), sc), combinaciones.end());
                        contador = 2;
                        activador = 1;
                    }
                }else if (lcp[temp] >= k1 && activador == 1){
                    string sc = str.substr(suffixArr[temp], k1);
                    if (sc.find_first_of("$BOUXZ")==std::string::npos){
                        contador = contador + 1;
                    }else{
                        arr[0].veces = contador;
                        if (contador > 0){
                            mypq.push(arr[0]);
                        }
                        ds = ds + 1;
                        activador = 0;
                    }
                }else if (lcp[temp] < k1 && activador == 0){
                    string sc = str.substr(suffixArr[temp], k1);
                    int scnumber = sc.size();
                    if (scnumber == k1){
                        if (sc.find_first_of("$BOUXZ")==std::string::npos){
                            ds = ds + 1;
                            arr[0].nombre = sc;
                            auxiliar.insert(sc);
                            arr[0].veces = 1;
                            mypq.push(arr[0]);
                        }
                    //arr[0].nombre = sc;
                    //arr[0].veces = 1;
                    //mypq.push(arr[0]);
                    }
                }else if (lcp[temp] < k1 && activador == 1){
                    arr[0].veces = contador;
                    if (contador > 0){
                        mypq.push(arr[0]);
                    }
                    ds = ds + 1;
                    activador = 0;
                }
            //}
        }

        for (int k2 = 0; k2 < combinaciones.size(); k2++){
            int antes = auxiliar.size();
            auxiliar.insert(combinaciones[k2]);
            int despues = auxiliar.size();
            if (despues > antes){
                arr[0].nombre = combinaciones[k2];
                arr[0].veces = 0;
                mypq.push(arr[0]);
            }
        }

        k_utilizados << "Diferentes substrings para " << k1 << " es " << ds << endl;
        k_utilizados << " \n";

        //for (int posicion = 0; posicion < 19; posicion++){
        while (!mypq.empty()){
            k_utilizados << mypq.top().nombre << " " << mypq.top().veces;
            //auxiliar.insert(mypq.top().nombre);
            mypq.pop();
            k_utilizados << endl;
        }

        /*while (!mypq.empty()){
            k_utilizados << mypq.top().nombre << " " << mypq.top().veces;
            mypq.pop();
            k_utilizados << endl;
        }*/

        t5 = clock();
        double t45 = (double(t5-t4)/CLOCKS_PER_SEC);
        k_utilizados << " \n";
        k_utilizados << "Tiempo utilizado: " << t45 << " segundos" << endl;
        k_utilizados << "--------------------------" << endl;
        k_utilizados << " \n";

        cout << "6- " << k1 << endl;
        //cout << nada.size() << endl;
        fflush(stdout);

        k_utilizados.close();


    }


    k_utilizados.open("resultados_k6to50.txt");
    for (int k1 = 6; k1 < 51; k1++){

        t4 = clock();
        subcadena arr[1];
        priority_queue<subcadena, vector<subcadena>, comparador> mypq;
        int activador = 0;
        int contador;
        int ds = 0;

        for (int temp = inicio; temp < n; temp++){
            //string sc = str.substr(suffixArr[temp], k1);
            //int scnumber = sc.size();
            //if (scnumber == k1 && sc.find_first_of("$BOUXZ")==std::string::npos){
                //cout << sc << " no hay pesos" << endl;
                if (lcp[temp] >= k1 && activador == 0){
                    string sc = str.substr(suffixArr[temp], k1);
                    if (sc.find_first_of("$BOUXZ")==std::string::npos){
                        arr[0].nombre = sc;
                        contador = 2;
                        activador = 1;
                    }
                }else if (lcp[temp] >= k1 && activador == 1){
                    string sc = str.substr(suffixArr[temp], k1);
                    if (sc.find_first_of("$BOUXZ")==std::string::npos){
                        contador = contador + 1;
                    }else{
                        arr[0].veces = contador;
                        if (contador > 0){
                            mypq.push(arr[0]);
                        }
                        ds = ds + 1;
                        activador = 0;
                    }
                }else if (lcp[temp] < k1 && activador == 0){
                    string sc = str.substr(suffixArr[temp], k1);
                    int scnumber = sc.size();
                    if (scnumber == k1){
                        if (sc.find_first_of("$BOUXZ")==std::string::npos){
                            ds = ds + 1;
                        }
                    //arr[0].nombre = sc;
                    //arr[0].veces = 1;
                    //mypq.push(arr[0]);
                    }
                }else if (lcp[temp] < k1 && activador == 1){
                    arr[0].veces = contador;
                    if (contador > 0){
                        mypq.push(arr[0]);
                    }
                    ds = ds + 1;
                    activador = 0;
                }
            //}
        }

        k_utilizados << "Diferentes substrings para " << k1 << " es " << ds << endl;
        k_utilizados << " \n";

        for (int posicion = 0; posicion < 20; posicion++){
        //while (!mypq.empty()){
            k_utilizados << mypq.top().nombre << " " << mypq.top().veces;
            mypq.pop();
            k_utilizados << endl;
        }
        t5 = clock();
        double t45 = (double(t5-t4)/CLOCKS_PER_SEC);
        k_utilizados << " \n";
        k_utilizados << "Tiempo utilizado: " << t45 << " segundos" << endl;
        k_utilizados << "--------------------------" << endl;
        k_utilizados << " \n";

        cout << "6- " << k1 << endl;
        fflush(stdout);
    }
    //tfinal = clock();
    //double tiempo_total = (double(tfinal-tinicial)/CLOCKS_PER_SEC);
    //resultados << "El tiempo total ocupado es de " << tiempo_total << " segundos" << endl;
    k_utilizados.close();
    cout << "7" << endl;
    fflush(stdout);

    /*while (!mypq.empty())
    {
        cout<<mypq.top().nombre<<" "<<mypq.top().veces;
        mypq.pop();
        cout<<endl;
    }*/

	return 0;
}
