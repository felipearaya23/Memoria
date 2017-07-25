#include <bits/stdc++.h>
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <algorithm>
#include <set>
#include <cstdio>
#include <ctime>
#include <iomanip>
#include <istream>
#include <cstdlib>
# include <limits.h>
# include <string.h>
# include <stdio.h>
# include <math.h>

# define NO_OF_CHARS 256

using namespace std;

void preKMP(string pattern, int f[]){

    int m = pattern.length(), k;
    f[0] = -1;

    for (int i = 1; i < m; i++){
        k = f[i - 1];
        while (k >= 0){
            if (pattern[k] == pattern[i - 1])
                break;
            else
                k = f[k];
        }

        f[i] = k + 1;
    }
}

//check whether target string contains pattern

bool KMP(string pattern, string target){

    int m = pattern.length();
    int n = target.length();
    int f[m];

    preKMP(pattern, f);
    int i = 0;
    int k = 0;

    while (i < n)    {
        if (k == -1){
            i++;
            k = 0;
        }
        else if (target[i] == pattern[k]){
            i++;
            k++;
            if (k == m)
                return 1;
        }
        else
            k = f[k];
    }
    return 0;
}




// A utility function to get maximum of two integers

int max(int a, int b){
    return (a > b) ? a : b;
}

// The preprocessing function for Boyer Moore's bad character heuristic

void badCharHeuristic(string str, int size, int* badchar){

    int i;

    // Initialize all occurrences as -1
    for (i = 0; i < 256; i++)
        badchar[i] = -1;

    // Fill the actual value of last occurrence of a character

    for (i = 0; i < size; i++)
        badchar[(int)str[i]] = i;
}

int search_bm(string txt, string pat){

    vector<int> retVal;
    int m = pat.length();
    int n = txt.length();

    int* badchar = new int[256];
    badCharHeuristic(pat, m, badchar);

    int s = 0; // s is shift of the pattern with respect to text
    while (s <= (n - m)){
        int j = m - 1;

        while (j >= 0 && pat[j] == txt[s + j])
            j--;

        if (j < 0){
            //printf("\npattern occurs at shift = %d", s);
            retVal.push_back(s);
            break;
            s += (s + m < n) ? m - badchar[(int)txt[s + m]] : 1;
        }
        else
            s += max(1, j - badchar[(int)txt[s + j]]);
    }
    delete[] badchar;
    int cantidad = retVal.size();
    return cantidad;
}

// Horspool algorithm
template <typename RAIter>
vector<size_t> hStrMatch(RAIter first, RAIter last,
                         const string & pattern)
{
    vector<size_t> result;        // Vector to hold return value
    size_t size = last-first;     // Size of text ("n"0
    size_t len = pattern.size();  // Length of pattern ("m")

    // Handle trivial case: zero-length pattern always matches
    if (len == 0)
    {
        result.push_back(0);
        return result;
    }

    // Make bad-symbol shift table
    vector<size_t> badsymbol(256, len);  // 256 possible char values
    for (size_t i = 0; i != len-1; ++i)
    {
        badsymbol[pattern[i]] = len - 1 - i;
    }

    // Do the search from the beginning of the text
    size_t loc = 0;  // loc is current search location in text
    while (loc+len <= size)
    {
        // Check match against pattern right-to-left
        size_t k = len;
        while (true)
        {
            --k;
            if (first[loc+k] != pattern[k])
                break;
            if (k == 0)
            {  // Found! Return location
                result.push_back(loc);
                return result;
            }
        }

        // Did not find yet; advance loc using bad-symbol shift table
        char c = first[loc+len-1];
        loc += badsymbol[c];
    }

    // Return not-found result
    return result;
}

/* Driver program to test above funtion */

int main(){

    ifstream fin("uniprot_sprot.fasta");
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
            outputfile << endl;
            cantidad_ss = cantidad_ss +1;
        }
        else{
            outputfile << line;
        }
    getline(fin, line);
    }

    outputfile.close();


    //metodo generador de substrings desde el 1
    /*vector <string> vocales;
    vector <string> substrings_1argo_1;
    ifstream bases("substringsbase.txt");
    string vocal;

    getline(bases,vocal);

    while (bases){

        vocales.push_back(vocal);
        getline(bases,vocal);
    }
    unsigned t0, t1;

    t0 = clock();


    int largo_abecedario = vocales.size();
    for (int i = 0; i<largo_abecedario; i++){
        ifstream proteinas("substrings.txt");
        string alternativa;
        while (proteinas){
            getline(proteinas, alternativa);
            //vector<size_t> hresult = hStrMatch(alternativa.begin(), alternativa.end(), vocales[i]);
            //if (!hresult.empty()){
            if (search_bm(alternativa, vocales[i])){
            //if (KMP(vocales[i], alternativa)){
                substrings_1argo_1.push_back(vocales[i]);
                break;
            }
        }
    }

    t1 = clock();
    double time = (double(t1-t0)/CLOCKS_PER_SEC);
    cout << "para 1 ds es " << substrings_1argo_1.size() << " en " << time << " segundos" << endl;

    vector <string> combinaciones = substrings_1argo_1;
    vector <string> auxiliar;
    vector <string> vacio;
    //vocales = vacio;
    int largo_ss1 = substrings_1argo_1.size();
    ifstream proteinas2("substrings.txt");
    string nueva_alt;
    //Ahora sedetermina el k a buscar
    for(int k = 2; k < 51; k++){
        t0 = clock();
        int cuenta = 0;
        int largo_actual = combinaciones.size();
        for (int a = 0 ; a < largo_actual; a++){
            for (int b = 0; b < largo_ss1; b++){
                auxiliar.push_back(combinaciones[a] + substrings_1argo_1[b]);
            }
        }
        combinaciones = vacio;
        int largo_auxiliar = auxiliar.size();
        for (int i2 = 0; i2<largo_auxiliar; i2++){
            while (getline(proteinas2, nueva_alt)){
                //vector<size_t> hresult = hStrMatch(nueva_alt.begin(), nueva_alt.end(), auxiliar[i2]);
                //if (!hresult.empty()){
                if (search_bm(nueva_alt, auxiliar[i2])){
                //if (KMP(auxiliar[i2], nueva_alt)){
                    combinaciones.push_back(auxiliar[i2]);
                    cuenta = cuenta + 1;
                    cout << cuenta << endl;
                    break;
                }
            }
            proteinas2.clear();
            proteinas2.seekg(0, std::ios_base::beg);

        }
        auxiliar = vacio;
        t1 = clock();
        double time = (double(t1-t0)/CLOCKS_PER_SEC);
        cout << "para " << k << " ds es " << combinaciones.size() << " en " << time << " segundos" << endl;
        //continuar
    }



    //metodo de lectura de archivo
    //ofstream outputfile;
    //outputfile.open("substrings.txt");
    */

    vector <string> prohibido {"X", "B", "O", "U", "Z"};
    int cant_prohibidos = prohibido.size();
    unsigned t0, t1;

    for(int k = 1; k < 51; k++){
        t0 = clock();
        //int cuenta = 0;
        set <string> unicos;
        ifstream proteinas("substrings.txt");
        //ofstream sec_diferentes;
        //sec_diferentes.open("diferentes_subs.txt");
        //sec_diferentes.close();
        string secuencia;
        for (int j = 0; j < cantidad_ss; j++){
            getline(proteinas, secuencia);
            int tam_sec = secuencia.size();
            int maximo = tam_sec-k+1;
            for (int i = 0; i < maximo; i++){
                int conta = 0;
                for (int i4 = 0; i4 < cant_prohibidos; i4++){
                    conta = secuencia.substr(i,k).find(prohibido[i4]);
                    if (conta > -1){
                        break;
                    }
                }
                if (conta == -1){
                    int j = unicos.size();
                    unicos.insert(secuencia.substr(i,k));
                    int j_new = unicos.size();
                    if (j_new > j){
                        cout << j_new << endl;
                    }
                }

            }
            if (unicos.size() == pow(20, k)){
                break;
            }
        }
        t1 = clock();
        double time = (double(t1-t0)/CLOCKS_PER_SEC);
        cout << "para " << k << " ds es " << unicos.size() << " en " << time << " segundos" << endl;

    }
    /*int k=30;
    ifstream proteinas("substrings.txt");
    ofstream sec_diferentes;
    sec_diferentes.open("diferentes_subs.txt");
    sec_diferentes.close();
    string secuencia;
    for (int j = 0; j < cantidad_ss; j++){
        getline(proteinas, secuencia);
        int tam_sec = secuencia.size();
        int maximo = tam_sec-k+1;
        for (int i = 0; i < maximo; i++){
            ifstream diferentes("diferentes_subs.txt");
            string palabra;
            int encontrado = 0;
            while (diferentes){
                getline(diferentes, palabra);
                if (KMP(secuencia.substr(i,k), palabra)){
                    encontrado = 1;
                    break;
                }
            }
            if (encontrado == 0){
                ofstream agregar;
                agregar.open("diferentes_subs.txt", std::ios_base::app | std::ios_base::out);
                agregar << secuencia.substr(i,k) << endl;
            }
        }
    }
    ifstream cuenta_final("diferentes_subs.txt");
    string substring_diferente;
    int cant_diff_subs = 0;
    getline(cuenta_final, substring_diferente);
    while (cuenta_final){
        cant_diff_subs = cant_diff_subs + 1;
        getline(cuenta_final, substring_diferente);
    }
    cout << cant_diff_subs << endl;*/

   //outputfile.close();

    /*string tar = "san and linux training";
    string pat = "lin";
    if (KMP(pat, tar))
        cout<<"'"<<pat<<"' found in string '"<<tar<<"'"<<endl;
    else
        cout<<"'"<<pat<<"' not found in string '"<<tar<<"'"<<endl;

    pat = "sanfoundry";
    if (KMP(pat, tar))
        cout<<"'"<<pat<<"' found in string '"<<tar<<"'"<<endl;
    else
        cout<<"'"<<pat<<"' not found in string '"<<tar<<"'"<<endl;*/

    //char txt[] = "mississippi";
    //char pat[] = "p";
    //string txt = "mississippi";
    //string pat = "p";
    //int valor = search_bm(txt, pat);
    //cout << valor << endl;
    //cout << cantidad_ss << endl;
    //cout << auxiliar.size() << endl;


    return 0;

}

