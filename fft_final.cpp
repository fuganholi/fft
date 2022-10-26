#include <iostream>
#include <string>
#include <iomanip>
#include <sstream>
#include <vector>
#include <fstream>
#include<math.h>
#include<complex>

using namespace std;

vector<complex<float>> fft(vector<complex<float>> a)
{
	//determina o numero de amostras
	int N = a.size();
	int K = N;

	//variavel interna para calculo
	complex<float> intSum;

	//vetor de saida
	vector<complex<float>> output;
	output.reserve(K);

	//calcula 
	for (int k=0; k<K; k++)
	{
		intSum = complex<float>(0,0);
		for(int n=0; n<N; n++)
		{
			float realPart = cos(((2*M_PI)/N) * k * n);
			float imagPart = sin(((2*M_PI)/N) * k * n);
			complex<float> w (realPart, -imagPart);
			intSum += a[n] * w;
		}
		output.push_back(intSum);
	}
	return output;

}

int main()
{

	//arquivo com os dados do acelerômetro
	string infile="1602245833-2715-NAO7856.txt";

	//arquivo de saida
    ofstream outfile ("output.txt");

	//valores lidos ao abrir o arquivo, por linha
	string X, Y, Z;

	//vetores com os todos valores do arquivo, ainda em string
	vector<string>xs;
	vector<string>ys;
	vector<string>zs;

	//vetores com os valores convertidos em float
	vector<complex<float>>xf;
	vector<complex<float>>yf;
	vector<complex<float>>zf;

	//vetores que receberao o resultado do calculo
	vector<float>x;
	vector<float>y;
	vector<float>z;
	vector<float>freq;

    //vetores com as strings de saída, para o arquivo output.txt
	vector<string>xOut;
	vector<string>yOut;
	vector<string>zOut;
	vector<string>freqOut;

	//contador de linhas de dados
	int i = 0;

	//abre o arquivo
	ifstream coeff(infile);

	//se abrir o arquivo
	if (coeff.is_open())
	{
		//ignora o cabecalho
		string line;
		getline(coeff, line);

		//enquanto nao chega ao fim do arquivo
		while (!coeff.eof())
		{
			//passa cada valor que fica entre virgulas para o vetor respectivo (x,y,z)
			getline(coeff, X, ',');
			xs.push_back(X);
			getline(coeff, Y, ',');
			ys.push_back(Y);
			getline(coeff, Z, '\n');
			zs.push_back(Z);
			i += 1;
		}

		//fecha o arquivo
		coeff.close();
	}

	//mensagem de erro de leitura do arquivo
	else cout << "Unable to open file";

	//converte os valores de entrada de string para float
	for(int k=0;k<(i-1);k++)
	{
		xf.push_back(stof(xs[k]));
		yf.push_back(stof(ys[k]));
		zf.push_back(stof(zs[k]));
	}

	//define a quantidade de amostras
	float N = i-1;

	//vetores das FFTs
    vector<complex<float>> Fx = fft(xf);
    vector<complex<float>> Fy = fft(yf);
    vector<complex<float>> Fz = fft(zf);

	//passa os valores da magnitude para o vetor dos resultados
    for (int q = 0; q < N; ++q)
    {
        x.push_back(sqrt((Fx[q].real()*Fx[q].real())+(Fx[q].imag()*Fx[q].imag())));
        y.push_back(sqrt((Fy[q].real()*Fy[q].real())+(Fy[q].imag()*Fy[q].imag())));
        z.push_back(sqrt((Fz[q].real()*Fz[q].real())+(Fz[q].imag()*Fz[q].imag())));

        x[q] = (2 * std::abs(x[q]/N));
        y[q] = (2 * std::abs(y[q]/N));
        z[q] = (2 * std::abs(z[q]/N));
    }

	//variaveis para gerar o vetor das frequencias
    float period = 2.715;
    float frequency = 1/period;
    float value = 0;
    int divisor = N/2;

	//passa os valores para o vetor
    for(int k = 0; k < divisor; k++)
    {
        freq.push_back(value);
        value += frequency;
    }

    for(int k = divisor+1; k < N; k++)
    {
        freq.push_back(-value);
        value -= frequency;
    }
        
    //converte os valores calculados para string
	for(int l=0;l<i-1;l++)
	{
		xOut.push_back(to_string(x[l]));
		yOut.push_back(to_string(y[l]));
		zOut.push_back(to_string(z[l]));
        freqOut.push_back(to_string(freq[l]));
	}

	//cria a string para o arquivo de saida
	string out = "";

	//acrescenta as strings dos vetores a string de saida
	for (int j = 0; j < N; j++) {
        if(freq[j] > 0)
        {
		    out += (xOut[j] + "," + yOut[j] + "," + zOut[j] + "," + freqOut[j] +  "\n");
        }
	}

	//escreve a variavel de saida no arquivo output.txt
	outfile << out << endl;

	//fecha o arquivo de saida
	outfile.close();

    return 0;

}
