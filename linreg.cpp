/* Linear Regression Model 
	by Ramon Sandoval, 2009
	
	When entering a text file into the program, use a tab delimited file, with x first, 
	sigma(x) or W(x) second, y third, and sigma(y) or W(y) last. Do not use a non-number,
	have not accounted for.
*/

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <string>
#include <cstdlib>
#include <ctime>
#define PI 3.141592

long ptnumber;
long double precision = 0.00001;
std::vector<long double> x;
std::vector<long double> sigx;
std::vector<long double> y;
std::vector<long double> sigy;
std::vector<long double> W;
std::vector<long double> Wx;
std::vector<long double> Wy;
std::vector<long double> U;
std::vector<long double> V;
long double bracketx;
long double brackety;
long double mydelta;
long double mygamma;
long double mybeta;
long double myalpha;
int kindoferror = 2;

//------------------------------------

bool Readfile(std::vector<long double> &x , std::vector<long double> &sigx, std::vector<long double> &y, std::vector<long double> &sigy) {
	
	std::cout << "File name: ";
	
	std::string textfile;
	std::cin >> textfile;

	std::string key (".txt");
	size_t foundinstr;

	foundinstr=textfile.rfind(key);
	if (foundinstr == std::string::npos)
		textfile += ".txt";
	
	
	char str[255];
	strcpy( str, textfile.c_str() );
	
	std::ifstream myfiletest;
	std::ifstream myfile;
	
	myfiletest.open(str);
	
	if (!myfiletest.is_open()) {
		std::cout << "Unable to open file\n" << std::endl;
		myfile.close();
		return false;
	}
	
	char ch;
	long count = 0;
	while (!myfiletest.eof()) {
        myfiletest.get(ch);
		if (ch == '	')
			count++;
		else if ( ch == '\n' ) {
			if (((count % 3) != 0) || (count == 0)) {
				//std::cout << "Number of tabs: " << count << std::endl;
				std::cout << "Unable to read file properly, is it formatted correctly?\n" << std::endl;
				myfiletest.close();
				return false;
			}
		}
    }
	
	myfiletest.close();

	myfile.open(str);
	
	std::vector<std::string> xstr (1, "");
	std::vector<std::string> sigxstr (1, "");
	std::vector<std::string> ystr (1, "");
	std::vector<std::string> sigystr (1, "");
	
	long i = 0;
	
	while (!myfile.eof()) {
		
		xstr.resize(i+1);
		sigxstr.resize(i+1);
		ystr.resize(i+1);
		sigystr.resize(i+1);
	
		std::getline(myfile, xstr[i], '	');
		std::getline(myfile, sigxstr[i], '	');   
		std::getline(myfile, ystr[i], '	');
		std::getline(myfile, sigystr[i]);
		
		x.push_back((long double)atof(xstr[i].c_str()));
		
		sigx.push_back((long double)atof(sigxstr[i].c_str()));
			if (sigx[i] == 0) {
				sigx[i] = 0.00000000001*x[i];
			}
			
		y.push_back((long double)atof(ystr[i].c_str()));
		
		sigy.push_back((long double)atof(sigystr[i].c_str()));
			if (sigy[i] == 0) {
				sigy[i] = 0.00000000001*y[i];
			}
		
		if (myfile.eof()) {
			x.resize(i);
			sigx.resize(i);
			y.resize(i);
			sigy.resize(i);
			break;
		}
	
		i++;
	}
	
	myfile.close();
	
	return true;
}

//------------------------------------

long double sum(std::vector<long double> arr) {
		long i;
		long double ans = 0;
		
		for( i = 0; i < arr.size(); i++ ) {
			ans += arr[i];
		}
		
		return ans;
}

//------------------------------------

std::vector<long double> operator + (std::vector<long double> a, std::vector<long double> b ) {
	long i;
	std::vector<long double> ans (a.size(),0);
	for ( i = 0; i < a.size(); i++) {
		ans[i] = a[i] + b[i];
	}
	
	return ans;
}

//------------------------------------

std::vector<long double> operator - (std::vector<long double> a, std::vector<long double> b ) {
	long i;
	std::vector<long double> ans (a.size(),0);
	for ( i = 0; i < a.size(); i++) {
		ans[i] = a[i] - b[i];
	}
	
	return ans;
}

//------------------------------------

std::vector<long double> operator * (std::vector<long double> a, std::vector<long double> b ) {
	long i;
	std::vector<long double> ans (a.size(),0);
	for ( i = 0; i < a.size(); i++) {
		ans[i] = a[i] * b[i];
	}
	
	return ans;
}

//------------------------------------

std::vector<long double> operator / (std::vector<long double> a, std::vector<long double> b ) {
	long i;
	std::vector<long double> ans (a.size(),0);
	for ( i = 0; i < a.size(); i++) {
		ans[i] = a[i] / b[i];
	}
	
	return ans;
}

//------------------------------------

std::vector<long double> operator * (std::vector<long double> a, long double b ) {
	long i;
	std::vector<long double> ans (a.size(),0);
	for ( i = 0; i < a.size(); i++) {
		ans[i] = a[i] * b;
	}
	
	return ans;
}

//------------------------------------

std::vector<long double> operator / (std::vector<long double> a, long double b ) {
	long i;
	std::vector<long double> ans (a.size(),0);
	for ( i = 0; i < a.size(); i++) {
		ans[i] = a[i] / b;
	}
	
	return ans;
}

//------------------------------------

std::vector<long double> operator * (long double a, std::vector<long double>  b ) {
	long i;
	std::vector<long double> ans (b.size(),0);
	for ( i = 0; i < b.size(); i++) {
		ans[i] = a * b[i];
	}
	
	return ans;
}

//------------------------------------

std::vector<long double> operator / (long double a, std::vector<long double>  b ) {
	long i;
	std::vector<long double> ans (b.size(),0);
	for ( i = 0; i < b.size(); i++) {
		ans[i] = a / b[i];
	}
	
	return ans;
}

//------------------------------------

std::vector<long double> operator + (long double a, std::vector<long double>  b ) {
	long i;
	std::vector<long double> ans (b.size(),0);
	for ( i = 0; i < b.size(); i++) {
		ans[i] = a + b[i];
	}
	
	return ans;
}

//------------------------------------

std::vector<long double> operator - (long double a, std::vector<long double>  b ) {
	long i;
	std::vector<long double> ans (b.size(),0);
	for ( i = 0; i < b.size(); i++) {
		ans[i] = a - b[i];
	}
	
	return ans;
}

//------------------------------------

std::vector<long double> operator + (std::vector<long double> a, long double b ) {
	long i;
	std::vector<long double> ans (a.size(),0);
	for ( i = 0; i < a.size(); i++) {
		ans[i] = a[i] + b;
	}
	
	return ans;
}

//------------------------------------

std::vector<long double> operator - (std::vector<long double> a, long double  b ) {
	long i;
	std::vector<long double> ans (a.size(),0);
	for ( i = 0; i < a.size(); i++) {
		ans[i] = a[i] - b;
	}
	
	return ans;
}

//------------------------------------

std::ofstream &operator << (std::ofstream &myfile, std::vector<long double> a) {
	long i;
	for ( i = 0; i < a.size(); i++) {
		myfile << a[i] << std::endl;
	}
	
	return myfile;
}

//------------------------------------

void writefile ( std::string str, std::vector<long double> a) {
	std::ofstream myfile;
	
	char textfile[255];
	strcpy( textfile , str.c_str() );
	myfile.open(textfile);
	
	myfile << a;
	
	myfile.close();
}

//------------------------------------

void writefile ( std::string str, double a) {
	std::ofstream myfile;
	
	char textfile[255];
	strcpy( textfile , str.c_str() );
	myfile.open(textfile);
	
	myfile << a;
	
	myfile.close();
}

//------------------------------------

bool InitVectors() {
	if (!Readfile(x, sigx, y, sigy)) {
		return false;
	}
	ptnumber = x.size();
	W.resize(ptnumber);
	Wx.resize(ptnumber);
	Wy.resize(ptnumber);
	U.resize(ptnumber);
	V.resize(ptnumber);
	
	return true;
}

//------------------------------------

void FindW(long double m) {
	long i;

	if (kindoferror == 1) {
		Wx = sigx;
		Wy = sigy;
	}
	
	//Use if reporting only weights, not with errors
	
	else {
		Wx = 1/(sigx*sigx);
		Wy = 1/(sigy*sigy);
	}

	W = (Wx*Wy)/((m*m*Wy)+Wx);
}

//------------------------------------

void Findbracket() {	
	bracketx = sum(W*x/sum(W));
	brackety = sum(W*y/sum(W));
	U = x - bracketx;
	V = y - brackety;
}

//------------------------------------

long double Findmydelta() {
	mydelta = sum(W*W*U*U/Wx);
}

long double Findmygamma() {
	mygamma = (-sum(W*U*V))/mydelta;
}

long double Findmybeta() {
	mybeta = (sum(W*W*V*V/Wx)-sum(W*U*U))/(3*mydelta);
}

long double Findmyalpha() {
	myalpha = (2*sum(W*W*U*V/Wx))/(3*mydelta);
}

//------------------------------------

long double RefindValues (long double m) {
	FindW(m);
	Findbracket();
	Findmydelta();
	Findmygamma();
	Findmybeta();
	Findmyalpha();
}

//------------------------------------

long double FindS(long double m) {
	long double S;

	RefindValues(m);
	S = sum(W*(V-(m*U))*(V-(m*U)));
	
	return S;
}

//------------------------------------

//Find vector element that holds the smallest S value
long Findsmallest(std::vector<long double> arr) {
	
	long i,j;
	bool found = false;
	
	for(i=0; i < arr.size(); i++) {
		for( j = 0; j < arr.size(); j++) {
			if (arr[i] > arr[j]) {
				found = false;
				break;
			}
			
			else {
				found = true;
			}
		}
		
		if ( found == true ) {
			return i;
		}
	}
}

//------------------------------------

long Findlargest(std::vector<long double> arr) {
	
	long i,j;
	bool found = false;
	
	for(i=0; i < arr.size(); i++) {
		for( j = 0; j < arr.size(); j++) {
			if (arr[i] < arr[j]) {
				found = false;
				break;
			}
			
			else {
				found = true;
			}
		}
		
		if ( found == true ) {
			return i;
		}
	}
}

//------------------------------------

long double FindJJ(long double m) {
	long double JJ;
	
	//writefile ("JJ.txt", -((2*m)/sum(W))*(W*W*U)/Wx);
	JJ = -((2*m)/sum(W))*sum((W*W*U)/Wx);
	return JJ;
}

//------------------------------------

long double FindHH(long double m) {
	long double HH;
	
	//writefile ("HH.txt", -((2*m)/sum(W))*((W*W*V)/Wx));
	HH = -((2*m)/sum(W))*sum((W*W*V)/Wx);
	return HH;
}

//------------------------------------

long double FindHHforAA(long double m) {
	long double HH;
	
	//writefile ("HH.txt", -((2*m)/sum(W))*((W*W*V)/Wx));
	HH = -((2)/sum(W))*sum((W*W*V)/Wx);
	return HH;
}

//------------------------------------

long double FindGGj(long j) {
	long double GGj;
	long i;
	long double krondel;
	
	std::vector<long double> arrayforkrondel (ptnumber,0);
	
	for ( i = 0; i < ptnumber; i++) {
		if ( i == j ) {
			krondel = 1;
		}
		else {
			krondel = 0;
		}
		
		arrayforkrondel[i] = (W[i]*W[i]*U[i]/Wx[i])*(krondel - W[j]/sum(W));
	}
	
	GGj = sum(arrayforkrondel);
	//writefile ("GGj.txt", arrayforkrondel);
	return GGj;
}

//------------------------------------

long double FindFFj(long j) {
	long double FFj;
	long i;
	long double krondel;
	
	std::vector<long double> arrayforkrondel (ptnumber,0);
	
	for ( i = 0; i < ptnumber; i++) {
		if ( i == j ) {
			krondel = 1;
		}
		else {
			krondel = 0;
		}
		
		arrayforkrondel[i] = (W[i]*W[i]*V[i]/Wy[i])*(krondel - W[j]/sum(W));
	}
	
	FFj = sum(arrayforkrondel);
	//writefile ("FFj.txt", arrayforkrondel);
	return FFj;
}

//------------------------------------

long double FindEEj(long j) {
	long double EEj;
	long i;
	long double krondel;
	
	std::vector<long double> arrayforkrondel (ptnumber,0);
	
	for ( i = 0; i < ptnumber; i++) {
		if ( i == j ) {
			krondel = 1;
		}
		else {
			krondel = 0;
		}
		
		arrayforkrondel[i] = (W[i]*W[i]*U[i]/Wy[i])*(krondel - W[j]/sum(W));
	}
	
	EEj = 2*sum(arrayforkrondel);
	//writefile ("EEj.txt", arrayforkrondel);
	return EEj;
}

//------------------------------------

long double FindDDj(long j) {
	long double DDj;
	long i;
	long double krondel;
	
	std::vector<long double> arrayforkrondel (ptnumber,0);
	
	for ( i = 0; i < ptnumber; i++) {
		if ( i == j ) {
			krondel = 1;
		}
		else {
			krondel = 0;
		}
		
		arrayforkrondel[i] = (W[i]*W[i]*V[i]/Wx[i])*(krondel - W[j]/sum(W));
	}
	
	DDj = sum(arrayforkrondel);
	//writefile ("DDj.txt", arrayforkrondel);
	return DDj;
}

//------------------------------------

long double FindCC(long double m) {
	long double CC;
	
	//writefile ("CC.txt", (-1)*(W*W/Wy)*(4*m*((W*U*V)/Wx) + V*FindJJ(m) + U*FindHH(m)));
	CC = -sum(((W*W)/Wy)*((4*m*((W*U*V)/Wx)) + V*FindJJ(m) + U*FindHH(m)));
	return CC;
}

//------------------------------------

long double FindBB(long double m) {
	long double BB;
	
	//writefile ("BB.txt", (-1)*(W*W)*(4*m*(W/Wx)*(U*U/Wy - V*V/Wx) - 2*V*FindHH(m)/Wx + 2*U*FindJJ(m)/Wy));
	BB = -sum((W*W)*((4*m*(W/Wx)*((U*U)/Wy - (V*V)/Wx)) - (2*V*FindHH(m))/Wx + (2*U*FindJJ(m))/Wy));
	return BB;
}

//------------------------------------

long double FindAA(long double m) {
	long double AA;
	
	//writefile ("AA.txt", 4*m*(W*W*W*U*V/(W*W) - sum(W)*FindHH(m)*FindJJ(m)));
	AA = (4*m)*sum((W*W*W*U*V)/(Wx*Wx)) - (sum(W)*FindHHforAA(m)*FindJJ(m));
	return AA;
}

//------------------------------------

long double Findpartialmxj(long double m, long j) {
	long double partialmxj;
	long double A, B;
	
	A = sum((W*W*U*V)/Wx);
	B = sum((W*W)*(((U*U)/Wy)-((V*V)/Wx)));
	
	//writefile ("A1.txt", (W*W*U*V)/Wx);
	//writefile ("B1.txt", (W*W)*(((U*U)/Wy)-((V*V)/Wx)));
	
	partialmxj = -((m*m*FindDDj(j) + m*FindEEj(j) - FindFFj(j))/(2*m*A + B - FindAA(m)*m*m + FindBB(m)*m - FindCC(m)));
	
	return partialmxj;
}

//------------------------------------

long double Findpartialmyj(long double m, long j) {
	long double partialmyj;
	long double A, B;
	
	A = sum((W*W*U*V)/Wx);
	B = sum((W*W)*(((U*U)/Wy)-((V*V)/Wx)));
	
	//writefile ("A2.txt", (W*W*U*V)/Wx);
	//writefile ("B2.txt", (W*W)*(((U*U)/Wy)-((V*V)/Wx)));
	
	partialmyj = -((m*m*FindGGj(j) - 2*m*FindDDj(j) - (FindEEj(j)/2))/(2*m*A + B - FindAA(m)*m*m + FindBB(m)*m - FindCC(m)));
	
	return partialmyj;
}

//------------------------------------

long double Findpartialcxj(long double m, long j) {
	long double partialcxj;

	partialcxj = (FindHH(m) - m*FindJJ(m) - bracketx)*Findpartialmxj(m,j) - m*W[j]/sum(W);
	
	return partialcxj;
}

//------------------------------------

long double Findpartialcyj(long double m, long j) {
	long double partialcyj;

	partialcyj = (FindHH(m) - m*FindJJ(m) - bracketx)*Findpartialmyj(m,j) + W[j]/sum(W);
	
	return partialcyj;
}

//------------------------------------

long double Finderror(long double m) {
	long double summed;
	long double other;
	long i;
	long double variance;
	long double deviation;
	
	std::vector<long double> partialmx;
	std::vector<long double> partialmy;

	RefindValues(m);	
	for ( i = 0; i < ptnumber; i++) {
		partialmx.push_back(Findpartialmxj(m,i));
		partialmy.push_back(Findpartialmyj(m,i));
	}
	
	/*
	//////////////////////////////////////////////////////////////////////////////////////////////////
	writefile("x.txt", x);
	writefile("sigx.txt", sigx);
	writefile("y.txt", y);
	writefile("sigy.txt", sigy);
	writefile("bracketx.txt", bracketx);
	writefile("brackety.txt", brackety);
	writefile("U.txt", U);
	writefile("V.txt", V);
	writefile("W.txt", W);
	writefile("Wy.txt", Wy);
	writefile("Wx.txt", Wx);
	writefile("partialmx.txt", partialmx);
	writefile("partialmy.txt", partialmy);
	////////////////////////////////////////////////////////////////////////////////////////////////////
	*/
	
	summed = sum((partialmx*partialmx/Wx) + (partialmy*partialmy/Wy));
	other = FindS(m)/(ptnumber - 2);
	variance = other*summed;
	deviation = pow(variance,0.5);
	
	return deviation;
}

//------------------------------------

long double Finderrorinc(long double m) {
	long double summed;
	long double other;
	long i;
	long double variance;
	long double deviation;
	
	std::vector<long double> partialcx;
	std::vector<long double> partialcy;

	RefindValues(m);
	
	for ( i = 0; i < ptnumber; i++) {
		partialcx.push_back(Findpartialcxj(m,i));
		partialcy.push_back(Findpartialcyj(m,i));
	}
	
	summed = sum((partialcx*partialcx/Wx) + (partialcy*partialcy/Wy));
	other = FindS(m)/(ptnumber - 2);
	variance = other*summed;
	deviation = pow(variance,0.5);

	return deviation;
}

//------------------------------------

long double tanmethod(long double &start, long double &end, long double step) {
	std::vector<long double> sumofsquares;
	std::vector<long double> m;
	long double i, begin;
	long smallest;
	begin = start;
	
	for (i = start; i < end; i +=step) {
		if (isnan(tan(i)) || isinf(tan(i)))
			continue;
		sumofsquares.push_back(FindS(tan(i)));
		m.push_back(tan(i));
	}

	smallest = Findsmallest(sumofsquares);
	start = begin + (smallest*step) - step;
	end = begin + (smallest*step) + step;
	
	return m[smallest];
}

//-------------------------------------------------------------------------------------------------------------------------------------------------------


int main () {

	long double mroot;
	long double dmroot = (PI/180);
	long double startsearch = ((-1)*PI/180);
	long double endsearch (181*PI/180);
	long double c;
	long double S;
	
	srand((unsigned)time(0));
	
	std::cout <<  "\n" <<  "\n";
	std::cout << "Least Squares Linear Regression Slope/Intercept Solver" << std::endl;
	std::cout << "by Ramon Sandoval, 2009" << std::endl << std::endl << std::endl << std::endl;
	
	std::cout << "Enter name of text file with data. \n** N.B. File must be tab delimited (no superfluous tabs), in this order **\n";
	std::cout << "** | x | sigma(x) or W(x) | y | sigma(y) or W(y) | **\n\n";
	
	//-----------------------------------------------
	
	while(!InitVectors()) {};
	std::cout << std::endl;
	
	//-----------------------------------------------
	
	std::cout << "What kind of reported error? In Weight or in sigma(x)? \n(type 1 for Weight, type 2 for sigma(x)): ";
	bool getout = false;
	std::string precstr5;
	char g[255];
	char gindex = 0;
	
	while (getout == false) {
		std::cin >> precstr5;
		strcpy( g, precstr5.c_str() );
		
		gindex = 0;
		
		while (g[gindex] == '\n') {
			gindex++;
		}
		
		if (g[gindex] == '1') {
			std::cout << "Weight it is then" << std::endl << std::endl;
			kindoferror = 1;
			getout = true;
		}
		else if (g[gindex] == '2') {
			std::cout << "sigma(x) it is then" << std::endl << std::endl;
			kindoferror = 2;
			getout = true; 
		}
		else {
			getout = false;
			std::cout << "Type properly" << std::endl << std::endl;
		}
	}
	
	//-----------------------------------------------
	
	std::cout << "Enter precision? (Y or N) ";
	
	std::string precstr;
	std::cin >> precstr;
	
	char a[255];
	strcpy( a, precstr.c_str() );
	
	char aindex = 0;
	
	while (a[aindex] == '\n') {
		aindex++;
	}
	
	if ( (a[aindex] == 'Y') ||  (a[aindex] == 'y') ) {
		std::cout << "Enter precision value: ";
		std::cin >> precision;
		std::cout << std::endl;
	}
	else if ( (a[aindex] == 'N') ||  (a[aindex] == 'n') ) {
		std::cout << std::endl;
	}
	else {
		std::cout << "I'm going to assume No then" << std::endl << std::endl;
	}
	
	//-----------------------------------------------
	
	std::cout << "Finding slope, please wait..." << std::endl << std::endl;
	
	long double difference = 2*precision;
	long stoploop = 0;
	while ( (difference > precision) && (stoploop < 20) ) {
		mroot = tanmethod(startsearch,endsearch,dmroot);
		difference = fabs(startsearch - endsearch);
		dmroot = difference/10;
		stoploop++;
	}
	
	if (difference > precision) {
		std::cout << "Too precise, did not want to loop that much" << "\n" << "\n";
	}
	
	RefindValues(mroot);
	c = brackety - mroot*bracketx;
	long double yorkmsquared = (FindS(mroot)/(ptnumber-2))/sum(W*U*U);
	long double yorkm = pow(yorkmsquared, 0.5);
	long double yorkcsquared = yorkm*yorkm*sum(W*x*x/sum(W));
	long double yorkc = pow(yorkcsquared, 0.5);

	//-----------------------------------------------
	
	std::cout << "Least Squares Regression slope: " << mroot << std::endl;
	std::cout << "Least Squares Regression intercept: " << c << std::endl;
	std::cout << "Least Squares Regression Sum of Weighted Squared Residuals: " << FindS(mroot) << std::endl;
	std::cout << "Least Squares Regression calculated slope error: " << Finderror(mroot) << std::endl;
	std::cout << "Least Squares Regression entered precision error: " << precision << std::endl;
	std::cout << "Least Squares Regression intercept error: " << Finderrorinc(mroot) << std::endl;
	
	//std::cout << "Least Squares Regression calculated slope error by York: ";
	//std::cout << yorkm << std::endl;
	//std::cout << "Least Squares Regression intercept error by York: ";
	//std::cout << yorkc << std::endl;
	
	putchar(7);
	
	return 0;
}
