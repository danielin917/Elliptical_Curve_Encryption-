#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <gmpxx.h>
#include <utility>
#include "ec_ops.h"
using namespace std;

Zp Zp::inverse() const{
//function inverse(a, n)
mpz_class  g, s, a, b;
mpz_t t;
mpz_init(t);
g = 1;
a = PRIME;
b = this -> value;
//use stl gcd function
   mpz_gcdext (g.get_mpz_t(), s.get_mpz_t(), t, a.get_mpz_t(), b.get_mpz_t());
   mpz_class output(t);
   Zp inv(output);
   return inv;

//Compute g, s, and t, such that as + bt = g = gcd (a, b). If t is NULL, that argument is not computed. 

//Set res to the greatest common divisor of operand1 and operand2. 
	// Implement the Extended Euclidean Algorithm to return the inverse mod PRIME		
}


ECpoint ECpoint::operator + (const ECpoint &a) const {
//account for infinity point
	if(this -> infinityPoint){
		return a;
	}
	else if(a.infinityPoint){
		return *this;
	} 
//first way of adding
	else if(!(*this == a) && (this -> x).getValue() != a.x.getValue()){
		mpz_class t,b,d;
	
		t = (a.y - this -> y).getValue() % PRIME;
		b = ( a.x - this -> x).getValue() % PRIME;
	
			Zp inv(b);
			//cout << b <<"????????\n" <<flush;
			Zp n(t);
			inv = inv.inverse();
			//cout << inv <<"+++++\n" <<flush;
			d = (inv*n).getValue();
			
		//cout << d <<"-------\n" <<flush;
		
		//cout << d <<endl;
		
		//account for negative mod
		mpz_class X,Y;
		X =(d*d)-(this -> x).getValue()-a.x.getValue();
			while(X < 0){
				X += PRIME;
			}
		Zp xx(X);
		Y =(this -> x - xx).getValue()*d - (this -> y).getValue();
			while(Y < 0){
				Y += PRIME;
			}
		
		Zp yy(Y);
		ECpoint ret(xx,yy);
		return ret;
		
	}
//second way of adding
	else if((*this == a) && (this -> y.getValue())*2 != 0){
		mpz_class t,b,d;
		
		t = (3*(this -> x * this -> x).getValue() + A) % PRIME;
		b = 2*(this -> y).getValue() % PRIME;
		
			Zp inv(b);
			Zp n(t);
			inv = inv.inverse();
			d = (inv*n).getValue();
			
			
		//account for negative mod
		mpz_class X,Y;
		X = (d*d) - 2*(this -> x).getValue();
			while(X < 0){
				X += PRIME;
			}
		Zp xx(X);
		Y = (this -> x - xx).getValue()*d - (a.y.getValue());
			while(Y < 0){
				Y += PRIME;
			}
		
		Zp yy(Y);
		ECpoint ret(xx,yy);
		return ret;
		
	
	}
//else infinity point
	else{
		ECpoint ret;
		ret.infinityPoint = 1;
		return ret;
	}
	// Implement  elliptic curve addition 		
}


ECpoint ECpoint::repeatSum(ECpoint p, mpz_class v) const {
	ECpoint x(1);
	mpz_t lsb;
	mpz_init(lsb);
	mpz_class i,h, b;
	h = 1;
	mpz_class a(1);
	long unsigned int shifter;
	size_t m = mpz_sizeinbase(v.get_mpz_t(),2);
	
	mpz_class M(m);
	
	//Double and Add method
	for(i = M; i >= 0; i--){

		x = x + x;
		mpz_t bin;
		mpz_init (bin);
		
		shifter = i.get_si();
		//cout << i<<flush;
		mpz_fdiv_q_2exp (bin, v.get_mpz_t(), shifter);
		
		//cout << bin<<"\n" << flush;
		 mpz_and (lsb, bin, h.get_mpz_t());
		// printf("%zu\n", m); 
		//cout << lsb << "     " << i <<"\n" <<flush;
		
		if(mpz_get_ui (lsb)){
			x = x + p;
		}
		
		//cout <<x.x << ","<< x.y <<"\n";
		
	}
	return x;
	//Find the sum of p+p+...+p (vtimes)		
}

Zp ECsystem::power(Zp val, mpz_class pow) {
	//repeated squaring
	//cout << val <<"\n";
	if(pow < 0){
	
	}
	else if(pow == 0){
		Zp ret(1); 
		return ret;
	}
	else if(pow == 1){
		Zp ret(val.getValue());
		return ret;
	}
	else if(pow % 2 == 0){
		
		return power(val*val, pow/2);
	}
	else{
		return val * power(val*val, (pow-1)/2);
	}
	//Find the sum of val*val+...+val (pow times)
}


mpz_class ECsystem::pointCompress(ECpoint e) {
	//It is the gamma function explained in the assignment.
	//Note: Here return type is mpz_class because the function may
	//map to a value greater than the defined PRIME number (i.e, range of Zp)
	//This function is fully defined.	
	mpz_class compressedPoint = e.x.getValue();
	compressedPoint = compressedPoint<<1;
	
	if(e.infinityPoint) {
		cout<<"Point cannot be compressed as its INF-POINT"<<flush;
		abort();
		}
	else {
		if (e.y.getValue()%2 == 1)
			compressedPoint = compressedPoint + 1;
		}
		//cout<<"For point  "<<e<<"  Compressed point is <<"<<compressedPoint<<"\n";
		return compressedPoint;

}

ECpoint ECsystem::pointDecompress(mpz_class compressedPoint){
	mpz_class bit = compressedPoint & 1;
	compressedPoint = compressedPoint & -2;
	compressedPoint = compressedPoint >> 1;
	Zp xr(compressedPoint);
	//compressed point now equals xr
	mpz_class z = compressedPoint*compressedPoint*compressedPoint + A*compressedPoint + B;
	Zp square(z);
	Zp residue = power(square, (PRIME + 1)/4);
	if((residue.getValue() & 1) == bit){
		ECpoint decompressed(xr, residue);
		return decompressed;
	
	}
	else{
		Zp prime(0);
		ECpoint decompressed(xr, prime - residue);
		return decompressed;
	
	}
	
	//Implement the delta function for decompressing the compressed point
}


pair<mpz_class, mpz_class> ECsystem::encrypt(ECpoint publicKey, mpz_class privateKey,mpz_class plaintext){
	// You must implement elliptic curve encryption
	//  Do not generate a random key. Use the private key that is passed from the main function
	ECpoint Q, R;
	Q= Q.repeatSum(G,privateKey);
	
	
	R= R.repeatSum(publicKey, privateKey); 
	
	
	mpz_class C1 = pointCompress(Q); 
	mpz_class C2 = pointCompress(R)^plaintext;
	
	return pair<mpz_class, mpz_class>(C1,C2);
	//Q=private * G R = Private * Public
	//set C1 = compress(Q) and C2 M xor Compress(R) 
}


mpz_class ECsystem::decrypt(pair<mpz_class, mpz_class> ciphertext){
	ECpoint R;
	R = pointDecompress(ciphertext.first);
	R = R.repeatSum(R,privateKey);
	
	return (ciphertext.second)^(pointCompress(R));
	// Implement EC Decryption
}


/*
 * main: Compute a pair of public key and private key
 *       Generate plaintext (m1, m2)
 *       Encrypt plaintext using elliptic curve encryption
 *       Decrypt ciphertext using elliptic curve decryption
 *       Should get the original plaintext
 *       Don't change anything in main.  We will use this to 
 *       evaluate the correctness of your program.
 */


int main(void){
	srand(time(0));
	ECsystem ec;
	mpz_class incrementVal;	
	
	/*
	cout << "addition test:" <<"\n";
	ECpoint R, Z;
	ECpoint Q(1);
	
	R.x = 2;
	R.y = 4;
	cout << "(" <<Q.x<<"," <<Q.y << ")+(" <<R.x <<"," << R.y <<")=\n";
	Z = Q + R;
	cout << "(" <<Z.x<<"," <<Z.y << ")";
	mpz_class n = 6;
	cout << "Repeated sum test:" <<"\n";
	cout << "(" <<Q.x<<"," <<Q.y << ")*"<<n<<"\n";
	cout << (Q.repeatSum(Q,n)).x <<"," << (Q.repeatSum(Q,n)).y <<"\n";
	Zp base(2);
	mpz_class pow(4);
	Zp power = ec.power(base, pow);
	cout << "--0-0-0-" <<  power << "-0-0-0-\n" <<flush; 
	
	*/
	
	cout<<"Public key is: "<<"\n";
	
	pair <ECpoint, mpz_class> keys = ec.generateKeys();
	
	mpz_class plaintext = MESSAGE;
	ECpoint publicKey = keys.first;
	cout<<"Public key is: "<<publicKey<<"\n";
	
	cout<<"Enter offset value for sender's private key"<<endl;
	cin>>incrementVal;
	mpz_class privateKey = XB + incrementVal;
	//cout << "Private key is: "<<privateKey<<"\n"; //remove this it is added
	pair<mpz_class, mpz_class> ciphertext = ec.encrypt(publicKey, privateKey, plaintext);	
	cout<<"Encrypted ciphertext is: ("<<ciphertext.first<<", "<<ciphertext.second<<")\n";
	mpz_class plaintext1 = ec.decrypt(ciphertext);
	
	cout << "Original plaintext is: " << plaintext << endl;
	cout << "Decrypted plaintext: " << plaintext1 << endl;


	if(plaintext == plaintext1)
		cout << "Correct!" << endl;
	else
		cout << "Plaintext different from original plaintext." << endl;	
	
	return 1;
	
}

