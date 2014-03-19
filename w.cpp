// MAPPING ASSEMBLY CODE
//  HONOURS PROJECT 
//  SMITH WATERMAAN AND INDEXING THE REFERENCE GENOME 


///////////////////////////////////////////////////////////////
#include<iostream>
#include<string>
#include<vector>
#include<map>
#include<stdlib.h>
#include<algorithm>
#include<omp.h>
using namespace std;
///////////////////////////////////////////////////////////////////
int align(string,string); 	// smith watermann algorithm
int func(char,char);		// matching
int max(int,int,int);		// returns maximum of three numbers
string reversecomp(string);	// reverse compliment of a string 
void mapread(int);		// MAIN MAPPING FUNCTION "i" IS THE READNO.
void calculatefreq();		// FOR CACULATING THE FREQUENCY OF THE BASES AT A POSITION
void printseq();		// PRINTING THE SEQUENCE
void *func1(void *);		// FOR ADDING THE PARALLELISM
////////////////////////////////////////////////////////
typedef struct reads{
	int len;			// read structure
	string r;		// contains read and starting postion of the read NEEDS TO MODIFIED FOR REAL DATA 
	int pos;		// gives the final position of the mapped read
	int rev;		// flag for the reversed=1 and forward=0
	string s;
}reads;
///////////////////////////////////////////////////////////////
bool bokka(const reads &a, const reads &b)
{	
	return a.pos < b.pos;
}
int gl; 		// genome length
int rl=100;			// read length
int sl=17;		// seed length
string genome;		// genome
int D=-1;
int maxpos=0;		// maximum position
int maxl=0;
vector<reads> pp;	// reads vector
map<string,vector<int> > mp;
map<string,vector<int> >::iterator it;
vector< vector <double> > freq;
int size;		// total number of reads 
//// introducing quality score ////////////
double qs[]={0.0 , 0.205671765276 , 0.36904265552 , 0.498812766373 , 0.601892829447 , 0.683772233983 , 0.748811356849 , 0.800473768503 , 0.841510680754 , 0.874107458821 , 0.9 , 0.920567176528 , 0.936904265552 , 0.949881276637 , 0.960189282945 , 0.968377223398 , 0.974881135685 , 0.98004737685 , 0.984151068075 , 0.987410745882 , 0.99 , 0.992056717653 , 0.993690426555 , 0.994988127664 , 0.996018928294 , 0.99683772234 , 0.997488113568 , 0.998004737685 , 0.998415106808 , 0.998741074588 , 0.999 , 0.999205671765 , 0.999369042656 , 0.999498812766 , 0.999601892829 , 0.999683772234 , 0.999748811357 , 0.999800473769 , 0.999841510681 , 0.999874107459 , 0.9999 , 0.999920567177 , 0.999936904266 , 0.999949881277 , 0.999960189283 , 0.999968377223 , 0.999974881136 , 0.999980047377 , 0.999984151068 , 0.999987410746 , 0.99999 , 0.999992056718 , 0.999993690427 , 0.999994988128 , 0.999996018928 , 0.999996837722 , 0.999997488114 , 0.999998004738 , 0.999998415107 , 0.999998741075 , 0.999999 , 0.999999205672 , 0.999999369043 , 0.999999498813 , 0.999999601893 , 0.999999683772 , 0.999999748811 , 0.999999800474 , 0.999999841511 , 0.999999874107 };
//////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
	int i,j,k;
	string s,rq;

	cin>>gl;
	reads re;
	cin>> genome;
	for(i=0;i<gl-sl;i++)
	{
		s.assign(genome,i,sl);
		mp[s].push_back(i);
	}
	cout<<mp.size()<<endl;
	cout<<"Started Reading\n";
	cin>>i;
	while(i!=-1)
	{
		cin>>s;
		cin>>rq;					// taking the read in 
		re= {i,s,-1,0,rq};
		pp.push_back(re);				// push it back
		cin>>i;
	}

        size=pp.size();
	cout<<"Reading done\n";
//////////////////////////////////////////////////////////////////////////////////
#pragma omp parallel for
	
	for(i=0;i<size;i++)
	{
		mapread(i);	// mapping the reads
	}
	
////////////////////////////////////////////////////////////////////////////////
/*	int x=size/32;
	int b[32][2];
	for(i=0;i<31;i++)
	{
		b[i][0]=i*x;
		b[i][1]=(i+1)*x;
	}
	b[31][0]=b[30][1];
	b[31][1]=size;
	int aa,ab;
	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&aa);
	MPI_Comm_rank(MPI_COMM_WORLD,&ab);
	for(i=b[ab][0];i<b[ab][1];i++)
		mapread(i);
	MPI_Finalize();

	cout<<"MAPPING DONE\n";*/
/////////////////////////////////////////////////////////////////////////////////	
//	sort( pp.begin(), pp.end() ,bokka);	
//	for(i=0;i<size;i++)
//		cout<<"real pos = "<<pp[i].b<<" map pos = "<<pp[i].pos<<endl;
	int notmap=0;
	for(i=0;i<size;i++)
	{
		if(pp[i].pos==-1)
		{
			notmap++;
		}
	if(maxpos<pp[i].pos) { maxpos=pp[i].pos;maxl=pp[i].len;}
	}
	cout<<" No of reads not mapped "<<notmap<<"\n";


	calculatefreq();
	printseq();
	return 0;
}	
////////////////////////////////////////////////////////////////////////////////////////////
int align(string s1, string s2)
{
	int mx=-98999;     
			// Smith Watermaan Algorithm
	int l=s1.length();
	int m=s2.length();
	int s[l+1][m+1];
	int i,j,k;
	for(i=0;i<l+1;i++)             // ONLY A
		s[i][0]=0;

	for(i=0;i<m+1;i++)
		s[0][i]=0;                          // ONLY B
	for(i=1;i<l+1;i++)
	{
		for(j=1;j<m+1;j++)
		{
			int m,in,de;
			m= s[i-1][j-1]+ func(s1[i],s2[j]);  // MATCH 
			de=s[i-1][j]+D;                     // DELETE
			in=s[i][j-1]+D;                    //  INSERT
			s[i][j]= max(m,de,in); 
			if(s[i][j]<0)
				s[i][j]=0;// MAXIMUM VALUE STORED
			if(mx<s[i][j])
				mx=s[i][j];
		}
	}
	return mx;   // MAXIMUM VALUE RETURNED
}
int func(char a, char b)
{
	if(a==b)
		return 1;
	else 
		return -1;
}
/////////////////////////////////////////////////////////////////////////////
int max(int a, int b, int c)
{
	int k=a;
	if(a<b)
		k=b;
		if(k<c)
		k=c;
	return k;
}
/////////////////////////////////////////////////////////////////////////////
void mapread(int no)
{
	int p;
	int n;
	int i,j,k=-1;
	int a1=-1,a2=-1;
	int mx;
	int px=-1;
	string s,s1,s2;
	int rev=0;
	int x1,x2,leng;
	rl=pp[no].len;
	while(rev!=2 && k==-1)
	{
		if(rev==1){
				pp[no].r=reversecomp(pp[no].r);
				pp[no].rev=1;
			}

		for(p=0;p<(rl-sl)&&(k==-1);p++)
		{
			mx=0;
			px=-1;
			s.assign(pp[no].r,p,sl);
			a1=a2=0;
			vector<int> v;
			v=mp[s];
			n=mp[s].size();
			for(i=0;i<n;i++)
			{
				
				if (p>0&&v[i]>0)
				{
					if(v[i]>=p)
					{
						x1=0;
						x2=v[i]-p;
						leng=p;

					}
					else
					{
						x1=p-v[i];
						x2=0;
						leng=v[i];

					}


					s1.assign(pp[no].r,x1,leng);
					s2.assign(genome,x2,leng);
					a1=align(s1,s2);
				}
				
				if(p<rl-sl)
				{
					if(v[i]+rl-p<gl)
					{
						x1=p+sl;
						x2=v[i]+sl;
						leng=rl-p-sl;
					}
					else
					{
						x1=p+sl;
						x2=v[i]+sl;
						leng=gl-v[i]-sl;
					}
					s1.assign(pp[no].r,x1,leng);			// all the checks 
					s2.assign(genome,x2,leng);
					a2=align(s1,s2);
				}
				
				
				
				if((a1+a2)>=58&&(mx<a1+a2))
				{
					px=v[i]-p;
					mx=a1+a2;				// taking both the sum to be 64 . SHould be changed
				}
			}
			k=px;							// k is a flag
		}
		rev++;
	}
	pp[no].pos=k;					//	k contains the final pos
	//if(maxpos<k) maxpos=k;


}
/////////////////////////////////////////////////////////////////////////////////
string reversecomp(string s)
{
	int i,j;
	string s1;
	i=s.length();
	for(j=i-1;j>=0;j--)
	{
		switch(s[j])
		{
			case 'A':
				s1.append("T");
				break;
			case 'T':
				s1.append("A");
				break;
			case 'G':
				s1.append("C");
				break;
			case 'C':
				s1.append("G");
				break;
			default:
				break;
		}
	}
	return s1;
}
///////////////////////////////////////////////////////////////////////////////////////////////
void calculatefreq()
{
	int i,j,k;
	//cout<<maxpos<<"\n";
	k=maxpos+maxl;
	//cout<<k<<"\n";
	freq.resize(k);
	for(i=0;i<k;i++)
		freq[i].resize(4);
	for(i=0;i<size;i++)
	{
		if(pp[i].pos!=-1)
		{
			for(j=0;j<rl;j++)
			{
				switch (pp[i].r[j])
				{
					case 'A':
						freq[pp[i].pos+j][0]+=qs[(int)pp[i].s[j]-33];
						break;
					case 'C':
						freq[pp[i].pos+j][1]+=qs[(int)pp[i].s[j]-33];
						break;
					case 'G':
						freq[pp[i].pos+j][2]+=qs[(int)pp[i].s[j]-33];
						break;
					case 'T':
						freq[pp[i].pos+j][3]+=qs[(int)pp[i].s[j]-33];
						break;
				}
			}
		}
	}

//	for(i=0;i<k;i++)
//		cout<<i<<" "<<"A= "<<freq[i][0]<<" C= "<<freq[i][1]<<" G= "<<freq[i][2]<<" T= "<<freq[i][3]<<endl;
}
///////////////////////////////////////////////////////////////////////////////////////////////
void printseq()
{
	int i,j,k=maxpos+rl;
	int max;
	int num;
	 string out;
	for(i=0;i<k;i++)
	{
		max=0;
		num=-1;
		for(j=0;j<4;j++)
			if(freq[i][j]>max)
			{
				max=freq[i][j];
				num=j;
			}
		
			switch(num){
				case 0:
					out.append("A");
					break;
				case 1:
					out.append("C");
					break;
				case 2:
					out.append("G");
					break;
				case 3:
					out.append("T");
					break;
				default:
					out.append("-");
					break;
			}
	}
	cout<<out<<endl;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////
void *func1(void *t)
{
	int *p=(int *)t;
	int i;
	for(i=p[0];i<p[1];i++)
		mapread(i);
}

		








		


