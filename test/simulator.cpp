/*
 * This program is a simple simulator of sequencing machines with a given error rate and genome reference.  
 * It can also inject a mutation rate into a reference before the sequencing simulation progress. 
 *
 * For example: 
 * simulator -mutate 0.05 -err 0.01 -orgRef ../human_ref.fa -mutateRef org.fa -errReads human_seq.fa
 *
 * The above command means that, the reference human_ref.fa will be mutated with a rate of 5%, and then 
 * the mutated reference will be stored in org.fa. After that, the mutated reference can be sequenced by the 
 * simulator with a coverage of 50X (default) and readlength of 100bp(default). Finally the simulated reads 
 * will be stored in human_seq.fa.  
 *
 *
 * Author: Jintao Meng
 * */

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <stdexcept>
#include <fstream>
#include <random>
#include <chrono>
#include <time.h>
#include <algorithm>
using namespace std;

class myParameters
{
public:
string refFileName;
string outFileName;
string errFileName;
int    readLength;
int    coverage;
double errorRate;
double mutateRate;
double deltaNormal;
bool   pairEndFlag;
double insertLength;
vector<double> qualityScore;
	
myParameters()
{
	refFileName=string("");
	outFileName=string("");
	errFileName=string("");
	readLength=100;
	coverage=50;
	errorRate=0.01;
	mutateRate=0.05;
	deltaNormal=0.4;
	pairEndFlag=false;
	insertLength=2048;

	srand((unsigned)time(NULL));

	double a[100] = {37.768, 37.359,37.414,37.431,37.437,37.484,37.438,37.418,37.437,37.416, 37.489, 37.462, 37.444, 37.417, 37.392, 37.434, 37.398, 37.367, 37.346, 37.320, 37.319, 37.284, 37.224, 37.200, 37.174, 37.174, 37.118, 37.066, 37.027, 36.984, 36.994, 36.924, 36.874, 36.810, 36.764, 36.768, 36.660, 36.542, 36.425, 36.445, 36.486, 36.394, 36.278, 36.210, 36.139, 36.174, 36.068, 35.970, 35.884, 35.759, 35.812, 35.696, 35.569, 35.426, 35.350, 35.375, 35.216, 35.102, 34.993, 34.845, 34.881, 34.728, 34.574, 34.428, 34.270, 34.313, 34.108, 33.949, 33.799, 33.672, 33.734, 33.558, 33.356, 33.171, 32.995, 33.095, 32.885, 32.640, 32.463, 32.258, 32.318, 32.107, 31.839, 31.592, 31.409, 31.466, 31.198, 30.911, 30.643, 30.386, 30.437, 30.127, 29.791, 29.471, 29.178, 29.156, 28.813, 28.440, 28.046, 27.589};
	for(int i=0;i<100;i++)		qualityScore.push_back(a[i]);
}


void setRefFileName(string fileName)
{
	refFileName = fileName;
}

void setMutatedRefFileName(string fileName)
{
	outFileName = fileName;
}

void setErrFileName(string fileName)
{
	errFileName = fileName;
}	

void setErrRate(double err)
{
	errorRate = err;
}
void setMutateRate(double rate)
{
	mutateRate = rate;	
}
void setCoverage(int cov)
{
	coverage=cov;
}
void setReadLen(int readLen)
{
	readLength=readLen;
}
void parametersInit(int argc, char *argv[])
{
	if(argc==1)
	{
		cout<<"Usage: readIO -mutate mutateRate -err errorRate -orgRef ref.fa -mutateRef org.fa -errReads err.fa"<<endl;
		exit(0);
	}

	int ErrFlag=0;
	for(int i=1;i<argc;i++)
	{
		if(strcmp(argv[i], "-mutate")==0)
		{
			double t = atof(argv[++i]);
			if(t<0||t>1)	
			{
				cout<<"mutation rate should be a number in [0,1]"<<endl;
				exit(0);
			}
			setMutateRate(t);	
		}
		else if(strcmp(argv[i],"-err")==0)	
		{
			double t = atof(argv[++i]);
			if(t<0||t>1)	
			{
				cout<<"error rate should be a number in [0,1]"<<endl;
				exit(0);
			}
			setErrRate(t);
		}
		else if(strcmp(argv[i], "-orgRef")==0)		setRefFileName(string(argv[++i]));
		else if(strcmp(argv[i], "-mutateRef")==0)	setMutatedRefFileName(string(argv[++i]));	
		else if(strcmp(argv[i], "-errReads")==0)	setErrFileName(string(argv[++i]));
		else if(strcmp(argv[i], "-cov")==0)		setCoverage(atoi(argv[++i]));
		else if(strcmp(argv[i], "-readLen")==0)		setReadLen(atoi(argv[++i]));
		else 
		{
			cout<<"Usage: readIO -mutate mutateRate -err errorRate -orgRef ref.fa -mutateRef org.fa -errReads err.fa"<<endl;
			exit(0);
		}
	}		

	cout<<argv[0]<<" MutateRate:"<<mutateRate<<" ErrRate:"<<errorRate<<" RefFileName:"<<refFileName<<" MutateRefFileName:"<<outFileName<<" ErrReadsFileName:"<<errFileName<<endl;

}
};

class readIO
{
public:
	vector<string> reads;
	vector<string> quality;
	bool           isFasta;            
	bool           isPairEnd;
	int            insertLength;

	void setInsertLength(bool isFasta, bool isPairEnd, int insertLength);
	void getReads(string FileName);
	void printReads(string FileName);
};

void readIO::printReads(string FileName)
{
	ofstream fs(FileName.c_str());
	for(int i=0;i<reads.size();i++)
	{
		fs<<">Read_"<<i<<endl;
		fs<<reads[i]<<endl;
	}
	fs.close();
}

void readIO::setInsertLength(bool isFasta, bool isPairEnd, int length)
{
	this->isFasta = isFasta;
	this->isPairEnd = isPairEnd;
	this->insertLength = length;	
}

void readIO::getReads(string FileName)
{
	ifstream fs(FileName.c_str());
      	if(FileName.length()==0 || !fs.is_open())  
       	{ 
		cout<<"Error opening file"<<endl; 
		exit(1); 
	}  

	string buf, str;
	reads.clear();

        while(fs.eof()==0)
        {
       	        getline(fs, buf, '\n');
		if(buf[0]=='>'||buf[0]=='+')	
		{
			if(str.length())	reads.push_back(str);
			str.clear();
		}
		else 	str += buf;	
        }
	if(str.length())	reads.push_back(str);
	fs.close();

	for(int i=0;i<reads.size();i++)
	{
		string buf;
		for(long long int j=0;j<(long long int) reads[i].length();j++)
			if(reads[i][j]=='A'||reads[i][j]=='C'||reads[i][j]=='G'||reads[i][j]=='T')	buf.push_back(reads[i][j]); 
		reads[i]=buf;		
	}	
	long long int tot=0;
	for(int i=0;i<reads.size();i++)	tot+=reads[i].length(); 
	cout<<"Reading reference data finished\n"<<"The total number of nucletides in the reference is "<<tot<<endl;
}

class readSimulator
{
	bool           isFasta;            
	bool           isPairEnd;
	int            insertLength;
	vector<string> orgReads;
	vector<string> errReads;
	vector<string> quality;	

public:
	string reverseComplement(string str);	
	vector<double> generateQuality(vector<double> qualityModel, double errorRate, double delta);
	void insertError(string& orgRead, double errorRate);
	string qualityToString(vector<double> quality);

	
	void mutateRef(vector<string>& refs, double mutateRate, string mutatedRefFileName);
	void generateReads(string orgFileNmae, string errFileNmae, vector<string>& ref, myParameters curParameters);
	void generatePairFastq(string orgFileNmae, string errFileName, vector<string> ref, myParameters curParameters);
//	void generateSingleFasta(string FileName, vector<string> ref, myParameters curParameters);
//	void generateFastq(string FileName);
};

string readSimulator::qualityToString(vector<double> quality)
{
	string ret;
	for(int i=0;i<quality.size();i++)
	{
		ret.push_back('!' + (int) (quality[i]+0.5) );
	}
	ret[0] = '!';	
	return ret;
}

string readSimulator::reverseComplement(string str)
{
	for(int i=0;i<str.length();i++)
	{
		if(str[i]=='A')		str[i]='T';
		else if(str[i]=='C')	str[i]='G';
		else if(str[i]=='G')	str[i]='C';
		else 			str[i]='A';
	}

	for(int i=0;i<str.length()/2;i++)
	{
		char tmp = str[i];
		str[i] = str[str.length()-1-i];
		str[str.length()-1-i] = tmp;
	}
	return str;
}


void readSimulator::insertError(string& read, double errorRate)
{
	char ch[4] = {'A', 'C', 'G', 'T'};
	for(unsigned int i=0;i<read.length();i++)
	{
		double number = rand()%10000; 
		if(number<errorRate*10000)	read[i] = ch[(rand()%4)];
	}
}

void readSimulator::mutateRef(vector<string>& refs, double mutateRate, string mutatedRefFileName)
{
	double totLen=0, totErr=0;
	for(int i=0;i<refs.size();i++)
	{
		totLen += refs[i].length();
		string ret = refs[i];

		insertError(refs[i], mutateRate);
		for(int j=0;j<refs[i].length();j++)	if(refs[i][j]!=ret[j])		totErr++;
	}
	cout<<"The total number of mutated bases in reference is: "<<totErr<<"("<<totErr*100.0/totLen<<"%)"<<endl;	
	if(mutatedRefFileName.length()==0)	return ;
	ofstream fs(mutatedRefFileName.c_str());
      	if(!fs.is_open())  
       	{ 
		cout<<"Error opening mutatedRefFile :"<<mutatedRefFileName<<endl; 
		exit(1); 
	}  

	for(int i=0;i<refs.size();i++)	fs<<"<"<<i<<endl<<refs[i]<<endl;
	fs.close();
}		

void readSimulator::generateReads(string orgFielNmae, string errFileName, vector<string>& ref, myParameters curParameters)
{
	if(errFileName.length()==0)	return ;
	ofstream fs(errFileName.c_str());
      	if(!fs.is_open())  
       	{ 
		cout<<"Error opening readFile :"<<errFileName<<endl; 
		exit(1); 
	}  


	long long int refLength = 0;
	vector<int> Pos(ref.size(), 0);
	for(int i=0;i<ref.size();i++)	{refLength += ref[i].length(); Pos[i]=refLength/10; }

	long long int readNum = refLength * curParameters.coverage / curParameters.readLength;  
	cout<<"The total number of reads generated is "<<readNum<<endl;

	long long int readID = 0;
	long long int errNums=0;
	while(readID<=readNum)
	{
		long long int localPos = rand()%(refLength/10); 
		long long int r = lower_bound(Pos.begin(), Pos.end(), localPos)-Pos.begin();

		if((long long int) ref[r].length()<=curParameters.readLength)	continue;	
		readID++;
		
		long long int pos = rand()%((long long int)ref[r].length()+1-curParameters.readLength);	
		fs<<">"<<readID<<"_refID:"<<r<<"_refLength:"<<ref[r].size()<<"_Read_pos:["<<pos<<"~"<<(pos+curParameters.readLength)<<"]"<<endl;
//		fs<<">"<<readID<<endl;	
	
		string read = ref[r].substr(pos, curParameters.readLength);
		string ret=read;
		insertError(ret, curParameters.errorRate);
		for(int j=0;j<ret.length();j++)	if(read[j]!=ret[j])	errNums++; 
		fs<<ret<<endl;
	}
	cout<<"The total number errous nucletides generated in the data is "<<errNums<<"("<<errNums*100.0/(refLength*curParameters.coverage)<<"%)"<<endl;
	fs.close();
}

int main(int argc, char *argv[])
{
	myParameters curParameters;
	curParameters.parametersInit(argc, argv);

	readIO  readsData;
	bool    fastaFlag = true;
	readsData.setInsertLength(fastaFlag, curParameters.pairEndFlag, curParameters.insertLength);
	
	readsData.getReads(curParameters.refFileName); 
	cout<<"Reading References Finished"<<endl;

	readSimulator readSim;
	readSim.mutateRef(readsData.reads, curParameters.mutateRate, curParameters.outFileName);		
	cout<<"Mutating References Finished"<<endl;
	readSim.generateReads(curParameters.outFileName, curParameters.errFileName, readsData.reads, curParameters);

        return 0;
}

