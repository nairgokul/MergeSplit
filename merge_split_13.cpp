#include<iostream>
#include<iomanip>
#include<vector>
#include<cstdlib>
#include<chrono>
#include<random>
#include<algorithm>
#include<ctime>
#include<fstream>
#include<string>

using namespace std;

// Random Number Generator Initialisation
unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
default_random_engine gen(seed1);

// Function Declarations
int get_index(vector<int> x,int key);
void clean_vector(vector<int>& x,int key);
void clean_vector(vector<uint64_t>& x,int key);
void clean_counter(vector<vector<int> >& x);
void clean_counter(vector<vector<uint64_t> >& x);
vector<int> seq(int start,int end,int step);
int hash_counter(int N,int N1,int i,int j,bool unhash);
void print_vec(vector<int>& x);
void print_vec(vector<double>& x);
void print_vec_of_vec(vector<vector<int> >& x);
void print_vec_of_vec(vector<vector<double> >& x);
void print_vec_of_vec(vector<vector<uint64_t> >& x);
void random_distribution(vector<vector<int> >& counter,int N,int N1,int s,bool print);
vector<vector<double> > vector_deepcopy(vector<vector<int> > x);
int custom_distr(vector<double> x,int start);
int custom_distr(vector<int> x,int start);
vector<double> row_sum(vector<vector<double> > x);
vector<int> row_sum(vector<vector<int> > x);
double het_split_rate(double Ps0,double d,int n,int k,double a);
double sum_vector(vector<double>& x,int start,int end);
int sum_vector(vector<int>& x,int start,int end);
int longest_group(vector<vector<int> >& x);
int longest_group(vector<vector<uint64_t> >& x);
int rand_int(int low,int high);
bool ber(double p);
vector<int> randomly_split(int n,int k);
const string currentDateTime();
void write_into_file(string file_name1,string file_name2,vector<vector<uint64_t> > group_counter);

// Function searches for element in an vector, and returns the first index of its occurence
int get_index(vector<int> x,int key)
{
	for(int i=0;i<x.size();i++)
	{
		if(x[i]==key)
		{
			return i;
		}
	}
	cout<<"Error: get_index(), "<<key<<" could not be found in x."<<endl;
	exit(0);
}

// Removes all trailing zeroes in x
void clean_vector(vector<int>& x,int key=0)
{
	while(x.back() == key && !x.empty())
	{
		x.pop_back();
	}
}

void clean_vector(vector<uint64_t>& x,int key=0)
{
	while(x.back() == key && !x.empty())
	{
		x.pop_back();
	}
}

// Removes empty entries and trailing zeroes in sub-vectors
void clean_counter(vector<vector<int> >& x)
{
	int i;
	for(i=0;i<x.size();i++)
	{
		if(!x[i].empty())
		{
			clean_vector(x[i]);
		}
	}
	while((x.back()).empty() && !x.empty())
	{
		x.pop_back();
	}
}

void clean_counter(vector<vector<uint64_t> >& x)
{
	int i;
	for(i=0;i<x.size();i++)
	{
		if(!x[i].empty())
		{
			clean_vector(x[i]);
		}
	}
	while((x.back()).empty() && !x.empty())
	{
		x.pop_back();
	}
}

// Simple sequence generator, only works for increasing sequences
vector<int> seq(int start,int end,int step=1)
{
	
	vector<int> x((float)(abs(end-start))/(float)(step),0);
	int j=0;
	for(int i=start;i<end;i+=step)
	{
		x[j]=i;
		j++;
	}
	return(x);
}

// Hash function that maps [group_size][no. of N1] to [group_size][j] indices in the counter
int hash_counter(int N,int N1,int i,int j,bool unhash = false)
{
	int lower = max(0,i-(N-N1)),higher = min(i,N1),k;
	if(unhash)
	{
		k = seq(lower,higher+1)[j];
		return(k);
	}
	else
	{
		k = get_index(seq(lower,higher+1),j);
		return(k);
	}
}

// Prints vectors
void print_vec(vector<int>& x)
{
	for(int i=0;i<x.size();i++)
	{
		cout<<x[i]<<endl;
	}
}		

void print_vec(vector<double>& x)
{
	for(int i=0;i<x.size();i++)
	{
		cout<<x[i]<<endl;
	}
}

// Prints vector of vectors
void print_vec_of_vec(vector<vector<int> >& x)
{
	int i,j;
	for(i=0;i<x.size();i++)
	{
		if(!x[i].empty())
		{
			for(j=0;j<x[i].size();j++)
			{
				cout<<x[i][j]<<" ";
			}
			cout<<endl;
		}
		else
		{
			cout<<0<<endl;
		}
	}
}

void print_vec_of_vec(vector<vector<double> >& x)
{
	int i,j;
	for(i=0;i<x.size();i++)
	{
		if(!x[i].empty())
		{
			for(j=0;j<x[i].size();j++)
			{
				cout<<x[i][j]<<" ";
			}
			cout<<endl;
		}
		else
		{
			cout<<0<<endl;
		}
	}
}

void print_vec_of_vec(vector<vector<uint64_t> >& x)
{
	int i,j;
	for(i=0;i<x.size();i++)
	{
		if(!x[i].empty())
		{
			for(j=0;j<x[i].size();j++)
			{
				cout<<x[i][j]<<" ";
			}
			cout<<endl;
		}
		else
		{
			cout<<0<<endl;
		}
	}
}

// Creates a random distribution by uniformly placing individuals in sites
void random_distribution(vector<vector<int> >& counter,int N,int N1,int s,bool print=false)
{
	vector<int> sub;
	int i,j;
	vector<vector<int> > y(s,vector<int>(2,0));
	uniform_int_distribution<int> dist(0,s-1);
	for(i = 0;i<N+1;i++)
	{
		counter.push_back(sub);
	}
	
	for(i=0;i<N1;i++)
	{
		y[dist(gen)][0]++;
	}
	
	for(i=0;i<(N-N1);i++)
	{
		y[dist(gen)][1]++;
	}
	
	for(i=0;i<y.size();i++)
	{
		if(y[i][0]+y[i][1] != 0)
		{
			j = hash_counter(N,N1,y[i][0]+y[i][1],y[i][0]);
			if(j>=counter[y[i][0]+y[i][1]].size())
			{
				counter[y[i][0]+y[i][1]].resize(j+1);
			}
			counter[y[i][0]+y[i][1]][j]++;
		}
	}
	clean_counter(counter);
	if(print)
	{
		print_vec_of_vec(counter);
	}
}

// Copies vector of vector of ints and returns vector of vector of doubles
vector<vector<double> > vector_deepcopy(vector<vector<int> > x)
{
	vector<vector<double> > y;
	vector<double> sub;
	for(int i=0;i<x.size();i++)
	{
		sub.clear();
		sub.assign(x[i].begin(),x[i].end());
		y.push_back(sub);
	}
	return(y);
}	

// Given a discrete probability mass function (may not be normalised) this function simulates it 
int custom_distr(vector<double> x,int start = 1)
{
	vector<double> y(x.begin(),x.end());
	int i;
	if(start<1)
	{
		start = 1;
	}
	
	for(i=1;i<y.size();i++)
	{
		if(i>=start)
		{
			y[i] = y[i]+y[i-1];
		}
		else
		{
			y[i]=0;
		}
	}
	double total = y.back();
	uniform_real_distribution<double> dist(0.0,1.0);
	double u = dist(gen);
	for(i=1;i<y.size();i++)
	{
		if(u<(y[i]/total))
		{
			return(i);
		}
	}
}

int custom_distr(vector<int> x,int start = 1)
{
	vector<double> y(x.begin(),x.end());
	int i;
	if(start<1)
	{
		start = 1;
	}
	
	for(i=1;i<y.size();i++)
	{
		if(i>=start)
		{
			y[i] = y[i]+y[i-1];
		}
		else
		{
			y[i]=0;
		}
	}
	double total = y.back();
	uniform_real_distribution<double> dist(0.0,1.0);
	double u = dist(gen);
	for(i=1;i<y.size();i++)
	{
		if(u<(y[i]/total))
		{
			return(i);
		}
	}
}

// Returns a vector, such that the i'th element is the sum of all the elements in the i'th row of the given vector 
vector<double> row_sum(vector<vector<double> > x)
{
	vector<double> y;
	for(int i=0;i<x.size();i++)
	{
		y.push_back(sum_vector(x[i],0,x[i].size()));
	}
	return(y);
}

vector<int> row_sum(vector<vector<int> > x)
{
	vector<int> y;
	for(int i=0;i<x.size();i++)
	{
		y.push_back(sum_vector(x[i],0,x[i].size()));
	}
	return(y);
}

// Returns the heterogenous split rate for a group of size n, with k individuals of type-1
double het_split_rate(double Ps0,double d,int n,int k,double a)
{
	return(Ps0 + (a*(k*(n-k)*d)/(n*n))-((1-a)*(k*(n-k)*d)/(n*n))+((1-a)*d/4));
}

// Returns the sum of elements in a vector between the indices 'start' and 'end'
double sum_vector(vector<double>& x,int start,int end)
{
	double total = 0;
	for(int i=start;i<end;i++)
	{
		total+=x[i];
	}
	return(total);
}

int sum_vector(vector<int>& x,int start,int end)
{
	int total = 0;
	for(int i=start;i<end;i++)
	{
		total+=x[i];
	}
	return(total);
}

int longest_group(vector<vector<int> >& x)
{
	int size = 0,largest;
	for(int i=0;i<x.size();i++)
	{
		if(x[i].size()>size)
		{
			size=x[i].size();
			largest=i;
		}
	}
	return(largest);
}

int longest_group(vector<vector<uint64_t> >& x)
{
	uint64_t size = 0;
	int largest;
	for(int i=0;i<x.size();i++)
	{
		if(x[i].size()>size)
		{
			size=x[i].size();
			largest=i;
		}
	}
	return(largest);
}	

int rand_int(int low,int high)
{
	uniform_int_distribution<int> dist(low,high);
	return(dist(gen));
}

bool ber(double p)
{
	bernoulli_distribution dist(p);
	return dist(gen);
}

vector<int> randomly_split(int n,int k)
{
	vector<int> permutation,return_vec;
	int j,swap,split_index,n11=0;
	for(int i=0;i<k;i++)
	{
		permutation.push_back(1);
	}
	for(int i=0;i<(n-k);i++)
	{
		permutation.push_back(0);
	}
	for(int i=0;i<n;i++)
	{
		j=rand_int(i,n-1);
		swap=permutation[i];
		permutation[i]=permutation[j];
		permutation[j]=swap;
	}
	split_index=rand_int(1,permutation.size()-1);
	for(int i=0;i<split_index;i++)
	{
		if(permutation[i]==1)
		{
			n11++;
		}
	}
	return_vec.push_back(split_index);
	return_vec.push_back(n11);
	return return_vec;
}

const string currentDateTime() 
{
    time_t now = time(0);
    struct tm tstruct;
    char buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
    return(buf);
}
// The main class involved in the simulation	
class MergeSplit
{
	public:
		// member variables 
		int N,N1,s;
		double Ps0,d,a,Pm;
		vector<vector<int> > counter;
		
		// member functions
		MergeSplit(int N,int N1,double Ps0,double d, double a, double Pm,int s);
		void next_event();
		int get_population();
		int get_group_count();
		
	private:
		void split();
		void merge();
};		

// Constructor
MergeSplit::MergeSplit(int N,int N1,double Ps0,double d, double a, double Pm,int s)
{
	this->N = N;
	this->N1 = N1;
	this->Ps0 = Ps0;
	this->d = d;
	this->a = a;
	this->Pm = Pm;
	this->s = s;
	random_distribution(this->counter,this->N,this->N1,this->s);
}

// Handles events
void MergeSplit::next_event()
{
	double split_rate=0,merge_rate=0,merge_weight=0;
	int k;
	for(int i=1;i<(this->counter.size()+1);i++)
	{
		if(i>1)
		{
			for(int j=0;j<this->counter[i].size();j++)
			{
				k=hash_counter(this->N,this->N1,i,j,true);
				split_rate+=(this->counter[i][j]*het_split_rate(this->Ps0,this->d,i,k,this->a));
			}
		}
		if(i<(this->counter.size()-1))
		{
			merge_rate+=sum_vector(this->counter[i],0,this->counter[i].size());
		}
	}
	merge_rate = merge_rate*(this->Pm);
	merge_weight = ((double) get_group_count())/((double) s);
	merge_weight=min(merge_weight,1.0);
	
	if(ber((split_rate/(split_rate+merge_rate))))
	{
		split();
	}
	else if(ber(merge_weight))
	{
		merge();
	}
	
	clean_counter(this->counter);
}

// Returns the population of the system, mostly for testing purposes
int MergeSplit::get_population()
{
	int n = 0,n_i;
	for(int i=0;i<this->counter.size();i++)
	{
		if(this->counter[i].size()>0)
		{
			n_i = 0;
			for(int j=0;j<this->counter[i].size();j++)
			{
				n_i+=this->counter[i][j];
			}
			n+=(i*n_i);
		}
	}
	return(n);
}	

// Returns the total number of groups
int MergeSplit::get_group_count()
{
	int n=0;
	for(int i=0;i<this->counter.size();i++)
	{
		for(int j=0;j<this->counter[i].size();j++)
		{
			n+=this->counter[i][j];
		}
	}
	return(n);
}

// Private function that handles the 'Split' event
void MergeSplit::split()
{
	vector<vector<double> > weighted_counter = vector_deepcopy(this->counter);
	int k;
	for(int i=0;i<weighted_counter.size();i++)
	{
		for(int j=0;j<weighted_counter[i].size();j++)
		{
			k = hash_counter(this->N,this->N1,i,j,true);
			weighted_counter[i][j] = weighted_counter[i][j]*het_split_rate(this->Ps0,this->d,i,k,this->a);
		}
	}
	
	int split_size_n = custom_distr(row_sum(weighted_counter),2);
	vector<double> temp = weighted_counter[split_size_n];

	temp.insert(temp.begin(),1,0);
	
	int split_size_k_hash = custom_distr(temp) - 1;
	
	int split_size_k = hash_counter(this->N,this->N1,split_size_n,split_size_k_hash,true);
	int new_pop1 = 0,new_pop2 = 0;
	uniform_int_distribution<int> dist1(0,split_size_k);
	uniform_int_distribution<int> dist2(0,split_size_n-split_size_k);
	
	while((new_pop1 == 0 && new_pop2 == 0) || (new_pop1 == split_size_k && new_pop2 == split_size_n - split_size_k))
	{
		if(split_size_k != 0)
		{
			new_pop1 = dist1(gen);
		}
		if((split_size_n-split_size_k) != 0)
		{
			new_pop2 = dist2(gen);
		}
	}
	int new_pop_hash1 = hash_counter(this->N,this->N1,new_pop1+new_pop2,new_pop1);
	if(this->counter[new_pop1+new_pop2].size()<=new_pop_hash1)
	{
		this->counter[new_pop1+new_pop2].resize(new_pop_hash1+1);
	}
	int new_pop_hash2 = hash_counter(this->N,this->N1,split_size_n - (new_pop1+new_pop2),split_size_k-new_pop1);
	if(this->counter[split_size_n - (new_pop1+new_pop2)].size()<=new_pop_hash2)
	{
		this->counter[split_size_n - (new_pop1+new_pop2)].resize(new_pop_hash2+1);
	}
	
	this->counter[split_size_n][split_size_k_hash]--;
	this->counter[new_pop1+new_pop2][new_pop_hash1]++;
	this->counter[split_size_n - (new_pop1+new_pop2)][new_pop_hash2]++;
}

// Private function that handles the 'Merge' event
void MergeSplit::merge()
{
	int merge_size_n1 = custom_distr(row_sum(this->counter));
	vector<int> temp = this->counter[merge_size_n1];
	temp.insert(temp.begin(),1,0);
	int merge_size_k1_hash = custom_distr(temp) - 1;
	this->counter[merge_size_n1][merge_size_k1_hash] --;
	int merge_size_k1 = hash_counter(this->N,this->N1,merge_size_n1,merge_size_k1_hash,true);
	
	int merge_size_n2 = custom_distr(row_sum(this->counter));
	temp.clear();
	temp.assign(this->counter[merge_size_n2].begin(),this->counter[merge_size_n2].end());
	temp.insert(temp.begin(),1,0);
	int merge_size_k2_hash = custom_distr(temp) - 1;
	this->counter[merge_size_n2][merge_size_k2_hash] --;
	int merge_size_k2 = hash_counter(this->N,this->N1,merge_size_n2,merge_size_k2_hash,true);
	
	int new_n=merge_size_n1+merge_size_n2,new_k=merge_size_k1+merge_size_k2;
	int new_k_hash = hash_counter(this->N,this->N1,new_n,new_k);
	
	if(this->counter.size()<=new_n)
	{
		this->counter.resize(new_n+1);
	}
	
	if(this->counter[new_n].size()<=new_k_hash)
	{
		this->counter[new_n].resize(new_k_hash+1);
	}
	
	this->counter[new_n][new_k_hash]++;
}

void write_into_file(string file_name1,vector<vector<uint64_t> > group_counter,int N,int N1,int s,double Ps0,double Pm,double d,int events_max)
{
	ofstream fp1;
	int i,j;
	
	fp1.open(file_name1.c_str());
	
	fp1<<"N,"<<N<<",N1,"<<N1<<",Ps0,"<<Ps0<<",d,"<<d<<",Pm,"<<Pm<<",s,"<<s<<",events,"<<events_max<<",\n";
	fp1<<",";
	
	for(i=0;i<longest_group(group_counter)+1;i++)
	{
		fp1<<i<<",";
	}
	fp1<<"\n";
	for(i=0;i<group_counter.size();i++)
	{
		fp1<<i<<",";
		if(!group_counter[i].empty())
		{
			for(j=0;j<group_counter[i].size();j++)
			{
				fp1<<group_counter[i][j]<<",";
			}
		}
		
		else
		{
			fp1<<0<<",";
			}
		fp1<<"\n";
	}
	fp1.close();
}

int main()
{	
	int N,N1,s,sampling_interval = 500,i,j,Z;
	uint64_t events_max;
	double Ps0,Pm,d,a;
	cout<<"Lattice Points (s)= ";
	cin>>s;
	cout<<"Population (N)= ";
	cin>>N;
	cout<<"Population of Type-1 (N1)= ";
	cin>>N1;
	cout<<"Homogeneous Split Rate (Ps0)= ";
	cin>>Ps0;
	cout<<"Merge Rate (Pm)= ";
	cin>>Pm;
	cout<<"d= ";
	cin>>d;
	cout<<"a= ";
	cin>>a;
	cout<<"Number of Events= ";
	cin>>events_max;
	
	clock_t start=clock();
	MergeSplit m(N,N1,Ps0,d,a,Pm,s);
	ofstream fp2;//,fp3;
	vector<vector<uint64_t> > group_counter;
	vector<uint64_t> v;
	string date_and_time(currentDateTime());
	string folder_name ("Output/");
	folder_name=folder_name+"10^"+to_string((int) log10((double) events_max))+"__Ps0="+to_string(Ps0)+"_Pm="+to_string(Pm)+"_d="+to_string(d)+"_a="+to_string(a);
	string mkdir="mkdir "+folder_name;
	system(mkdir.c_str());
	
	string file_name1="group_distr_"+date_and_time+".csv";
	string temp_file_name1;
	string file_name2=folder_name+"/freq_distr_"+date_and_time+".csv";
	//string file_name3=folder_name+"/Zvst_";;
	temp_file_name1=folder_name+"/."+file_name1;
	file_name1=folder_name+"/"+file_name1;
	file_name1=file_name1+"_"+".csv";
	
	//file_name2=file_name2+"_"+to_string(Ps0)+"_"+to_string(Pm)+"_"+to_string(d)+".csv";	
	//file_name3=file_name3+date_and_time+".csv";
	//fp3.open(file_name3.c_str());
	cout<<"Temporary files are written into "<<temp_file_name1<<endl;
	for(i=0;i<200;i++)
	{
		v.clear();
		v.assign(i+1,0);
		group_counter.push_back(v);
	}
	
	for(i=0;i<events_max;i++)
	{
		if(i%(events_max / 10) == 0)
		{
			write_into_file(temp_file_name1,group_counter,N,N1,s,Ps0,Pm,d,events_max);
			cout<<endl;
			print_vec_of_vec(m.counter);
			cout<<endl;
		}
		if(i%(events_max/1000)==0)
		{
			Z=m.get_group_count();
			cout<<"\r"<<"Progress="<<setw(log10(events_max)+1)<<i;
			cout<<'/'<<events_max<<"  (N,N1,Ps0,Pm,d)=("<<N<<" "<<N1<<" "<<Ps0<<" "<<Pm<<" "<<d<<" "<<a<<")"<<"  n=";
			cout<<setw(log10(N)+1)<<Z<<" (t="<<(clock() - start)/((double) CLOCKS_PER_SEC)<<")"<<flush;
			//fp3<<i<<","<<Z<<"\n";
		}
		if(i%sampling_interval == 0 && ((100*i)/events_max) >= 10)
		{
			for(j=0;j<m.counter.size();j++)
			{
				for(int k=0;k<m.counter[j].size();k++)
				{
					if(group_counter.size()<=j)
					{
						group_counter.resize(j+1);
					}
					if(group_counter[j].size()<=hash_counter(N,N1,j,k,true))
					{
						group_counter[j].resize(hash_counter(N,N1,j,k,true)+1);
					}
					group_counter[j][hash_counter(N,N1,j,k,true)] += m.counter[j][k];
					clean_counter(group_counter);
				}
			}
		}
		
		m.next_event();
	}
	//fp3.close();
	cout<<endl;
	clean_counter(group_counter);
	print_vec_of_vec(group_counter);
	
	vector<vector<double> > freq_counter;
	vector<double> sub;
	
	for(i=0;i<group_counter.size();i++)
	{
		sub.clear();
		sub.assign(group_counter[i].begin(),group_counter[i].end());
		freq_counter.push_back(sub);
	}
	
	double total_count=0;
	
	for(i=0;i<freq_counter.size();i++)
	{
		for(j=0;j<freq_counter[i].size();j++)
		{
			total_count=total_count+freq_counter[i][j];
		}
	}
	
	for(i=0;i<freq_counter.size();i++)
	{
		for(j=0;j<freq_counter[i].size();j++)
		{
			freq_counter[i][j]=freq_counter[i][j]/total_count;
		}
	}
	cout<<endl;
	write_into_file(file_name1,group_counter,N,N1,s,Ps0,Pm,d,events_max);
	fp2.open(file_name2.c_str());
	fp2<<"N,"<<N<<",N1,"<<N1<<",Ps0,"<<Ps0<<",d,"<<d<<",Pm,"<<Pm<<",events,"<<events_max<<",\n";
	fp2<<",";
	
	for(i=0;i<longest_group(group_counter)+1;i++)
	{
		fp2<<i<<",";
	}
	fp2<<"\n";
	for(i=0;i<group_counter.size();i++)
	{
		fp2<<i<<",";
		if(!group_counter[i].empty())
		{
			for(j=0;j<group_counter[i].size();j++)
			{
				fp2<<freq_counter[i][j]<<",";
			}
		}
		
		else
		{
			fp2<<0<<",";
		}
		fp2<<"\n";
	}
	fp2.close();
	
	cout<<"\nData written into "<<file_name1<<", frequencies written into "<<file_name2<<endl;
	cout<<"time elapsed="<<(clock() - start)/((double) CLOCKS_PER_SEC)<<endl;
	return(0);
}
