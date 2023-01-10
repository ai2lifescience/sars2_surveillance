#include <cstdio>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <thread>
#include <pthread.h>
#include <map>
#include <unordered_map>
#include <time.h>
#include <math.h>
using namespace std;

class Nucleotide//碱基类
{
public:
     string Reference;//该位置的参考基因组碱基
     int TotalDepth;//该位置的总测序深度 等于A+T+C+G
     vector<int>ATCGDepth;//ATCG各自的测序深度
     vector<int>ATCGForwardReads;//ATCG各自的正向深度
     vector<int>ATCGBackwardReads;//ATCG各自的反向深度
     string Position;
     int DELDepth;//deletion
     int DELForwardReads;
     int DELBackwardReads;
     string DEL;
     int INSDepth;//insertion
     int INSForwardReads;
     int INSBackwardReads;
     string INS;

     int MainnucCount;//主碱基reads数
};

static vector<Nucleotide> RawReadcounts;//原始Varscan Readcounts
static vector<Nucleotide> BetaCovReadcounts;//使用Kraken结果Beta冠状病毒的Varscan Readcounts
unordered_map<string, int> ATCGMap;
bool UseBetaCov = false;//是否使用Kraken结果矫正

//string split函数
vector<string> splitStr(string str, char delimiter)
{
    vector<string> r;
    string tmpstr;
    while (!str.empty()) {
        int ind = str.find_first_of(delimiter);
        if (ind == -1) {
            r.push_back(str);
            str.clear();
        }
        else {
            r.push_back(str.substr(0, ind));
            str = str.substr(ind + 1, str.size() - ind - 1);
        }
    }
    return r;
}

static string arg_i = "NA";
static string arg_k = "NA";
static double arg_t = 0.7;
static int arg_m = 5;
static double arg_sbias = 0;
static string arg_o = "TBD";
static bool arg_base_ambiguous = false;

//判断主碱基
string BasenucCaller(Nucleotide Base, double t)
{
    string ATCG = "ATCG";
    unordered_map<string,string> Nucbase_Ambinuc;
    Nucbase_Ambinuc.insert(pair<string,string>("A","A"));
    Nucbase_Ambinuc.insert(pair<string, string>("T","T"));
    Nucbase_Ambinuc.insert(pair<string, string>("C","C"));
    Nucbase_Ambinuc.insert(pair<string, string>("G","G"));
    Nucbase_Ambinuc.insert(pair<string, string>("AG","R"));
    Nucbase_Ambinuc.insert(pair<string, string>("TC","Y"));
    Nucbase_Ambinuc.insert(pair<string, string>("AC","M"));
    Nucbase_Ambinuc.insert(pair<string, string>("TG","K"));
    Nucbase_Ambinuc.insert(pair<string, string>("CG","S"));
    Nucbase_Ambinuc.insert(pair<string, string>("AT","W"));
    Nucbase_Ambinuc.insert(pair<string, string>("ATC","H"));
    Nucbase_Ambinuc.insert(pair<string, string>("TCG","B"));
    Nucbase_Ambinuc.insert(pair<string, string>("ACG","V"));
    Nucbase_Ambinuc.insert(pair<string, string>("ATG","D"));
    Nucbase_Ambinuc.insert(pair<string, string>("ATCG","N"));

    vector<int> FourNuccount;
    vector<string> FourNuc;
    if (Base.ATCGDepth[0] != 0) { FourNuc.push_back("A"); FourNuccount.push_back(Base.ATCGDepth[0]); }
    if (Base.ATCGDepth[1] != 0) { FourNuc.push_back("T"); FourNuccount.push_back(Base.ATCGDepth[1]); }
    if (Base.ATCGDepth[2] != 0) { FourNuc.push_back("C"); FourNuccount.push_back(Base.ATCGDepth[2]); }
    if (Base.ATCGDepth[3] != 0) { FourNuc.push_back("G"); FourNuccount.push_back(Base.ATCGDepth[3]); }
    int i, j, k;
    string tmp;
    for (i = 0; i < FourNuc.size(); i++)//降序排序
        for (j = i + 1; j < FourNuc.size(); j++)
        {
            if (FourNuccount[i] < FourNuccount[j])
            {
                k = FourNuccount[i];
                FourNuccount[i] = FourNuccount[j];
                FourNuccount[j] = k;

                tmp = FourNuc[i];
                FourNuc[i] = FourNuc[j];
                FourNuc[j] = tmp;
            }
        }
    double th = 0;
    string thseq = "";
    int totalnuccount = Base.TotalDepth;
    for (i = 0; i < FourNuc.size(); i++)
        if ((double)FourNuccount[i] / totalnuccount + th >= t)
        {
            thseq += FourNuc[i];
            if (FourNuc.size() - i > 1)
                for (k = i + 1; k < FourNuc.size(); k++)
                    if (FourNuccount[k] == FourNuccount[i])
                        thseq += FourNuc[k];
            break;
        }
        else { thseq += FourNuc[i]; th += (double)FourNuccount[i] / totalnuccount; }

    string ATCGseq = "";
    for (i = 0; i < ATCG.size(); i++)
        if (thseq.find(ATCG[i])!=-1) ATCGseq += ATCG[i];

    //if (ATCGseq.Length > 1)
    //    Console.Write("gacha");
    return Nucbase_Ambinuc[ATCGseq];
}
//call共识序列
void CallCSS()//mut5col cssseq
{
    string seqout = arg_o;
    string mutout = arg_o + ".mut5col";
    ofstream writeseq(seqout);
    ofstream writemut(mutout);
    ATCGMap.insert(pair<string, int>("A", 0)); ATCGMap.insert(pair<string, int>("T", 1));
    ATCGMap.insert(pair<string, int>("C", 2)); ATCGMap.insert(pair<string, int>("G", 3));
    vector<string> name1 = splitStr(arg_i, '/');
    string name = name1[name1.size() - 1];
    int i, j, k;
    string outputseq = "";
    for (i = 1; i < RawReadcounts.size(); i++)
    {
        if (RawReadcounts[i].TotalDepth == 0)outputseq += "N";
        else
        {
            if (RawReadcounts[i].ATCGDepth[ATCGMap[RawReadcounts[i].Reference]] == RawReadcounts[i].TotalDepth)//是否存在多态性
                outputseq += RawReadcounts[i].Reference;
            else
            {
                if (UseBetaCov)
                {
                    if (BetaCovReadcounts[i].TotalDepth < arg_m)
                        outputseq += "N";
                    else
                    {
                        if (BetaCovReadcounts[i].DELDepth >= BetaCovReadcounts[i].TotalDepth * arg_t)//Deletion
                        {
                            vector<string> linz = splitStr(BetaCovReadcounts[i].DEL, '-');
                            string mut5colout = name;
                            mut5colout += "\t";
                            mut5colout += to_string(i);
                            mut5colout += "\t";
                            mut5colout += to_string(i+stoi(linz[1]));
                            mut5colout += "\t";
                            mut5colout += RawReadcounts[i].Reference + linz[2];
                            mut5colout += "\t";
                            mut5colout += RawReadcounts[i].Reference;
                            mut5colout += "\t";
                            mut5colout += name;
                            writemut << mut5colout << endl;
                            outputseq += RawReadcounts[i].Reference;
                            i += stoi(linz[1]);
                        }
                        else if (BetaCovReadcounts[i].INSDepth >= BetaCovReadcounts[i].TotalDepth * arg_t)//Insertion
                        {
                            vector<string> linz = splitStr(BetaCovReadcounts[i].INS, '-');
                            string mut5colout = name;
                            mut5colout += "\t";
                            mut5colout += to_string(i);
                            mut5colout += "\t";
                            mut5colout += to_string(i);
                            mut5colout += "\t";
                            mut5colout += RawReadcounts[i].Reference;
                            mut5colout += "\t";
                            mut5colout += RawReadcounts[i].Reference + linz[2];
                            mut5colout += "\t";
                            mut5colout += name;
                            writemut << mut5colout << endl;
                            outputseq += RawReadcounts[i].Reference + linz[2];
                        }
                        else//SNP
                        {
                            string RYMKSWHBVDN = "RYMKSWHBVDN";
                            string Mainnuc = BasenucCaller(BetaCovReadcounts[i], arg_t);
                            if (arg_base_ambiguous == false && RYMKSWHBVDN.find(Mainnuc) != -1)
                                outputseq += "N";
                            else
                            {
                                if (Mainnuc == BetaCovReadcounts[i].Reference)
                                    outputseq += Mainnuc;
                                else
                                {
                                    double sbias = (double)BetaCovReadcounts[i].ATCGForwardReads[ATCGMap[Mainnuc]] / (BetaCovReadcounts[i].ATCGBackwardReads[ATCGMap[Mainnuc]] + BetaCovReadcounts[i].ATCGForwardReads[ATCGMap[Mainnuc]]);
                                    if (sbias >= arg_sbias && sbias <= (1 - arg_sbias))
                                    {
                                        outputseq += Mainnuc;
                                        string mut5colout = name;
                                        mut5colout += "\t";
                                        mut5colout += to_string(i);
                                        mut5colout += "\t";
                                        mut5colout += to_string(i);
                                        mut5colout += "\t";
                                        mut5colout += RawReadcounts[i].Reference;
                                        mut5colout += "\t";
                                        mut5colout += Mainnuc;
                                        mut5colout += "\t";
                                        mut5colout += name;
                                        writemut << mut5colout << endl;
                                    }
                                    else
                                        outputseq += "N";
                                }
                            }
                        }
                    }
                }
                else
                {
                    
                    if (RawReadcounts[i].TotalDepth < arg_m)
                        outputseq += "N";
                    else
                    {
                        if (RawReadcounts[i].DELDepth >= RawReadcounts[i].TotalDepth * arg_t)//Deletion
                        {
                            vector<string> linz = splitStr(RawReadcounts[i].DEL, '-');
                            string mut5colout = name;
                            mut5colout += "\t";
                            mut5colout += to_string(i);
                            mut5colout += "\t";
                            mut5colout += to_string(i + stoi(linz[1]));
                            mut5colout += "\t";
                            mut5colout += RawReadcounts[i].Reference + linz[2];
                            mut5colout += "\t";
                            mut5colout += RawReadcounts[i].Reference;
                            mut5colout += "\t";
                            mut5colout += name;
                            writemut << mut5colout << endl;
                            outputseq += RawReadcounts[i].Reference;
                            i += stoi(linz[1]);
                        }
                        else if (RawReadcounts[i].INSDepth >= RawReadcounts[i].TotalDepth * arg_t)//Insertion
                        {
                            vector<string> linz = splitStr(RawReadcounts[i].INS, '-');
                            string mut5colout = name;
                            mut5colout += "\t";
                            mut5colout += to_string(i);
                            mut5colout += "\t";
                            mut5colout += to_string(i);
                            mut5colout += "\t";
                            mut5colout += RawReadcounts[i].Reference;
                            mut5colout += "\t";
                            mut5colout += RawReadcounts[i].Reference + linz[2];
                            mut5colout += "\t";
                            mut5colout += name;
                            writemut << mut5colout << endl;
                            outputseq += RawReadcounts[i].Reference + linz[2];
                        }
                        else
                        {
                            string RYMKSWHBVDN = "RYMKSWHBVDN";
                            string Mainnuc = BasenucCaller(RawReadcounts[i], arg_t);
                            if (arg_base_ambiguous == false && RYMKSWHBVDN.find(Mainnuc) != -1)
                                outputseq += "N";
                            else
                            {
                                if (Mainnuc == RawReadcounts[i].Reference)
                                    outputseq += Mainnuc;
                                else
                                {
                                    double sbias = (double)RawReadcounts[i].ATCGForwardReads[ATCGMap[Mainnuc]] / (RawReadcounts[i].ATCGBackwardReads[ATCGMap[Mainnuc]] + RawReadcounts[i].ATCGForwardReads[ATCGMap[Mainnuc]]);
                                    if (sbias >= arg_sbias && sbias <= (1 - arg_sbias))
                                    {
                                        outputseq += Mainnuc;
                                        string mut5colout = name;
                                        mut5colout += "\t";
                                        mut5colout += to_string(i);
                                        mut5colout += "\t";
                                        mut5colout += to_string(i);
                                        mut5colout += "\t";
                                        mut5colout += RawReadcounts[i].Reference;
                                        mut5colout += "\t";
                                        mut5colout += Mainnuc;
                                        mut5colout += "\t";
                                        mut5colout += name;
                                        writemut << mut5colout << endl;
                                    }
                                    else
                                        outputseq += "N";
                                }
                            }
                        }
                    }
                }
            }
        }
        
    }
    writeseq << ">" << name << endl;
    writeseq << outputseq << endl;
    writemut.close();
    writeseq.close();
}


//读入数据
void ReadData()
{
    int i, j, k;
    ifstream readraw(arg_i);
    if (!readraw.is_open())
    {
        printf("ERROR /// CCS: Can Not Open File: %s /// \n", arg_i.c_str());
        return;
    }
    Nucleotide q;
    RawReadcounts.push_back(q);
    BetaCovReadcounts.push_back(q);
    string line;
    getline(readraw, line);
    while (getline(readraw, line))
    {
        Nucleotide a;
        vector<string> line1 = splitStr(line, '\t');
        a.Position = line1[1];
        a.Reference = line1[2];
        a.TotalDepth = 0;
        a.ATCGDepth.push_back(0); a.ATCGDepth.push_back(0); a.ATCGDepth.push_back(0); a.ATCGDepth.push_back(0);
        a.ATCGForwardReads.push_back(0); a.ATCGForwardReads.push_back(0); a.ATCGForwardReads.push_back(0); a.ATCGForwardReads.push_back(0);
        a.ATCGBackwardReads.push_back(0); a.ATCGBackwardReads.push_back(0); a.ATCGBackwardReads.push_back(0); a.ATCGBackwardReads.push_back(0);
        a.DELDepth = 0; a.DELForwardReads = 0; a.DELBackwardReads = 0;
        a.INSDepth = 0; a.INSForwardReads = 0; a.INSBackwardReads = 0;
        for (j = 5; j < line1.size(); j++)
        {
            vector<string> linemc = splitStr(line1[j], ':');
            if (j == 5)a.MainnucCount = stoi(linemc[1]);
            if (linemc[0] == "A")
            {
                a.ATCGDepth[0] = stoi(linemc[1]); a.ATCGForwardReads[0] = stoi(linemc[5]); a.ATCGBackwardReads[0] = stoi(linemc[6]);
            }
            if (linemc[0] == "T")
            {
                a.ATCGDepth[1] = stoi(linemc[1]); a.ATCGForwardReads[1] = stoi(linemc[5]); a.ATCGBackwardReads[1] = stoi(linemc[6]);
            }
            if (linemc[0] == "C")
            {
                a.ATCGDepth[2] = stoi(linemc[1]); a.ATCGForwardReads[2] = stoi(linemc[5]); a.ATCGBackwardReads[2] = stoi(linemc[6]);
            }
            if (linemc[0] == "G")
            {
                a.ATCGDepth[3] = stoi(linemc[1]); a.ATCGForwardReads[3] = stoi(linemc[5]); a.ATCGBackwardReads[3] = stoi(linemc[6]);
            }
            if (linemc[0][0] == 'D' && a.DELDepth < stoi(linemc[1]))
            {
                a.DEL = linemc[0]; a.DELDepth = stoi(linemc[1]); a.DELForwardReads = stoi(linemc[5]); a.DELBackwardReads = stoi(linemc[6]);
            }
            if (linemc[0][0] == 'I' && a.INSDepth < stoi(linemc[1]))
            {
                a.INS = linemc[0]; a.INSDepth = stoi(linemc[1]); a.INSForwardReads = stoi(linemc[5]); a.INSBackwardReads = stoi(linemc[6]);
            }
        }
        a.TotalDepth = a.ATCGDepth[0] + a.ATCGDepth[1] + a.ATCGDepth[2] + a.ATCGDepth[3];
        a.TotalDepth += a.DELDepth + a.INSDepth;
        RawReadcounts.push_back(a);
    }
    readraw.close();
    if (arg_k != "NA")
    {
        ifstream readbeta(arg_k);
        if (!readbeta.is_open())
        {
            printf("ERROR /// CCS: Can Not Open File: %s /// \n", arg_k.c_str());
            return;
        }
        UseBetaCov = true;
        getline(readbeta, line);
        while (getline(readbeta, line))
        {
            Nucleotide a;
            vector<string> line1 = splitStr(line, '\t');
            a.Reference = line1[2];
            a.TotalDepth = 0;
            a.ATCGDepth.push_back(0); a.ATCGDepth.push_back(0); a.ATCGDepth.push_back(0); a.ATCGDepth.push_back(0);
            a.ATCGForwardReads.push_back(0); a.ATCGForwardReads.push_back(0); a.ATCGForwardReads.push_back(0); a.ATCGForwardReads.push_back(0);
            a.ATCGBackwardReads.push_back(0); a.ATCGBackwardReads.push_back(0); a.ATCGBackwardReads.push_back(0); a.ATCGBackwardReads.push_back(0);
            a.DELDepth = 0; a.DELForwardReads = 0; a.DELBackwardReads = 0;
            a.INSDepth = 0; a.INSForwardReads = 0; a.INSBackwardReads = 0;
            for (j = 5; j < line1.size(); j++)
            {
                vector<string> linemc = splitStr(line1[j], ':');
                if(j==5)a.MainnucCount = stoi(linemc[1]);
                if (linemc[0] == "A")
                {
                    a.ATCGDepth[0] = stoi(linemc[1]); a.ATCGForwardReads[0] = stoi(linemc[5]); a.ATCGBackwardReads[0] = stoi(linemc[6]);
                }
                if (linemc[0] == "T")
                {
                    a.ATCGDepth[1] = stoi(linemc[1]); a.ATCGForwardReads[1] = stoi(linemc[5]); a.ATCGBackwardReads[1] = stoi(linemc[6]);
                }
                if (linemc[0] == "C")
                {
                    a.ATCGDepth[2] = stoi(linemc[1]); a.ATCGForwardReads[2] = stoi(linemc[5]); a.ATCGBackwardReads[2] = stoi(linemc[6]);
                }
                if (linemc[0] == "G")
                {
                    a.ATCGDepth[3] = stoi(linemc[1]); a.ATCGForwardReads[3] = stoi(linemc[5]); a.ATCGBackwardReads[3] = stoi(linemc[6]);
                }
                if (linemc[0][0] == 'D' && a.DELDepth == 0)
                {
                    a.DEL = linemc[0]; a.DELDepth = stoi(linemc[1]); a.DELForwardReads = stoi(linemc[5]); a.DELBackwardReads = stoi(linemc[6]);
                }
                if (linemc[0][0] == 'I' && a.INSDepth == 0)
                {
                    a.INS = linemc[0]; a.INSDepth = stoi(linemc[1]); a.INSForwardReads = stoi(linemc[5]); a.INSBackwardReads = stoi(linemc[6]);
                }
            }
            a.TotalDepth = a.ATCGDepth[0] + a.ATCGDepth[1] + a.ATCGDepth[2] + a.ATCGDepth[3];
            a.TotalDepth += a.DELDepth + a.INSDepth;
            a.MainnucCount;
            BetaCovReadcounts.push_back(a);
        }
        readbeta.close();
    }
}
//5个输入 
/*
1.-i readcounts文件 
2.-o output文件，默认【input.fa】
3.-k kraken readcounts文件（可选）
4.-t Basenuc 主碱基需要达到的最低频率，默认【0.7】
5.-m Make a call的碱基最低覆盖reads数，默认【5】
6.-Sbias 链偏好性过滤参数，默认【0.1】
7.-BaseAmbiguous 是否保留兼并碱基，默认【false N替换】
*/
int main(int argc, char* argv[])
{
    printf("CallConsensusSequence v0.1.1:\n");
    int i, j, k;
    if (argc > 2)
    {
        for(i=1;i<argc;i+=2)
        {
            string a = argv[i];
            if (a == "-i")arg_i = argv[i + 1];
            if (a == "-o")arg_o = argv[i + 1];
            if (a == "-k")arg_k = argv[i + 1];
            if (a == "-t")arg_t = stof(argv[i + 1]);
            if (a == "-m")arg_m = stoi(argv[i + 1]);
            if (a == "-sbias")arg_sbias = stof(argv[i + 1]);
            if (a == "-baseAmbiguous" && argv[i+1] == "true")arg_base_ambiguous = true;
        }
        if (arg_o == "TBD")arg_o = arg_i + ".CSS.fa";

        printf("START\n");
        ReadData();//读入数据
        printf("Readin Finished\n");
        CallCSS();//call共识序列

        printf("Done.\n");
        time_t now = time(0);
        cout << ctime(&now) << endl;
    }
    else
    {
        printf("1 Input Is Needed and 5 Parameter Can Be Specified\n");
        printf("--------------------------------------------------------\n");
        printf("1. -i  Target readcounts file \n");
        printf("2. -o  Outfile dictionary/file (optional, same as input)\n");
        printf("3. -k  Kraken readcounts file (optional)\n");
        printf("4. -t  Basenuc mainnuc minimum freq (optional, 0.7)\n");
        printf("5. -m  Minimum position depth to make a call (optional, 5)\n");
        printf("6. -sbias  (optional, 0)\n");
        printf("7. -baseAmbiguous  Whether keep ambiguous bases: true / false (optional, false, replaced by N)\n");
        printf("--------------------------------------------------------\n\n");
    }
    return 0;
}