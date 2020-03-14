#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <vector>
#include <stack>
#include <queue>
#include <fstream>
#include <sstream>
#include <iostream>
#include <math.h>
#include <random>

using namespace std;

typedef chrono::high_resolution_clock myclock;
typedef struct node
{
	int id;
	int parentID;
	int branchID;
	double sol;
	double lb;
	int solStatus;		// Optimal: 1, Relax: 2, Infeasible: 3
	int depth;
	node(int d, int parent, int branch, int status, double s, double b, int dep) :
		id(d), parentID(parent), branchID(branch), solStatus(status), sol(s), lb(b), depth(dep) {}
} node;

typedef struct tree
{
    int depth;
    vector<double> p;
    vector<double> q;
    vector<double> r;
    tree (int d) : depth(d) {}
} tree;

void insertNewNode(queue<node*>& nodeQ, int d, int parent, int branch, int status, double s, double b, int dep)
{
	node* newNode = new node(d, parent, branch, status, s, b, dep);
	nodeQ.push(newNode);
}

vector<vector<double>> getProbs(char* probsname)
{
	vector<vector<double>> res;
	ifstream inFile(probsname);
	double p,q,r;
	while (inFile >> p >> q >> r) {
		res[0].push_back(p);
		res[1].push_back(q);
		res[2].push_back(r);
	}
	return res;
}

void genProb(char* pathname, char* filename)
{
	queue<node*> nodesToBranch;
	//myclock::time_point beginning = myclock::now();
	// Since depth starts at 0 here, the depth should add one when counting
	int depth = 14;
	int optShowDepth = depth;
	int smallEndDepth = depth;
	int bigShowDepth = 1;

	double globLB = 15;
	double optVal = 20;
	double globUB = 25;
	double optProbBase = 0.05;
	double optProbConst = 1.1;
	double smallProbBase = 0;
	double smallProbConst = 0.8;
	double bigProbBase = 0;
	double bigProbConst = 0.3;
	//double tol = 1/pow(2,depth);

	vector<double> optProbVec;
	vector<double> smallProbVec;
	vector<double> bigProbVec;

	//int probBase = 1000;
	int nodeCount = 0;
	double status;
	double probOpt, probSmall, probBig;
	int optCount = 0;
	//myclock::duration d;

	// Fill probability array P
	optProbVec.resize(depth+1);
	for (int i = 0; i <= depth; i++) {
        if (i < optShowDepth)
            optProbVec[i] = 0;
        else {
            optProbVec[i] = optProbBase * pow(optProbConst,i-optShowDepth);
        }
	}
	// Fill probability array Q
	smallProbVec.resize(depth+1);
	for (int i = 0; i <= depth; i++) {
        if (i >= smallEndDepth)
            smallProbVec[i] = 0;
        else {
            smallProbVec[i] = smallProbBase + pow(smallProbConst,i);
            if (smallProbVec[i] > 1)
                smallProbVec[i] = 1;
            //else if (smallProbVec[i] < tol)
                //smallProbVec[i] = 0;
        }
	}
	// Fill probability array R
	bigProbVec.resize(depth+1);
	for (int i = 0; i <= depth; i++) {
        if (i < bigShowDepth)
            bigProbVec[i] = 0;
        else {
            bigProbVec[i] = bigProbBase + pow(bigProbConst,depth - i);
            if (bigProbVec[i] > 1)
                bigProbVec[i] = 1;
            //else if (bigProbVec[i] < tol)
                //bigProbVec[i] = 0;
        }
	}

	char path[255];
	sprintf(path, "%s/%s",pathname,filename);
	FILE* probFile = fopen(path, "w");
    fprintf(probFile,"%s\n",filename);
    for (int i = 0; i <= depth; i++) {
        fprintf(probFile,"%f %f %f\n", optProbVec[i], smallProbVec[i], bigProbVec[i]);
    }
    //fprintf(probFile,"%d %d %f %f\n",depth,optShowDepth,optProbBase / 1000,optProbConst / 1000);

	node* curNode = new node(0, 0, 0, 2, globLB, globLB, 0);
	nodesToBranch.push(curNode);

    //d = myclock::now() - beginning;
    //mt19937 gen(rd());
    uniform_real_distribution<> dist(0, 1);
    default_random_engine gen{static_cast<long unsigned int>(myclock::now().time_since_epoch().count())};
    //srand(myclock::now().time_since_epoch().count());

	while (!nodesToBranch.empty()) {
		curNode = nodesToBranch.front();
		nodesToBranch.pop();

		// Get probability of different kind.
		probOpt = optProbVec[curNode->depth];
		probSmall = smallProbVec[curNode->depth];
		probBig = bigProbVec[curNode->depth];

		if (curNode->depth == depth) {
            if (curNode->lb < optVal) {
                //status = rand() % probBase;
                curNode->solStatus = 1;
                curNode->sol = globUB + 1;      // Make sure no pruning can be done with this
            } else if (curNode->lb == optVal) {
                //status = rand() % probBase;
                status = dist(gen);
                curNode->solStatus = 1;
                curNode->sol = (status < probOpt) ? optVal : globUB + 1;
                if (curNode->sol == optVal)
                    optCount++;
            } else {
                curNode->solStatus = 1;
                curNode->sol = globUB + 1;
            }
		} else if (curNode->depth == 0) {
            insertNewNode(nodesToBranch, ++nodeCount, curNode->id, 1, -1, -1, curNode->sol, curNode->depth + 1);
            insertNewNode(nodesToBranch, ++nodeCount, curNode->id, 2, -1, -1, curNode->sol, curNode->depth + 1);
		} else if (curNode->lb < optVal) {
            // First the results when the node
            // No Optimal solution can come from node with H < Z^*
            curNode->solStatus = 2;
            //status = rand() % probBase;
            status = dist(gen);
            if (probSmall > 0 && status < probSmall) {
                // Current node still has
                curNode->sol = globLB;
                insertNewNode(nodesToBranch, ++nodeCount, curNode->id, 1, -1, -1, curNode->sol, curNode->depth + 1);
                insertNewNode(nodesToBranch, ++nodeCount, curNode->id, 2, -1, -1, curNode->sol, curNode->depth + 1);
            } else if (probBig > 0 && status < (probSmall + probBig)) {
                curNode->sol = globUB;
                insertNewNode(nodesToBranch, ++nodeCount, curNode->id, 1, -1, -1, curNode->sol, curNode->depth + 1);
                insertNewNode(nodesToBranch, ++nodeCount, curNode->id, 2, -1, -1, curNode->sol, curNode->depth + 1);
            } else {
                curNode->sol = optVal;
                insertNewNode(nodesToBranch, ++nodeCount, curNode->id, 1, -1, -1, curNode->sol, curNode->depth + 1);
                insertNewNode(nodesToBranch, ++nodeCount, curNode->id, 2, -1, -1, curNode->sol, curNode->depth + 1);
            }
		} else if (curNode->lb == optVal) {
		    //status = rand() % probBase;
		    status = dist(gen);
		    if (probOpt > 0 && status < probOpt) {
                curNode->solStatus = 1;
                curNode->sol = optVal;
		    } else if (probBig > 0 && status < (probOpt + probBig)) {
		        curNode->solStatus = 2;
		        curNode->sol = globUB;
                insertNewNode(nodesToBranch, ++nodeCount, curNode->id, 1, -1, -1, curNode->sol, curNode->depth + 1);
                insertNewNode(nodesToBranch, ++nodeCount, curNode->id, 2, -1, -1, curNode->sol, curNode->depth + 1);
		    } else {
		        curNode->solStatus = 2;
		        curNode->sol = optVal;
                insertNewNode(nodesToBranch, ++nodeCount, curNode->id, 1, -1, -1, curNode->sol, curNode->depth + 1);
                insertNewNode(nodesToBranch, ++nodeCount, curNode->id, 2, -1, -1, curNode->sol, curNode->depth + 1);
		    }
		} else {
		    curNode->solStatus = 2;
		    curNode->sol = globUB;
            insertNewNode(nodesToBranch, ++nodeCount, curNode->id, 1, -1, -1, curNode->sol, curNode->depth + 1);
            insertNewNode(nodesToBranch, ++nodeCount, curNode->id, 2, -1, -1, curNode->sol, curNode->depth + 1);
		}

		fprintf(probFile, "%d %d %d %d %f\n",
			curNode->id, curNode->parentID, curNode->branchID, curNode->solStatus, curNode->sol);

		delete curNode;
	}
    if (optCount == 0)
        cout << "Opt sol not generated." << endl;
	fclose(probFile);
}

void genProb(tree* t, char* pathname, char* filename)
{
	queue<node*> nodesToBranch;
	// Since depth starts at 0 here, the depth should add one when counting
	int depth = t->depth;

	double globLB = 15;
	double optVal = 20;
	double globUB = 25;
    int nodeCount = 0;
	double status;
	double probOpt, probSmall, probBig;
	int optCount = 0;

	char path[255];
	sprintf(path, "%s/%s",pathname,filename);
	FILE* probFile = fopen(path, "w");
    fprintf(probFile,"%s\n",filename);
    fprintf(probFile,"%d\n",depth);
    for (int i = 0; i <= depth; i++) {
        fprintf(probFile,"%f %f %f\n", t->p[i], t->q[i], t->r[i]);
    }
    //fprintf(probFile,"%d %d %f %f\n",depth,optShowDepth,optProbBase / 1000,optProbConst / 1000);

	node* curNode = new node(0, 0, 0, 2, globLB, globLB, 0);
	nodesToBranch.push(curNode);

    uniform_real_distribution<> dist(0, 1);
    default_random_engine gen{static_cast<long unsigned int>(myclock::now().time_since_epoch().count())};

	while (!nodesToBranch.empty()) {
		curNode = nodesToBranch.front();
		nodesToBranch.pop();

		// Get probability of relaxation solution
		probOpt = t->p[curNode->depth];
		probSmall = t->q[curNode->depth];
		probBig = t->r[curNode->depth];

		if (curNode->depth == depth) {
            if (curNode->lb < optVal) {
                //status = rand() % probBase;
                curNode->solStatus = 1;
                curNode->sol = globUB + 1;      // Make sure no pruning can be done with this
            } else if (curNode->lb == optVal) {
                //status = rand() % probBase;
                status = dist(gen);
                curNode->solStatus = 1;
                curNode->sol = (status < probOpt) ? optVal : globUB + 1;
                if (curNode->sol == optVal)
                    optCount++;
            } else {
                curNode->solStatus = 1;
                curNode->sol = globUB + 1;
            }
		} else if (curNode->depth == 0) {
            insertNewNode(nodesToBranch, ++nodeCount, curNode->id, 1, -1, -1, curNode->sol, curNode->depth + 1);
            insertNewNode(nodesToBranch, ++nodeCount, curNode->id, 2, -1, -1, curNode->sol, curNode->depth + 1);
		} else if (curNode->lb < optVal) {
            // First the results when the node
            // No Optimal solution can come from node with H < Z^*
            curNode->solStatus = 2;
            //status = rand() % probBase;
            status = dist(gen);
            if (probSmall > 0 && status < probSmall) {
                // Current node still has
                curNode->sol = globLB;
                insertNewNode(nodesToBranch, ++nodeCount, curNode->id, 1, -1, -1, curNode->sol, curNode->depth + 1);
                insertNewNode(nodesToBranch, ++nodeCount, curNode->id, 2, -1, -1, curNode->sol, curNode->depth + 1);
            } else if (probBig > 0 && status < (probSmall + probBig)) {
                curNode->sol = globUB;
                insertNewNode(nodesToBranch, ++nodeCount, curNode->id, 1, -1, -1, curNode->sol, curNode->depth + 1);
                insertNewNode(nodesToBranch, ++nodeCount, curNode->id, 2, -1, -1, curNode->sol, curNode->depth + 1);
            } else {
                curNode->sol = optVal;
                insertNewNode(nodesToBranch, ++nodeCount, curNode->id, 1, -1, -1, curNode->sol, curNode->depth + 1);
                insertNewNode(nodesToBranch, ++nodeCount, curNode->id, 2, -1, -1, curNode->sol, curNode->depth + 1);
            }
		} else if (curNode->lb == optVal) {
		    //status = rand() % probBase;
		    status = dist(gen);
		    if (probOpt > 0 && status < probOpt) {
                curNode->solStatus = 1;
                curNode->sol = optVal;
		    } else if (probBig > 0 && status < (probOpt + probBig)) {
		        curNode->solStatus = 2;
		        curNode->sol = globUB;
                insertNewNode(nodesToBranch, ++nodeCount, curNode->id, 1, -1, -1, curNode->sol, curNode->depth + 1);
                insertNewNode(nodesToBranch, ++nodeCount, curNode->id, 2, -1, -1, curNode->sol, curNode->depth + 1);
		    } else {
		        curNode->solStatus = 2;
		        curNode->sol = optVal;
                insertNewNode(nodesToBranch, ++nodeCount, curNode->id, 1, -1, -1, curNode->sol, curNode->depth + 1);
                insertNewNode(nodesToBranch, ++nodeCount, curNode->id, 2, -1, -1, curNode->sol, curNode->depth + 1);
		    }
		} else {
		    curNode->solStatus = 2;
		    curNode->sol = globUB;
            insertNewNode(nodesToBranch, ++nodeCount, curNode->id, 1, -1, -1, curNode->sol, curNode->depth + 1);
            insertNewNode(nodesToBranch, ++nodeCount, curNode->id, 2, -1, -1, curNode->sol, curNode->depth + 1);
		}

		fprintf(probFile, "%d %d %d %d %f\n",
			curNode->id, curNode->parentID, curNode->branchID, curNode->solStatus, curNode->sol);

		delete curNode;
	}
    if (optCount == 0)
        cout << "Opt sol not generated." << endl;
	fclose(probFile);
}

int main()
{
    //char* path = "D:/ResearchTests/Theory/Simulation/D15P005S08B03";
    char pname[255];
    char fname[255];
    //char file[255];
    int numInst = 1000;

    int depth = 11; // Actual depth - 1
	int optShowDepth = depth;
	int smallEndDepth = depth;
	int bigShowDepth = 1;
	double optProb, smallProb, bigProb;
	double optProbMulti = 1.1;
	double smallProbBase = 0;
	double bigProbBase = 0;
	double optProbConst = 0.02;
	double smallProbConst = 0.7;
	double bigProbConst = 0.3;

    tree* tprob = new tree(depth);
    tprob->p.resize(depth+1);
    tprob->q.resize(depth+1);
    tprob->r.resize(depth+1);
    for (int op = 5; op <= 5; op++) {
        for (int sp = 50; sp <= 50; sp+=5) {
            for (int bp = 72; bp <= 78; bp+=2) {
                optProb = double(op)/100;
                smallProb = double(sp)/100;
                if (bp >= 10)
                    bigProb = double(bp)/100;
                else
                    bigProb = double(bp)/10;
                // Optimal solution probability
                for (int i = 0; i <= depth; i++) {
                    if (i < optShowDepth)
                        tprob->p[i] = 0;
                    else
                        tprob->p[i] = optProb * pow(optProbMulti,i-optShowDepth);
                }
                // Relaxation less than optimal solution probability
                for (int i = 0; i <= depth; i++) {
                    if (i >= smallEndDepth)
                        tprob->q[i] = 0;
                    else {
                        tprob->q[i] = smallProbBase + pow(smallProb,i);
                        if (tprob->q[i] > 1)
                            tprob->q[i] = 1;
                    }
                }
                // Relaxation greater than optimal solution probability
                for (int i = 0; i <= depth; i++) {
                    if (i < bigShowDepth)
                        tprob->r[i] = 0;
                    else {
                        tprob->r[i] = bigProbBase + pow(bigProb,depth - i);
                        if (tprob->r[i] > 1)
                            tprob->r[i] = 1;
                    }
                }
                sprintf(pname,"D:/ResearchTests/Theory/Simulation/D%dB0%d/D%dP00%dS0%dB0%d", depth+1, bp, depth+1,op,sp,bp);
				//sprintf(pname, "D:/ResearchTests/Theory/Simulation/D%dS0%d/D%dP00%dS0%dB0%d", depth + 1, sp, depth + 1, op, sp, bp);
                cout << pname << endl;
                for (int i = 0; i < numInst; i++) {
                    sprintf(fname,"inst%d",i);
                    genProb(tprob, pname, fname);
                }
            }
        }
    }
    delete tprob;

//    for (int i = 0; i < numInst; i++) {
//        sprintf(fname,"inst%d",i);
//        //sprintf(file,"%s/%s",path,fname);
//        genProb(path,fname);
//    }
	return 0;
}
