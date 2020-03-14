#pragma once
// Header for the branch-and-bound algorithm

#include "util.h"
#include <map>
#include <stack>
#include <queue>
#include <list>

using namespace std;

class bnbProbData;
class bnbData;
class bnbNode;
class bnbSolve;
class bnbSearch;
class bnbBranch;
class BranchAndBound;

typedef chrono::high_resolution_clock myclock;
typedef multimap<double, bnbNode*> NodesMap;
typedef map<int, NodesMap> ContourMap;

typedef struct node
{
	int id;
	int parentID;
	int branchID;
	double sol;
	int solStatus;		// Optimal: 1, Feasible: 2, Infeasible: 3
	node(int d, int parent, int branch, int status, double s) :
		id(d), parentID(parent), branchID(branch), solStatus(status), sol(s) {}
} node;

typedef struct options
{
	Mode md;
	TBMode tb;
	int mesrBest;
	int solveMode;
	double timeLimit;
	int iterLimit;
	int lweight, rweight;
	options(double time, int iter, Mode m, TBMode t, int mesr, int sm, int lw, int rw) : 
		timeLimit(time), iterLimit(iter), md(m), tb(t), mesrBest(mesr), solveMode(sm), lweight(lw), rweight(rw) {}
} options;

/* Store entire tree for problem instance */
class bnbProbData
{
public:
	bnbProbData() {}
	~bnbProbData();
	void buildTree(const char* filename);
	node* getNode(int id) { return mNodes[id]; }
	vector<node*> getChildNodes(int parentID);
	string mInstName;
	int Depth, OptStartDepth;
	double OptProbBase, OptProbConst;
private:
	vector<node*> mNodes;
	map<int, vector<node*>> mTree;
};

/* Store search tree and general algorithm info */
class bnbData
{
public:
	bnbData() {}
	~bnbData();
	// create search tree and store relavent data
	void initialize(options opt);
	// operations on search tree
	void addNode(bnbNode* node);
	void delNode(bnbNode* node);
	bnbNode* getNode(int nodeID);
	void dumpAllNodes();
	// operations on global info
	double updateLb();
	double updateUb(double sol);

	void updateTraceSearch(double lb);
	vector<int>& getTraceSearch() { return mTraceSearch; }
	// operations for accessing private variables
	int getNumUnexpNodes() { return mNumUnexpNodes; }
	bool hasUnexpNodes() { return (mNumUnexpNodes > 0); }
	Mode getSearchMode() { return mSearchMode; }
	TBMode getContSearchMode() { return mContSearchMode; }
	int getMesrBest() { return mMesrBest; }
	int getSolveMode() { return mSolveMode; }
	double getTimeLimit() { return mTimeLimit; }
	int getIterLimit() { return mIterLimit; }
	int getLweight() { return mLweight; }
	int getRweight() { return mRweight; }

	double globLB, globUB;
	int mNumIter, mNumNodes;
	long mElapsTime;
	bnbProbData* mBnBProbData;
	bnbSearch* mBnBSearch;
	bnbBranch* mBnBBranch;
	bnbSolve* mBnBSolve;
private:
	map<int, bnbNode*> mUnexpNodesHeap;
	map<double, int> mLbCount;
	vector<int> mTraceSearch;
	int mNumUnexpNodes;
	Mode mSearchMode;
	TBMode mContSearchMode;
	int mMesrBest;
	int mSolveMode;
	double mTimeLimit;
	int mIterLimit;
	int mLweight, mRweight;
	int mCurCycle;
};

/* Determine which node to explore next */
class bnbSearch
{
public:
	bnbSearch() {}
	void initialize(bnbData* data);
	// choices of search strategy
	int getNextNode();
	//int getNextNodeBFS();
	int getNextNodeDFS();
	int getNextNodeCBFS();
	int calContour(bnbNode* node);
	// operations on auxilary tree
	void addNode(bnbNode* node);
	void delNode(bnbNode* node);
	void setIterator() { mCurContour = mContours.begin(); }
	
	int getNodeInLevel(int d, int i) { return mNodeInLevels[d][i]; }
private:
	int mLweight, mRweight;
	int mMesrBest;
	Mode mMode;
	TBMode mTBMode;
	bnbData* mData;
	ContourMap mContours;
	ContourMap::iterator mCurContour;
	stack<bnbNode*> mNodesStack;
	vector<vector<int>> mNodeInLevels;
};

/* Generate the child nodes */
class bnbBranch
{
public:
	bnbBranch() {}
	void initialize(bnbData* data);
	// operations on generating new nodes
	int createNewNodes(bnbNode* node);			// create nodes from curNode
	bnbNode* getNewNode();						// get the created nodes one at a time to insert in tree
	int numNodes;
private:
	int mSolveMode;
	bnbData* mData;
	queue<bnbNode*> mNewNodes;
};

class bnbSolve
{
public:
	bnbSolve() {}
	void initialize(bnbData* data);
	int solveNode(bnbNode* node);
private:
	int mSolveMode;
	bnbData* mData;
};

/* Store info on individual nodes */
class bnbNode
{
public:
	bnbNode() {}
	bnbNode(bnbNode* parentNode);			// Creating child nodes
	bnbNode(bnbData* data, node* nd);
	bnbNode(bnbNode* parentNode, node* nd);						// Creating node from a node struct
	//bnbNode(bnbData* problem);			// Creating root node
	~bnbNode();
	// node info
	int mNodeID, mParentID;
	int mContour, mDepth;
	int mNumLBrch, mNumRBrch;
	double mLBound, mRexSol;
	int mSolStatus;
	bool mIsOptimal;
private:
	bnbData* mData;
};

/* Control solving process */
class BranchAndBound
{
public:
	BranchAndBound(const char* filename, char* jsonFile, char* levelFile, options opt);
	
	void solve();
	void addNode(bnbNode* node);
	void delNode(bnbNode* node);
	void printSolToJson();
	void cleanup();

	bnbData* mData;
private:
	// readFile(const char* filename);
	FILE* mJson;
	FILE* mLevels;
	bnbProbData* mProbData;
	bnbBranch* mBranch;
	bnbSearch* mSearch;
	bnbSolve* mSolve;
};

class RunStats
{
public:
	RunStats(const char* jsonFile);

	void updateRunStats(vector<int>& traceSearch);
	void printRunStats();
	void cleanup();
private:
	FILE* mCycleFile;
	vector<int> lNodesByCycle;
	vector<int> eNodesByCycle;
	vector<int> gNodesByCycle;
	int mCurMaxCycle;
};