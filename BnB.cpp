#include "BnB.h"

void bnbProbData::buildTree(const char* filename)
{
	ifstream inFile(filename);
	int id, parentID, branchID, solStatus, numNodes;
	//int depth;
	double sol;
	double p, q, r;
	getline(inFile, mInstName);
	inFile >> Depth;
	for (int i = 0; i <= Depth; i++) {
		inFile >> p >> q >> r;
	}
	while (inFile >> id >> parentID >> branchID >> solStatus >> sol) {
		mNodes.push_back(new node(id, parentID, branchID, solStatus, sol));
		if (id != 0) {
			if (mTree.count(parentID) == 0) {
				mTree[parentID] = { mNodes[id] };
			}
			else {
				mTree[parentID].push_back(mNodes[id]);
			}
		}
	}
}

vector<node*> bnbProbData::getChildNodes(int parentID)
{
	vector<node*> re;
	if (mTree.count(parentID) != 0) {
		re = mTree[parentID];
	}
	return re;
}

bnbProbData::~bnbProbData()
{
	mTree.clear();
	for (auto nodeptr : mNodes) {
		delete nodeptr;
	}
	mNodes.clear();
}

void bnbData::initialize(options opt)
{
	mNumUnexpNodes = 0;
	mNumNodes = 0;
	mNumIter = 0;
	mCurCycle = 0;
	globLB = MinInt;
	globUB = MaxInt;

	// store options info
	mSearchMode = opt.md;
	mContSearchMode = opt.tb;
	mMesrBest = opt.mesrBest;
	mSolveMode = opt.solveMode;
	mTimeLimit = opt.timeLimit;
	mIterLimit = opt.iterLimit;
	mLweight = opt.lweight;
	mRweight = opt.rweight;
}

void bnbData::addNode(bnbNode* node)
{
	mUnexpNodesHeap.insert({ node->mNodeID, node });
	mNumUnexpNodes++;

	if (mLbCount.count(node->mLBound) == 0) {
		mLbCount[node->mLBound] = 1;
	} else {
		mLbCount[node->mLBound]++;
	}
}

void bnbData::delNode(bnbNode* node)
{
	mUnexpNodesHeap.erase(node->mNodeID);
	mNumUnexpNodes--;

	if (mLbCount.count(node->mLBound) == 1 && mLbCount[node->mLBound] > 0) {
		mLbCount[node->mLBound]--;
		if (mLbCount[node->mLBound] == 0)
			mLbCount.erase(node->mLBound);
	}
	updateLb();
}

double bnbData::updateLb()
{
	if (!mLbCount.empty()) {
		auto lbIter = mLbCount.begin();
		double lb = lbIter->first;
		if (lbIter->second <= 0)
			throw ERROR << "Retrive global lower bound error.";
		if (lb <= globUB)
			globLB = lb;
		return globLB;
	} else {
		return MinInt;
	}
}

double bnbData::updateUb(double sol)
{
	if (sol < globUB)
		globUB = sol;
	return globUB;
}

void bnbData::updateTraceSearch(double lb)
{
	if (globUB == 20)
		return;
	if (lb == 20)
		mTraceSearch.push_back(1);
	else if (lb < 20)
		mTraceSearch.push_back(0);
	else
		mTraceSearch.push_back(2);
	mCurCycle++;
}

bnbNode* bnbData::getNode(int nodeID)
{
	return mUnexpNodesHeap[nodeID];
}

void bnbData::dumpAllNodes()
{
	for (auto bnbNodeSet : mUnexpNodesHeap) {
		delete bnbNodeSet.second;
	}
}

bnbData::~bnbData()
{
	if (!mUnexpNodesHeap.empty())
		dumpAllNodes();
}

void bnbSearch::initialize(bnbData* data)
{
	mMode = data->getSearchMode();
	mTBMode = data->getContSearchMode();
	mMesrBest = data->getMesrBest();
	mLweight = data->getLweight();
	mRweight = data->getRweight();
	mData = data;

	// initialize vectors
	mNodeInLevels.resize(data->mBnBProbData->Depth + 1);
	for (int i = 0; i < mNodeInLevels.size(); i++)
	{
		mNodeInLevels[i].resize(3, 0);
	}
}

int bnbSearch::getNextNode()
{
	int id = -1;
	switch (mMode) {
	case DFS:
		id = getNextNodeDFS();
		break;
	default:
		id = getNextNodeCBFS();
		break;
	}
	return id;
}

int bnbSearch::getNextNodeDFS()
{
	if (mNodesStack.empty()) 
		return -1;
	bnbNode* node = mNodesStack.top();
	mNodesStack.pop();
	return node->mNodeID;
}

int bnbSearch::getNextNodeCBFS()
{
	int n;
	bnbNode* node;
	vector<bnbNode*> candidates;
	mCurContour++;
	if (mCurContour == mContours.end())
		mCurContour = mContours.begin();

	double search;
	search = ((mCurContour->second).begin())->first;
	auto range = mCurContour->second.equal_range(search);

	switch (mTBMode) {
	case FIFO:
		// FIFO tie breaking
		node = (range.first)->second;
		break;
	case LIFO:
		// LIFO tie breaking
		node = (--range.second)->second;
		break;
	case ARB:
		// Arbitrary Tie Breaking
		for (auto iter = range.first; iter != range.second; iter++)
			candidates.push_back(iter->second);
		n = rand() % candidates.size();
		node = candidates[n];
		break;
	default:
		// FIFO tie breaking
		node = (range.first)->second;
		break;
	}
	return node->mNodeID;
}

void bnbSearch::addNode(bnbNode* node)
{
	double best;
	switch (mMesrBest) {
	case 1:
		best = node->mLBound;
		break;
	default:
		best = node->mLBound;
		break;
	}
	if (mMode == DFS) {
		mNodesStack.push(node);
		node->mContour = -1;
	} else {
		node->mContour = calContour(node);
		mContours[node->mContour].insert({ best, node });
	}
	if (node->mLBound < 20)
		mNodeInLevels[node->mDepth][0]++;
	else if (node->mLBound > 20)
		mNodeInLevels[node->mDepth][2]++;
	else
		mNodeInLevels[node->mDepth][1]++;
}

void bnbSearch::delNode(bnbNode* node)
{
	double search;
	switch (mMesrBest) {
	case 1: 
		search = node->mLBound;
		break;
	default:
		search = node->mLBound;
		break;
	}
	if (mMode != DFS) {
		int contID = node->mContour;
		auto range = mContours[contID].equal_range(search);

		for (auto iter = range.first; iter != range.second; iter++) {
			if ((iter->second)->mNodeID == node->mNodeID) {
				mContours[contID].erase(iter);
				break;
			}
		}

		if (mContours[contID].empty()) {
			if (mCurContour == mContours.find(contID)) {
				if (mCurContour == mContours.begin())
					mCurContour = mContours.end();
				mCurContour--;
			}
			mContours.erase(contID);
		}
	}
	if (node->mLBound < 20)
		mNodeInLevels[node->mDepth][0]--;
	else if (node->mLBound > 20)
		mNodeInLevels[node->mDepth][2]--;
	else
		mNodeInLevels[node->mDepth][1]--;
}

int bnbSearch::calContour(bnbNode* node)
{
	int contour;
	switch (mMode) {
	case BFS:
		contour = 0;
		break;
	case DepthCont:
		contour = node->mDepth;
		break;
	case WeightCont:
		contour = node->mNumLBrch * mLweight + node->mNumRBrch * mRweight;
		break;
	default:
		contour = 0;
		break;
	}
	return contour;
}

void bnbBranch::initialize(bnbData* data)
{
	numNodes = 0;
	mSolveMode = data->getSolveMode();
	mData = data;
}

int bnbBranch::createNewNodes(bnbNode* curNode)
{
	int count = 0;
	switch (mSolveMode) {
	case 1:
		while (!mNewNodes.empty())
			mNewNodes.pop();
		auto allChildNodes = mData->mBnBProbData->getChildNodes(curNode->mNodeID);
		if (!allChildNodes.empty()) {
			for (auto newNode : allChildNodes) {
				bnbNode* newBnBNode = new bnbNode(curNode, newNode);
				mNewNodes.push(newBnBNode);
				count++;
			}
		}
		break;
	}
	numNodes += count;
	mData->mNumNodes = numNodes;
	return count;
}

bnbNode* bnbBranch::getNewNode()
{
	if (mNewNodes.empty())
		return nullptr;
	bnbNode* node = mNewNodes.front();
	mNewNodes.pop();
	return node;
}

void bnbSolve::initialize(bnbData* data)
{
	mSolveMode = data->getSolveMode();
	mData = data;
}

int bnbSolve::solveNode(bnbNode* node)
{
	int status = -1;
	switch (mSolveMode) {
	case 1:
		status = node->mSolStatus;
		break;
	}
	return status;
}

bnbNode::bnbNode(bnbNode* parentNode)
{
	mNodeID = parentNode->mNodeID; mParentID = parentNode->mNodeID;
	mDepth = parentNode->mDepth + 1;
	mNumLBrch = parentNode->mNumLBrch;
	mNumRBrch = parentNode->mNumRBrch;
	mLBound = parentNode->mRexSol;
	mData = parentNode->mData;
	mContour = -1;
}

bnbNode::bnbNode(bnbData* data, node* nd)
{
	mNodeID = nd->id; mParentID = -1;
	mDepth = 0;
	mNumLBrch = mNumRBrch = 0;
	mLBound = MinInt;
	mRexSol = nd->sol;
	mSolStatus = nd->solStatus;
	mContour = -1;
	mData = data;
}

bnbNode::bnbNode(bnbNode* parentNode, node* nd)
{
	mNodeID = nd->id; mParentID = parentNode->mNodeID;
	mDepth = parentNode->mDepth + 1;
	if (nd->branchID == 1) {
		mNumLBrch = parentNode->mNumLBrch + 1;
		mNumRBrch = parentNode->mNumRBrch;
	} else {
		mNumLBrch = parentNode->mNumLBrch;
		mNumRBrch = parentNode->mNumRBrch + 1;
	}
	mLBound = parentNode->mRexSol;
	mRexSol = nd->sol;
	mSolStatus = nd->solStatus;
	mContour = -1;
	mData = parentNode->mData;
}

bnbNode::~bnbNode()
{
	mData->delNode(this);
	mData->mBnBSearch->delNode(this);
}

BranchAndBound::BranchAndBound(const char* filename, char* jsonFile, char* levelFile, options opt)
{
	mProbData = new bnbProbData();
	mProbData->buildTree(filename);
	if (jsonFile != nullptr)
		mJson = fopen(jsonFile, "a");
	if (levelFile != nullptr)
		mLevels = fopen(levelFile, "w");
	mData = new bnbData();
	mData->initialize(opt);
	mData->mBnBProbData = mProbData;
	mBranch = new bnbBranch();
	mBranch->initialize(mData);
	mSearch = new bnbSearch();
	mSearch->initialize(mData);
	mSolve = new bnbSolve();
	mSolve->initialize(mData);

	mData->mBnBBranch = mBranch;
	mData->mBnBSearch = mSearch;
	mData->mBnBSolve = mSolve;
}

void BranchAndBound::addNode(bnbNode* node)
{
	mData->addNode(node);
	mSearch->addNode(node);
}

void BranchAndBound::delNode(bnbNode* node)
{
	mData->addNode(node);
	mSearch->addNode(node);
}

void BranchAndBound::solve()
{
	int status, count;
	bnbNode* curNode = new bnbNode(mData, mProbData->getNode(0));
	addNode(curNode);
	mSearch->setIterator();
	while (mData->hasUnexpNodes()) {
		if (mData->mNumIter % 10 == 0)
			printf("NID    Iter    NLB     rexSol     gLB    gUB    contr  \n");
		curNode = mData->getNode(mSearch->getNextNode());
		
		if (curNode->mLBound < mData->globUB) {
			if (mLevels != nullptr && curNode->mDepth == mProbData->Depth)
			{
				fprintf(mLevels, "%d %d\n", mData->mBnBSearch->getNodeInLevel(curNode->mDepth, 0),
					mData->mBnBSearch->getNodeInLevel(curNode->mDepth, 2));
			}
			status = mSolve->solveNode(curNode);
			switch (status) {
			case 1:
				// optimal, update global upper bound if possible
				mData->updateUb(curNode->mRexSol);
				printf("%d    %d    %.2f    %.2f    %.2f    %.2f    %d	\n",
					curNode->mNodeID, mData->mNumIter, curNode->mLBound,
					curNode->mRexSol, mData->globLB, mData->globUB, curNode->mContour);
				printf("At ite %d, feasible solution %.2f found.\n", mData->mNumIter, curNode->mRexSol);
				mData->updateTraceSearch(curNode->mLBound);
				break;
			case 2:
				// feasible, branching
				count = mBranch->createNewNodes(curNode);
				printf("%d    %d    %.2f    %.2f    %.2f    %.2f    %d	\n",
					curNode->mNodeID, mData->mNumIter, curNode->mLBound,
					curNode->mRexSol, mData->globLB, mData->globUB, curNode->mContour);
				while (count > 0) {
					bnbNode* newNode = mBranch->getNewNode();
					if (newNode != nullptr) {
						addNode(newNode);
						printf("At iter %d, node %d is created.\n", mData->mNumIter, newNode->mNodeID);
					}
					count--;
				}
				break;
			case 3:
				// infeasible, do nothing
				printf("%d    %d    %.2f    ---    %.2f    %.2f    %d	\n",
					curNode->mNodeID, mData->mNumIter, curNode->mLBound,
					mData->globLB, mData->globUB, curNode->mContour);
				printf("At iter %d, node %d is infeasible.", mData->mNumIter, curNode->mNodeID);
				break;
			case -1:
				printf("Solve node %d error.", curNode->mNodeID);
				throw ERROR << "Node cannot solve.";
				break;
			}
			mData->mNumIter++;
		} else {
			printf("%d    %d    %.2f    ---    %.2f    %.2f    %d	\n",
				curNode->mNodeID, mData->mNumIter, curNode->mLBound,
				mData->globLB, mData->globUB, curNode->mContour);
			printf("At iter %d, node %d is in pruned by LB.\n", mData->mNumIter, curNode->mNodeID);
		}
		delete curNode;
	}
	
	if (mData->globUB < MaxInt) {
		printf("The optimal solution is: %.2f.\n", mData->globUB);
	}
}

void BranchAndBound::printSolToJson()
{
	if (mJson == nullptr) return;
	fprintf(mJson, "{");

	fprintf(mJson, "\"name\": \"%s\", \"iter\": %d, \"obj_value\": %.2f, \"obj_lb\": %.2f, \"total_nodes\": %d",
		mProbData->mInstName.c_str(), mData->mNumIter, mData->globUB, mData->globLB, mData->mNumNodes);

	fprintf(mJson, "}\n");
}

void BranchAndBound::cleanup()
{
	delete mSolve;
	delete mBranch;
	delete mSearch;
	delete mData;
	delete mProbData;

	if (mJson != nullptr)
		fclose(mJson);
	if (mLevels != nullptr)
		fclose(mLevels);
}

RunStats::RunStats(const char* jsonFile)
{
	if (jsonFile != nullptr)
		mCycleFile = fopen(jsonFile, "w");
	mCurMaxCycle == 0;
}

void RunStats::updateRunStats(vector<int>& traceSearch)
{
	int n = traceSearch.size();
	if (n > mCurMaxCycle)
	{
		mCurMaxCycle = n;
		lNodesByCycle.resize(n, 0);
		eNodesByCycle.resize(n, 0);
		gNodesByCycle.resize(n, 0);
	}
	for (int i = 0; i < n; i++)
	{
		switch (traceSearch[i])
		{
		case 0:
			lNodesByCycle[i]++;
			break;
		case 1:
			eNodesByCycle[i]++;
			break;
		case 2:
			gNodesByCycle[i]++;
		}
	}
}

void RunStats::printRunStats()
{
	if (mCycleFile == nullptr) return;
	for (int i = 0; i < mCurMaxCycle; i++)
	{
		fprintf(mCycleFile, "%d %d %d\n", lNodesByCycle[i], eNodesByCycle[i], gNodesByCycle[i]);
	}
}

void RunStats::cleanup()
{
	fclose(mCycleFile);
}