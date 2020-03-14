#include "BnB.h"

int main()
{
	const char* probPath = "D:/ResearchTests/Theory/Simulation";
	const char* outPath = "D:/ResearchTests/Theory/Simulation";
	const char* probNamePre = "inst";
	char json[255];
	char level[255];
	char runStatsFile[255];
	char prob[255];
	char probName[255];
	int numInsts = 1000;

	Mode md = DepthCont;
	TBMode tb = ARB;
	int mesrBest = 1;
	int solveMode = 1;
	double timeLimit = 3600;
	int iterLimit = 1000000;
	int lweight = 1;
	int rweight = 1;
	int depth = 11;

	options opt = options(timeLimit, iterLimit, md, tb, mesrBest, solveMode, lweight, rweight);
	for (int op = 5; op <= 5; op++)
	{
		for (int sp = 50; sp <= 50; sp+=3)
		{
			for (int bp = 7; bp <= 7; bp+=2)
			{
				if (md == DepthCont)
				{
					sprintf(json, "%s/D%dB0%d/D%dP00%dS0%dB0%d/outDepthContArb", outPath, depth + 1, bp, depth + 1, op, sp, bp);
					sprintf(runStatsFile, "%s/D%dB0%d/D%dP00%dS0%dB0%d/runStatsContArb", outPath, depth + 1, bp, depth + 1, op, sp, bp);
				}
				else if (md == BFS)
				{
					sprintf(json, "%s/D%dB0%d/D%dP00%dS0%dB0%d/outBFSArb", outPath, depth + 1, bp, depth + 1, op, sp, bp);
					sprintf(runStatsFile, "%s/D%dB0%d/D%dP00%dS0%dB0%d/runStatsBFSArb", outPath, depth + 1, bp, depth + 1, op, sp, bp);
				}
				//auto f = fopen(json, "w");
				//fclose(f);
				RunStats* run = new RunStats(runStatsFile);
				for (int i = 0; i < numInsts; i++) 
				{
					sprintf(probName, "%s%d", probNamePre, i);
					sprintf(prob, "%s/D%dB0%d/D%dP00%dS0%dB0%d/%s", probPath, depth + 1, bp, depth + 1, op, sp, bp, probName);
					sprintf(level, "%s/D%dB0%d/D%dP00%dS0%dB0%d/%d", probPath, depth + 1, bp, depth + 1, op, sp, bp, i);

					srand(static_cast<long unsigned int>(myclock::now().time_since_epoch().count()));
					BranchAndBound* test = new BranchAndBound(prob, json, nullptr, opt);
					test->solve();
					run->updateRunStats(test->mData->getTraceSearch());
					//test->printSolToJson();
					test->cleanup();
					delete test;
				}
				run->printRunStats();
				run->cleanup();
				delete run;
			}
		}
	}

	//cin.get();
	return 0;
}