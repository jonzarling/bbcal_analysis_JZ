#ifndef _DEventProcessor_BBCAL_tree_
#define _DEventProcessor_BBCAL_tree_

#include <JANA/JEventProcessor.h>
using namespace jana;


#include "ANALYSIS/DTreeInterface.h"
#include <TH2.h>

class DEventProcessor_BBCAL_tree:public JEventProcessor{
	public:
		DEventProcessor_BBCAL_tree(){};
		~DEventProcessor_BBCAL_tree(){};
		const char* className(void){return "DEventProcessor_BBCAL_tree";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(JEventLoop *eventLoop, int32_t  runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
		
		// Convert crate, slot, channel mapping goes to internal numbering
		pair<int,int> getBBCALIndex(int crate, int slot, int channel);

		// Use cases, other inputs
		Bool_t isCosmicSetup = false; 	// Flag to switch to cosmic configuration
		Bool_t isRawMode     = false; 	// Flag to switch to raw waveform input configuration
		int nRawArrSize = 300;          // Max array size
		int NChannels = 16; // Number of channels per side, 4x4
		// For reading in gain constants
		void ParseGainFile(string fname);
		FILE* myfile;
		Double_t C_N[16]; Double_t C_S[16];  		

		// Hall D classes for ROOT I/O
		DTreeInterface* dTreeInterface;
		DTreeFillData   dTreeFillData;

};

#endif // _DEventProcessor_BBCAL_tree_

